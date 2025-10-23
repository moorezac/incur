"""
Batch TIFF Frame Extraction and Cellpose Processing

This script processes a folder of multi-frame TIFF files by:
1. Extracting frames from each TIFF to a temporary directory
2. Running Cellpose segmentation on the extracted frames
3. Saving Cellpose outputs for each original TIFF file

Usage:
    python batch_cellpose.py input_folder/ [--output output_folder/]
"""

import argparse
import gc
import os
import sys
import tempfile
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from pathlib import Path
from datetime import datetime
import io
import warnings
import logging

import numpy as np
import tifffile
import torch
from tqdm import tqdm


@contextmanager
def suppress_stdout():
    with open(os.devnull, "w") as devnull:
        old_stdout = sys.stdout
        sys.stdout = devnull
        try:
            yield
        finally:
            sys.stdout = old_stdout


class DualOutput:
    """Write to both file and terminal (for tqdm only)."""
    def __init__(self, file, original_stderr):
        self.file = file
        self.terminal = original_stderr
        self._buffer = ""
        
    def write(self, message):
        # Check if this is from tqdm (crude but effective)
        if any(marker in str(message) for marker in ['\r', '\x1b[', '%|', '█', '░']):
            self.terminal.write(message)
            self.terminal.flush()
        else:
            # Buffer and write complete lines to file
            self._buffer += str(message)
            if '\n' in self._buffer:
                lines = self._buffer.split('\n')
                for line in lines[:-1]:
                    if line.strip():  # Only write non-empty lines
                        self.file.write(line + '\n')
                self._buffer = lines[-1]
            self.file.flush()
            
    def flush(self):
        if self._buffer.strip():
            self.file.write(self._buffer + '\n')
            self._buffer = ""
        self.file.flush()
        self.terminal.flush()
        
    def fileno(self):
        return self.terminal.fileno()
        
    def isatty(self):
        return self.terminal.isatty()
        
    def __getattr__(self, attr):
        # Forward any other attributes to terminal
        return getattr(self.terminal, attr)


with suppress_stdout():
    from cellpose import io, models, transforms
    
# Set environment variables to reduce cellpose verbosity
os.environ['CELLPOSE_VERBOSE'] = '0'
os.environ['CELLPOSE_LOG_LEVEL'] = 'ERROR'


def extract_frames_to_temp(tiff_path: Path, temp_dir: Path, error_log_path: Path | None = None) -> list[Path]:
    """
    Extract all frames from a multi-frame TIFF to a temporary directory.
    Each extracted frame is saved as a multi-channel TIFF (shape: Y, X, C).

    Handles arbitrary TIFF axis orders (e.g., TYXC, TCYX, CYXT, ZYXC, etc.).
    """
    frame_paths = []

    try:
        with tifffile.TiffFile(tiff_path) as tif:
            series = tif.series[0]
            axes = series.axes
            data = series.asarray()
            ndim = data.ndim

            # Identify axis indices
            axis_map = {ax: i for i, ax in enumerate(axes)}

            # Find the "frame" dimension (prefer T, else Z, else default 0)
            frame_dim = axis_map.get("T", axis_map.get("Z", 0))

            # Find the channel dimension (default to None if missing)
            channel_dim = axis_map.get("C", None)

            # Move frame dimension to front for easy iteration
            if frame_dim != 0:
                data = np.moveaxis(data, frame_dim, 0)

            # Handle single-frame files
            if data.ndim == 3 and "C" in axes and data.shape[channel_dim] < 10:
                # Single timepoint but multi-channel — wrap in new axis
                data = data[np.newaxis, ...]

            # Iterate over frames (now always data[frame_i, ...])
            for i, frame in enumerate(data):
                # Move channels to last dimension if not already
                if frame.ndim == 3 and channel_dim is not None:
                    # Ensure C is last
                    if np.argmin(frame.shape) == channel_dim:
                        frame = np.moveaxis(frame, channel_dim, -1)
                elif frame.ndim == 2:
                    # Add singleton channel
                    frame = frame[..., np.newaxis]

                frame_path = temp_dir / f"frame_{i:04d}.tiff"
                tifffile.imwrite(frame_path, frame, imagej=True)
                frame_paths.append(frame_path)

        return frame_paths

    except Exception as e:
        error_msg = f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] Error extracting frames from {tiff_path.name}: {e}\n"
        if error_log_path:
            with open(error_log_path, "a") as f:
                f.write(error_msg)
        return []


def run_cellpose(
    frame_dir: Path,
    output_dir: Path | None,
    tiff_name: str,
    original_tiff_path: Path,
    channels: list[int] = [0],
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    niter: int | None = None,
    tile_norm_blocksize: int | None = None,
    batch_size: int = 8,
    suppress_diameter_errors: bool = False,
    error_log_path: Path | None = None,
):
    """
    Run Cellpose on extracted frames.
    """
    def log_error(msg: str):
        if error_log_path:
            with open(error_log_path, "a") as f:
                f.write(f"[{datetime.now().strftime('%Y-%m-%d %H:%M:%S')}] {msg}\n")
    
    # Prepare output directory
    tiff_output_dir = None
    if output_dir is not None:
        tiff_output_dir = output_dir / tiff_name
        tiff_output_dir.mkdir(parents=True, exist_ok=True)

    frame_files = sorted(frame_dir.glob("*.tiff"))
    if not frame_files:
        log_error(f"No frames found in {frame_dir}")
        return

    try:
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        model = models.CellposeModel(gpu=True)
        all_masks = []
        failed_frames = []

        for frame_path in tqdm(frame_files, desc=f"  {tiff_name}", leave=False):
            try:
                img = io.imread(frame_path)
                if img.ndim == 2:
                    img = np.expand_dims(img, axis=2)

                max_channel = img.shape[-1] - 1
                for ch in channels:
                    if ch > max_channel:
                        raise ValueError(
                            f"Channel index {ch} exceeds available channels (0-{max_channel})"
                        )

                img = transforms.move_min_dim(img, force=True)
                img = transforms.normalize_img(img, normalize=True, percentile=(1, 99))
                img_processed = img[..., channels]

                if img_processed.ndim == 2:
                    img_processed = img_processed[..., None]

                normalize_dict = {}
                if tile_norm_blocksize is not None:
                    normalize_dict["tile_norm_blocksize"] = tile_norm_blocksize

                eval_params = {
                    "batch_size": batch_size,
                    "flow_threshold": flow_threshold,
                    "cellprob_threshold": cellprob_threshold,
                }
                if diameter is not None:
                    eval_params["diameter"] = diameter
                if normalize_dict:
                    eval_params["normalize"] = normalize_dict
                if niter is not None:
                    eval_params["niter"] = niter

                try:
                    masks, flows, styles = model.eval(img_processed, **eval_params)
                except Exception as e:
                    # If diameter causes an error, try without it
                    if diameter is not None and "diameter" in str(e):
                        if not suppress_diameter_errors:
                            log_error(f"Diameter {diameter} caused error for {frame_path.name}, using automatic detection")
                        eval_params.pop("diameter", None)
                        masks, flows, styles = model.eval(img_processed, **eval_params)
                    else:
                        raise e
                
                all_masks.append(masks)

                if output_dir is not None:
                    mask_path = tiff_output_dir / f"{frame_path.stem}_mask.tiff"
                    io.imsave(mask_path, masks)

                    roi_path = tiff_output_dir / f"{frame_path.stem}_rois.zip"
                    io.save_rois(masks, roi_path)
                    
            except Exception as e:
                log_error(f"Error processing frame {frame_path.name} in {tiff_name}: {e}")
                failed_frames.append(frame_path.name)
                # Append a blank mask to maintain frame alignment
                if len(all_masks) > 0:
                    blank_mask = np.zeros_like(all_masks[0])
                else:
                    # Try to create a reasonable blank mask
                    try:
                        img = io.imread(frame_path)
                        blank_mask = np.zeros(img.shape[:2], dtype=np.uint8)
                    except:
                        blank_mask = np.zeros((512, 512), dtype=np.uint8)  # Default size
                all_masks.append(blank_mask)

        if all_masks:
            # Combine masks
            all_masks_concat = np.stack(all_masks, axis=0)
            all_masks_concat_path = original_tiff_path.parent / f"{tiff_name}_masks.tiff"
            tifffile.imwrite(all_masks_concat_path, all_masks_concat.astype("uint8"))
            
            if failed_frames:
                log_error(f"Completed {tiff_name} with {len(failed_frames)} failed frames: {', '.join(failed_frames)}")

        del model

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        gc.collect()

        if output_dir is not None and not failed_frames:
            pass  # Silent success

    except Exception as e:
        log_error(f"Critical error running Cellpose on {tiff_name}: {e}")
        return


def process_tiff_folder(
    input_folder: Path,
    output_folder: Path | None = None,
    channels: list[int] = [0],
    diameter: float | None = None,
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    niter: int | None = None,
    tile_norm_blocksize: int | None = None,
    batch_size: int = 8,
    skip_existing: bool = True,
    suppress_diameter_errors: bool = False,
    error_log: str | None = None,
    quiet: bool = False,
):
    """
    Process all TIFF files in a folder.
    """
    input_folder = Path(input_folder)

    if output_folder is not None:
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)
        
    # Set up error logging
    error_log_path = None
    log_file = None
    original_stdout = sys.stdout
    original_stderr = sys.stderr
    
    try:
        if error_log:
            error_log_path = Path(error_log)
            # Create header for new log session
            with open(error_log_path, "a") as f:
                f.write(f"\n{'='*60}\n")
                f.write(f"Batch processing started at {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"Input folder: {input_folder}\n")
                f.write(f"{'='*60}\n\n")
            
            # Open log file for the entire session
            log_file = open(error_log_path, "a", buffering=1)
            
            # Create dual output that preserves tqdm
            dual_out = DualOutput(log_file, original_stderr)
            
            # Redirect both stdout and stderr to our dual output
            sys.stdout = dual_out
            sys.stderr = dual_out
            
            # Capture warnings
            def custom_warning(message, category=UserWarning, filename='', lineno=-1, file=None, line=None):
                log_file.write(f"[WARNING] {category.__name__}: {message}\n")
                log_file.flush()
            warnings.showwarning = custom_warning
            
            # Configure logging to go to our log file
            logging.basicConfig(
                level=logging.WARNING,
                format='[%(levelname)s] %(message)s',
                force=True,
                handlers=[logging.FileHandler(error_log_path, mode='a')]
            )
            
            # Suppress specific noisy warnings if needed
            warnings.filterwarnings('ignore', category=DeprecationWarning)
            warnings.filterwarnings('ignore', category=FutureWarning)

        # Find TIFF files
        tiff_files = list(input_folder.glob("*.tiff")) + list(input_folder.glob("*.tif"))
        tiff_files = [
            f
            for f in tiff_files
            if not f.name.endswith(("masks.tiff", "merged.tiff", "combined.tiff"))
        ]

        if not tiff_files:
            print(f"No files found in {input_folder}")
            return

        skipped = 0
        if skip_existing:
            pre_filtered = len(tiff_files)
            tiff_files = [
                f for f in tiff_files if not (f.parent / f"{f.stem}_masks.tiff").exists()
            ]
            skipped = pre_filtered - len(tiff_files)
            if skipped:
                print(f"Skipping {skipped} file(s) with existing _masks.tiff\n")

        if not tiff_files:
            print("All files already processed — nothing to do.")
            return

        print(f"Processing {len(tiff_files)} new file(s)\n")
        
        failed_files = []

        # Use tqdm 
        for tiff_path in tqdm(tiff_files, desc="Processing", total=len(tiff_files)):
            mask_path = tiff_path.parent / f"{tiff_path.stem}_masks.tiff"

            # Each file gets its own temporary workspace
            with tempfile.TemporaryDirectory() as temp_dir:
                temp_dir = Path(temp_dir)

                # Extract individual frames
                frame_paths = extract_frames_to_temp(tiff_path, temp_dir, error_log_path)
                if not frame_paths:
                    failed_files.append(tiff_path.name)
                    continue

                # Run Cellpose segmentation
                run_cellpose(
                    frame_dir=temp_dir,
                    output_dir=output_folder,
                    tiff_name=tiff_path.stem,
                    original_tiff_path=tiff_path,
                    channels=channels,
                    diameter=diameter,
                    flow_threshold=flow_threshold,
                    cellprob_threshold=cellprob_threshold,
                    niter=niter,
                    tile_norm_blocksize=tile_norm_blocksize,
                    batch_size=batch_size,
                    suppress_diameter_errors=suppress_diameter_errors,
                    error_log_path=error_log_path,
                )

        print("\nComplete!")
        
        if error_log_path and failed_files:
            print(f"\nEncountered errors with {len(failed_files)} file(s). Check {error_log_path} for details.")
        elif error_log_path:
            # Check if any errors were logged
            with open(error_log_path, "r") as f:
                content = f.read()
                if "Error" in content:
                    print(f"\nSome errors were encountered. Check {error_log_path} for details.")
                    
    finally:
        # Restore original stdout/stderr
        sys.stdout = original_stdout
        sys.stderr = original_stderr
        
        if log_file:
            log_file.close()
    
    # Print completion message to terminal (stderr) so it's always visible
    if not quiet:
        if error_log_path:
            print("\nProcessing complete! Check log for details.", file=sys.stderr)
        else:
            print("\nProcessing complete!", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description="Extract frames from .tiff/.tiff files and run segmentation",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python batch_cellpose.py tiff_folder/
  python batch_cellpose.py tiff_folder/ --output results/
  python batch_cellpose.py data/ --channels 0 1 2
  python batch_cellpose.py data/ --flow_threshold 0.6 --cellprob_threshold 0.5
  python batch_cellpose.py data/ --diameter 30
  python batch_cellpose.py data/ --error-log processing.log
  python batch_cellpose.py data/ --error-log errors.log --quiet
  python batch_cellpose.py data/ --diameter 30 --suppress-diameter-errors --error-log errors.log
        """,
    )

    parser.add_argument("input_folder", help="Folder containing TIFF files to process")
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Optional output folder for individual frame results",
    )
    parser.add_argument(
        "--channels",
        "-c",
        nargs="+",
        type=int,
        default=[0],
        help="Channel indices to use for segmentation (default: 0)",
    )
    parser.add_argument(
        "--diameter",
        "-d",
        type=float,
        default=None,
        help="Expected cell diameter in pixels (default: automatic detection)",
    )
    parser.add_argument(
        "--flow_threshold",
        "--flow",
        type=float,
        default=0.4,
        help="Flow error threshold for masks (default: 0.4)",
    )
    parser.add_argument(
        "--cellprob_threshold",
        "--cellprob",
        type=float,
        default=0.0,
        help="Cell probability threshold (default: 0.0)",
    )
    parser.add_argument(
        "--niter",
        type=int,
        default=None,
        help="Number of iterations for dynamics (default: None)",
    )
    parser.add_argument(
        "--tile_norm_blocksize",
        "--tile",
        type=int,
        default=None,
        help="Block size for tile normalization (default: None)",
    )
    parser.add_argument(
        "--batch_size",
        "-b",
        type=int,
        default=8,
        help="Batch size for processing (default: 8)",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force reprocessing even if _masks.tif files already exist (default: skip existing)",
    )
    parser.add_argument(
        "--suppress-diameter-errors",
        action="store_true",
        help="Suppress warnings when diameter parameter causes errors",
    )
    parser.add_argument(
        "--error-log",
        type=str,
        default=None,
        help="Path to save error log file (default: no error logging)",
    )
    parser.add_argument(
        "--quiet",
        "-q",
        action="store_true",
        help="Quiet mode - suppress all terminal output except progress bars",
    )
    args = parser.parse_args()

    try:
        process_tiff_folder(
            args.input_folder,
            args.output,
            channels=args.channels,
            diameter=args.diameter,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
            niter=args.niter,
            tile_norm_blocksize=args.tile_norm_blocksize,
            batch_size=args.batch_size,
            skip_existing=not args.force,
            suppress_diameter_errors=args.suppress_diameter_errors,
            error_log=args.error_log,
            quiet=args.quiet,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
