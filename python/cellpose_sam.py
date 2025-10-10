"""
Batch TIFF Frame Extraction and Cellpose Processing

This script processes a folder of multi-frame TIFF files by:
1. Extracting frames from each TIFF to a temporary directory
2. Running Cellpose segmentation on the extracted frames
3. Saving Cellpose outputs for each original TIFF file

Requirements:
    pip install tifffile cellpose tqdm

Usage:
    python batch_cellpose.py input_folder/ [--output output_folder/]
"""

import argparse
import gc
import os
import sys
import tempfile
from contextlib import contextmanager
from pathlib import Path

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


with suppress_stdout():
    from cellpose import io, models, transforms


def extract_frames_to_temp(tiff_path: Path, temp_dir: Path) -> list[Path]:
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

        # tqdm.write(
        #     f"  Extracted {len(frame_paths)} multi-channel frames from {tiff_path.name}"
        # )
        return frame_paths

    except Exception as e:
        tqdm.write(f"  Error extracting frames from {tiff_path.name}: {e}")
        return []


def run_cellpose(
    frame_dir: Path,
    output_dir: Path | None,
    tiff_name: str,
    original_tiff_path: Path,
    channels: list[int] = [0],
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    niter: int | None = None,
    tile_norm_blocksize: int | None = None,
    batch_size: int = 8,
):
    """
    Run Cellpose on extracted frames.
    """
    # Prepare output directory
    tiff_output_dir = None
    if output_dir is not None:
        tiff_output_dir = output_dir / tiff_name
        tiff_output_dir.mkdir(parents=True, exist_ok=True)

    frame_files = sorted(frame_dir.glob("*.tiff"))
    if not frame_files:
        tqdm.write(f"  No frames found in {frame_dir}")
        return

    # tqdm.write(f"  Running Cellpose on {len(frame_files)} frames...")

    try:
        if torch.cuda.is_available():
            torch.cuda.empty_cache()

        model = models.CellposeModel(gpu=True)
        all_masks = []

        for frame_path in tqdm(frame_files, desc=f"  {tiff_name}", leave=False):
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
            if normalize_dict:
                eval_params["normalize"] = normalize_dict
            if niter is not None:
                eval_params["niter"] = niter

            masks, flows, styles = model.eval(img_processed, **eval_params)
            all_masks.append(masks)

            if output_dir is not None:
                mask_path = tiff_output_dir / f"{frame_path.stem}_mask.tiff"
                io.imsave(mask_path, masks)

                roi_path = tiff_output_dir / f"{frame_path.stem}_rois.zip"
                io.save_rois(masks, roi_path)

        # Combine masks
        all_masks_concat = np.stack(all_masks, axis=0)
        all_masks_concat_path = original_tiff_path.parent / f"{tiff_name}_masks.tiff"
        tifffile.imwrite(all_masks_concat_path, all_masks_concat.astype("uint8"))
        # tqdm.write(f"  Combined masks saved to: {all_masks_concat_path}")

        # Merge masks and original
        # orig = io.imread(original_tiff_path)
        # mask = io.imread(all_masks_concat_path)
        # mask = mask.transpose(1, 0, 2)

        # if orig.ndim == 3:
        #    orig = np.expand_dims(orig, axis=1)
        # mask = np.expand_dims(mask, axis=1)
        # combined = np.concatenate([mask, orig], axis=1)

        # combined_path = original_tiff_path.parent / f"{tiff_name}_combined.tiff"
        # tifffile.imwrite(combined_path, combined, imagej=True, metadata={"axes": "TCYX"})
        # tqdm.write(f"  Combined image saved to: {combined_path}")

        del model

        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        gc.collect()

        if output_dir is not None:
            tqdm.write(f"  Individual frame results saved to {tiff_output_dir}")

    except Exception as e:
        tqdm.write(f"  Error running Cellpose on {tiff_name}: {e}")


def process_tiff_folder(
    input_folder: Path,
    output_folder: Path | None = None,
    channels: list[int] = [0],
    flow_threshold: float = 0.4,
    cellprob_threshold: float = 0.0,
    niter: int | None = None,
    tile_norm_blocksize: int | None = None,
    batch_size: int = 8,
    skip_existing: bool = True,
):
    """
    Process all TIFF files in a folder.
    """
    input_folder = Path(input_folder)

    if output_folder is not None:
        output_folder = Path(output_folder)
        output_folder.mkdir(parents=True, exist_ok=True)

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

    for tiff_path in tqdm(tiff_files, desc="Processing", total=len(tiff_files)):
        mask_path = tiff_path.parent / f"{tiff_path.stem}_masks.tiff"

        # Each file gets its own temporary workspace
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir = Path(temp_dir)

            # Extract individual frames
            frame_paths = extract_frames_to_temp(tiff_path, temp_dir)
            if not frame_paths:
                tqdm.write(f"  Skipping {tiff_path.name} (no frames found)")
                continue

            # Run Cellpose segmentation
            run_cellpose(
                frame_dir=temp_dir,
                output_dir=output_folder,
                tiff_name=tiff_path.stem,
                original_tiff_path=tiff_path,
                channels=channels,
                flow_threshold=flow_threshold,
                cellprob_threshold=cellprob_threshold,
                niter=niter,
                tile_norm_blocksize=tile_norm_blocksize,
                batch_size=batch_size,
            )

    print("\nComplete!")


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
    args = parser.parse_args()

    try:
        process_tiff_folder(
            args.input_folder,
            args.output,
            channels=args.channels,
            flow_threshold=args.flow_threshold,
            cellprob_threshold=args.cellprob_threshold,
            niter=args.niter,
            tile_norm_blocksize=args.tile_norm_blocksize,
            batch_size=args.batch_size,
            skip_existing=not args.force,
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
