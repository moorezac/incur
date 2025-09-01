import os
import sys
import pandas as pd
import numpy as np
from typing import List, Dict, Tuple, Optional, Union

import tifffile
from tifffile import imread, imwrite, TiffFile

import xml
import xml.etree.ElementTree as ET

from datetime import datetime

def check_folders_contain_same_names(directory_paths):
    """
    Check if all folders contain exactly the same set of names (files and subfolders).
    
    Args:
    folder_paths (list): A list of folder paths to compare
    
    Returns:
    tuple: (boolean indicating if all folders have the same names, set of common names)
    """
    if not directory_paths:
        return False, set()
    
    # Get the list of names (files and subfolders) for the first folder
    first_directory_names = set(os.listdir(directory_paths[0]))
    
    # Compare names with each subsequent folder
    for directory in directory_paths[1:]:
        # Get names for the current folder
        current_directory_names = set(os.listdir(directory))
        
        # If the sets of names are different, return False
        if current_directory_names != first_directory_names:
            return False, first_directory_names
    
    # If we've made it through all folders, they have the same names
    return True, first_directory_names

def extract_metadata_from_multipage_tiff(tiff):
    """
    Extract metadata from each page of a multi-page TIFF file
    """
    try:
        # List to store metadata for each page
        page_metadata_list = []
        
        # Iterate through each page in the TIFF
        for page_index, page in enumerate(tiff.pages):
            # Base metadata dictionary for this page
            metadata_dict = {
                'page_index': page_index,
                'image_width': page.shape[1] if len(page.shape) > 1 else None,
                'image_height': page.shape[0] if len(page.shape) > 0 else None,
                'dtype': str(page.dtype)
            }
            
            # Try to extract XML metadata from page tags
            for tag in page.tags.values():
                if tag.name == 'ImageDescription':
                    try:
                        # Try to parse XML metadata
                        xml_str = tag.value
                        root = ET.fromstring(xml_str)
                        
                        # Extract props from MetaData section
                        for prop in root.findall('.//prop'):
                            prop_name = prop.get('id', 'Unknown')
                            prop_value = prop.get('value', 'Unknown')
                            prop_type = prop.get('type', 'string')
                            
                            # Convert value based on type
                            if prop_type == 'float':
                                try:
                                    prop_value = float(prop_value)
                                except:
                                    pass
                            elif prop_type == 'int':
                                try:
                                    prop_value = int(prop_value)
                                except:
                                    pass
                            elif prop_type == 'bool':
                                prop_value = prop_value.lower() == 'true'
                            
                            metadata_dict[prop_name] = prop_value
                        
                        # Extract custom props
                        for prop in root.findall('.//custom-prop'):
                            prop_name = prop.get('id', 'Unknown')
                            prop_value = prop.get('value', 'Unknown')
                            metadata_dict[prop_name] = prop_value
                    
                    except ET.ParseError:
                        # If XML parsing fails, store raw description
                        metadata_dict['image_description'] = xml_str
                    except Exception as e:
                        metadata_dict['metadata_extraction_error'] = str(e)
            
            # Add timing information if possible
            if hasattr(page, 'time'):
                metadata_dict['page_time'] = page.time
            
            # Try to extract acquisition time from metadata
            if 'acquisition-time-local' in metadata_dict:
                metadata_dict['acquisition_time'] = metadata_dict['acquisition-time-local']
            
            page_metadata_list.append(metadata_dict)
        
        return page_metadata_list
    
    except Exception as e:
        return [{
            'error': str(e)
        }]

def convert_16bit_to_8bit(
    input_data: Union[str, np.ndarray, List[np.ndarray]], 
    output_path: Union[str, List[str]] = None
) -> Union[np.ndarray, List[np.ndarray]]:
    """
    Convert 16-bit image(s) to 8-bit.
    
    Args:
        input_data: Can be:
            - A string path to a single TIFF file
            - A numpy array representing a single image
            - A list/series of numpy arrays
        output_path: Optional. Can be:
            - A string path for single image output
            - A list of paths for multiple images
            - None to skip saving to disk
    
    Returns:
        Converted 8-bit image(s) as numpy array(s)
    """
    
    def process_single_image(img_array: np.ndarray) -> np.ndarray:
       
        min_val = np.min(img_array)
        max_val = np.max(img_array)
        
        # avoid division by zero
        if max_val == min_val:
            return np.zeros_like(img_array, dtype=np.uint8)
        
        img_8bit = ((img_array - min_val) / (max_val - min_val) * 255).astype(np.uint8)
        return img_8bit

    # handle different input types
    if isinstance(input_data, str):
        # Single image path
        img_16bit = np.array(Image.open(input_data))
        result = process_single_image(img_16bit)
        
        if output_path:
            Image.fromarray(result).save(output_path)
        
        return result
    
    elif isinstance(input_data, np.ndarray):
        # Single numpy array
        result = process_single_image(input_data)
        
        if output_path:
            Image.fromarray(result).save(output_path)
        
        return result
    
    elif isinstance(input_data, (list, tuple, np.ndarray)) and len(input_data) > 0:
        # list/series of numpy arrays
        results = []
        
        for i, img_array in enumerate(input_data):
            result = process_single_image(img_array)
            results.append(result)
            
            if output_path and isinstance(output_path, (list, tuple)):
                if i < len(output_path):
                    Image.fromarray(result).save(output_path[i])
        
        return results
    
    else:
        raise ValueError("invalid input type")

def process_individual_well(file_name, directory_list):
    """
    For each well, extract metadata from each channel (as multi-page TIFF file), filter on common times, merge, and save

    Args:
    file_name: Name of the file across all directories to merge. Must be consistent across directories
    directory_list (list): A list of folder paths to compare
    """

    # Read in metadata information
    tiff_list = [TiffFile(os.path.join(directory, file_name)) for directory in directory_list]

    # Extract metadata for all images across all channels
    metadata_list = [extract_metadata_from_multipage_tiff(tiff=tiff) for tiff in tiff_list]

    # Extract the datetime 
    datetime_list = [pd.DataFrame(metadata)["acquisition-time-local"] for metadata in metadata_list]

    # Find the common datetimes
    common_datetime = set(datetime_list[0]).intersection(*[set(dates) for dates in datetime_list[1:]])

    # Find which indexes to keep for each channel
    index_list = [datetime[datetime.isin(common_datetime)].index for datetime in datetime_list]

    # Read in the actual data
    tiff_np_list = [tifffile.imread(os.path.join(directory, file_name)) for directory in directory_list]

    # Filter down to the common datetime
    tiff_np_list_filtered = [np.take(a = tiff, indices = index, axis = 0) for tiff, index in zip(tiff_np_list, index_list)]

    # Make 8-bit
    tiff_np_list_filtered = [convert_16bit_to_8bit(tiff_np) for tiff_np in tiff_np_list_filtered]
    
    # Merge them into a multi-channel time series
    time_series = np.stack(tiff_np_list_filtered, axis = 1)

    # Save the multi-channel time series as a TIFF
    output_path = os.path.join(output_directory, file_name)
    tifffile.imwrite(output_path, time_series, metadata = {'axes': 'TCYX'}, imagej = True)

    # Save both full metadata and common datetime
    file_name_sans_ext = os.path.splitext(file_name)[0]
    
    metadata_df = [pd.DataFrame(metadata) for metadata in metadata_list]
    metadata_df = pd.concat(metadata_df)
    metadata_df.to_csv(os.path.join(output_directory, file_name_sans_ext + "_metadata.csv"))
    
    common_datetime_df = pd.DataFrame({'datetime': sorted(common_datetime)})
    common_datetime_df.to_csv(os.path.join(output_directory, file_name_sans_ext + "_datetime.csv"), header = False)

    print(f"processed {file_name}")

if __name__ == "__main__":
    if len(sys.argv) == 0:
        print("usage: provide a list of directories to merge, with the last being the directory to save merged images")
        sys.exit(1)

    directory_list = sys.argv[1:-1]
    output_directory = sys.argv[-1]

    # Create output directory if it doesn't exist
    os.makedirs(output_directory, exist_ok = True)

    # Check that directories all contain the same files
    all_equal, file_names = check_folders_contain_same_names(directory_paths = directory_list)

    if(all_equal == False):
        print("mismatch in names of files between directories")
        sys.exit(1)

    # Process
    for file_name in file_names:
        process_individual_well(file_name, directory_list)
