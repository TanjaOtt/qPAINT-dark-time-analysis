######################################################################################
# TdCDF Analysis Script
# @author: Soohyen Jang
# Research group of Mike Heilemann, Goethe University Frankfurt a.M.

# Input Files Required:
#     - HDF5 files (.hdf5) containing localization data after DBSCAN clustering and filtering
#     - Corresponding YAML files (.yaml) with metadata
    
# Output Files Generated:
#     1. Intermediate files (for data without dark column):
#         - *_dark.hdf5 & *_dark.yaml: Contains calculated dark times
#     2. Final results:
#         - *_TdCDF.hdf5 & *_TdCDF.yaml: Contains final analysis results
#         - *_sample_fits.png: Visualization of sample fits

# Usage:
#     1. Set your input directory in the main() function below
#     2. (Optional) Adjust analysis parameters in the Config class if needed
#     3. Run the script

# Parameters (adjustable in Config class):
#     - exposure_time: Camera exposure time in seconds (default: 0.15s)
#     - min_points_for_fit: Minimum points required for fitting (default: 6)
#     - max_td_value: Maximum allowed Td value in seconds (default: 5000s)
#     - num_sample_fits: Number of sample fits to plot (default: 10)
#     - num_histogram_bins: Number of bins for histogram (default: 100)
######################################################################################

INPUT_DIR = r"Folder path"



import os
from pathlib import Path
from typing import List, Tuple
from dataclasses import dataclass
import h5py
import yaml
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


"""
Please set input_dir where you have your HDF5 and YAML files
all the files in the directory will be processed
"""


"""
If you want to set your own parameters, please modify the Config class
"""
@dataclass
class Config:
    """Configuration parameters for the analysis"""
    input_dir: Path
    output_dir: Path = None # Make output_dir optional
    exposure_time: float = 0.15  # s
    min_points_for_fit: int = 10
    min_histogram_bins: int = 100
    max_td_value: float = 10e3  # s
    num_sample_fits: int = 10
    num_histogram_bins: int = 100
    bin_min: float = 0.0 # Minimum bin value for cumulative frequency
    bin_max: float = 5000.0 # Maximum bin value for cumulative frequency    
    bin_step: float = 10.0 # Bin step size for cumulative frequency
    
    def __post_init__(self):
        """Convert paths and create output directory"""
        self.input_dir = Path(self.input_dir)

        # If output_dir is not specified, create 'output' in input_dir
        if self.output_dir is None:
            self.output_dir = self.input_dir / 'dark_time_analysis'
        else:
            self.output_dir = Path(self.output_dir)

        # Create output directory if it doesn't exist
        self.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nOutput directory created/verified: {self.output_dir}")


def read_input_files(hdf5_file: Path, yaml_file: Path, config:Config) -> Tuple[pd.DataFrame, List]:
    """
    Read input HDF5 and YAML files and ensure 'dark' column is present
    """
    # Read HDF5 file
    with h5py.File(hdf5_file, 'r') as f:
        # Get the structured array and convert directly to DataFrame
        df = pd.DataFrame(f['locs'][()])
        
    # Read YAML file
    with open(yaml_file, 'r') as yaml_file:
        yaml_docs = list(yaml.safe_load_all(yaml_file))
    
    if 'dark' not in df.columns:
        original_name = hdf5_file.stem
        df = calculate_and_save_dark(df, yaml_docs, config.output_dir, original_name)
        print("\nDark column calculated from frame differences")
    else:
        print("\nUsing existing dark column")

    return df, yaml_docs

def calculate_and_save_dark(df: pd.DataFrame, yaml_docs: List, output_dir: Path, original_name: str) -> pd.DataFrame:
    """
    Calculate dark time for each group from 'frame' column and add as a new column 'dark'
    save intermediate file with _dark.hdf5
    Args:
        df: Input Dataframe
        outptu_dir: Directory to save the intermediate file
        original_name: Original filename without extension

    Returns:
        pd.DataFrame: Dataframe with dark time added    
    """
    
    # Ensure output directory exists
    output_dir.mkdir(parents=True, exist_ok=True)
    # Sort by group
    df = df.sort_values(by=['group', 'frame'])
    # Calculate dark time
    df['dark'] = df.groupby('group')['frame'].diff().fillna(0)
    
    # Save intermediate file
    intermediate_hdf5 = output_dir / f'{original_name}_dark.hdf5'
    intermediate_yaml = output_dir / f'{original_name}_dark.yaml'

    # Save intermediate HDF5 file
    with h5py.File(intermediate_hdf5, 'w') as f:
        dtype = [(name, df[name].dtype) for name in df.columns]
        structured_array = np.empty(len(df), dtype=dtype)
        for name in df.columns:
            structured_array[name] = df[name]
        f.create_dataset('locs', data=structured_array) 
    
    # Create footer for intermediate YAML
    footer_template = {
        'dark_time_calculation': {
            'dark_time_calculateion_method': 'frame_difference within groups',
        }
    }
    
    # Save intermediate YAML file
    output_docs = yaml_docs + [footer_template]
    with open(intermediate_yaml, 'w') as f:
        yaml.safe_dump_all(output_docs, f, default_flow_style=False)

    print(f"\nIntermediate file with calcuated dark times saved to:")
    print(f"Data: {intermediate_hdf5}")
    print(f"Data: {intermediate_yaml}")

    return df

def exp_func(x: np.ndarray, b: float, A: float, Td: float) -> np.ndarray:
    return b + A * (1 - np.exp(-x / Td))

def fit_cumulative_frequency(df: pd.DataFrame, config: Config) -> pd.DataFrame:
    """Fit cumulative frequency using configuration parameters"""
    df_with_td = df.copy()
    unique_groups = df['group'].unique()
    total_groups = len(unique_groups)
    successful_fits = 0
    failed_fits = 0

    for group in unique_groups:
        group_data = df[df['group'] == group]
        
        # Use config parameters instead of hard-coded values
        if len(group_data) < config.min_points_for_fit:
            df_with_td.loc[df_with_td['group'] == group, 'Td'] = np.nan
            failed_fits += 1
            continue
            
        dark_values = group_data['dark'].values * config.exposure_time

        # Use fixed binning strategy with parameters from Config
        hist, bin_edges = np.histogram(
            dark_values, 
            bins=np.arange(config.bin_min, config.bin_max, config.bin_step), 
            density=True
        )
        
        # Calculate bin centers
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

        # Calculate cumulative frequency 
        cum_freq = np.cumsum(hist) * np.diff(bin_edges)
        
        try:
            # Initial parameter guesses
            p0 = [0.01, 0.3, 100]
            popt, _ = curve_fit(exp_func, bin_centers, cum_freq, p0=p0)
            
            # Add Td to dataframe for this group
            df_with_td.loc[df_with_td['group'] == group, 'Td'] = popt[2]
            successful_fits += 1
            #print(f"Successfully fit group {group}, Td = {Td:.2f}")
            
        except (RuntimeError, ValueError):
            #print(f"Warning: Could not fit exponential function for group {group}: {str(e)}")
            df_with_td.loc[df_with_td['group'] == group, 'Td'] = np.nan
            failed_fits += 1

    # print summary of fits
    print("\nFitting Summary:")
    print(f"Total groups: {total_groups}")
    print(f"Successful fits: {successful_fits} ({successful_fits/total_groups*100:.1f}%)")
    print(f"Failed fits: {failed_fits} ({failed_fits/total_groups*100:.1f}%)")
    
    return df_with_td


def plot_sample_fits(df: pd.DataFrame, config: Config, output_path: Path):
    """Plot sample fits using configuration parameters"""
    fig, axs = plt.subplots(2, 5, figsize=(20, 8))
    axs = axs.ravel()
    
    valid_groups = df[
        (~df['Td'].isna()) & 
        (df['Td'] > 0) & 
        (df['Td'] < config.max_td_value)
    ]['group'].unique()
    
    num_plots = min(config.num_sample_fits, len(valid_groups))
    sample_groups = np.random.choice(valid_groups, num_plots, replace=False)
    
    def exp_func(x, b, A, Td):
        return b + A * (1 - np.exp(-x / Td))
    
    # Plot each sample
    for idx, group in enumerate(sample_groups):
        group_data = df[df['group'] == group]
        dark_values = group_data['dark'].values * config.exposure_time
        
        # Calculate histogram and cumulative frequency
        hist, bin_edges = np.histogram(dark_values, bins=min(50, len(dark_values)), density=True)
        cum_freq = np.cumsum(hist) * np.diff(bin_edges)
        x_data = bin_edges[1:]
        
        # Get Td and calculate fit
        Td = group_data['Td'].iloc[0]
        x_fit = np.linspace(min(x_data), max(x_data), 100)
        p0 = [np.min(cum_freq), np.max(cum_freq)-np.min(cum_freq), Td]
        popt, _ = curve_fit(exp_func, x_data, cum_freq, p0=p0)
        y_fit = exp_func(x_fit, *popt)
        
        # Create plot
        axs[idx].plot(x_data, cum_freq, 'b.', label='Data', alpha=0.5)
        axs[idx].plot(x_fit, y_fit, 'r-', label=f'Fit (Td={Td:.2f}s)')
        axs[idx].set_title(f'Group {group}')
        axs[idx].set_xlabel('Dark Time [s]')
        axs[idx].set_ylabel('Cumulative Frequency')
        axs[idx].legend()
    
    # Adjust layout and save
    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

def save_all_results(df: pd.DataFrame, yaml_docs: List, config: Config, original_name: str, num_sample_fits: int = 10) -> Tuple[Path, Path, Path]:
    """
    Save all results: HDF5 data, YAML metadata, and sample fits plot
    
    Args:
        df: DataFrame with the processed data
        yaml_docs: List of YAML documents
        config: Configuration parameters
        original_name: Original filename without '.hdf5'
        num_sample_fits: Number of sample fits to plot (default=10)
    """
    
    # Ensure ouput directory exists
    config.output_dir.mkdir(parents=True, exist_ok=True)
    # Create filenames using the original name
    hdf5_filename = config.output_dir / f'{original_name}_TdCDF.hdf5'
    yaml_filename = config.output_dir / f'{original_name}_TdCDF.yaml'
    plot_filename = config.output_dir / f'{original_name}_sample_fits.png'
    
    # Save HDF5
    with h5py.File(hdf5_filename, 'w') as f:
        dtype = [(name, df[name].dtype) for name in df.columns]
        structured_array = np.empty(len(df), dtype=dtype)
        for name in df.columns:
            structured_array[name] = df[name]
        f.create_dataset('locs', data=structured_array)
    
    # Save YAML with footer
    footer_template = {
        'TdCDF calculated using TdCDF_of_single_clusters_v4.py': {
            'Td [s]': float(df['Td'].mean()),
            'fit fuction used': 'b + A * (1 - np.exp(-x / Td))'
        }
    }
    
    output_docs = yaml_docs + [footer_template]
    with open(yaml_filename, 'w') as f:
        yaml.safe_dump_all(output_docs, f, default_flow_style=False)
    
    # Save sample fits plot
    plot_sample_fits(df, config, plot_filename)
    
    # Create a single-row dataframe for each group
    single_row_df = df.groupby('group').first().reset_index()

    # Create fitlenames for the single-row fiels
    single_row_hdf5_filename = config.output_dir / f'{original_name}_TdCDF_single_row.hdf5'
    single_row_yaml_filename = config.output_dir / f'{original_name}_TdCDF_single_row.yaml'    

    # Save single-row hdf5
    with h5py.File(single_row_hdf5_filename, 'w') as f:
        dtype = [(name, single_row_df[name].dtype) for name in single_row_df.columns]
        structured_array = np.empty(len(single_row_df), dtype=dtype)
        for name in single_row_df.columns:
            structured_array[name] = single_row_df[name].to_numpy()
        f.create_dataset('locs', data=structured_array)

    # Save single-row yaml
    single_footer_template = {
        'TdCDF calculated using TdCDF_of_single_clusters_v4.py': {
            'Td [s]': float(df['Td'].mean()),
            'fit fuction used': 'b + A * (1 - np.exp(-x / Td))'
        }
    }

    single_output_docs = yaml_docs + [single_footer_template]
    with open(single_row_yaml_filename, 'w') as f:
        yaml.safe_dump_all(single_output_docs, f, default_flow_style=False)

    print(f"\nResults saved to:")
    print(f"Data: {hdf5_filename}")
    print(f"Metadata: {yaml_filename}")
    print(f"Sample fits: {plot_filename}")
    print(f"Single-row data: {single_row_hdf5_filename}")
    print(f"Single-row metadata: {single_row_yaml_filename}")
    
    return hdf5_filename, yaml_filename, plot_filename, single_row_hdf5_filename, single_row_yaml_filename

def process_all_files(config: Config):
    """Process all files using configuration parameters"""
    try:
        # Ensure output directory exists at the start
        config.output_dir.mkdir(parents=True, exist_ok=True)
        print(f"\nOuput directory verified: {config.output_dir}")

        print(f"\nSearching for HDF5 files in: {config.input_dir}")
        all_hdf5_files = list(config.input_dir.glob('*.hdf5'))
        hdf5_files = [f for f in all_hdf5_files]
        # If you want to read hdf5 files that ends with '_filtered_in.hdf5'
        # hdf5_files = [f for f in all_hdf5_files if f.name.endswith('_filtered_in.hdf5')]    
        # If you want to read hdf5 files that contains 'dbscan' and 'in'
        #hdf5_files = [f for f in all_hdf5_files if 'dbscan' in f.name.lower() and 'in' in f.name.lower()]

        print(f"\nSearching for HDF5 files in: {config.input_dir}")
        print(f"Found {len(hdf5_files)} HDF5 files")
        for f in hdf5_files:
            print(f"  {f.name}") 
 
    
        if not hdf5_files:
            print("No .hdf5 files found in input directory")
            return
        
        results = []
        for hdf5_file in hdf5_files:
            print(f"\nChecking file: {hdf5_file.name}")
            # Skip files that are already processed
            if '_TdCDF.hdf5' in hdf5_file.name:
                print(f"Skipping processed file: {hdf5_file.name}")
                continue
            # Find corresponding YAML file
            yaml_file = hdf5_file.with_suffix('.yaml')
            if not yaml_file.exists():
                print(f"Warning: No matching YAML file for {hdf5_file.name}")
                continue
            
            try:
                print(f"\nProcessing file pair:")
                print(f"HDF5: {hdf5_file.name}")
                print(f"YAML: {yaml_file.name}")
                
                # Read input files
                df, yaml_docs = read_input_files(hdf5_file, yaml_file, config)
            
                # Process data
                df = fit_cumulative_frequency(df, config)
                
                # Print summary
                print("\nTd values for each group:")
                td_summary = df.groupby('group')['Td'].first().reset_index()
                print(td_summary)
            
                # Save all results
                if '_dark.hdf5' in hdf5_file.name:
                    original_name = hdf5_file.stem.replace('_dark.hdf5','_TdCDF.hdf5')
                else:
                    original_name = hdf5_file.stem  

                result = save_all_results(df, yaml_docs, config, original_name)
                results.append(result)
            
            except Exception as e:
                print(f"Error processing {hdf5_file.name}: {str(e)}")
                continue
    
        # Print summary of all processed files
        print("\nProcessing Summary:")
        print(f"Total files processed: {len(results)}")
        print("\nOutput files:")
        for hdf5, yaml, plot, single_row_hdf5, single_row_yaml in results:
            print(f"\nFile set:")   
            print(f"Data: {hdf5}")
            print(f"Metadata: {yaml}")
            print(f"Plot: {plot}")
            print(f"Single-row data: {single_row_hdf5}")
            print(f"Single-row metadata: {single_row_yaml}")
    except Exception as e:
        print(f"Error processing all files: {str(e)}")
        raise

def main():
    """Main function to run the analysis"""
    config = Config(
        input_dir=INPUT_DIR
        # If you do not specify output_dir, it will create 'output' in input_dir
        #output_dir=r"C:\Users\pcoffice79\Desktop\test\dark_time\output"
    )
    process_all_files(config)

if __name__ == "__main__":
    main()
