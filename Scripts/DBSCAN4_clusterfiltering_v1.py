######################################################################################
# DBSCAN cluster filtering script with 4 input files
# @author: Soohyen Jang 
# Research group of Mike Heilemann, Goethe University Frankfurt a.M.
#
# Input Files Required:
#     - dbscan.hdf5 or clustered_hdf5
#     - dbscan_centers.hdf5 or clustered_centers.hdf5
#     - dbscan.yaml or clustered.yaml   
#     - dbscan_centers.yaml or clustered_centers.yaml
# Output Files Generated:
#     1. *_filtered_in.hdf5 & *_filtered_in.yaml: Contains filtered-in data
#     2. *_filtered_out.hdf5 & *_filtered_out.yaml: Contains filtered-out data
#     3. *_centers_filtered_in.hdf5 & *_centers_filtered_in.yaml: Contains centers filtered-in data
#     4. *_centers_filtered_out.hdf5 & *_centers_filtered_out.yaml: Contains centers filtered-out data
# Usage:
#     1. Set your input directory in the main() function below
# Parameters (adjustable in Config class):
#     - apply_area_filter: Whether to apply area filter (default: True)
#     - max_area: Maximum area for filtering (default: 0.25)
######################################################################################


from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, List, Dict, Optional
import h5py
import yaml
import pandas as pd
import numpy as np
from scipy.stats import norm


INPUT_DIR = r"E:\SMLM_FAMP\20251020_Vidhya\251020_FAM134B_SiR_Hy5_1nM_4"


@dataclass
class Config:
    input_dir: Path
    output_dir: Path | None = None # Make output_dir optional
    apply_area_filter: bool = True # This need to be defined by the user
    max_area: float = 1
    min_std_frame: Optional[float] = 3000 # If you want to set the min value for std_frame, set it here
    max_std_frame: Optional[float] = 9000 # If you want to set the max value for std_frame, set it here
    std_frame_filter_percentage: float = 20.0 # If you do not set the min and max for std_frame, you can remove small values of std_frame 
    max_localizations_per_group: int = 50 # Set the value to remove fiducial markers

@dataclass
class FilSet:
    base_name: str
    dbscan_hdf5: Path
    centers_hdf5: Path
    dbscan_yaml: Path
    centers_yaml: Path

class DataProcessor:
    def __init__(self, config: Config):
        self.config = config
        # If you want to save output files in an input directory, use the line below
        if self.config.output_dir is None:
            self.config.output_dir = self.config.input_dir 
        
        # If you want to save output files in a specific directory, use the line below
        # if self.config.output_dir is None:
        #     self.config.output_dir = self.config.input_dir / 'cluster_filtered'
        #self.config.output_dir.mkdir(parents=True, exist_ok=True)

    def find_file_sets(self) -> List[FilSet]:
        """Find all sets of files with the same base name"""
        file_sets = []
        dbscan_files = list(self.config.input_dir.glob('*_dbscan.hdf5'))
        
        for dbscan_file in dbscan_files:
            base_name = dbscan_file.stem.replace('_dbscan', '')
            
            # Skip filtered files
            if 'filter' in base_name.lower():
                print(f"Skipping filtered file: {base_name}")
                continue
            
            # Construct expected file paths
            file_set = FilSet(
                base_name = base_name,
                dbscan_hdf5 = dbscan_file,
                centers_hdf5 = self.config.input_dir / f"{base_name}_dbscan_centers.hdf5",
                dbscan_yaml = self.config.input_dir / f"{base_name}_dbscan.yaml",
                centers_yaml = self.config.input_dir / f"{base_name}_dbscan_centers.yaml"
            )

            # Verity all files exist
            if self._validate_file_set(file_set):
                file_sets.append(file_set)
        return file_sets

    def process_file_set(self, file_set: FilSet) -> bool:
        print(f"\nProcessing {file_set.base_name}...")

        try:
            # Read input files
            df_dbscan, df_centers, yaml_docs_dbscan, yaml_docs_centers = self._read_input_files(file_set)   

            # Process data
            min_frame, max_frame, min_std_frame, max_std_frame = self._calculate_bounds(df_centers)
            # Count group sizes in the original dbscan DataFrame
            group_counts_dbscan = df_dbscan['group'].value_counts()
            valid_groups_by_size = group_counts_dbscan[group_counts_dbscan <= self.config.max_localizations_per_group].index
            filtered_centers = self._filter_centers(
                df_centers, min_frame, max_frame, min_std_frame, max_std_frame, valid_groups_by_size
            )
            valid_groups = filtered_centers['group'].unique()
            dbscan_filter_in, dbscan_filter_out = self._filter_dbscan_groups(df_dbscan, valid_groups)
            centers_filter_in, centers_filter_out = self._filter_dbscan_groups(df_centers, valid_groups)

            # Save the results
            self._save_filtered_data(
                dbscan_filter_in, dbscan_filter_out, centers_filter_in, centers_filter_out, yaml_docs_dbscan, yaml_docs_centers,
                file_set.base_name, min_frame, max_frame, min_std_frame, max_std_frame
            )

            print(f"Successfully processed {file_set.base_name}:")
            print(f"  Filtered in: {len(dbscan_filter_in)} entries")
            print(f"  Filtered out: {len(dbscan_filter_out)} entries")
            print(f"  Centers Filtered in: {len(centers_filter_in)} entries")
            print(f"  Centers Filtered out: {len(centers_filter_out)} entries")
            return True
        
        except Exception as e:
            print(f"Error processing {file_set.base_name}:")
            print(f"  {str(e)}")
            return False

    def _read_input_files(self, file_set: FilSet) -> Tuple[pd.DataFrame, pd.DataFrame, List, List]:
        """
        Read HDF5 and YAML files for both dbscan and dbscan_centers
        """
        print(f"\nReading files:")
        print(f"DBSCAN HDF5: {file_set.dbscan_hdf5.name}")
        print(f"Centers HDF5: {file_set.centers_hdf5.name}")
        print(f"DBSCAN YAML: {file_set.dbscan_yaml.name}")
        print(f"Centers YAML: {file_set.centers_yaml.name}")
    

        # Read the files
        with h5py.File(file_set.dbscan_hdf5, 'r') as dbscan_file:
            df_dbscan = pd.DataFrame(dbscan_file['locs'][()])
    
        with h5py.File(file_set.centers_hdf5, 'r') as centers_file :
            df_centers = pd.DataFrame(centers_file['locs'][()])
    

        with open(file_set.dbscan_yaml, 'r') as dbscan_yaml_file:
            yaml_docs_dbscan = list(yaml.safe_load_all(dbscan_yaml_file))
    
        with open(file_set.centers_yaml, 'r') as centers_yaml_file:
            yaml_docs_centers = list(yaml.safe_load_all(centers_yaml_file))

        return df_dbscan, df_centers, yaml_docs_dbscan, yaml_docs_centers
    
    def _validate_file_set(self, file_set: FilSet) -> bool:
        """Verify that all files in the set exist"""
        required_files = {
            'DBSCAN HDF5': file_set.dbscan_hdf5,
            'Centers HDF5': file_set.centers_hdf5,
            'DBSCAN YAML': file_set.dbscan_yaml,
            'Centers YAML': file_set.centers_yaml
        }
        missing_files = [name for name, path in required_files.items() if not path.exists()]
        
        if missing_files:
            print(f"\nSkipping {file_set.base_name} due to missing files:")
            for file_type in missing_files:
                print(f"  - Missing {file_type}")
            return False
        return True



    def _calculate_bounds(self, df: pd.DataFrame) -> Tuple[float, float, float, float]:
        """
    Calculate the min and max values for the frame column
    """
        
        # Fit a Gaussian to the 'frame' column from centers.hdf5 and calculate min/max
        mu_frame, std_frame = norm.fit(df['frame'])
        min_frame = mu_frame - 2 * std_frame
        max_frame = mu_frame + 2 * std_frame

        if self.config.min_std_frame is not None and self.config.max_std_frame is not None:
            # Use user-defined boundaries
            min_std_frame = self.config.min_std_frame
            max_std_frame = self.config.max_std_frame
        else:
            # Remove the beginning --% of std_frame
            percentage_to_remove = self.config.std_frame_filter_percentage / 100.0
            ten_percent_index = int(len(df)*percentage_to_remove)
            filtered_std_frame = df['std_frame'].iloc[ten_percent_index:]
        
            # Fit a Gaussian to the 'std_frame' column from centers.hdf5 and calculate min/max
            mu_std_frame, std_std_frame = norm.fit(filtered_std_frame)
            min_std_frame = mu_std_frame - 2 * std_std_frame
            max_std_frame = mu_std_frame + 2 * std_std_frame

        return min_frame, max_frame, min_std_frame, max_std_frame


    def _filter_centers(self, df: pd.DataFrame, min_frame: float, max_frame: float, min_std_frame: float, max_std_frame: float, valid_groups_by_size) -> pd.DataFrame:
        """
        Filter the DBSCAN groups based on frame and std_frame ranges
        """
        # Filter the Dataframe by frame and std_frame
        frame_filter = (df['frame'] >= min_frame) & (df['frame'] <= max_frame)
        std_frame_filter = (df['std_frame'] >= min_std_frame) & (df['std_frame'] <= max_std_frame)
    
        # Combine filters
        combined_filter = frame_filter & std_frame_filter
        filter_centers_in = df[combined_filter]
    
        # Filter using the area 
        if self.config.apply_area_filter and 'area' in filter_centers_in.columns:
            area_filter = (filter_centers_in['area'] <= self.config.max_area )
            filter_centers_in = filter_centers_in[area_filter]

        # Only keep centers from groups with <= max_localizations_per_group in the original dbscan
        filter_centers_in = filter_centers_in[filter_centers_in['group'].isin(valid_groups_by_size)]
        
        return filter_centers_in

    def _filter_dbscan_groups(self, df: pd.DataFrame, valid_groups: np.ndarray) -> Tuple[pd.DataFrame, pd.DataFrame]:
        """
        Filter the DBSCAN groups based on the valid_groups series
        """    
        # Filter by valid groups
        dbscan_filter_in = df[df['group'].isin(valid_groups)]
        dbscan_filter_out = df[~df['group'].isin(valid_groups)]
    
        return dbscan_filter_in, dbscan_filter_out


    #yaml.add_representer(OrderedDict, lambda dumper, data: dumper.represent_mapping('tag:yaml.org,2002:map', data.items()))

    def _save_filtered_data(self, df_in: pd.DataFrame, df_out: pd.DataFrame, centers_in: pd.DataFrame, centers_out: pd.DataFrame, yaml_docs_dbscan: List, yaml_docs_centers: List, base_name: str, 
        min_frame: float, max_frame: float, min_std_frame: float, max_std_frame: float) -> None:
        """
        Save the filtered data to HDF5 and YAML files
        
        Args:
            df_in: DataFrame containing filtered-in data
            df_out: DataFrame containing filtered-out data
        yaml_docs: Original YAML documents to include
            base_name: Base name for the output files
        """
    # Prepare file paths
        files = {
            'in': {
                'dbscan_hdf5': self.config.output_dir / f"{base_name}_dbscan_filtered_in.hdf5",
                'dbscan_yaml': self.config.output_dir / f"{base_name}_dbscan_filtered_in.yaml",
                'centers_hdf5': self.config.output_dir / f"{base_name}_dbscan_centers_filtered_in.hdf5",
                'centers_yaml': self.config.output_dir / f"{base_name}_dbscan_centers_filtered_in.yaml"
                },
            'out': {
                'dbscan_hdf5': self.config.output_dir / f"{base_name}_dbscan_filtered_out.hdf5",
                'dbscan_yaml': self.config.output_dir / f"{base_name}_dbscan_filtered_out.yaml",
                'centers_hdf5': self.config.output_dir / f"{base_name}_dbscan_centers_filtered_out.hdf5",
                'centers_yaml': self.config.output_dir / f"{base_name}_dbscan_centers_filtered_out.yaml"
            }
        }

      # Footer template with statistics in specific order
        def create_in_footer() -> dict:
            return {
                'Cluster filtered - in': {
                    'Min frame': float(min_frame),
                    'Max frame': float(max_frame),
                    'Min std_frame': float(min_std_frame),
                    'Max std_frame': float(max_std_frame),
                    'Area filter applied': self.config.apply_area_filter,
                    'Max area': float(self.config.max_area) 
                }
        }
    
        # Footer template with statistics
        def create_out_footer() -> dict:
            return  {
                'Cluster filtered - out': {
                    'Min frame': float(min_frame),
                    'Max frame': float(max_frame),
                    'Min std_frame': float(min_std_frame),
                    'Max std_frame': float(max_std_frame),
                    'Area filter applied': self.config.apply_area_filter,
                    'Max area': float(self.config.max_area) 
                }
            }

        # Save both filtered-in and filtered-out data
        for filter_type, data in [
            ('in', (df_in, centers_in, yaml_docs_dbscan, yaml_docs_centers)),
            ('out', (df_out, centers_out, yaml_docs_dbscan, yaml_docs_centers))
        ]:
            dbscan_df, centers_df, dbscan_yaml, centers_yaml = data

            # Save HDF5
            with h5py.File(files[filter_type]['dbscan_hdf5'], 'w') as f:
                dtype = [(name, dbscan_df[name].dtype) for name in dbscan_df.columns]
                structured_array = np.empty(len(dbscan_df), dtype=dtype)
                for name in dbscan_df.columns:
                    structured_array[name] = dbscan_df[name]
                f.create_dataset('locs', data=structured_array)
        
            # Save Centers HDF5 and YAML
            with h5py.File(files[filter_type]['centers_hdf5'], 'w') as f:
                dtype = [(name, centers_df[name].dtype) for name in centers_df.columns]
                structured_array = np.empty(len(centers_df), dtype=dtype)
                for name in centers_df.columns:
                    structured_array[name] = centers_df[name]
                f.create_dataset('locs', data=structured_array)

            # Save YAML with footer
            footer = create_in_footer() if filter_type == 'in' else create_out_footer()

            output_docs_dbscan = dbscan_yaml + [footer]
            with open(files[filter_type]['dbscan_yaml'], 'w') as f:
                yaml.dump_all(output_docs_dbscan, f, default_flow_style=False, sort_keys=False)

            output_docs_centers = centers_yaml + [footer]
            with open(files[filter_type]['centers_yaml'], 'w') as f:
                yaml.dump_all(output_docs_centers, f, default_flow_style=False, sort_keys=False)
        
            print(f"Saved {filter_type} files:")
            print(f"  HDF5: {files[filter_type]['dbscan_hdf5'].name}")
            print(f"  YAML: {files[filter_type]['dbscan_yaml'].name}")
            print(f"  Centers HDF5: {files[filter_type]['centers_hdf5'].name}")
            print(f"  Centers YAML: {files[filter_type]['centers_yaml'].name}")            


def main():
    """Main function to run the analysis"""
    config = Config(
        input_dir = Path(INPUT_DIR),
        apply_area_filter = True     
    )
    
    processor = DataProcessor(config)
    file_sets = processor.find_file_sets()

    if not file_sets:
        print(f"No valid files found to process.")
        return
    
    processed_count = 0
    total_files = len(file_sets)

    print(f"\nFound {total_files} file sets to process.")

    for i, file_set in enumerate(file_sets, 1):
        print(f"\nProcessing file set {i} of {total_files}:")
        if processor.process_file_set(file_set):
            processed_count += 1

    print(f"\nProcessed {processed_count} out of {total_files} files.")

if __name__ == "__main__":
    main()



