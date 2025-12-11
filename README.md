# qPAINT-dark-time-analysis
Scripts that process DBSCANed hdf5 files to extract dark time for quantitative PAINT


**DBSCAN_cluster_filtering.py**
Filters out clusters based on ther binding events distribution using mean frame, standard deviation of mean frame, and area
Input Files Required:
    - dbscan.hdf5 or clustered_hdf5
    - dbscan_centers.hdf5 or clustered_centers.hdf5
    - dbscan.yaml or clustered.yaml   
    - dbscan_centers.yaml or clustered_centers.yaml

Output Files Generated:
    1. *_filtered_in.hdf5 & *_filtered_in.yaml: Contains filtered-in data
    2. *_filtered_out.hdf5 & *_filtered_out.yaml: Contains filtered-out data
    3. *_centers_filtered_in.hdf5 & *_centers_filtered_in.yaml: Contains centers filtered-in data
    4. *_centers_filtered_out.hdf5 & *_centers_filtered_out.yaml: Contains centers filtered-out data

Usage:
    1. Set your input directory in the main() function below and run

Parameters (adjustable in Config class):
    - apply_area_filter: Whether to apply area filter (default: True)
    - max_area: Maximum area for filtering (default: 0.25)


**TdCDF_of_single_clusters.py**
Extract the dark times of individual clusters. Then, fit the dark time distribution of each cluster with a mono-exponential function to determine a mean dark time per cluser.
Input Files Required:
    - HDF5 files (.hdf5) containing localization data after DBSCAN clustering and filtering
    - Corresponding YAML files (.yaml) with metadata
    
Output Files Generated:
    1. Intermediate files (for data without dark column):
        - *_dark.hdf5 & *_dark.yaml: Contains calculated dark times
    2. Final results:
        - *_TdCDF.hdf5 & *_TdCDF.yaml: Contains final analysis results
        - *_sample_fits.png: Visualization of sample fits

Usage:
    1. Set your input directory in the main() function below
    2. (Optional) Adjust analysis parameters in the Config class if needed
    3. Run the script

Parameters (adjustable in Config class):
    - exposure_time: Camera exposure time in seconds (default: 0.15s)
    - min_points_for_fit: Minimum points required for fitting (default: 6)
    - max_td_value: Maximum allowed Td value in seconds (default: 5000s)
    - num_sample_fits: Number of sample fits to plot (default: 10)
    - num_histogram_bins: Number of bins for histogram (default: 100)
