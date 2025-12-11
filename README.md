# qPAINT-dark-time-analysis
Scripts that process DBSCANed hdf5 files to extract dark time for quantitative PAINT

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
