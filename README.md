# qPAINT-dark-time-analysis

Scripts that process DBSCANed SMLM data to extract dark times for quantitative PAINT
- DBSCAN_clusterfiltering: cluster filtering of dbscaned clusters using Picasso software
- TdCDF_of_single_clusters: extract dark times of a single cluster and determine the mean dark time

## DBSCAN_cluster_filtering.py
Filters out clusters based on the binding events distribution using the mean frame, the standard deviation of the mean frame, and the area
**Input Files Required:**  
- dbscan.hdf5 or clustered_hdf5      
- dbscan_centers.hdf5 or clustered_centers.hdf5  
- dbscan.yaml or clustered.yaml 
- dbscan_centers.yaml or clustered_centers.yaml

**Output Files Generated:**
1. *_filtered_in.hdf5 & *_filtered_in.yaml: Contains filtered-in data  
2. *_filtered_out.hdf5 & *_filtered_out.yaml: Contains filtered-out data  
3. *_centers_filtered_in.hdf5 & *_centers_filtered_in.yaml: Contains centers filtered-in data  
4. *_centers_filtered_out.hdf5 & *_centers_filtered_out.yaml: Contains centers filtered-out data  

**Usage:**  
1. Set your input directory in the main() function below and run  

**Parameters (adjustable in Config class):**  
- apply_area_filter: Whether to apply area filter (default: True)  
- max_area: Maximum area for filtering (default: 0.25)  

## TdCDF_of_single_clusters.py  
Extract the dark times of individual clusters. Then, fit the dark time distribution of each cluster with a mono-exponential function to determine a mean dark time per cluster.  

**Input Files Required:**  
- HDF5 files (.hdf5) containing localization data after DBSCAN clustering and filtering  
- Corresponding YAML files (.yaml) with metadata  
 
**Output Files Generated:**  
1. Intermediate files (for data without dark column):  
   - *_dark.hdf5 & *_dark.yaml: Contains calculated dark times  
2. Final results  
   - *_TdCDF.hdf5 & *_TdCDF.yaml: Contains final analysis results  
   - *_sample_fits.png: Visualization of sample fits  

**Usage:**  
1. Set your input directory in the main() function below  
2. (Optional) Adjust analysis parameters in the Config class if needed  
3. Run the script  

**Parameters (adjustable in Config class):**  
- exposure_time: Camera exposure time in seconds (default: 0.15s)  
- min_points_for_fit: Minimum points required for fitting (default: 6)  
- max_td_value: Maximum allowed Td value in seconds (default: 5000s)  
- num_sample_fits: Number of sample fits to plot (default: 10)  
- num_histogram_bins: Number of bins for histogram (default: 100)  

## Install  
```
conda create --name dark_time_analysis python=3.10
conda activate dark_time_analysis
conda install numpy pandas scipy matplotlib h5py pyyaml pathlib
```  
## Run  
Set the folder where you have dbscan.hdf5 & DBSCAN.yaml & dbscan_centers.hdf5 & dbscsan_centers.yaml 
```
python TdCDF_of_single_clusters.py
```





