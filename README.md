# qPAINT-dark-time-analysis batch processing

This script extends the functionality of the TdCDF_of_single_clusters_v01.py script from https://github.com/SoohyenJang/qPAINT-dark-time-analysis/tree/main/Scripts
by enabeling batch processing of input files. Multiple HDF5 files are recursively processed within a specified folder and all its subfolders, and an output folder is created automatically in the input directory and all its subfolders. Only files with a user defined file suffix are processed within the specified folder and subfolders. The existing analysis procedure for individual files was retained.

Script that process DBSCANed SMLM data to extract dark times for quantitative PAINT
- TdCDF_of_single_clusters: extract dark times of a single cluster and determine the mean dark time

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
1. Set your input directory (INPUT_DIR) in the the user input section below (multiple HDF5 files are recursively processed within the specified folder and its subfolders)   
2. Set the desired file suffix (HDF5_SUFFIX) of the HDF5 files containing your clustered localization data in the user input section below 
   (only HDF5 files with this file suffix will be processed) 
3. (Optional) Adjust analysis parameters in the Config class if needed
4. Run the script  

**Parameters (adjustable in Config class):**  
- exposure_time: Camera exposure time in seconds (default: 0.15s)  
- min_points_for_fit: Minimum points required for fitting (default: 6)
- min_histogram_bins: int (default: 5)
- max_td_value: Maximum allowed Td value in seconds (default: 5000s)  
- num_sample_fits: Number of sample fits to plot (default: 10)  
- num_histogram_bins: Number of bins for histogram (default: 100)
- bin_min: float: Minimum bin value for cumulative frequency (default: 0)
- bin_max: float: Maximum bin value for cumulative frequency (default: 5000)
- bin_step: float: Bin step size for cumulative frequency (default: 10)

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





