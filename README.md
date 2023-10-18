This repository contains replication code for ERA5-Land and GMFD Uncover The Effect of Daily Temperature Extremes on Agricultural Yields, which is an unpublished manuscript currently in review. The code in this repository analyzes climate and agricultural yield data from various sources. These data are available for download from the [Harvard Dataverse](https://doi.org/10.7910/DVN/XRMDBW).

## System requirements and installation guide
The code was written and tested in R, version 4.2.1. Packages required to run the code are in `R_analysis/init.R`. This script is set up to install packages to your current R environment if they are not already installed. You can also install these packages manually with `install.packages()` or a conda environment. In addition, users must have [jupyter](https://jupyter.org/) installed to run the analysis, as the functions written in R are run in notebooks to produce the final outputs (i.e., figures and tables) for the paper. 

## Instructions for use

To run the code, download the compressed files from the link above and unzip all files in the same directory. Modify the `data_path.txt` file in the root folder of this repository with the directory containing the unzipped data files. The notebooks and functions will look for the input data sets in this specified path. The code files are organized as follows:

- `R_analysis/`
    - `analysis/` - Contains functions for loading data running each step of the econometric analysis in the paper, i.e., weather-yield relationships, cross validation, and climate projections.
    - `figures/` - Contains functions for plotting results from the analysis, including figures in the main text and Supplementary Information.
    - `notebooks/` - Contains Jupyter notebooks, which import functions from the other folders to generate and store output for the paper.
- `outputs/` - stores intermediate and final outputs, including all tables and figures for the paper.

The file names in `R_analysis/notebooks` correspond to the output they generate. For example, the notebook `R_analysis/notebooks/Fig1_response_functions.ipynb` produces Figure 1 in the paper.
