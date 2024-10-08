# GUV-symmetry-breaking-analysis
This repository contains the code used for symmetry breaking analysis of giant unilamellar vesicles (GUVs) in the article: 

Shiva Razavi, Felix Wong, Bedri Abubaker-Sharif, Hideaki T. Matsubayashi, Hideki Nakamura, Nhung Thi Hong Nguyen, Douglas N. Robinson, Baoyu Chen, Pablo A. Iglesias, and Takanari Inoue. Synthetic control of actin polymerization and symmetry breaking in active protocells. _Science Advances_ **10**, eadk9731 (2024). [DOI:10.1126/sciadv.adk9731](https://doi.org/10.1126/sciadv.adk9731)

### Code 1: Generating masks to segment GUVs. 
The _G01_input_ folder contains raw tiff files separated by channels from experimental data. Run Code_1_segmentation.m to generate a mask for the target. Output is in the _G01_acs_v0_ folder. There are tif files for the smoothed images (“_final.tif”), and tif files for the masks (the main one being “_mask_inclusion.tif”). To check the mask, open “_mask_fused.tif” which has the mask overlayed on the channel used to generate it (here, a membrane marker channel). Search the code file for "\*\*SETUP\*\*" which indicate sections of the code with variables that can be adjusted as needed. 

### Code 2: Generating kymographs for biochemical properties (membrane  marker, ActA, and actin signal) and physical properties (local velocities, cumulative displacements, and curvature) of GUVs. 
The "target_settings.csv" file contains parameters for each guv that need to be adjusted (naming, number of frames to use, time of rapamycin addition, inner-outer mask size parameter, and center angle for kymographs). Search the code file for "\*\*SETUP\*\*" which indicate sections of the code with variables can be adjusted as needed. Output for data from a single target's relevant tif files (e.g. tifs in _G01_acs_v0_ folder) is in the _G01_singlerun_output_ folder. If multiple targets are identified and segmented in separate folders (_G01_input/G01_acs_v0_, _G02_input/G02_acs_v0_ etc,), then the settings for each target can be listed in the “target_settings.csv” document and Code_2_kymographs.m can be setup to read data from a table (set the “runtype” variable to 2).  

### Code 3: Correlation analysis and principal component analysis to study symmetry breaking. 
It also contains code for analysis of radial actin polymerization rate for locally stimulated GUVs. 

Code for generating simulations of ActA and Actin in supplemental figure 15 is found in: simActA_Actin_NOmarkernormalization.m

### Output files
The following output file types are found in the _G01_output_kymos_tseries_ folder: 

- ".m" files: extra copy of the MATLAB code used to generate the output
- ".png" and ".fig" files: image files (and corresponding MATLAB figure files) for biophysical kymographs. File names containing "Var" refer to biophysical kymographs with varying number of boundary points (to track changes in GUV size). File names containing "Fix" refer to biophysical kymographs with fixed number of boundary points. 
- ".tif" files: time-series of inner/outer masks used to segment each GUV and compute boundary quantities. Masks are overlayed on the grayscale GUV images for different fluorescent channels, with the following color-coding: 1) green zones indicate the region outside the GUV; 2) red zones indicate the region inside the GUV; 3) the space between the dashed cyan and dashed red lines indicates the and red zones indicates the boundary region of the GUV from which biochemical quantities are computed; 4) the dashed white line indicates the boundary from which physical quantities are computed.
- ".mat" files: saved MATLAB variables for biophysical quantities

## Requirements
MATLAB 2022a (or later versions) is required to run these scripts. 

## Additional Information

Additional details regarding the approaches used can be found in the Supplementary Information section of the paper. Please read the citation details below if you would like to use a part of this repository in your research.

The GUV TIFF and kymograph dataset containing all input files and outputs based on this code can be found on Dryad: https://doi.org/10.5061/dryad.tqjq2bw66

## Authors
The code was primarily developed by Bedri Abubaker-Sharif and Pablo A. Iglesias (Johns Hopkins University, Baltimore, MD, USA).

## Citation
This program is free software (please see License for details). However, if you are using this code in your work, please cite our work as:

Shiva Razavi et al., Synthetic control of actin polymerization and symmetry breaking in active protocells. _Science Advances_ **10**, eadk9731 (2024). [DOI:10.1126/sciadv.adk9731](https://doi.org/10.1126/sciadv.adk9731)

## License
Copyright © 2023 Shiva Razavi, Felix Wong, Bedri Abubaker-Sharif, Hideaki T. Matsubayashi, Hideki Nakamura, Eduardo Sandoval, Douglas N. Robinson, Baoyu Chen, Jian Liu, Pablo A. Iglesias, and Takanari Inoue.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
