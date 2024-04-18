# GUV-symmetry-breaking-analysis
This repository contains the code used for symmetry breaking analysis of giant unilamellar vesicles (GUVs) in the article: 

Shiva Razavi, Felix Wong, Bedri Abubaker-Sharif, Hideaki T. Matsubayashi, Hideki Nakamura, Eduardo Sandoval, Douglas N. Robinson, Baoyu Chen, Jian Liu, Pablo A. Iglesias, and Takanari Inoue. “Synthetic control of actin polymerization and symmetry breaking in active protocells”. 
DOI:   

Code 1 generates masks to segment GUVs. 
G01_input folder contains raw tiff files separated by channels from experimental data. Run Code_1_segmentation.m to generate a mask for the target. Output is in the “G01_acs_v0” folder. There are tif files for the smoothed images (“_final.tif”), and tif files for the masks (the main one being “_mask_inclusion.tif”). To check the mask, open “_mask_fused.tif” which has the mask overlayed on the channel used to generate it (here, a membrane marker channel). Use “ctrl+f” to search for {{"**SETUP**"}} which indicate sections of the code with variables that can be adjusted as needed. 

Code 2 generate kymographs for biochemical properties (membrane  marker, ActA, and actin signal) and physical properties (curvature and eccentricity) of GUVs. The "target_settings.csv" file contains parameters for each guv that need to be adjusted (naming, number of frames to use, time of rapamycin addition, inner-outer mask size parameter, and center angle for kymographs). Use “ctrl+f” to search for <"**SETUP**"> which indicate sections of the code with variables can be adjusted as needed. Output for data from a single target's relevant tif files (e.g. tifs in “G01_acs_v0” folder) is in the G01_singlerun_output folder. If multiple targets are identified and segmented in separate folders (“G01_input”, “G02_input” etc,), then the settings for each target can be listed in the “target_settings.csv” document and Code_2_kymographs.m can be setup to read data from a table (set the “runtype” variable = 2).  

Code 3 applies correlation analysis and principal component analysis to study symmetry breaking. It also contains code for analysis of radial actin polymerization rate for locally stimulated GUVs. 

Code for generating simulations of ActA and Actin in supplemental figure 15 is found in: simActA_Actin_NOmarkernormalization.m

Additional details regarding the approaches used can be found in the Supplementary Information section of the paper. 

Please read the citation details below if you would like to use a part of this repository in your research.

## Requirements
MATLAB 2022a (or later versions) is required to run these programs. Additional code dependencies:

## Authors
The code was primarily developed by Bedri Abubaker-Sharif and Pablo A. Iglesias (Johns Hopkins University, Baltimore, MD, USA).

## Citation
This program is free software (please see License for details). However, if you are using this code in your work, please cite our work as:

Razavi, S. et al. Synthetic control of actin polymerization and symmetry breaking in active protocells. (JOURNAL INFO).

## License
Copyright © 2023 Shiva Razavi, Felix Wong, Bedri Abubaker-Sharif, Hideaki T. Matsubayashi, Hideki Nakamura, Eduardo Sandoval, Douglas N. Robinson, Baoyu Chen, Jian Liu, Pablo A. Iglesias, and Takanari Inoue.

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see https://www.gnu.org/licenses/.
