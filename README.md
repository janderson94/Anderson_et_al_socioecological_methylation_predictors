### Anderson et al: DNA methylation signatures of early life adversity are exposure-dependent in wild baboons 

![image](./misc/baboon.jpg)

#### This repo is setup to be reproducible from intermediate files. Intermediate files should be downloaded from a zenodo dataset (link forthcoming) and dropped into the directory where this repository is cloned. Set your working to the respository, and you should be set up.

These scripts break down into 5 parts:

#### 1. Main models.
Reads in the MACAU results and calculates standardized betas and p-values for each variable, stores them for downstream analyses. 
#### 2. Figures
Takes the output from Main Models and produces all of the main and supplementary figures in the manuscript.
#### 3. Misc_for_mansucript
Reproduces results in text, as they appear.
#### 4. Baboon_mstar_scripts_final_ms.txt 
Demonstrates the initial processing workflow of raw mstarr data.
#### 5. Baboon_mstarr_scripts_final_ms.R
Code for generating permutation-based p-values and associated q-values for identifying regulatory and methylation-dependent regulatory activity.

Code for trimming, mapping, and filtering methylation data can be found https://github.com/janderson94/BaboonEpigeneticAging/tree/master/Trim_map_filter
