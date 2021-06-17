## DMP microbiome description

This folder contains scripts used to analyze basic statistics for microbiome data and produce Fig.1/B, Fig.S1 and Fig.S2 (panels A-D)

### Dependencies:

- codes depend on following R packages: ggplot2,iNEXT, coin, vegan, pheatmap, viridis. Running *install.packages(c("ggplot2","iNEXT","coin","vegan","pheatmap","viridis"))* is usually enough to install these. 
- <name>_Mockdata.R are implementations using *mock data* distributed with this git repo. These are intended for code testing without having to download the real DMP microbiome data from the EGA, but produce results different from the results of DMP manuscript. Make sure to unzip the zipped mock data in 
- Mock data (distributed in this github repo folder *Mock_data* and *heritability_analysis_v2/mock_data*) is required to execute <name>_Mockdata.R scripts. NOTE: mock data in *heritability_analysis_v2/mock_data* has to be unzipped
- EGA data (can be downloaded from EGA (study EGAS00001005027, see manuscript for details), data should be unzipped into *DMP_data_EGA* folder (empty folder is part of git repo) is required to execute <name>_EGAdata.R scripts

### Notes:

- make sure all dependancies are installed / ready
- by default, scripts are intended to be executed from the DMP folder (root of this repo), for example using: "Rscript microbiome_description/DMP_rarefaction_Fig1B_S1_Mockdata.R" command. If codes are run from R studio or from different path, codes should be edited to set working directory to DMP folder
- scripts were tested in R version <xyz> (see sessioninfo.txt for details)
- scripts normally produce some warning messages, this is intended behavior

### Folders:

- *plots* : folder contains expected output of scripts for Mock data and EGA data
- *misc* : folder contains details of environment (OS, R, R libraries) used to test scripts
- *data* : folder contains MetaCyc pathway grouping file (tested_pathway_Name_class_cleaned_category_new.txt) required to run DMP_taxa_pathways_variance_FigS2_CD.R and DMP_taxa_pathways_variance_Mockdata_FigS2_CD.R and tables which are expected output of these scripts

### Files:

- DMP_rarefaction_Fig1B_S1_EGAdata.R : script that performs resampling and rarefaction analysis of numbers of expected micobial features in the data, implementation on EGA data (requires EGA data to be downloaded and unzipped into *DMP_data_EGA* folder). Produces Figures 1/B and Supplementary Figure 1
- DMP_rarefaction_Fig1B_S1_Mockdata.R : as *DMP_rarefaction_Fig1B_S1_EGAdata.R*, but uses Mock Data included in this repo. 
- DMP_microbiome_FigS2A_1C_EGAdata_MockPhenotypes.R : script that generate microbiome ordination biplot (Fig 1/C) and microbiome ordination plot colored by *P. copri* abundance (Supplementary Figure 2/A), implementation on EGA data. NOTE: please note that this script uses DMP microbiome data (has to be downloaded from EGA) and *Mock phenotype data* as the distribution of DMP phenotype data is under restricted access to maintain the privacy of participants and cannot be made availabe at EGA or public domain. DMP phenotype data can be requested from Lifelines biobank (see manuscript for details)
- DMP_microbiome_FigS2A_1C_Mockdata_MockPhenotypes.R : as *DMP_microbiome_FigS2A_1C_EGAdata_MockPhenotypes.R*, but uses Mock microbiome data avaliable with this repo.
- DMP_microbiome_FigS2B_EGAdata.R : generates average phylum composition piechart (Fig S2/B) using real data (avaliable at EGA)
- DMP_microbiome_FigS2B_Mockdata.R : as *DMP_microbiome_FigS2B_EGAdata.R*, but implemented on mock data distributed with this repo
- DMP_taxa_pathways_variance_FigS2_CD_EGAdata.R : plots phyla and pathway class distribution across the cohort (Fig S2/C, Fig S2/D) and quantifies variances of taxa vs pathways. Implementation on real DMP data
- DMP_taxa_pathways_variance_FigS2_CD_Mockdata.R : as *DMP_taxa_pathways_variance_FigS2_CD_EGAdata*, but uses Mock Data avaliable on this github (Note: make sure that mock data in heritability_analysis_v2/mock_data is unzipped)