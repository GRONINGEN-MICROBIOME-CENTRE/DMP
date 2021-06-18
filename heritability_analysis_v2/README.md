## DMP heritability analysis, heritability analysis

The folder contains R scripts used to execute heritability analysis for DMP project, mock data for testing of codes, and expected results from running codes on the mock data.

### Notes:

- mock_data folder contains *mock* data, this data is intended to be used only for testing of codes. It was designed to produce results similar to real data, but it is *not* real DMP data. 
- real microbiome data used for DMP project can be obtained from EGA study EGAS00001005027, real kinship matrices, family information and cohousing data can be obtained from Lifelines Biobank
- example results generated from mock data are provided (zipped) in mock_data_heritability_results. These results are expected output from DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R
- results permutation runs are NOT provided due to github space constrains, these can be generated using DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R (see Instructions)
- example consolidated results produced from mock data are provided in root folder. These are expected results from running DMP_heritability_v10_mockdata_collectdata.R on mock data heritability models
- example plots generated from mock data results are provided in Plots folder. These are expected results from running DMP_heritability_v10_mockdata_plotresults.R
- codes were developed & tested on R 3.6.0 (see R_session_info.txt for list of used packages)

### Dependancies:

- R libraries compositions,coxme, lme4qtl, RLRsim, lmerTest, foreach, plyr, ggplot2, tidyr
- the majority of R libraries can be installed from CRAN with *install.packages(c("compositions","coxme", "RLRsim", "lmerTest", "foreach", "plyr", "ggplot2", "tidyr"))*, but lme4qtl does require additional steps (see https://github.com/variani/lme4qtl for details) 
- Mock Data (distributed with this repo, see *mock_data* folder, data has to be unzipped)
- Note: heritability analysis performed in DMP project used private data of participants which could not be shared on this repo or EGA to maintain privacy of participants. This data can be requested from Lifelines biobank (see manuscript for details)

### Instructions for code testing:

- download the codes and data
- unzip the mock data (mock_data/*.zip)
- Run these scripts form root folder of this repo (example: *Rscript heritability_analysis_v2/DMP_heritability_v10_mockdata_taxa.R*) or set working directories to the root of this repo (DMP folder) - following lines of code should be edited: DMP_heritability_v10_mockdata_taxa.R [line 70]; DMP_heritability_v10_mockdata_pwy.R [line 70]; DMP_heritability_v10_mockdata_plotresults.R [line 42]; DMP_heritability_v10_mockdata_collectdata.R [line 15]
- run DMP_heritability_v10_mockdata_taxa.R to calculate heritability models for taxa and DMP_heritability_v10_mockdata_pwy.R for pathways
- run DMP_heritability_v10_mockdata_collectdata.R to consolidate data and calculate study-wide FDRs from permutation runs
- run DMP_heritability_v10_mockdata_plotresults.R to generate heritability plots
- NOTE: DMP_heritability_v10_mockdata_plotresults.R can be run directly on the provided consolidated data, example input data is provided
- NOTE: DMP_heritability_v10_mockdata_collectdata.R requires permutation runs produced by MP_heritability_v10_mockdata_taxa.R and DMP_heritability_v10_mockdata_pwy.R 
- NOTE: these scripts are intended for code testing and as such have number of permutations set to 10 for quick execution. "Real world" calculation should be 1000 - 10,000 permutations to obtain more precise empyrical p-values and FDR estimates. 
- NOTE: confidence interval profiling can fail if random effect(s) of the model approach 0, this is "intended" behavior of profiling functions and would require major changes in R stats libraries to bypass. As such, in DMP we also constructed reduced models with fewer random effects (no cohousing, no family effect) to work around this "feature"

### Files:

- README.md : this readme file
- bacpaths.txt : list of bacterial pathways used in DMP project (required by DMP_heritability_v10_mockdata_collectdata.R script)
- DMP_heritability_v10_mockdata_taxa.R : script for heritability analysis of microbiome taxa, requires unzipped mock data
- DMP_heritability_v10_mockdata_pwys.R : script for heritability analysis of microbiome pathways, requires unzipped mock data
- DMP_heritability_v10_mockdata_collectdata.R : script for collecting results and study-wide FDR analysis using permutation runs. Requires results of DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwys.R
- DMP_heritability_v10_mockdata_plotresults.R : script for plotting heritability results. Requires results of DMP_heritability_v10_mockdata_collectdata.R
- results_mockdata_taxa.csv : results of heritability analysis of taxa
- results_mockdata_pwys.csv : results of heritability analysis of pathways
- results_mockdata_withFDRs_and_CIs_taxa.csv : results of heritability analysis of taxa with FDRs calculated from permutation runs
- results_mockdata_withFDRs_and_CIs_pwys.csv : results of heritability analysis of pathways with FDRs calculated from permutation runs
- R_session_info.txt: list of R packages and versions used in development & testing
- ./mock_data : mock data for code testing (zipped)
- ./mock_data_heritability_results : results of heritability analysis of mock data (zipped)
- ./mock_data_permutation_runs : example results of permutation runs of one taxon and one pathway
- ./Plots : heritability plots generated from heritability analysis of mock data