## DMP microbiome-phenotype association analysis

The folder contains scripts used to calculate associations between microbial features (taxa, pwys, virulence factors, antibiotic resistance genes) and phenotypes studied in DMP project. 

### Notes:

- real microbiome data used for DMP project can be obtained from EGA study EGAS00001005027, phenotypes used in DMP can be obtained from Lifelines biobank (see DMP manuscript for details)
- DMP_association_analysis_plots_esizes.R uses DMP summary statistics and re-produces DMP manuscript figures
- as DMP phenotypes could not be shared on EGA or github due to privacy regulations, the codes distributed here implement the workflow and analysis on the *Mock data* that is distributed with this repo for demonstration / testing purposes. Results of the analysis on the DMP data is provided in DMP supplementary material and results used in DMP figures are also provided in csv format in *data* folder (this data is sufficient to reproduce the figures).
- codes were developed & tested on R 3.6.x (see R_session_info.txt for list of used packages)

### Dependancies:

- R libraries: *plyr, pheatmap, ggplot2, coin, vegan, foreach, data.table, ggtree*
- these libs usually install correctly using: *install.packages(c('plyr','pheatmap', 'coin','ggplot2','pheatmap','vegan','foreach','data.table'))*
- *ggtree* package is part of bioconductor and can be installed using *BiocManager::install("ggtree")*; *BiocManager* can be installed from CRAN using *install.packages("BiocManager")*.
- DMP Associations summary statistics and phenotype grouping (distributed in *data* folder) - these files are used to re-produce DMP Figures 3/b, 4/a and Supplementary Figures S7 - S10
- Mock data (distributed in ../Mock_data) - these files are used for code testing / demonstration of association analysis (performed by DMP_runAssociations.R). 

### Files:

- README.md : this readme file
- DMP_association_analysis_plots_esizes.R : script for plotting microbiome-association heatmaps (DMP Figures 3/b, 4/a and Supplementary Figures S7 - S10)
- R_sessionInfo_DMP_association_analysis_plots.txt : writedown of R environment used for testing of DMP_association_analysis_plots_esizes.R
- DMP_Mockdata_runAssociations.R : script that performs microbiome-phenotype association analysis on *Mock data*. 
- DMP_ScriptsAssociationWithCorrections.R : functions implementing microbiome-phenotype association analysis performed in the DMP project. For demonstration / code testing purposes, analysis is implemented on *Mock data* distributed with this github repo. NOTE: this file contains functions implementating the analysis and is not intended for direct execution - *DMP_Mockdata_runAssociations.R* is wrapped around these scripts that performs the analysis
- DMP_Fig3a.R : script that generates graphic elements used to construct DMP Figures 3/A and S5.
- *data* : summary statistics of DMP taxa-phenotype associations analysis and phenotype grouping files these files are used by *DMP_association_analysis_plots_esizes.R* and *DMP_Fig3a.R*
- *plots* : plots produced by *DMP_association_analysis_plots_esizes.R* and *DMP_Fig3a.R*
- *MockData.output* : expected results of *DMP_Mockdata_runAssociations.R*
