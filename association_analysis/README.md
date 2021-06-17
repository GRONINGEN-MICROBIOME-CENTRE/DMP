## DMP microbiome-phenotype association analysis

The folder contains scripts used to calculate associations between microbial features (taxa, pwys, virulence factors, antibiotic resistance genes) and phenotypes studied in DMP project

### Notes:

- mock_data folder contains *mock* data, this data is intended to be used only for testing of codes. It was designed to produce results similar to real data, but it is *not* real DMP data. 
- real microbiome data used for DMP project can be obtained from EGA study EGAS00001005027, real kinship matrices, family information and cohousing data can be obtained from Lifelines Biobank
- example results generated from mock data are provided (zipped) in mock_data_heritability_results. These results are expected output from DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R
- results permutation runs are NOT provided due to github space constrains, these can be generated using DMP_heritability_v10_mockdata_taxa.R & DMP_heritability_v10_mockdata_pwy.R (see Instructions)
- example consolidated results produced from mock data are provided in root folder. These are expected results from running DMP_heritability_v10_mockdata_collectdata.R on mock data heritability models
- example plots generated from mock data results are provided in Plots folder. These are expected results from running DMP_heritability_v10_mockdata_plotresults.R
- codes were developed & tested on R 3.6.0 (see R_session_info.txt for list of used packages)

### Dependancies:

- R libraries: *plyr, pheatmap, ggplot2, coin, vegan*
- these libs usually install correctly using: *install.packages('plyr','pheatmap', 'coin', 'ggplot2','pheatmap','vegan')*
- Mock Data (distributed with this repo, see *mock_data* folder, data has to be unzipped)
- Note: heritability analysis performed in DMP project used private data of participants which could not be shared on this repo or EGA to maintain privacy of participants. This data can be requested from Lifelines biobank (see manuscript for details)

### Files:

- README.md : this readme file
