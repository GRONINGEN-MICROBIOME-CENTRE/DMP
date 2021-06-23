## Prediction of health and diseases from microbiome data

The folder contains scripts used to train and test models for prediction of health and diseases on DMP data.

### Notes:

- real microbiome data used for DMP project can be obtained from EGA study EGAS00001005027, phenotypes used in DMP can be obtained from Lifelines biobank (see DMP manuscript for details)
- DMP_association_analysis_plots_esizes.R uses DMP summary statistics and re-produces DMP manuscript figures
- as DMP phenotypes could not be shared on EGA or github due to privacy regulations, the codes distributed here implement the workflow and analysis on the *Mock data* that is distributed with this repo for demonstration / testing purposes. Results of the analysis on the DMP data is provided in DMP supplementary material
- codes were developed & tested on R 3.6.x (see R_session_info.txt for list of used packages)

### Dependancies:

- R libraries: *glmnet, pROC, doSNOW, vegan, foreach, corrplot*
- these libs usually install correctly using: *install.packages(c(glmnet, pROC, doSNOW, vegan, foreach, corrplot))*
- Mock data (distributed in ../Mock_data) - these files are used for code testing / demonstration of association analysis (performed by DMP_health_disease_prediction.R). 

### Files:

- README.md : this readme file
- R_sessionInfo.txt : writedown of R environment used for script testing
- DMP_health_disease_prediction.R : script for constructing predictive models for diseases and health (lack of diseases), implemented on mock data included in this repo for code testing and evaluation; also produces Fig4b-equivalent with *mock data* results
- *Mockdata.out* : expected results for *DMP_health_disease_prediction.R* - results for *mock data* prediction of health (lack of diseases) and individual diseases using microbiome, clinical phenotypes and combination of microbiome.  phenotypes.
- diseases_signature_correlation.txt : summary statistics produced by running DMP_health_disease_prediction.R on real DMP data, used by DMP_Fig4b.R
- DMP_Fig4b.R : script for reproducing DMP plot Fig4b using real DMP data
- DMP_Fig4b.pdf : plot produced by running DMP_Fig4b.R
- Mockdata_Fig4b.pdf : plot produced by running DMP_health_disease_prediction.R