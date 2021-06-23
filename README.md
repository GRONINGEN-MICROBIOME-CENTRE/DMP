By: Weersma Group, Fu Group and Zhernakova Group, UMCG, The Netherlands 

Last update: 17/06/2021

## Dutch Microbiome Project

This github repo describes workflow and codes used in The Dutch Microbiome Project (DMP) study:

### Contents:

- Microbiome profiling
- Heritability analysis
- Identification of core and keystone microbes 
- Microbiome clustering
- Calculation of microbiome variance explained by phenotypes
- Microbiome-phenotypes association analyses
- Calculation of microbiome signatures predictive of diseases and health
- Calculation of General Microbiome Health Index (GMHI)
- Miscellaneous scripts
- Supporting R scripts
- Supporting data
- Installation instructions
- Code use instructions

### Microbiome profiling and description of microbiome

Metagenomes were profiled consistent with previous data analysis of 1000IBD and Lifelines-DEEP10 cohorts, as follows. KneadData tools (v0.5.1) were used to process metagenomic reads (in fastq format) by trimming the reads to PHRED quality 30 and removing Illumina adapters. Following trimming, the KneadData integrated Bowtie2 tool (v2.3.4.1) was used to remove reads that aligned to the human genome (GRCh37/hg19).
Taxonomic composition of metagenomes was profiled by MetaPhlAn2 tool (v2.7.2) using the MetaPhlAn database of marker genes mpa_v20_m200. Profiling of genes encoding microbial biochemical pathways was performed using the HUMAnN2 pipeline (v0.11.1) integrated with the DIAMOND alignment tool (v0.8.22), UniRef90 protein database (v0.1.1) and ChocoPhlAn pan-genome database (v0.1.1). As a final quality control step, samples with unrealistic microbiome composition (eukaryotic or viral abundance > 25% of total microbiome content or total read depth < 10 million) were excluded, leaving 8,208 samples for further analyses. Analyses were performed using locally installed tools and databases on CentOS (release 6.9) on the high-performance computing infrastructure available at our institution and using the MOLGENIS data platform2.

Rarefaction and extrapolation (R/E) sampling curves for estimation of total richness of species and genera in the population were constructed using a sample size-based interpolation/ extrapolation algorithm implemented in the iNEXT package for R. Variance of phyla in the population was compared to variance of pathways using F Test for comparison of variances and total variance of microbial taxa vs pathways were compared using permutation tests. 

- An example of metagenome processing jobs is provided in *microbiome_profiling* folder
- Rarefaction curves, variance comparison and other microbiome description codes are in *microbiome_description* folder

### Heritability analysis

We estimated the heritability of bacterial taxa and pathways using linear mixed models. In particular, we fitted the following model using the function relmatLmer from the lme4qtl package for R:

  Y ~ age + age squared + sex + read depth + stool frequency + 1|ID + 1|FAM + 1|cohousing 
  
where Y is the relative abundance of the bacterial taxon or pathway, transformed using the centred additive log-ratio (CLR) transformation, with the geometric mean calculated on the taxonomic level of species used as denominator. Age, age squared, sex, read depth and stool frequency are participant-specific factors modelled as fixed effect covariates, while the remaining three terms are random effects representing a polygenic additive effect (1|ID, equivalent to twice the kinship matrix), history of a familial shared environment (1|FAM, implemented as unique family identifier) and current co-housing (1|cohousing, implemented as unique housing location identifier). Narrow sense heritability was estimated as the proportion of variance explained by the polygenic additive effect over total variance, using the profile function of the same R package. We restricted the analysis of heritability to the relative abundances taxa and pathways present in at least 5% of individuals and focused on families in which at least two individuals had available microbiome data. Significance of fixed effects was calculated using type-II analysis of variance using Wald tests implemented in anova function for R package car, while the significance of random effects was calculated using likelihood ratio tests implemented in ranova function for R package lmerTest. Confidence intervals were estimated using the function profile from R stats package.

- Codes used for heritability are provided in *heritability_analysis_v2* 
- Codes used for pair-wise comparison of microbiome similarity for cohabiting and non-cohabiting participants is provided in *microbiome_cohousing_pairwise*

### Analysis of core and keystone microbiome features

To identify core microbial species and pathways, we used a bootstrapping-based selection approach. We randomly sampled 1% to 100% of the samples of the cohort a hundred times and calculated the standard deviation of the presence rate of each microbial species/pathway at different sampling percentages. Microbial features with a presence rate of more than 95% of samples were defined as the core microbiome. 
To analyse microbiome community structure, we constructed microbial species and pathway co-abundance networks using SparCC tool for network inference. Relative abundances of taxa were converted to estimated read counts by multiplying abundance percentages by total sequenced reads per sample after quality control. For pathway analysis, the read counts (RPKM) from HUMAnN2 were directly used for SparCC. Significant co-abundance was controlled at FDR 0.05 level using 100 permutations. In each permutation, the abundance of each microbial feature was randomly shuffled across samples.

- Scripts for identification of core microbiome and keystone feature are in *core_keystone_microbes* folder


### Microbiome clustering

To identify microbial clusters and assess the presence of gut enterotypes in our cohort, we performed the partitioning around the medoid method on the relative abundances of microbial species and used the Calinski-Harabasz index to select the optimal number of clusters, as previously published in a study of gut enterotypes. Enrichment of phenotypes in each cluster was assessed by logistic regression in R

- Codes used for clustering, plotting of clusters and enrichment analysis are in *microbiome_clustering* folder

### Calculation of microbiome variance explained by phenotypes

The microbiome composition variance explained by phenotypes was calculated by permutational multivariate analysis of variance using distance matrices, implemented in the adonis function for R package vegan (v.2.4-6), using 20,000 permutations and a Bray-Curtis distance matrix calculated using relative abundances of microbial species. A separate analysis was performed to calculate the microbiome functional potential explained by phenotypes using equivalent methodology. The functional dissimilarity matrix was calculated using the Bray-Curtis dissimilarity index calculated on the relative abundances of MetaCyc microbial biochemical pathways.

- Scripts used for calculation are in *microbiome_variance*

### Microbiome-phenotypes association analyses

Prior to the association analysis of phenotypes and microbiome features, the microbiome data was transformed using the clr transformation. The geometric mean for clr transformation of relative abundances of taxa was calculated on species-level and applied to higher levels. The associations between phenotypes and microbial features (microbial taxa, MetaCyc functional pathways, CARD and VFDB entities) were calculated using linear regression, adjusting for age, sex and BMI of the individual along with Bristol stool scale of the faecal sample and technical factors (DNA concentration, sequencing read depth, sequencing batch and sampling season). Benjamini-Hochberg correction was used to control for multiple testing with the number of tests equal to the number of tested feature - phenotype pairs. 

- *association_analysis* folder contains scripts used for association analysis

### Calculation of microbiome signatures predictive of diseases and health

We calculated the microbial signatures predictive of the 36 most common (Ncases > 100) diseases in our dataset. In addition, we defined a healthy phenotype as an absence of any self-reported disease. To build prediction models for common diseases, the dataset was randomly split into training (90%) and test (10%) sets. Next, we performed elastic net L1/L2 regularized regression (R package glmnet v.4.0) on the training set, using Shannon diversity, clr-transformed microbial taxa, clr-transformed MetaCyc bacterial pathways and age, sex and BMI as fixed covariates (not penalized in the models). The model for each disease was calculated independently using five-fold cross-validation to select the optimal lambda penalization factor (at L1/L2 mixing parameter alpha fixed at 0.5). The lambda with minimal cross-validation error was used in the downstream analysis. In total, we defined three probabilistic models: a null signature that only includes effects of general covariates (age, sex and BMI), a microbiome signature that includes all selected microbiome features and a combined signature that includes both the effects of microbiome features and general covariates.

- *health_disease_prediction* lists scripts used for training and testing of prediction models

### Calculation of General Microbiome Health Index (GMHI)

We calculated the recently developed Gut Microbiome Health Index (GMHI) for DMP data, using the parameters identified in GMHI study (Gupta et al., Nature Comm. 2020) on our data. 

- *ghmi* folder contains scripts used for calcuation of GMHI

### Miscellaneous scripts

- *misc_scripts* folder contains miscellaneous used in data analysis not covered by other groups. For example, scripts for comparison of diet questionnaires over time and assigment of gastrointestinal disorders using ROMEIII criteria

### Supporting R scripts

*r_scripts_library* folder contains R functions used for various tasks (such as generation of supplementary data, plots, data parsing...). These functions were developed by Weersma group for internal use and are provided as-is, without comprehensive documentation

### Supporting data

- Mock phenotypes: *Mock_data* folder contains mock data (phenotypes and microbiome data) intended for testing of codes avaliable on this github repo. *heritability_analysis_v2/mock_data* contains mock data intended for testing of heritability analysis codes. This data is zipped to reduce file size and should be unzipped into *heritability_analysis_v2/mock_data*
-	DMP metagenomics: The raw metagenomic sequencing data and processed taxonomy and pathway tables are avaliable on European Genome-Phenome Archive (EGA) under accession number EGAS00001005027 (https://ega-archive.org/studies/EGAS00001005027) 
- Phenotypes: A core set of clinical characteristics has is avaliable at EGA (study EGAS00001005027). In compliance with the informed consent and privacy regulations for the Lifelines biobank, extensive data for >200 clinical characteristics can only be requested via the LifeLines directly through https://www.lifelines.nl/researcher. For a detailed description of data available within the Lifelines biobank see: https://www.lifelines.nl/researcher/data-and-biobank/wiki.

### Installation and code use instructions

- *git clone*  will pull the data, alternatively download and unzip the repo
- Dependencies: scripts were tested on R (v.3.6.x) running on Ubuntu 20.04.2 LTS, CentOS Linux 7 and Windows 10. Older version of R might lack some of required libraries, while R 4.x might require installing older versions of libraries if these are not avaliable yet. README files in individual folders list dependencies and *R_session_info* files list detailed list of used R libraries. The majority of dependencies should install correctly using install.packages(<dependency list>), with the exception of *lme4qtl* package used in heritability analysis which requires additional steps (see https://github.com/variani/lme4qtl for details)
- Running scripts: scripts are intended to be executed from github root folder (DMP folder), for example: *Rscript ./microbiome_clustering/ClusterAnalysis_FigS3_MockData.R* for microbiome cluster analysis script. If using R studio, please make sure working directory is the location of DMP github root folder before running the script or replace *setwd* line in the script to point to appropriate location

### Code use instructions / reproducing DMP manuscript results
- Figure 1/B and Supplementary Fig. S1 : run *microbiome_description/DMP_rarefaction_Fig1B_S1_EGAdata.R*. Note: DMP microbiome data deposited on EGA is required to run the script. *microbiome_description/DMP_rarefaction_Fig1B_S1_Mockdata.R* is implementation using mock data avaliable with this repo intended for demonstration/testing without the need to download EGA data.
- Figure 1/C: reproducing this figure requires phenotype data avaliable through Lifelines and microbiome data deposited at EGA. For code demonstration we provide: *microbiome_description/DMP_microbiome_FigS2A_1C_EGAdata_MockPhenotypes.R* which will correctly reproduce ordination plot and fit it with *Mock data* phenotypes, and *DMP_microbiome_FigS2A_1C_Mockdata_MockPhenotypes.R* for implementation using mock microbome data and mock phenotypes which does not require downloading data from EGA.
- Supplementary Fig. S2/A : run *microbiome_description/DMP_microbiome_FigS2A_1C_EGAdata_MockPhenotypes.R* to fully reproduce the figure or *DMP_microbiome_FigS2A_1C_Mockdata_MockPhenotypes.R* for demonstration code which does not require downloading any extra data.
- Figure 2 and Fig S4: run *heritability_analysis_v2/DMP_heritability_v10_realdata_plotresults.R* to re-produce Figure 2/A-B, *microbiome_cohousing_pairwise/DMP_microbiome_distance_Fig2CDE_FigS4.R* to re-produce Figure 2/C-E and Fig. S4. Note: these codes plot results, see *heritability_analysis_v2* for codes that construct heritability models and implementation on mock data for demonstration purposes (full data can be obtained from EGA and Lifelines biobank)
- Figure 3/a and Fig S5: run *association_analysis/DMP_Fig3a.R* to generate elements used to construct Fig 3/a and Fig S5
- Figure 3/b : *association_analysis/DMP_association_analysis_plots_esizes.R* produces the plot from summary statistics of association analysis (distrubuted in this github repo and with DMP manuscript supplements). Association analysis requires DMP phenotypes and microbiome data (avaliable through EGA & Lifelines), 
- Figure 3/c : *microbiome_variance/DMP_variance_explained_plots_Fig3C.R* . Codes that perform variance explained analysis requires phenotype data (avaliable from Lifelines), demonstration codes implemented on mock data are: *microbiome_variance/Mockdata_DMP_adonis_taxa.R*, *microbiome_variance/Mockdata_DMP_adonis_pwys.R*
- Figure 4/a : *association_analysis/DMP_association_analysis_plots_esizes.R* produces the plot
- Figure 4/b : *health_disease_prediction/DMP_Fig4b.R* generates the plot from summary statistics included in this github repo, *health_disease_prediction/DMP_health_disease_prediction.R* generates predictive models and plots the results using mock data
- Supplementary Fig S3 : *microbiome_clustering/ClusterAnalysis_FigS3_EGAdata.R* performs clustering of microbiome (Fig S3/a-b requires microbiome data to be downloaded from EGA, Fig S3/c requires DMP phenotypes available from Lifelines). *microbiome_clustering/ClusterAnalysis_FigS3_MockData.R* is a demonstration code that performs analysis (Fig S3/a-c) on mock data bundled with this repo.
- Supplementary Fig S4 : run *microbiome_cohousing_pairwise/DMP_microbiome_distance_Fig2CDE_FigS4.R*
- Supplementary Fig S6 : run *gmhi/DMP_GMHI.Rmd* 
- Supplementary Fig S7-S10 : *association_analysis/DMP_association_analysis_plots_esizes.R* generates these plots
