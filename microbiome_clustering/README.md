## DMP microbiome clustering

The folder contains scripts used to execute microbiome clustering and enrichment analysis (Supplementary Figure 3)

### Dependencies: 

- R libraries ggplot2, ggridges, ggsci, randomcolorR, cluster, MASS, clusterSim, oddsratio, foreach. These packages should install correctly from CRAN by executing *install.packages(c("ggplot2", "ggridges", "ggsci", "randomcoloR", "cluster", "MASS", "clusterSim", "oddsratio","foreach"))* on Windows 10 system. On Ubuntu (tested on 18.04.5 LTS), installation of *randomcoloR* package requires *v8* and *curl* libs to build. These can be installed from command line using *sudo apt-get install libcurl4-openssl-dev libv8-dev*. 
- Mock Data implementation: Mock data (*Mock_data/taxa.txt*)
- EGA distributed DMP data: microbiome profiles, should be downloaded and unzipped into *DMP_data_EGA* folder
- Note: cluster enrichment analysis (Fig S3/C) uses private data of participants which could not be shared on this repo or EGA to maintain privacy of participants and codes here implement this analysis on the mock phenotype data. The DMP phenotype data can be requested from Lifelines biobank (see manuscript for details). 


### Folders: 

- Microbiome_clustering.R : script for calculation of microbiome clusters, optimal number of clusters and cluster composition
- Microbiome_clustering_plots.R : script for plotting of microbiome clustering results (Supplementary Figure 2/a)
- R_sessioninfo.txt : details of R environment used to test the codes
- *plots* : example results (plots for Fig S3/a-c using mock data, Fig S3/a-b for DMP data)
