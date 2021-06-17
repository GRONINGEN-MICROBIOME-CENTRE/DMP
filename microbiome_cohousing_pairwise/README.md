## DMP microbiome clustering

The folder contains scripts used to test and plot pair-wise distances between microbiomes of cohousing, related and unrelated DMP participants (Fig 2C,D,E and Supplementary Figure 4). Script uses summary of pairwise bray-curtis distances (data file *microbiome_pairs_braycurtis.csv*)

### Dependencies: 

- R libraries ggplot2, dplyr, plyr. These packages should install correctly from CRAN by executing *install.packages(c("ggplot2", "dplyr", "plyr"))*

### Notes: 

- script uses summary statistics data provided in this github repo. Raw data used to generate summary statistics is not shared on EGA or github to protect privacy of study participants and is avaliable from Lifelines biobank (see manuscript for details)

### Files: 

- DMP_microbiome_distance_Fig2CDE_FigS4.R : script
- R_sessioninfo.txt : details of R environment used to test the codes
- *data* : input file (summary of pair-wise bray-curtis distances per tested group)
- *plots* : example results (plots for Fig.S2/c,d,e and Fig.S4)
