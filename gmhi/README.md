## GMHI analysis

This folder contains scripts used to calculate Gut Microbiome Health Index (GMHI) on the DMP data. GMHI was introduced in study by Gupta et al. (Nature Comm. 2020; https://www.nature.com/articles/s41467-020-18476-8 & https://github.com/jaeyunsung/GMHI_2020)

Briefly, the GMHI classifier is first applied on the discovery cohort to find the parameters which optimizes the balanced accuracy, i.e. the prevalence fold change (_theta<sub>f</sub>_) and the prevalence difference (_theta<sub>d</sub>_). Then the same classifier is run on the validation cohort, and from those results, the balanced accuracy corresponding to the optimized parameters (i.e. _theta<sub>f</sub>_ and _theta<sub>d</sub>_) found on the discovery cohort is extracted.

**The .RMarkdown file contains:**

**ANALYSIS 1:** Run the GMHI classifier to predict disease status ("healthy" vs "unhealthy") based on our samples. We used a 90-10 training-testing split to generate a discovery (training data) and validation (testing data) cohort.
* **_Testing step 1a:_** Calculate GMHI using the optimized parameters found in the training step for each stool metagenome in the validation cohort to investigate whether the distributions of GMHI differ between "healthy" and "unhealthy" disease status
* **_Testing step 1b:_** Find the balance accuracy in the validation cohort that corresponds to the optimized _theta<sub>f</sub>_ and _theta<sub>d</sub>_ found in the training step  

**ANALYSIS 2:** Run the GMHI classifier on all samples for which we have species-level information using the optimized parameters found in Gupta et al. (2020).
