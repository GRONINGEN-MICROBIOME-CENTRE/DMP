# ============================================================
# ============================================================
#  >> DMP pair-wise microbiome distance plots (Fig 2.c,d,e)
# ============================================================
# ============================================================
library(ggplot2)
library(plyr)
library(dplyr)

# > set WD: 
# NOTE: Make sure the working directory is 
#       set to root of git repo (DMP folder)
# =====================================
# example: 
#setwd('D:/Vbox/shared/dag/git_14_05/DMP/')
setwd('.')

# load helper scripts
source('r_scripts_library/R_Microbiome_scripts.R')

# LOAD DATA
allResAnnot <- read.table('microbiome_cohousing_pairwise/data/microbiome_pairs_braycurtis.csv',sep=',',header = T,stringsAsFactors = F)

# prep variable names for plots
toPlotVar <- c("BC_Spec","BC_PWY","BC_VF","BC_CARD") 
toPlotYs <- c("Bray-Curtis distance of Species",
              "Bray-Curtis distance of Pathways",
              "Bray-Curtis distance of Virulence factors",
              "Bray-Curtis distance of Antibiotic resistance")
# prep pallete
cbPalette <- c("#E69F00", "#CC79A7", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")
# text size
txtSize = 16

# >>> Fig 2/c, Supplementary Figure 4
# >>> children vs siblings vs rnd vs partners, cohabitating only
# ====================================================================
dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("PARTNERS","RND.PAIR","PARENT_CHILD","SIBLINGS"),]
dfToPlot <- dfToPlot[dfToPlot$COHAB | dfToPlot$RELATIONSHIP.0=="RND.PAIR",]
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("RND.PAIR","PARTNERS","PARENT_CHILD","SIBLINGS"))
dfToPlot$RELATIONSHIP.0 <- mapvalues(dfToPlot$RELATIONSHIP.0 , from = c("PARENT_CHILD", "SIBLINGS"), to = c("PAR_CH", "SIBL"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display = '',
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.005,ylim = c(0.0,1.2))
  # save plot dataframe
  write.table(t[[1]],paste0('microbiome_cohousing_pairwise/plots/microbiome_similarity_relationship_0_cohab_',tpv,'_datatable.csv'),sep=',',row.names = F)
  # make plot
  g <- t[[2]] + xlab("Relationship") + ggtitle('') + ylab(toPlotYs[c]) + theme_classic() +
     scale_color_manual(values = cbPalette) + theme(legend.position="none") + theme(text = element_text(size=txtSize)) + 
    ylim(0.0,1.3)
  print(g)
  ggsave(plot=g,filename = paste0('microbiome_cohousing_pairwise/plots/microbiome_similarity_relationship_0_cohab_',tpv,'.png'),width =5.5,height=6,dpi = 320)
}

# Figure 2/d, Supplementary Figure 4
# >>> children vs siblings vs rnd vs partners, non-cohabitating only
# ====================================================================
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#999999","#0072B2","#D55E00")
dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("RND.PAIR","PARENT_CHILD","SIBLINGS"),]
dfToPlot <- dfToPlot[!dfToPlot$COHAB | dfToPlot$RELATIONSHIP.0=="RND.PAIR" | dfToPlot$RELATIONSHIP.0=="PARTNERS",]
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("RND.PAIR","PARENT_CHILD","SIBLINGS"))
dfToPlot$RELATIONSHIP.0 <- mapvalues(dfToPlot$RELATIONSHIP.0 , from = c("PARENT_CHILD", "SIBLINGS"), to = c("PAR_CH", "SIBL"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display = '',ylim = c(0.0,1),
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.005)
  
  # save plot dataframe
  write.table(t[[1]],paste0('microbiome_cohousing_pairwise/plots/microbiome_similarity_relationship_0_noncohab_',tpv,'_datatable.csv'),sep=',',row.names = F)
  # make plot
  g <- t[[2]] + xlab("Relationship (cohabiting pairs)") + ggtitle('') + ylab(toPlotYs[c]) + theme_classic() +
     scale_color_manual(values = cbPalette) + theme(legend.position="none") + theme(text = element_text(size=txtSize)) + ylim(c(0.0,1.0))
  print(g)
  ggsave(plot=g,filename = paste0('microbiome_cohousing_pairwise/plots/microbiome_similarity_relationship_0_noncohab_',tpv,'.png'),width =4.5,height=6,dpi = 320)
}

# Fig 2/e, Supplementary Figure 4
# >>> 1st-degree (cohab) vs 1st-degree (sep) vs random
# ======================================================================
cbPalette2 <- c("#E69F00", "#555555","#D55E00","#56B4E9", "#009E73", "#CC79A7","#F0E442","#0072B2")

dfToPlot <- allResAnnot[allResAnnot$RELATIONSHIP.0 %in% c("RND.PAIR","SIBLINGS","PARENT_CHILD"),]
dfToPlot$RELATIONSHIP.0[(dfToPlot$RELATIONSHIP.0=="SIBLINGS" | dfToPlot$RELATIONSHIP.0=="PARENT_CHILD") & dfToPlot$COHAB] <- "1stDEG.COH"
dfToPlot$RELATIONSHIP.0[(dfToPlot$RELATIONSHIP.0=="SIBLINGS" | dfToPlot$RELATIONSHIP.0=="PARENT_CHILD") & !dfToPlot$COHAB] <- "1stDEG.SEP"
dfToPlot$RELATIONSHIP.0 <- factor(dfToPlot$RELATIONSHIP.0,levels=c("RND.PAIR","1stDEG.SEP","1stDEG.COH"))
c = 0
for (tpv in toPlotVar) {
  c <- c + 1
  t <- testOneFeature(dataIn=dfToPlot,feature = tpv,responseVar = "RELATIONSHIP.0",display='',
                      saveFolder = F,doPlots = T,doSave = F,retPlot = T,cutoff = 0.05,ylim = c(0.0,1.25))
  # save plot dataframe
  write.table(t[[1]],paste0('microbiome_cohousing_pairwise/plots/microbiome_similarity_relationship_1stDeg_',tpv,'_datatable.csv'),sep=',',row.names = F)
  # make plot
  g <- t[[2]] + xlab("Relationship (cohabiting vs non-cohabitating pairs)") + 
    ggtitle('') + ylab(toPlotYs[c]) + theme_classic() + scale_color_manual(values = cbPalette2) + theme(legend.position="none") + theme(text = element_text(size=txtSize)) + 
    ylim(0.0,1.05)
  print(g)
  ggsave(plot=g,filename = paste0('microbiome_cohousing_pairwise/plots/microbiome_similarity_relationship_1stDeg_',tpv,'.png'),width =4.5,height=6,dpi = 320)
}

#sessionInfo()