# =============================================================
# =============================================================
# Fig 1/B & Fig S1: Rarefaction plots [implementation on Mock Data]
# =============================================================
# =============================================================

# ============== MAIN ==========================
# > set Working Directory
# NOTE: Make sure this path is correct
#       (should be (root) "DMP" folder location downloaded from github)
# =====================================
# example: 
#setwd('D:/Vbox/shared/dag/git_14_05/DMP/')
setwd('.')
# > choose between using mock data (default behavior)
#   and real DMP data:
useMockData = F

# > install packages if required
#install.packages(c("ggplot2","iNEXT","coin","vegan","pheatmap"))

# > load helper functions and libraries
# ===========================
library(ggplot2)
source('r_scripts_library/R_Microbiome_scripts.R')
source('r_scripts_library/R_Misc.R')

# LOAD DATA
# ===========================================
# if using mock data (default behavior)
if (useMockData) {
  mockS = 'MockData_'
  inDFmeta <- read.table('Mock_data/taxa.txt')
  inPWYs <- read.table('Mock_data/pathways.txt')
  inVFs <- read.table('Mock_data/VFs.txt')
  inCs <- read.table('Mock_data/CARDs.txt')
} else {
  # > load taxa
  mockS=''
  inDFmeta <- read.table('DMP_data_EGA/DMP_metaphlan2_raw.csv',sep=',',header=T)
  rownames(inDFmeta) <- inDFmeta$ID
  inDFmeta$ID <- NULL
  # > load pathways
  inPWYs <- read.table('DMP_data_EGA/DMP_humann2_metacyc_raw.csv',sep=',',header=T)
  rownames(inPWYs) <- inPWYs$ID
  # > load CARDs
  inCs <- read.table('DMP_data_EGA/DMP_CARD_raw.csv',sep=',',header=T)
  # > load VFs
  inVFs <- read.table('DMP_data_EGA/DMP_VFDB_raw.csv',sep=',',header=T)
}

# >> TAXA 
# ==================================================================
inDF <- inDFmeta
inDFmm <- filterMetaGenomeDF(inDF,presPerc = -1, minMedRelAb = -1, minMRelAb = -1, keepDomains = "All",keepLevels = c("S","G","F","O","C","P","K"))

nrBoots = 10
rcurves <- NULL
lvlN <- c("Genera","Species","Families")
cnt = 0
for (lvl in c("G","S","F")) {
  cnt = cnt + 1
  print (paste0('>> preparing curves for ',lvl))
  inDFlvl <- filterMetaGenomeDF(inDFmm,keepLevels = lvl,presPerc = -1,minMRelAb = -1,minMedRelAb = -1)
  inDFlvl$Cohort <- "DAG3"
  rcurve <- doRarefaction(inDF = inDFlvl,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F)
  rcurve$Taxon <- as.character(lvlN[cnt])
  rcurves <- rbind.data.frame(rcurves,rcurve)
}

# >> Pathways
# ==================================================================
# > prep for rarefaction
inPWYsCoh <- inPWYs
inPWYsCoh$Cohort <- "DAG3"
# > rarefy
print (' >> RAREFING PWYS')
rcurveP <- doRarefaction(inDF = inPWYsCoh,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F,doTaxa = F)
rcurveP$Taxon <- "Pathways"
rcurves <- rbind.data.frame(rcurves,rcurveP)

# >> VFs
# ==================================================================

# > prep for rarefaction & rarefy
inVFs$Cohort <- "DAG3"
print (' >> RAREFING VFs')
rcurveV <- doRarefaction(inDF = inVFs,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F,doTaxa = F)
rcurveV$Taxon <- "VFs"
rcurves <- rbind.data.frame(rcurves,rcurveV)

# >> CARDs
# ==================================================================
# > prep for rarefaction & rarefy
inCs$Cohort <- "DAG3"
print (' >> RAREFING CARDs')
rcurveC <- doRarefaction(inDF = inCs,replacements = F,bootstraps = nrBoots,extrapolate=F,doAll=F,doTaxa = F)
rcurveC$Taxon <- "CARDs"
rcurves <- rbind.data.frame(rcurves,rcurveC)

# extra annotation for plots
rcurves$linetype = "solid"
rcurves$linetype[rcurves$Taxon %in% c('CARDs','VFs','Pathways')] = "dashed"
cbPalette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#E69F00","#F0E442", "#999999")
rcurves$Taxon <- factor(rcurves$Taxon,levels=c("Species","Genera","Families","Pathways","VFs","CARDs"))

# save plot data table
# ======================
write.table(rcurves,paste0('microbiome_description/plots/',mockS,'Fig_1B_microbiome_rarefaction_datatable.csv'),sep=',',row.names = F)

# draw the plot
# ===================
g <- ggplot(rcurves,aes(x=nr,y=spec.nr.mn,col=Taxon)) + geom_line(size=1.25,linetype=rcurves$linetype) + 
  geom_errorbar(aes(ymin=spec.nr.mn-spec.nr.sd, ymax=spec.nr.mn+spec.nr.sd),size=1.05, colour="black") + geom_point(shape=21,size=1.75,fill="white") + theme_bw() + 
  xlab("Number of samples") + ylab("Number of features") + theme(legend.position="right") + 
  theme(text = element_text(size = 17)) + scale_color_manual(values = cbPalette)
print(g)
ggsave(plot = g,filename = paste0('microbiome_description/plots/',mockS,'Fig_1B_microbiome_rarefaction.png'),width = 6.5*1.2*1.2,height = 6*0.8*1.2)

# ========================================================
# ========================================================
# iNEXT library implementation of extrapolated rarefaction 
#    curves (supplementary Figure 1)
# ========================================================
# ========================================================
library (iNEXT)
if (useMockData) {
  dag3S <- filterMetaGenomeDF(inDFmm,keepLevels = "S",presPerc = -1,minMRelAb = 0.0,minMedRelAb = -1)
  dag3S.t <- t.data.frame(dag3S)
  dag3S.t.pa <- dag3S.t
  dag3S.t.pa[dag3S.t.pa > 0] <- 1
  dag3S.t.pa.rs <- rowSums(dag3S.t.pa)
  D_abund <- iNEXT(dag3S.t.pa.rs, datatype = 'abundance',knots = 250,endpoint = sum(dag3S.t.pa.rs)*1.25)
  D_abund$DataInfo$n <- 2000
  D_abund$iNextEst$m <- D_abund$iNextEst$m/sum(dag3S.t.pa.rs)*2000
} else {
  dag3S <- filterMetaGenomeDF(inDFmm,keepLevels = "S",presPerc = -1,minMRelAb = 0.000001,minMedRelAb = -1)
  dag3S.t <- t.data.frame(dag3S)
  dag3S.t.pa <- dag3S.t
  dag3S.t.pa[dag3S.t.pa > 0] <- 1
  dag3S.t.pa.rs <- rowSums(dag3S.t.pa)
  D_abund <- iNEXT (dag3S.t.pa.rs, datatype = 'abundance',knots = 100,endpoint = sum(dag3S.t.pa.rs)*3)
  D_abund$DataInfo$n <- 8208
  D_abund$iNextEst$m <- D_abund$iNextEst$m/sum(dag3S.t.pa.rs)*8208
}
gg <- ggiNEXT(D_abund, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE) + theme_classic() + ylab("Number of Species") + 
  xlab("Sample size") + theme(text = element_text(size = 18))
print(gg)
ggsave(plot=gg,filename = paste0('microbiome_description/plots/',mockS,'Fig_S1_microbiome_rarefaction_species_extrapolated.png'),width=6,height=6)

# GENERA
if (useMockData) {
  dag3G <- filterMetaGenomeDF(inDFmm,keepLevels = "G",presPerc = -1,minMRelAb = 0.0000,minMedRelAb = -1)
  dag3G.t <- t.data.frame(dag3G)
  dag3G.t.pa <- dag3G.t
  dag3G.t.pa[dag3G.t.pa > 0] <- 1
  dag3G.t.pa.rs <- rowSums(dag3G.t.pa)
  D_abundG <- iNEXT (dag3G.t.pa.rs, datatype = 'abundance',knots = 250,endpoint = sum(dag3G.t.pa.rs)*1.25)
  D_abundG$DataInfo$n <- 2000
  D_abundG$iNextEst$m <- D_abundG$iNextEst$m/sum(dag3G.t.pa.rs)*2000
} else {
  dag3G <- filterMetaGenomeDF(inDFmm,keepLevels = "G",presPerc = -1,minMRelAb = 0.000001,minMedRelAb = -1)
  dag3G.t <- t.data.frame(dag3G)
  dag3G.t.pa <- dag3G.t
  dag3G.t.pa[dag3G.t.pa > 0] <- 1
  dag3G.t.pa.rs <- rowSums(dag3G.t.pa)
  D_abundG <- iNEXT (dag3G.t.pa.rs, datatype = 'abundance',knots = 100,endpoint = sum(dag3G.t.pa.rs)*3)
  D_abundG$DataInfo$n <- 8208
  D_abundG$iNextEst$m <- D_abundG$iNextEst$m/sum(dag3G.t.pa.rs)*8208
}
gg <- ggiNEXT(D_abundG, type=1, se=TRUE, facet.var="none", color.var="site", grey=FALSE) + theme_classic() + ylab("Number of Genera") + 
  xlab("Sample size") + theme(text = element_text(size = 18))
print(gg)
ggsave(plot=gg,filename = paste0('microbiome_description/plots/',mockS,'Fig_S1_microbiome_rarefaction_genera_extrapolated.png'),width=6,height=6)
