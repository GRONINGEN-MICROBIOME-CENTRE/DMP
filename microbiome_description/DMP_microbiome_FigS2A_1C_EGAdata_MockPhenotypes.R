# =====================================================================
# DMP, codes for bray-curtis ordination analysis and plotting Fig S2/A
# =====================================================================

# NOTES:
# ====================================
# > this is implementation of script that generates
#   Fig S2/a (prevotella-colored beta diversity plot) using real data 
#   (real data is avaliable on EGA and is not shared on github)

# ============== MAIN ==========================
# ==============================================
# > set WD: 
# NOTE: Make sure this path is correct
#       (should be (root) "DMP" folder location downloaded from github)
# NOTE2: make sure the data is downloaded from EGA
#       and stored in DMP_data_EGA sub-folder
# =====================================
# example: 
#setwd('D:/Vbox/shared/dag/git_14_05/DMP/')
setwd('.')

# ==== load libraries ====
library(vegan)
library(ggplot2)
library(viridis)
source('r_scripts_library/R_Microbiome_scripts.R')

# data location
inMicrobiome <- 'DMP_data_EGA/DMP_metaphlan2_raw.csv'

# ======================================================================================================================================
# > CALCULATE B/C matrices
# ======================================================================================================================================
print('============================================================')
print(' >>>>         PREPARING B/C MATRIX                     <<<< ')
print('============================================================')

print('   > loading microbiome ... ')
inMBss <- read.table(inMicrobiome,header=T,sep=',')#[1:2000,]
rownames(inMBss) <- inMBss$ID
inMBss$ID <- NULL

# work with species
print('    >> subsetting taxa')
inMBsss <- filterMetaGenomeDF(inDF=inMBss,presPerc=-1,minMRelAb=-1,minMedRelAb=-1,keepLevels=c("S"),rescaleTaxa = F,
                                keepDomains = "All")

print('   > calculating B/C distance matrix')  
distMatrix <- vegdist(inMBsss,method = "bray")
print('    >> done ')

# ======================================================================================================================================
# > CALCULATE PCoA
# ======================================================================================================================================
# prep ordination (PcoA)
print('   > calculating PCoA (this takes a while)... ')
pCoa <- cmdscale( distMatrix, eig = T,k = 2 )
# calculate variance explained
varExp <- (eigenvals(pCoa)/sum(eigenvals(pCoa)))[1:2]
xVar <- as.numeric(varExp[1]*100)
yVar <- as.numeric(varExp[2]*100)
  # get axes
pCoaVecs <- as.data.frame(pCoa$points)
pCoaVecsDF <- cbind.data.frame(data.frame(PCoA=c(1:2)),data.frame(VarExplained=varExp))
colnames(pCoaVecs) <- paste0("PCo",c(1:2))
pCoaVecs$ID <- row.names(pCoaVecs)
print(' >> DONE')

# make plot (PCo1 vs PCo2)
# ==================================================
inMBssID <- inMBss
inMBssID$ID <- row.names(inMBssID)
inMBssPrev <- inMBssID[,c("ID","k__Bacteria.p__Bacteroidetes.c__Bacteroidia.o__Bacteroidales.f__Prevotellaceae.g__Prevotella.s__Prevotella_copri")]
colnames(inMBssPrev) <- c("ID","P.copri")
pCoaVecsTaxa <- merge(pCoaVecs,inMBssPrev,by="ID")
g <- ggplot(pCoaVecsTaxa,aes(x=PCo1,y=PCo2,color=P.copri)) + geom_point(size=1.25) + theme_classic() + scale_color_gradientn(colours = rainbow(5)) + 
  xlab(paste0('PCo1 ',round(xVar,2),' %')) + ylab(paste0('PCo2 ',round(yVar,2),' %')) + theme(text = element_text(size = 16))
# save
ggsave(plot=g,filename = 'microbiome_description/plots/FigS2a_PCoA.pdf')
ggsave(plot=g,filename = 'microbiome_description/plots/FigS2a_PCoA.png',height = 5,width = 5.5, dpi = 600)

# save PCoA (used in biplot)
write.table(pCoaVecsDF,'temp_files/PCoA_varExp.csv',sep=',')
#write.table(pCoaVecsTaxa,'temp_files/PCoA_ForPlot.csv',sep=',')
saveRDS(pCoa,'temp_files/PCoA_ForPlot.RDS')

# =========== ADD Phenotype VECTORS FOR BIPLOT (Fig 1/c) ===============================
# NOTE: manuscript plot requires full phenotype data to generate this plot, 
#       this data can be requested from Lifelines biobank
# This code is functional variant which uses MOCK data with similar distributions
# to produce results approximately equal to manuscript plot 1/c
# =====================================================================================
# REAL DATA
#phe <- read.table('../../phenos/v27/DAG3_metadata_merged_ready_v27.csv',sep=',',header=T)

# MOCK DATA
phe <- read.table('Mock_data/MOCK_phenotypes_Fig1c.csv',sep=',',header=T)

# read microbiome (EGA)
inMicrobiome <- 'DMP_data_EGA/DMP_metaphlan2_raw.csv'
print('   > loading microbiome ... ')
inMBss <- read.table(inMicrobiome,header=T,sep=',')#[1:2000,]
rownames(inMBss) <- inMBss$ID
inMBss$ID <- NULL

# read PCoA (used in biplot)
pCoaVecsDF <- read.table('temp_files/PCoA_varExp.csv',sep=',')
pCoa <- readRDS('temp_files/PCoA_ForPlot.RDS')
#pCoaVecsTaxa <- read.table('temp_files/PCoA_ForPlot.csv',sep=',')
xVar <- pCoaVecsDF$VarExplained[1]*100
yVar <- pCoaVecsDF$VarExplained[2]*100

rownames(phe) <- phe$DAG3_sampleID
sp <- inMBss

# # align data between tables
phe3=phe[row.names(phe)%in%row.names(sp),]
phe3$SEX=1
phe3$SEX[phe3$ANTHRO.Sex=="F"]=2
phe3$HEALTH2 <- 1
phe3$HEALTH2[ phe3$MED.HEALTH.RAND.Health.General=="poor" ]=2
phe3$HEALTH2[ phe3$MED.HEALTH.RAND.Health.General=="mediocre" ]=3
phe3$HEALTH2[ phe3$MED.HEALTH.RAND.Health.General=="good" ]=4
phe3$HEALTH2[ phe3$MED.HEALTH.RAND.Health.General=="very good" ]=5
phe3$HEALTH2[ phe3$MED.HEALTH.RAND.Health.General=="excellent" ]=6
phe3=phe3[match(rownames(sp),rownames(phe3)),]
phe3$BMI <- phe3$ANTHRO.BMI
phe3$AGE <- phe3$ANTHRO.AGE
phe3$BRISTOL <- phe3$META.POOP.BristolMean
phe4=phe3[,c("AGE","BMI","BRISTOL","SEX","HEALTH2")]

# prep PCoA
myPCOAv2=pCoa

# fit phenotype vectors on the PCOA data
en2=envfit(myPCOAv2,phe4, na.rm=T, permutations=9999)
en_coord_cont = as.data.frame(scores(en2, "vectors")) * ordiArrowMul(en2)
#en_coord_sp = as.data.frame(scores(en2, display = "species"))
coord2=en_coord_cont/2
row.names(coord2)=c("AGE","BMI","BRISTOL STOOL SCALE","SEX","HEALTH STATUS")

# plot PCOA with phenotype vectors
# NOTE: in the manuscript plot, positions of labels were slightly adjusted in post-processing for visibility
g <- ggplot(pCoaVecsTaxa,aes(x=PCo1,y=PCo2,color=P.copri)) + geom_point(size=1.25,alpha=0.5) + theme_classic() + 
  scale_color_viridis(direction = -1) +
  #scale_color_gradient2(low="yellow",mid="green",high="darkblue",midpoint = 0.25) + 
  xlab(paste0('PCo1 ',round(xVar,2),' %')) + ylab(paste0('PCo2 ',round(yVar,2),' %')) +  
  geom_segment(aes(x = 0, y = 0, xend = Dim1*0.5, yend = Dim2*0.5), data = coord2, size =1, colour = "black",
               arrow = arrow(length = unit(0.03, "npc")) ) +
  geom_label(data = coord2, aes(x = Dim1*0.55+(sign(Dim1)*0.035), y = Dim2*0.55+(sign(Dim1)*0.035)), colour = "black", 
             fontface = "bold", label = row.names(coord2)) + 
  xlim(-0.70,0.25) + ylim(-0.35,0.5) + theme(text = element_text(size = 12))
print(g)

# save
ggsave(plot=g,filename = 'microbiome_description/plots/FigS1c_PCoA_MockData.png',height = 5,width = 5.5, dpi = 600)


