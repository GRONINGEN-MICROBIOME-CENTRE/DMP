# ==========================================
# By: Weersma Group, UMCG, 2021
#
# Script used for microbiome clustering
# on DMP microbiome data
# 
# NOTES: 
# - this is implementation of scripts for DMP Microbiome data avaliable at EGA
# (DMP manuscript Figures S3/a, S3/b)

# - reproduction of Fig S3/c on real data requires
# full set of phenotypes used in DMP data, this data
# can be requested from the Lifelines biobank
# (see manuscript for details)
# ==========================================

# load libraries
library(ggplot2)
library(ggridges)
library(ggsci)
library(randomcoloR)
library(cluster)
library(MASS)
library(clusterSim)
library(oddsratio)

# SETTINGS:
# =======================
# PLEASE SET wd to appropriate path (DMP folder of github repo)
# example:
#setwd('D:/Vbox/shared/dag/git_10_06/DMP/')
setwd('.')
# ============================================

# change the metaphlan result to a composition table, select top n most abundant features
CompositionTable <- function(x,n) { 
  require(foreach)
  #  x[is.na(x)]=0
  mean_value <- data.frame(Taxa=colnames(x), Mean_abundance=colSums(x)/nrow(x))
  most <- as.character(mean_value[order(mean_value$Mean_abundance,decreasing = T),]$Taxa[1:n])
  print(paste("Most abundant taxa is",most,sep = " "))
  
  composition_table <- foreach(i=1:length(most),.combine = rbind) %do%  {
    return.string = data.frame(ID = rownames(x), Relative=x[,most[i]],Level=colnames(x[,most[i],drop=F]))
  }
  
  first <- composition_table[grep(most[1],composition_table$Level),]
  first <- first[order(first$Relative,decreasing = T),]
  level <- as.factor(first$ID)
  composition_table$ID <- factor(composition_table$ID,levels = level)
  
  return(composition_table)
}

# calculator for JSD distance
dist.JSD <- function(inMatrix, pseudocount=0.000001, ...) {
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(inMatrix))
  matrixRowSize <- length(rownames(inMatrix))
  colnames <- colnames(inMatrix)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  inMatrix = apply(inMatrix,1:2,function(x) ifelse (x==0,pseudocount,x))
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(inMatrix[,i]),
                             as.vector(inMatrix[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

# partition around medoid (PAM) clustering
pam.clustering=function(x,k) { # x is a distance matrix and k the number of clusters
  cluster = as.vector(pam(as.dist(x), k, diss=TRUE)$clustering)
}

# ==================================
# species level distribution [Figure S3/A]
# ==================================
dag3=read.table("DMP_data_EGA/DMP_metaphlan2_filtered.csv",sep = ",",header = T,check.names = F,stringsAsFactors = F,row.names = 1)

dag3_species=dag3[,grep("s__",colnames(dag3))]
dag3_species=dag3_species[,grep("t__",colnames(dag3_species),invert = T)]
dag3_species <- dag3_species/100.0
colnames(dag3_species)=lapply(colnames(dag3_species),function(x){
  strsplit(x,"s__")[[1]][2]
})
dag3_species_plot=CompositionTable(dag3_species,10)
dag3_species_plot$Level=factor(dag3_species_plot$Level,levels = c("Bacteroides_vulgatus","Bacteroidales_bacterium_ph8","Alistipes_shahii",
                                                                  "Sutterella_wadsworthensis","Alistipes_onderdonkii","Bacteroides_uniformis",
                                                                  "Subdoligranulum_unclassified","Prevotella_copri",
                                                                  "Faecalibacterium_prausnitzii","Alistipes_putredinis"))
dag3_species_plot[dag3_species_plot==0]=NA
dag3_species_plot=na.omit(dag3_species_plot)
dag3_species_plot$Relative=-log2(dag3_species_plot$Relative)
set.seed(10)
n <- 10
palette <- distinctColorPalette(n)

g <- ggplot(dag3_species_plot, aes(x = Relative, y = Level,fill=Level,color=Level)) + theme_classic() +
  geom_density_ridges(alpha=0.8,scale = 2)+scale_fill_manual(values = palette) +
  scale_color_manual(values = palette)+ylab("")+xlab("Norm. Rel. Abundance (-log2)")

print(g)
ggsave(plot = g, filename = 'microbiome_clustering/plots/FigS3a.png')

# ==================================
# species clustering using PAM
# ==================================
dag3_species=as.data.frame(t(dag3_species))
dag3_species.dist=dist.JSD(dag3_species)
dag3_species.cluster=pam.clustering(dag3_species.dist, k=3)

# > perform clustering
nclusters = index.G1(t(dag3_species), dag3_species.cluster, d = dag3_species.dist, centrotypes = "medoids")
nclusters=NULL

# > calculate CH index to identify optimal number of clusters
for (k in 1:10) {
  if (k==1) {
    nclusters[k]=NA
  } else {
    dag3_species.cluster_temp=pam.clustering(dag3_species.dist, k)
    nclusters[k]=index.G1(t(dag3_species),dag3_species.cluster_temp,  d = dag3_species.dist,
                          centrotypes = "medoids")
  }
}
plot(nclusters, type="h", xlab="k clusters", ylab="CH index")
cluster=data.frame(row.names = colnames(dag3_species),Cluster=dag3_species.cluster)
write.table(cluster,"microbiome_clustering/plots/DMP_species_PAM.txt",row.names = T,quote = F,sep = "\t")
dag3_species=as.data.frame(t(dag3_species))

# ==================================
# P.copri abundance per cluster
# ==================================

clusters=read.table("microbiome_clustering/plots/DMP_species_PAM.txt",sep = "\t",header = T,stringsAsFactors = F)
clusters$sampleID=rownames(clusters)

cluster_plot=merge(clusters,dag3_species[,"Prevotella_copri",drop=F],by.x="sampleID",by.y = "row.names",all=F)
cluster_plot$Cluster=factor(cluster_plot$Cluster,levels = c("2","1"))
cluster_plot[cluster_plot==0]=NA
cluster_plot=na.omit(cluster_plot)
cluster_plot$Prevotella_copri=-log2(cluster_plot$Prevotella_copri)
cluster_plot[is.na(cluster_plot)]=0
g <- ggplot(cluster_plot, aes(x=Cluster, y=Prevotella_copri,fill=Cluster)) + 
  geom_violin()+
  geom_boxplot(width=0.1, fill="white")+scale_fill_npg()+theme_classic()+ylab("Prevotella Copri, Norm. Rel. Abundance (-log2)")
print(g)
ggsave(plot = g, filename = 'microbiome_clustering/plots/FigS3b.png')


# ===========================================
# associations between diseases and clusters
# ===========================================
# NOTE: calculating these associations requires
#  DMP phenotypes, these can be required from Lifelines biobank

# see ClusterAnalysis_FigS3_MockData.R for mock data variant of the code, 















