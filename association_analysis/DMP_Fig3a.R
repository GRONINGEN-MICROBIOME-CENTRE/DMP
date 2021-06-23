# By: A. Kurilshikov, R. Gacesa (UMCG, 2021)
#
# Script plots DMP Fig 3/a, Fig S5
# (summary of results of association analysis)
# 
# ===========================================

library(foreach)
library(ggtree)
library(ggplot2)

# set wd:
# NOTE: make sure it is set to root of github rep (DMP folder)
#       (or leave it as '.' if executing script from command line from the DMP folder)
# example:
#setwd('C:/Users/Ranko/Documents/ub_shared/dag/gits/git_23_06/DMP')
setwd('.')

print(' > loading and prepping data ...')
# input file (summary statistics of microbiome associations, without corrections for anthropometrics)
v27_out = read.table("association_analysis/data/DMP_taxa_microbiome_associations.txt",
                     header=T,as.is= T,sep="\t")

species = v27_out[grep("s__",v27_out$taxon),]
all_layers = list()
for(i in c("^EXP.DIET.",
           "^ANTHRO",
           "^EXP.GREENSPACE|EXP.PETS|EXP.POLLUTION|EXP.SMOKING",
           "^EXP.BIRTH.|EXP.EARLYLIFE",
           "^DEMOGRAPHICS|SOCIOEC",
           "^MED.BLOOD|MED.INDICES|MED.URINE",
           "^MED.MEDS.|META.Antibiotics")) {
  #print(i)
  subset = species[grep(i,species[,1]),]
  
  include1=data.frame(taxon = subset$taxon, i = as.integer(subset$FDR<0.05), stringsAsFactors = F)
  
  all_layers[[i]] = data.frame(taxon = rownames(table(include1)), i = table(include1)[,2])
}

#proper direction for effect sizes
subset_diseases = species[grep("^MED.DISEASES|MED.SURGERY|MED.HEALTH.RAND",species[,1]),]
subset_diseases[grep("^MED.DISEASES|MED.SURGERY",subset_diseases$phenotype),"effect.size.asInteger"] = 
  as.numeric(sub(".*:","",subset_diseases[grep("^MED.DISEASES|MED.SURGERY",subset_diseases$phenotype),"effect.size"]))
subset_diseases$effect.size.asInteger = as.numeric(subset_diseases$effect.size.asInteger)
subset_diseases[grep("MED.DISEASES.None.No.Diseases",subset_diseases$phenotype),"effect.size.asInteger"] = 
  -subset_diseases[grep("MED.DISEASES.None.No.Diseases",subset_diseases$phenotype),"effect.size.asInteger"]
subset_diseases[grep("^MED.HEALTH.RAND.",subset_diseases$phenotype),"effect.size.asInteger"] = 
  -subset_diseases[grep("^MED.HEALTH.RAND.",subset_diseases$phenotype),"effect.size.asInteger"]

subset = subset_diseases
subset$FDR[sign(subset$effect.size.asInteger)>0] = 1
include1 =data.frame(taxon = subset$taxon, i = as.integer(subset$FDR<0.05),stringsAsFactors = F)
all_layers[["DISEASES.NEG"]] = data.frame(taxon = rownames(table(include1)), i = table(include1)[,2])

subset = subset_diseases
subset$FDR[sign(subset$effect.size.asInteger)<0] = 1
include1 =data.frame(taxon = subset$taxon, i = as.integer(subset$FDR<0.05),stringsAsFactors = F)
all_layers[["DISEASES.POS"]] = data.frame(taxon = rownames(table(include1)), i = table(include1)[,2])

table1 = Reduce(function(x,y)merge(x=x,y=y ,by="taxon",all = T), all_layers)
colnames(table1) = c("TAXON","Diet","Antropometrics","Exposome.current","Exposome.earlylife",
                     "Geography.Socioeconomics","Medical.measurements","Medication","Health/Disease.beneficial",
                     "Health/Disease.harmful")
rownames(table1) = table1[,1]
table1 = table1[,-1]

max_perGroup = apply(table1,2,max)
table2 = sweep(table1,2,max_perGroup,"/")
table2 = table2[,9:1]

print(' > generating plot elements ...')
#drawing
heatBar = function(x,multiplier.x = 0.95, multiplier.y = 0.9,grid.lwd = 0.3, 
                   colors = rev(c(
                     "green3",
                     "maroon3",
                     "forestgreen",
                     "maroon4",
                     "peachpuff4",
                     "brown3",
                     "lightblue4",
                     "seagreen3",
                     "red4"
                   ))){
  xlimit = nrow(x)
  ylimit = ncol(x)
  par(mar = c(0,0,0,0))
  plot(0,xlim = c(0,xlimit),ylim=c(0,ylimit),type = 'n')
#lines
  for (i in 0:xlimit) segments(x0 = i,x1 = i,y0=0,y1=ylimit,lwd = grid.lwd,col = "gray60")
  for (i in 0:ylimit) segments(x0 = 0, x1 = xlimit, y0 = i, y1 = i,lwd = grid.lwd,col = "gray60")
  
#rects
  for (i in 1:nrow(x)){
    for(j in 1:ncol(x)){
      if (x[i,j] > 0) {
        rect(
          xleft = i-1 + (1-multiplier.x)/2,
          ybottom = j-1 + (1-multiplier.y)/2,
          xright = i - (1-multiplier.x)/2,
          ytop = j-1 + (1-multiplier.y)/2 + multiplier.y * x[i,j],
          col = colors[j],
          border = NA)   
      }
    }
  }
}

#filtering taxa
names.matrix = do.call(rbind,strsplit(rownames(table1),split = "[.]"))
rownames(names.matrix) = rownames(table1)

names.matrix2 = names.matrix


dist.mat = foreach(i = 1:nrow(names.matrix2),.combine = rbind)  %:%
  foreach(j = 1:nrow(names.matrix2),.combine = cbind) %do% {
    sum(names.matrix2[i,] != names.matrix2[j,])
  }
rownames(dist.mat) = rownames(names.matrix2)
colnames(dist.mat) = rownames(names.matrix2)
dist.mat = dist.mat/7

tree = hclust(as.dist(dist.mat))
d = fortify(tree)
d = subset(d,isTip)

# save figure elements
# > matrix of numbers of associations
pdf("association_analysis/plots/Fig3a_matrix.pdf",height = 1.5,width = 6.5)
heatBar(table2[74:1,],multiplier.x = 0.75,multiplier.y = 0.8)
dev.off()
# > clustering tree
pdf("association_analysis/plots/Fig3a_tree.pdf",height = 1,width = 6)
ggtree(tree)+scale_x_reverse() + coord_flip()
dev.off()
# > labels for tree
write.table(tree$labels,file = "association_analysis/plots/Fig3a_labels.txt",row.names = F,col.names = F)
print(' > done!')
