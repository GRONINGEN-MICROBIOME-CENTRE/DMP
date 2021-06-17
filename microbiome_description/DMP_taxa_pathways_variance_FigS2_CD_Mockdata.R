# =====================================================================
# By: Weersma Group (UMCG, 2021)
# DMP calculation of plots of microbiome taxonomy and pathway variance
# =====================================================================
#
# NOTES:
# - Mock data implementation of microbiome phyla and pathway variance

#setwd('D:/Vbox/shared/dag/git_10_06/DMP/')
setwd('.')

source('r_scripts_library/R_Microbiome_scripts.R')

library(ggplot2)
library(reshape2)
library(stringi)   
library(stringr)

# load MOCK DATA
# ========================
inDFtaxa <- read.table('heritability_analysis_v2/mock_data/DMP_mock_microbiome_taxa_filtered.csv',sep=',',header=T)
inPh <- filterMetaGenomeDF(inDFtaxa,keepLevels = "P",keepDomains = c("Bacteria"),
                           presPerc = 0.01,minMRelAb = 0.0002,rescaleTaxa = T)
inDFpwys <- read.table('heritability_analysis_v2/mock_data/DMP_mock_microbiome_pwys_filtered.csv',sep=',',header=T)
inDFpwysNoID <- inDFpwys
inDFpwysNoID$ID <- NULL

# TAXA plot
# ==============
taxa <- inPh
colnames(taxa) <- gsub('k__Bacteria\\.','',colnames(taxa) )
colnames(taxa) <- gsub('p__','',colnames(taxa) )
taxaNoID <- taxa; taxaNoID$ID <- NULL
taxa <- taxaNoID
taxa$ID <- row.names(taxa)
taxa$ID <- as.factor(taxa$ID)
taxa$ID <- factor(taxa$ID,levels = taxa$ID[order(taxa$Bacteroidetes,decreasing = T)])

# make it long
taxaL <- melt(taxa,id.vars = c("ID"))
taxaL$variable <- factor(taxaL$variable,levels = colnames(taxaNoID)[order(colMeans(taxaNoID),decreasing = T)])
taxaL <- taxaL[order(taxaL$variable,decreasing = T),]
taxaL <- taxaL[order(taxaL$ID,decreasing = T),]
#set colors
colPal <- c("#e64b35","#4dbbd5","#00a087","#3c5488","#f39b7f","#8491b4","#91d1c2","#dc0000","#7e6148","#b09c85")
# make plot
g <- ggplot(taxaL,aes(x=ID,y=value,fill=variable,col=variable)) + geom_col() + theme_classic() + xlab('') +
  scale_fill_manual(values=colPal) + scale_color_manual(values=colPal) + ylim(0,1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab('')
print(g)
ggsave(plot = g,filename = 'microbiome_description/plots/MockData_FigS2c.png',dpi = 600,width = 12,height = 6)

# pathways plot
# =================
# =================
# load pathway merger scheme
pathCats <- read.table('microbiome_description/data/tested_pathway_Name_class_cleaned_category_new.txt',quote='',header=T,check.names = F,sep='\t')
pathCats$BP2 <- gsub('-','.',pathCats$`Bacterial pathway`)

# align IDs
inDFpwys$ID <- rownames(inDFpwys)
testPwys <- inDFpwys[inDFpwys$ID %in% taxa$ID,]
testPwys$ID <- rownames(testPwys)
testPwys <- inDFpwys[inDFpwys$ID %in% taxa$ID,]
testPwys <- testPwys[taxa$ID,]

testPwys$ID <- rownames(testPwys)
testPwys$ID <- as.factor(testPwys$ID)

# merge pathways into pathway classes
colnames(testPwys) <- str_split_fixed(colnames(testPwys), pattern = "\\.\\.", 2)[,1]
pathCats2=pathCats[which(pathCats$BP2 %in% colnames(testPwys)),]
pathMrg=matrix(NA,nrow = nrow(testPwys),ncol = length(unique(pathCats2$class_name)))
# init matrix
colnames(pathMrg) <- as.character(unique(pathCats2$class_name))
rownames(pathMrg) <- rownames(testPwys)
for(i in 1:ncol(pathMrg)){
  # init pwy class
  pcN <- colnames(pathMrg)[i]
  # find matching pwys
  mpwys <- pathCats2$BP2[pathCats2$class_name==pcN]
  # grab and sum
  if (length(mpwys) >= 2) {
    pathMrg[,i] <- rowSums(testPwys[,mpwys])
  } else if (length(mpwys) == 1) {
    pathMrg[,i] <- testPwys[,mpwys]
  } else {
    pathMrg[,i] <- 0
  }
}
# grab top-9 and sum rest to others
pathMrg2 <- pathMrg
toKeep <- colMeans(pathMrg2)[order(colMeans(pathMrg2),decreasing = T)][1:9]
pathMrg2 <- pathMrg2[,colnames(pathMrg2) %in% names(toKeep)]
toKeep <- c(toKeep,"Other"=1-sum(toKeep))
pwyLvls <- names(toKeep)
others <- 1 - rowSums(pathMrg2)
pathMrg2 <- as.data.frame(pathMrg2)
pathMrg2$Other <- others

# make it long
pathMrg2$ID <- rownames(pathMrg2)
pwysL <- melt(pathMrg2,id.vars = c("ID"))
pwysL$variable <- factor(pwysL$variable,levels=pwyLvls)
pwysL <- pwysL[order(pwysL$variable,decreasing = T),]
pwysL <- pwysL[order(pwysL$ID,decreasing = T),]
g <- ggplot(pwysL,aes(x=ID,y=value,fill=variable,col=variable)) + geom_col() + theme_classic() +
  scale_fill_manual(values=colPal) + scale_color_manual(values=colPal) + xlab('') + ylim(0,1) +
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ylab('')
print(g)
ggsave(plot = g,filename = 'microbiome_description/plots/MockData_FigS2d.png',dpi = 600,width = 12,height = 6)

# ========== variance tests ===================
# > test variance of each pwy group VS each phylum
taxaN <- colnames(taxaNoID)
pwysN <- colnames(pathMrg2)[colnames(pathMrg2)!="ID"]

res <- NULL
for (tN in taxaN) {
  varT <- var(taxaNoID[[tN]])
  for (pN in pwysN) {
    varP <- var(pathMrg2[[pN]])

    tst <- var.test(taxaNoID[[tN]],
                    pathMrg2[[pN]], 
                    alternative = "greater")
    
    oneRow <- data.frame(taxon=tN,
                         pathway_class=pN,
                         taxon.var=varT,
                         pathway_class.var=varP,
                         var.ratio=format(tst$estimate,digits = 3),
                         var.ratio.CI=paste0(format(tst$conf.int[1],digits = 3),'-',format(tst$conf.int[2],digits=3)), 
                         var.F.stat=tst$statistic,
                         var.pV=tst$p.value)
    res <- rbind.data.frame(res,oneRow)
  }
}
rownames(res) <- NULL
res$var.FDR < p.adjust(res$var.pV)
res <- res[order(res$var.pV,decreasing = F),]
write.table(res,'microbiome_description/data/Mockdata_taxa_vs_pwys_variance_tests.csv',sep=',',row.names = F)
