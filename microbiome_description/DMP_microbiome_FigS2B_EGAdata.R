# ============================================
# ============================================
# Fig S2/b, implementation with REAL data
# ============================================
# ============================================

# NOTES:
# ====================================
# > this is implementation of script that generates
#   Fig S2/b (piechart with phyla) using real data 
#   (real data is avaliable on EGA and is not shared on github)

# ============== MAIN ==========================
# ==============================================
# > set WD: 
# NOTE: Make sure this path is correct
#       (should be (root) "DMP" folder location downloaded from github)
# NOTE2: make sure data is downloaded from EGA
#       and stored in DMP_data_EGA sub-folder
# =====================================
# example: 
#setwd('D:/Vbox/shared/dag/git_14_05/DMP/')
setwd('.')

# > load helper functions and libraries
# ===========================
library(ggplot2)
source('r_scripts_library/R_Microbiome_scripts.R')
shortenNames2 <- function(names) {
  outnames <- c()
  names <- as.character(names)
  cns <- (strsplit(names,'\\.?__'))
  for (c in cns) {
    outnames <- c(outnames,c[length(c)])
  }
  outnames
}
percent <- function(x, digits = 2, format = "f", ...) {     
  paste0(formatC(x * 100, format = format, digits = digits, ...), "%")
}
# pie chart maker
# ============================================
makePieChart <- function(inDFt,cutLabel=0.01,doLabels=F) {
  avg_DF <- data.frame()
  blank_theme <- theme_minimal()+
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.border = element_blank(),
      panel.grid=element_blank(),
      axis.ticks = element_blank(),
      plot.title=element_text(size=14, face="bold")
    )
  for (c in colnames(inDFt)[grep('__',colnames(inDFt))]) {
    avg_DF <- rbind.data.frame(avg_DF, data.frame(Taxon=c,mean=mean(inDFt[[c]]),sd=sd(inDFt[[c]]) ))
  }
  avg_DF$Taxon <- shortenNames2(avg_DF$Taxon)
  avg_DF$lbl <- avg_DF$mean
  avg_DF$NR <- 0
  avg_DF$lbl[avg_DF$lbl < cutLabel] <- 0
  avg_DF$lbl <- percent(avg_DF$lbl)
  avg_DF$lbl[avg_DF$lbl == "0.0%"] <- ""
  avg_DF <- avg_DF[order(avg_DF$mean,decreasing = T),]
  avg_DF$Taxon <- factor(avg_DF$Taxon, levels = avg_DF$Taxon[order(avg_DF$mean,decreasing = T)])
  c = 1
  for (r in c(1:nrow(avg_DF))) {
    if (avg_DF$lbl[r] != "") {
      avg_DF$NR[r] = c
      c = c + 1
    }
  }
  pie <- ggplot(avg_DF,aes(x="",y=mean,col=Taxon,fill=Taxon)) + geom_col() + coord_polar("y", start=0)
  pie <- pie + blank_theme + theme(axis.text.x=element_blank()) 
  if (doLabels) pie <- pie + geom_text(aes(y = mean,x=0.9 + NR %% 2 * 0.3,label = (lbl)), size=5,col="black")
  pie
}

# LOADDATA
# =========================
inDFmeta <- read.table('DMP_data_EGA/DMP_metaphlan2_raw.csv',sep=',',header = T)

# >> Pie chart with Taxa
# NOTE: numbers are saved to csv file, 
#       manuscript plot had the numbers added in post-processing
# ==================================================================
inDFmeta$ID <- NULL
inDFt = filterMetaGenomeDF(inDFmeta,keepLevels = "P",keepDomains = c("Bacteria","Archaea"),
                           presPerc = 0.001,minMRelAb = 0,rescaleTaxa = T)
pieTaxa <- makePieChart(inDFt,cutLabel=0.5,doLabels = F)

meanTaxa <- (apply(inDFt,MARGIN= 2,FUN= mean ))
sdTaxa <- (apply(inDFt,MARGIN = 2,FUN = sd ))
tblTaxa <- cbind.data.frame(meanTaxa,sdTaxa)
print(pieTaxa)
ggsave(plot = pieTaxa,filename = 'microbiome_description/plots/FigS2b_Piechart_Phyla.png',width = 10*1.1,height = 10)
write.table(tblTaxa,file = 'microbiome_description/plots/FigS2b_Piechart_Phyla_numbers.csv',sep=',')

#print(sessionInfo())