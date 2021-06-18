# ========================================================================
# By: R.Gacesa, Weersma Group, UMCG (2020)
#
# Script plots results of heritability analysis
# on real data of DMP project (see real_data sub folder / supplementary tables S5)
# (Panels A & B in Figure 2 in main DMP manuscript )

# =========================================================================

# parse confidence intervals
parseCIrange <- function(varN,toget="range") {
  ret = c()
  for (cc in c(1:nrow(inDFs))) {
    if (!is.na(inDFs[[varN]][cc]) & ! grepl('Inf',inDFs[[varN]][cc])) {
      ccc <- unlist(strsplit(inDFs[[varN]][cc],'-') )
      if (toget=="range") {
        ret <- c(ret,abs(as.numeric(ccc[2])-as.numeric(ccc[1])) )
      } else if (toget=="high") {
        ret <- c(ret,as.numeric(ccc[2]))
      } else if (toget=="low") {
        ret <- c(ret,as.numeric(ccc[1]))
      }
    } else {
      varms <- paste0('MS.',varN)
      ccc <- unlist(strsplit(inDFs[[varms]][cc],'-') )
      if (toget=="range") {
        ret <- c(ret,abs(as.numeric(ccc[2])-as.numeric(ccc[1])) )
      } else if (toget=="high") {
        ret <- c(ret,as.numeric(ccc[2]))
      } else if (toget=="low") {
        ret <- c(ret,as.numeric(ccc[1]))
      }
    }
  }
  ret
}

library(ggplot2)
library(tidyr)
# NOTES:
# ======================
# load log & results file
# set working directory: NOTE this has to be set to appropriate path!
#setwd('# placeholder path DAG3_heritability - set this one to path with this file! #')

#setwd('D:/Vbox/shared/dag/git_14_05/DMP/')
setwd('.')

if (!dir.exists('heritability_analysis_v2/real_data_plots')) {
  dir.create('heritability_analysis_v2/real_data_plots')
}

# input: heritability results (taxa)
inDFm <- read.table('heritability_analysis_v2/real_data/DMP_heritability_withFDRs_and_CIs_taxa.csv',sep=',',header=T,quote = '"',fill = T,stringsAsFactors = F)

# make plots (taxa)
inDFm$FTYPE <- NA
inDFm$FTYPE[grep('^s_',inDFm$Trait.short)] <-"Taxon.S"
inDFm$FTYPE[grep('^g_',inDFm$Trait.short)] <-"Taxon.G"
inDFm$FTYPE[grep('^f_',inDFm$Trait.short)] <-"Taxon.F"
inDFm$FTYPE[grep('^c_',inDFm$Trait.short)] <-"Taxon.C"
inDFm$FTYPE[grep('^p_',inDFm$Trait.short)] <-"Taxon.P"
inDFm$FTYPE[grep('^k_',inDFm$Trait.short)] <-"Taxon.K"
inDFm$FTYPE[grep('^o_',inDFm$Trait.short)] <-"Taxon.O"

# debug output
#print(inDFm$FEATURE[is.na(inDFm$FTYPE)])

fTypes <- c("Taxon.S","Taxon.G","Taxon.F","Taxon.C","Taxon.O","Taxon.P")
#cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73")#, "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbPalette <- c("#E69F00","#999999", "#009E73","#56B4E9")#, "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

for (oneType in fTypes) {
#  oneType <- "Taxon.S"
  
  inDFs <- inDFm[inDFm$FTYPE==oneType,]
  inDFs$Taxon <- factor(as.character(inDFs$Trait.short),levels=inDFs$Trait.short[order(inDFs$VE_ID)])
  inDFs <- inDFs[order(inDFs$VE_ID),]
  inDFs$LBL = ""
  inDFs$LBL[inDFs$SW_PV_ID <= 0.05] <- "*"
  inDFs$LBL[inDFs$SW_FDR_ID <= 0.1] <- "**"
  inDFs$H2 <- inDFs$VE_ID
  inDFs$VE_Cohousing <- inDFs$VE_COHOUSING.ID_DMP
  inDFs$VE_Family <- inDFs$VE_famID
  inDFs$VE_Environment <- inDFs$VE_Residual
  if (grepl('Taxon',oneType)) {
    inDFs$Taxon_Shortname <- as.character(inDFs$Taxon)
    inDFs$Taxon_Shortname <- factor(as.character(inDFs$Taxon_Shortname),levels=inDFs$Taxon_Shortname[order(inDFs$VE_ID)])
  }
  # parse confidence intervals
  inDFs$CI_ID_low <- parseCIrange("MS.CI_ID","low")
  inDFs$CI_ID_low[inDFs$CI_ID_low < 0] <- 0
  inDFs$CI_ID_high <- parseCIrange("MS.CI_ID","high")
  inDFs$CI_ID_high[inDFs$CI_ID_high > 1] <- 1
  
  inDFtoPlot <- inDFs[,c("Taxon_Shortname","LBL","H2","VE_Cohousing","VE_Family","VE_Environment","CI_ID_low","CI_ID_high")]
  #inDFtoPlot$VE_Environment <- 1-inDFtoPlot$H2-inDFtoPlot$VE_Cohousing-inDFtoPlot$VE_Family
  #inDFtoPlot$SEH2 <- inDFs$SEH2
  inDFtoPlot$LBL <- inDFs$LBL
  inDFtoPlotL <- gather(inDFtoPlot,"Var.Exp","Var.Exp.NR", H2:VE_Environment,factor_key = T)
  inDFtoPlotL$Var.Exp <- as.character(inDFtoPlotL$Var.Exp)
  inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "H2"] <- "Additive genetics"
  inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VE_Environment"] <- "Environment"
  inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VE_Cohousing"] <- "Cohousing"
  inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VE_Family"] <- "Family"
  
  inDFtoPlotL$Var.Exp <- factor(as.character(inDFtoPlotL$Var.Exp),level = c("Environment","Cohousing","Family","Additive genetics"))
  inDFtoPlotL$Var.Exp <- factor(as.character(inDFtoPlotL$Var.Exp),level = c("Cohousing","Family","Environment","Additive genetics"))
  
  # save data frame used to make plot
  write.table(inDFtoPlotL,paste0('heritability_analysis_v2/real_data_plots/plot_datatable_heritability_v2_taxa_',oneType,'_all.csv'),sep=',',row.names = F)
  
  # make plot
  g <- ggplot(inDFtoPlotL,aes(x=Taxon_Shortname,y=Var.Exp.NR,fill=Var.Exp)) + 
    scale_fill_manual(values = cbPalette) +
    geom_col(col="black", width=1,size=0.75) + 
    theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(-0.01,1.01) +
    theme(axis.text.y = element_text(face="bold")) +
    geom_errorbar(ymin=inDFtoPlotL$CI_ID_low,ymax=inDFtoPlotL$CI_ID_high,width=0.25, linetype='solid') + 
    geom_text(data = inDFtoPlotL,
              aes(x = Taxon_Shortname, y=CI_ID_high, 
                  label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
              color="black", vjust=+0.75, angle = 0, hjust=-1,size=6) + ylim(-0.01,1.01) + 
    ylab('Microbiome variance explained') + xlab('') + 
    theme(legend.position="bottom") + 
    theme(text = element_text(size = 14)) + 
    coord_flip() + theme_classic()
  print(g)
  ggsave(plot = g, filename = paste0('heritability_analysis_v2/real_data_plots/heritability_Rev1_',oneType,'.png'),height = 1.25+8/50*nrow(inDFs),width = 9,limitsize = F)
  
  # smaller plot (top 20 signals), high resolution for manuscript
  topN <- 20
  if (nrow(inDFs) > topN) {
    topFeatures <- inDFs[order(inDFs$VE_ID,decreasing = T),]$Trait.short[1:topN]
    inDFtoPlotLs <- inDFtoPlotL[inDFtoPlotL$Taxon_Shortname %in% topFeatures,]
    # save data
    write.table(inDFtoPlotLs,paste0('real_data_plots/plot_data_heritability_v2_taxa_',oneType,'_top.csv'),sep=',',row.names = F)
    # make plot
    g <- ggplot(inDFtoPlotLs,aes(x=Taxon_Shortname,y=Var.Exp.NR,fill=Var.Exp)) + 
      scale_fill_manual(values = cbPalette) +
      geom_col(col="black", width=1,size=0.75) + 
      theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(-0.01,1.01) +
      theme(axis.text.y = element_text(face="bold")) +
      geom_errorbar(ymin=inDFtoPlotLs$CI_ID_low,ymax=inDFtoPlotLs$CI_ID_high,width=0.25, linetype='solid') + 
      geom_text(data = inDFtoPlotLs,
                aes(x = Taxon_Shortname, y=CI_ID_high,
                    label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
                color="black", vjust=+0.75, angle = 0, hjust=-1,size=6) + ylim(-0.01,1.01) +
      ylab('Microbiome variance explained') + xlab('') + 
      theme(legend.position="bottom") + 
      coord_flip() + theme_classic() +
      theme(text = element_text(size = 30))
    print(g)
    ggsave(plot = g, filename = paste0('heritability_analysis_v2/real_data_plots/heritability_Rev1_',oneType,'_top.png'),height = 2.50+16/50*topN,width = 18,limitsize = F)
  }
}

# ================================
# PATHWAY PLOTS
# ================================
# input: heritability results (taxa)
inDFm <- read.table('heritability_analysis_v2/real_data/DMP_heritability_withFDRs_and_CIs_pwys.csv',sep=',',header=T,quote = '"',fill = T,stringsAsFactors = F)

# make plots (pathways)
inDFm$FTYPE <- "PWYS"
cbPalette <- c("#E69F00","#999999", "#009E73","#56B4E9")#, "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

oneType <- "PWYS"

inDFs <- inDFm[inDFm$FTYPE==oneType,]
inDFs$Taxon <- factor(as.character(inDFs$Trait),levels=inDFs$Trait[order(inDFs$VE_ID)])
inDFs <- inDFs[order(inDFs$VE_ID),]
inDFs$LBL = ""
inDFs$LBL[inDFs$SW_PV_ID <= 0.05] <- "*"
inDFs$LBL[inDFs$SW_FDR_ID <= 0.1] <- "**"
inDFs$H2 <- inDFs$VE_ID
inDFs$VE_Cohousing <- inDFs$VE_COHOUSING.ID_DMP
inDFs$VE_Family <- inDFs$VE_famID
inDFs$VE_Environment <- inDFs$VE_Residual
inDFs$Taxon_Shortname <- as.character(inDFs$Taxon)
inDFs$Taxon_Shortname <- factor(as.character(inDFs$Taxon_Shortname),levels=inDFs$Taxon_Shortname[order(inDFs$VE_ID)])

inDFs$CI_ID_low <- parseCIrange(varN = "CI_ID","low")
inDFs$CI_ID_low[inDFs$CI_ID_low < 0] <- 0
inDFs$CI_ID_high <- parseCIrange("CI_ID","high")
inDFs$CI_ID_high[inDFs$CI_ID_high > 1] <- 1

inDFtoPlot <- inDFs[,c("Taxon_Shortname","LBL","H2","VE_Cohousing","VE_Family","VE_Environment","CI_ID_low","CI_ID_high")]
#inDFtoPlot$VE_Environment <- 1-inDFtoPlot$H2-inDFtoPlot$VE_Cohousing-inDFtoPlot$VE_Family
#inDFtoPlot$SEH2 <- inDFs$SEH2
inDFtoPlot$LBL <- inDFs$LBL
inDFtoPlotL <- gather(inDFtoPlot,"Var.Exp","Var.Exp.NR", H2:VE_Environment,factor_key = T)
inDFtoPlotL$Var.Exp <- as.character(inDFtoPlotL$Var.Exp)
inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "H2"] <- "Additive genetics"
inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VE_Environment"] <- "Environment"
inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VE_Cohousing"] <- "Cohousing"
inDFtoPlotL$Var.Exp[inDFtoPlotL$Var.Exp == "VE_Family"] <- "Family"

inDFtoPlotL$Var.Exp <- factor(as.character(inDFtoPlotL$Var.Exp),level = c("Environment","Cohousing","Family","Additive genetics"))
inDFtoPlotL$Var.Exp <- factor(as.character(inDFtoPlotL$Var.Exp),level = c("Cohousing","Family","Environment","Additive genetics"))


# save data frame used to make plot
write.table(inDFtoPlotL,paste0('heritability_analysis_v2/real_data_plots/plot_data_heritability_v2_pwys_all.csv'),sep=',',row.names = F)

# make plot
g <- ggplot(inDFtoPlotL,aes(x=Taxon_Shortname,y=Var.Exp.NR,fill=Var.Exp)) + 
  scale_fill_manual(values = cbPalette) +
  geom_col(col="black", width=1,size=0.75) + 
  theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(-0.01,1.01) +
  theme(axis.text.y = element_text(face="bold")) +
  geom_errorbar(ymin=inDFtoPlotL$CI_ID_low,ymax=inDFtoPlotL$CI_ID_high,width=0.25, linetype='solid') + 
  geom_text(data = inDFtoPlotL,
            aes(x = Taxon_Shortname, y=CI_ID_high, 
                label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
            color="black", vjust=+0.75, angle = 0, hjust=-1,size=6) + ylim(-0.01,1.01) + 
  ylab('Microbiome variance explained') + xlab('') + 
  theme(legend.position="bottom") + theme_classic() +
  theme(text = element_text(size = 14)) + 
  coord_flip()
print(g)
ggsave(plot = g, filename = paste0('heritability_analysis_v2/real_data_plots/heritability_Rev1_PWYS.png'),height = 1.25+8/50*nrow(inDFs),width = 20,limitsize = F)

# smaller plot (top 20)
topN <- 20
topFeatures <- inDFm[order(inDFm$VE_ID,decreasing = T),]$Trait[1:topN]
inDFtoPlotLs <- inDFtoPlotL[inDFtoPlotL$Taxon_Shortname %in% topFeatures,]

write.table(inDFtoPlotLs,paste0('real_data_plots/plot_data_heritability_v2_pwys_top.csv'),sep=',',row.names = F)
g <- ggplot(inDFtoPlotLs,aes(x=Taxon_Shortname,y=Var.Exp.NR,fill=Var.Exp)) + 
  scale_fill_manual(values = cbPalette) +
  geom_col(col="black", width=1,size=0.75) + 
  theme(axis.text.x = element_text(angle = 0,face="bold")) + ylim(-0.01,1.01) +
  theme(axis.text.y = element_text(face="bold")) +
  geom_errorbar(ymin=inDFtoPlotLs$CI_ID_low,ymax=inDFtoPlotLs$CI_ID_high,width=0.25, linetype='solid') + 
  geom_text(data = inDFtoPlotLs,
            aes(x = Taxon_Shortname, y=CI_ID_high,
                label = format(LBL, nsmall = 0, digits=1, scientific = FALSE)),
            color="black", vjust=+0.75, angle = 0, hjust=-1,size=6) + ylim(-0.01,1.01) +
  ylab('Microbiome variance explained') + xlab('') + 
  theme(legend.position="bottom") + theme_classic() +
  theme(text = element_text(size = 26)) + 
  coord_flip() + theme(legend.position = "none")
print(g)
ggsave(plot = g, filename = 'heritability_analysis_v2/real_data_plots/heritability_Rev1_PWYS_top20.png',height = 2.50+16/50*topN,width = 18,limitsize = F)

