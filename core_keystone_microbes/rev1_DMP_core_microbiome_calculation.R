# ========================================================
# By: Weersma Group, UMCG (2020)
#
# codes for identification of core microbiome
# ========================================================
# load libraries
# ========================
library(matrixStats)
library(stringr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(gridExtra)

# =======================
# set path:
# =======================
setwd('D:/Vbox/shared/dag/git_14_05/DMP/')

# load data
taxa=read.table("DMP_data_EGA/DMP_metaphlan2_raw.csv",header = T,sep = ",")
rownames(taxa) <- taxa$ID
taxa$ID <- NULL
taxa <- taxa[,grep('s__',colnames(taxa))]
taxa <- taxa[,-grep('t__',colnames(taxa))]
#taxa <- taxa[,grep('k__Bacteria|k__Archaea',colnames(taxa))]

pwys=read.table("DMP_data_EGA/DMP_humann2_metacyc_raw.csv",header = T,sep = ",")
rownames(pwys) <- pwys$ID
pwys$ID <- NULL

sp <- merge(taxa,pwys,by="row.names")
sp$Row.names <- NULL

## 1) calculate core microbes
# filter by prevalence
sp=sp[rowSums(sp)>0,colSums(sp>0)>nrow(sp)*0.05] #140 species 347 pathways

### 100x sampling based core microbiome selection
result=NULL
sd=NULL
percent=seq(0.01,1,by=0.01)
for (i in percent) {
  print(i)
  tmp_percent=NULL
  for(j in 1:100){
    permutation=sample(row.names(sp),round(nrow(sp)*i,0))
    per_sp=sp[permutation,]
    tmp_percent=rbind(tmp_percent,colSums(per_sp>0)/round(nrow(sp)*i,0))
  }
  result=cbind(result,colMeans(tmp_percent))
  sd=cbind(sd,colSds(tmp_percent))
}
colnames(result)=1:100
colnames(sd)=1:100

res2 <- result
row.names(res2)=str_split_fixed(row.names(result),".s__",2)[,2]
row.names(res2)[row.names(res2)==""] <- row.names(result)[row.names(res2)==""]
res2=cbind(res2,rowMeans(res2))
colnames(res2)[101]="mean_rate"

# plot prevalence VS resampling for each microbiome feature
pdf("core_keystone_microbes/plots/core_microbiome_path_sampling_threshold.pdf",width = 10,height = 10,useDingbats = F)
par(mfrow=c(5,2),mgp=c(1,2,0.2))
for(i in 1:nrow(res2)){
  plot(as.numeric(colnames(res2)[1:100]),res2[i,1:100],type = "b",cex=0.4,frame = FALSE, pch = 19, col = "lightpink", xlab = "Sampling percentage (%)", ylab = "Presence rate",main = row.names(res2)[i],cex.main=0.8,xlim = c(0,100),ylim = c(min(res2[i,1:100])-max(sd[i,1:100]),max(res2[i,1:100])+max(sd[i,1:100])))
  segments(as.numeric(colnames(res2)[1:100])-0.02,res2[i,1:100]-sd[i,1:100],as.numeric(colnames(res2)[1:100])+0.02,res2[i,1:100]+sd[i,1:100],col = "gray60",lty = 1,cex=0.02)
}
dev.off()

data=data.frame(t(res2[,1:100]))
data$id=row.names(data)
data=melt(data,id="id")
data$value=data$value*100
data$id=as.numeric(as.character(data$id))
p1=ggplot(data, aes(id, value, group=variable)) +geom_line(color="black", size=0.05,alpha=0.8)+xlab(label = "Sampling percentage (%)") + ylab(label = "Presence rate (%)")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank())+ theme(legend.position = "none")+scale_x_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100)) 

data=data.frame(res2)
data=data[order(data$mean_rate,decreasing = T),]
data$order=1:nrow(data)
data$mean_rate=data$mean_rate*100
data_tmp=data
data_tmp$name=row.names(data_tmp)
data_tmp$name[which(data_tmp$mean_rate<80)]=NA
p2=ggplot(data_tmp,aes(order,mean_rate))+geom_point(shape=20,alpha=1,cex=0.2)+geom_hline(yintercept = 80,linetype=3,colour="mediumpurple1",alpha=0.5)+labs(x="Rank of core microbiome",y="Presence rate (%)")+theme_bw()+theme(panel.grid=element_blank(),axis.title=element_text(size=8))+theme(legend.position="none")+geom_text_repel(aes(order,mean_rate,label=factor(data_tmp$name)),size=1,alpha=0.85,colour="black",direction = "both",segment.size=0.1,segment.alpha = 0.7,segment.color="red")

res2=data.frame(res2)
res2$mean_rate=res2$mean_rate*100
data=data.frame(percent=1:100,number=NA)
for(i in 1:100){
  data$number[i]=length(which(res2$mean_rate>=i))
}

p3=ggplot(data, aes(percent, number)) +geom_line(color="black", size=0.5,alpha=0.8)+xlab(label = "Presence rate (%)") + ylab(label = "Number of core microbiome")+guides(linetype=FALSE,shape=FALSE)+theme_bw()+theme(panel.grid=element_blank())+ theme(legend.position = "none")+scale_x_continuous(limits = c(0,100),breaks = c(0,20,40,60,80,100)) 
write.table(res2,file = "core_keystone_microbes/tables/res2_dag3_core_microbiome_path.txt",quote = F,sep = "\t")

pdf(file="core_keystone_microbes/plots/core_microbiome_path.pdf",useDingbats=F)
grid.arrange(p1,p3,p2,ncol=1,nrow = 3)
dev.off()