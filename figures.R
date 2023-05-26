#Code for generating figures

#Load dependencies
library(Cairo)
library(data.table)
library(extrafont)
library(ggbreak)
library(ggplot2)
library(ggrepel)
library(ggridges)
library(gridExtra)
library(Gviz)
library(qvalue)
library(reshape2)
library(rtracklayer)
library(Rsamtools)
library(scales)
library(stringr)
library(tidyverse)
library(Unicode)
library(UpSetR)
library(viridis)


################################################
# Figure 1
################################################
#Effect sizes between cumulative early adversity in the low and high quality environments.

#Exclude loci with too high of SE beta
#Set pvalues to NA if se_beta >.6
for (f in 1:3){
  results_model2[which(results_model2[,f+3]>.6),f+6]<-NA
}

#Focus on significant cumulative EA sites in low or high quality: 20% FDR
keep<-which(qvalue(results_model2$cumulative_high_quality_pvalue)$qvalues<.2 |
            qvalue(results_model2$cumulative_low_quality_pvalue)$qvalues<.2 )
data_1b<-results_model2[keep,2:3]
rownames(data_1b)<-rownames(results_model2)[keep]
colnames(data_1b)<-c("High Quality","Low Quality")
data_1b$site<-rownames(data_1b)
data_1b<-melt(data_1b,id.vars = "site")

colors=c("#5F4B8BFF","#E69A8DFF")
p1b<-ggplot(data=data_1b,aes(abs(value),fill=variable))+
  geom_density(alpha=0.7,position="identity",bw=0.02)+theme_minimal()+
  geom_histogram(aes(y=..density..),col="Black", binwidth=0.04,alpha=0.4,position="identity")+
  xlab("|Standardized betas|")+
  scale_fill_manual(values = colors)+
  geom_vline(xintercept = median(abs(data_1b$value)[data_1b$variable==levels(data_1b$variable)[1]]),col=colors[1])+
  geom_vline(xintercept = quantile(abs(data_1b$value)[data_1b$variable==levels(data_1b$variable)[1]],probs=.025),col=colors[1],lty=2)+
  geom_vline(xintercept = quantile(abs(data_1b$value)[data_1b$variable==levels(data_1b$variable)[1]],probs=.975),col=colors[1],lty=2)+
  geom_vline(xintercept = median(abs(data_1b$value)[data_1b$variable==levels(data_1b$variable)[2]]),col=colors[2])+
  geom_vline(xintercept = quantile(abs(data_1b$value)[data_1b$variable==levels(data_1b$variable)[2]],probs=.025),col=colors[2],lty=2)+
  geom_vline(xintercept = quantile(abs(data_1b$value)[data_1b$variable==levels(data_1b$variable)[2]],probs=.975),col=colors[2],lty=2)+
  xlim(c(0,2))+
  theme(text=element_text(color="black"),axis.text=element_text(size=16,color="black"),axis.title = element_text(size=20))+ylab("Density")+
  guides(fill="none")


p1c<-ggplot(data=data_1b)+
  geom_point(aes(variable,value,group=site),alpha=0.2,size=2)+
  geom_line(aes(variable,value,group=site,),lwd=0.05,alpha=0.1)+
  theme_minimal()+
  scale_color_manual(values = colors)+
  ylab("Standardized beta (cumulative early adversity effect)")+xlab("")+
  geom_hline(yintercept = 0,lty=2)+
  scale_alpha_continuous(range = c(.0125,1))+guides(alpha="none",color="none")+
  theme(text=element_text(color="black"),axis.text=element_text(size=16,color="black"),axis.title = element_text(size=20))+
  scale_x_discrete(labels=c("High-quality habitat","Low-quality habitat"))


grid.arrange(p1b,p1c,layout_matrix=matrix(c(1,2),nrow=1))

rm(list = grep(ls(),invert=T,pattern="results",value = T))




################################################
# Figure 2
################################################

#Here we want to visualize the results of rank from model1 age, habitat quality, and nested early adversity from model 3.
results<-cbind.data.frame(results_model1[13:14],results_model3[,25:36])
results2<-cbind.data.frame(results_model1[8:9],results_model3[,13:24])

sig<-as.data.frame(colnames(results))
colnames(sig)<-"variable"
sig$num_sig<-NA

colnames(results2)<-colnames(results)
for(f in sig$variable){
  sig$num_sig[sig$variable==f]<-length(which(
    qvalue(results[which(results2[,colnames(results2)==f]<.6),colnames(results)==f])$qvalues<.1
    
  ))
  
}

sig$type<-c("adult","adult",rep("early",11),"adult")
sig$variable[c(3,5,7,9,11,13)]
sig$orders<-c(2,1,13,9,10,11,12,7,8,5,6,3,4,1)
sig$spacing<-c(1.75,1.25,13,9.75,10.25,11.75,12.25,7.75,8.25,5.75,6.25,3.75,4.25,0.75)
sig$variable<-c("Male rank","Female rank","Habitat quality",
                "Maternal loss (high-quality habitat)","Maternal loss (low-quality habitat)",
                "Drought (high-quality habitat)","Drought (low-quality habitat)",
                "Group size (high-quality habitat)","Group size (low-quality habitat)",
                "Maternal rank (high-quality habitat)","Maternal rank (low-quality habitat)",
                "Close-in-age younger sibling (high-quality habitat)","Close-in-age younger sibling (low-quality habitat)",
                "Age")

sig$type[c(3,5,7,9,11,13)]<-c("early_low")

sig$spacing<-c(1.75,1.25,10.25,7.25,7.75,8.75,9.25,5.75,6.25,4.25,4.75,2.75,3.25,0.75)
p2a <-sig[,] %>% mutate(variable = fct_reorder(variable, orders)) %>%
  ggplot()+
  ylab("")+xlab("")+
  geom_hline(yintercept = 0)+
  geom_segment(aes(x=spacing,xend=spacing,y=0,yend=num_sig),lty=2)+
  geom_point(aes(spacing,num_sig,fill=as.factor(type)),pch=21,size=5)+
  coord_flip()+theme_bw()+
  scale_fill_manual(values=c("Maroon","#736795","#DCA69C"),
                    labels=c("Adult","Early life (high-quality)","Early life (low-quality)"))+
  scale_x_continuous(breaks=sig$spacing,labels=sig$variable)+
  xlab("")+ylab("")+
  theme(axis.text=element_text(size=16,color="black"),axis.text.x = element_text(size=12,angle = 45,hjust=1),axis.title.y=element_blank(),
        axis.title.x=element_blank(),panel.grid = element_blank(),text=element_text(size=20,color="black"),
        axis.text.x.top = element_blank(),
        legend.position = "top",legend.title = element_blank(),legend.text = element_text(size=18))+
  scale_y_log10(breaks=c(10,1000,100000),
                   labels=c(10,1000,100000))+xlab("Number of significant CpG sites (FDR <10%)")+guides(fill=FALSE)




#2B Example CpG site that shows an interaction for drought (with raw methylation ratios plotted)
a1<-fread("./Data/mratios/all_mratio.txt")
b1<-fread("./Data/info_all.txt")
b1<-b1[grep(b1$V1,pattern="^chr"),]
rownames(a1)<-b1$V1
early_mat_n295<-read.table("./Data/early_mat_n295.txt",header=T)

###############
#81025
top<-rownames(results_model3)[as.numeric(81025)]
par(mfrow=c(1,1),pty="s")

unlist(a1[rownames(a1)==top,])

f1=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==0],na.rm=T)
s1=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==0],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==0])))

f2=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==1],na.rm=T)
s2=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==1],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==1])))

f3=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==0],na.rm=T)
s3=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==0],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==0])))

f4=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==1],na.rm=T)
s4=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==1],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==1])))

colors<-c(c("#51427B","#D48F84"),c("#51427B","#D48F84"))


p2b<-ggplot()+
  geom_segment(aes(x=0,xend=1,y=f1,yend=f3,col="High quality"))+
  geom_segment(aes(x=0,xend=1,y=f2,yend=f4,col="Low quality"))+
  geom_point(aes(0,f1,fill="High quality"),size=7,pch=21)+
  geom_errorbar(aes(x=0,ymin=f1-s1,ymax=f1+s1,col="High quality"),width=0.1)+
  geom_point(aes(0,f2,fill="Low quality"),size=7,pch=21)+
  geom_errorbar(aes(x=0,ymin=f2-s2,ymax=f2+s2,col="Low quality"),width=0.1)+
  geom_point(aes(1,f3,fill="High quality"),size=7,pch=21)+
  geom_errorbar(aes(x=1,ymin=f3-s3,ymax=f3+s3,col="High quality"),width=0.1)+
  geom_point(aes(1,f4,fill="Low quality"),size=7,pch=21)+scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  geom_errorbar(aes(x=1,ymin=f4-s4,ymax=f4+s4,col="Low quality"),width=0.1)+
  scale_x_continuous(breaks = c(0,1),labels = c("No","Yes"))+ylab("Methylation level")+theme_bw()+
  theme(legend.title = element_blank())+xlab("Experienced drought?")+
  theme(text=element_text(size=16,color="black"),axis.text = element_text(size=16,color="black"),axis.title = element_text(size=20),
        legend.text = element_text(size=20))+guides(fill="none",colour="none")



top<-rownames(results_model3)[as.numeric(69289)]
par(mfrow=c(1,1),pty="s")

unlist(a1[rownames(a1)==top,])


f12=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==0],na.rm=T)
s12=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==0],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==0])))

f22=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==1],na.rm=T)
s22=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==1],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==0 & early_mat_n295[,15]==1])))

f32=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==0],na.rm=T)
s32=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==0],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==0])))

f42=mean(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==1],na.rm=T)
s42=sd(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==1],na.rm=T)/sqrt(length(na.exclude(unlist(a1[rownames(a1)==top,])[early_mat_n295[,11]==1 & early_mat_n295[,15]==1])))


p2c<-ggplot()+
  geom_segment(aes(x=0,xend=1,y=f12,yend=f32,col="High quality"))+
  geom_segment(aes(x=0,xend=1,y=f22,yend=f42,col="Low quality"))+
  geom_point(aes(0,f12,fill="High quality"),size=7,pch=21)+
  geom_errorbar(aes(x=0,ymin=f12-s12,ymax=f12+s12,col="High quality"),width=0.1)+
  geom_point(aes(0,f22,fill="Low quality"),size=7,pch=21)+
  geom_errorbar(aes(x=0,ymin=f22-s22,ymax=f22+s22,col="Low quality"),width=0.1)+
  geom_point(aes(1,f32,fill="High quality"),size=7,pch=21)+
  geom_errorbar(aes(x=1,ymin=f32-s32,ymax=f32+s32,col="High quality"),width=0.1)+
  geom_point(aes(1,f42,fill="Low quality"),size=7,pch=21)+
  scale_fill_manual(values=colors)+scale_color_manual(values=colors)+
  geom_errorbar(aes(x=1,ymin=f42-s42,ymax=f42+s42,col="Low quality"),width=0.1)+
  scale_x_continuous(breaks = c(0,1),labels = c("No","Yes"))+ylab("Methylation level")+theme_bw()+
  theme(legend.title = element_blank())+xlab("Experienced drought?")+
  theme(text = element_text(size=16,color="black"),axis.title = element_text(size=20),
        legend.text = element_text(size=20),axis.text = element_text(color="black"))+guides(fill="none",colour="none")



#Plot 2C 
b2<-results_model3[,2:11]
se_beta<-results_model3[,14:23]

for(f in 1:10){
  b2[which(se_beta[,f]>.6),f]<-NA
}

b_melt<-melt(b2)
b_melt$variable<-as.character(b_melt$variable)
b_melt$variable<-gsub(b_melt$variable,pattern="_bhat",replacement="")
b_melt$qual<-"Low-quality"
b_melt$qual[grep(b_melt$variable,pattern="high_quality")]<-"High-quality"
b_melt$Var2<-NA
b_melt$Var2<-gsub(gsub(gsub(b_melt$variable,pattern="_high_quality",replacement = ""),pattern="_low_quality",replacement = ""),pattern="_",replacement="")

b_melt$Var2[b_melt$Var2=="biggrp"]<-"Group size"
b_melt$Var2[b_melt$Var2=="competingsib"]<-"Close-in-age younger sibling"
b_melt$Var2[b_melt$Var2=="drought"]<-"Drought"
b_melt$Var2[b_melt$Var2=="lowrankmom"]<-"Maternal rank"
b_melt$Var2[b_melt$Var2=="matloss"]<-"Maternal loss"


p2d<-ggplot(data=b_melt,aes(abs(value),fill=qual))+
  geom_density(alpha=0.6)+
  geom_histogram(aes(y=..density..),col="Black", binwidth=0.04,alpha=0.4,position="identity")+
  facet_wrap(~Var2,ncol=1)+
  theme_minimal()+
  xlab("|Standardized beta|")+
  scale_fill_manual(values = colors)+
  theme(legend.title = element_blank())+ylab("Density")+
  theme(text = element_text(size=24),axis.text = element_text(size=16,color="black"),axis.title = element_text(size=20),legend.position = c(1,4),
        legend.text = element_text(size=20))+xlim(c(0,2))


grid.arrange(p2a,p2b,p2c,p2d,layout_matrix=rbind(c(1,1,4),c(1,1,4),c(2,3,4)))


#2E
significance<-matrix(NA,nrow=length(rownames(results_model3)),ncol=12)
rownames(significance)<-rownames(results_model3)
for(f in 1:12){
  significance[,f]<-as.numeric(qvalue(results_model3[,f+24])$qvalues<.1)
  significance[which(results_model3[,f+12]>.6),f]<-0
}

#Lower FDR threshold to 20% if one of the the variables of interest is significant
for(f in c(1,3,5,7,9,11)){
  sig<-which(apply(significance[,c(1,3,5,7,9,11)],1,function(x){length(which(x==1))})>0)
  
  significance[sig,f]<-as.numeric(qvalue(results_model3[,f+24])$qvalues<.2)[sig]
  significance[which(results_model3[,f+12]>.6),f]<-0
}

significance<-as.data.frame(significance)
colnames(significance)<-colnames(pvalue)
colnames(significance)

colnames(significance)[c(1,3,5,7,9,11)]<-c("Habitat Quality","Maternal Loss","Drought","Group Size","Low Maternal Rank","Close-in-age Younger Sibling")


upset(significance[,c(1,3,5,7,9,11)],
          nsets = 6,order.by = "freq",nintersects = 20,
          main.bar.color = hcl.colors(30,"Dark Mint")[1:20],
          mainbar.y.label = "Number of CpG sites",text.scale = 1.5,sets.x.label = "")



rm(list = grep(ls(),invert=T,pattern="results",value = T))



################################################################################
#3A heatmap of top effects in gene promoter
################################################################################

#Focus on (Male rank, age, habitat quality, drought, and a batch effect (e.g. new batch))
overlap<-read.table("./results/sites_in_each_basic_compartment.bed")
overlap$V7<-as.character(overlap$V7)
overlap$V7[overlap$V7=="."]<-"Unannotated"
table(overlap$V7)

OR<-matrix(NA,nrow=6,ncol=4)
p<-matrix(NA,nrow=6,ncol=4)

results_temp<-cbind.data.frame(results_model1[,13],results_model3[,c(25,29,36)])
results2_temp<-cbind.data.frame(results_model1[,8],results_model3[,c(13,17,24)])

colnames(OR)<-colnames(p)<-c("Male rank","Habitat quality","Drought (low-quality habitat)","Age")
rownames(OR)<-rownames(p)<-unique(overlap$V7)

overlap$site<-paste(overlap$V1,overlap$V2,sep="_")
u<-unique(overlap$V7)

for(x in 1:4){
  sigs<-cbind.data.frame(rownames(results_model1)[results2_temp[,x]<.6],as.numeric(qvalue(results_temp[results2_temp[,x]<.6,x])$qvalues<.1))
  names(sigs)<-c("site","sig")
  overlap2<-merge(overlap,sigs,by="site")
  
  #Ratio of sig to non-sig in each bin
  
  for(f in 1:6){
    OR[f,x]<-fisher.test(matrix(c(length(which(sigs$sig==0 & (!sigs$site %in% overlap2$site[overlap2$V7==u[f]]))),
                                  length(which(sigs$sig==0 & (sigs$site %in% overlap2$site[overlap2$V7==u[f]]))),
                                  length(which(sigs$sig==1 & (!sigs$site %in% overlap2$site[overlap2$V7==u[f]]))),
                                  length(which(sigs$sig==1 & (sigs$site %in% overlap2$site[overlap2$V7==u[f]])))),
                                nrow=2))$estimate
    
    p[f,x]<-fisher.test(matrix(c(length(which(sigs$sig==0 & (!sigs$site %in% overlap2$site[overlap2$V7==u[f]]))),
                                 length(which(sigs$sig==0 & (sigs$site %in% overlap2$site[overlap2$V7==u[f]]))),
                                 length(which(sigs$sig==1 & (!sigs$site %in% overlap2$site[overlap2$V7==u[f]]))),
                                 length(which(sigs$sig==1 & (sigs$site %in% overlap2$site[overlap2$V7==u[f]])))),
                               nrow=2))$p.value
  }
}

OR_melt<-melt(OR)

#Number of significant overlapping sites
p_melt<-melt(p)
OR_melt3<-merge(OR_melt,p_melt,by=c("Var1","Var2"))

OR_melt$regulatory<-NA
OR_melt$regulatory[OR_melt$Var1=="enhancer"]<-1
OR_melt$regulatory[OR_melt$Var1=="gene"]<-2
OR_melt$regulatory[OR_melt$Var1=="promoter"]<-3
OR_melt$regulatory[OR_melt$Var1=="cpg_island"]<-4
OR_melt$regulatory[OR_melt$Var1=="cpg_shores"]<-5
OR_melt$regulatory[OR_melt$Var1=="Unannotated"]<-6

OR_melt$Var1<-as.character(OR_melt$Var1)
OR_melt$Var1[OR_melt$Var1=="enhancer"]<-"Enhancer"
OR_melt$Var1[OR_melt$Var1=="gene"]<-"Gene body"
OR_melt$Var1[OR_melt$Var1=="promoter"]<-"Promoter"
OR_melt$Var1[OR_melt$Var1=="cpg_island"]<-"CpG island"
OR_melt$Var1[OR_melt$Var1=="cpg_shores"]<-"CpG shore"
OR_melt$Var1[OR_melt$Var1=="Unannotated"]<-"Unannotated"

OR_melt$order<-NA
OR_melt$order[OR_melt$Var2=="Male rank"]<-1
OR_melt$order[OR_melt$Var2=="Drought (low-quality habitat)"]<-2
OR_melt$order[OR_melt$Var2=="Habitat quality"]<-3
OR_melt$order[OR_melt$Var2=="Age"]<-4.5


OR_melt %>% mutate(Var1=fct_reorder(Var1,regulatory)) %>%
  mutate(Var2=fct_reorder(Var2,order)) %>%
  ggplot()+geom_tile(aes(x = Var1,order,fill=log2(value)),col="Black",width=.95,height=.95)+
  scale_y_continuous(breaks=c(1,2,3,4.5),labels=unique(as.character(OR_melt$Var2[order(OR_melt$order)])))+
  theme_minimal()+
  xlab("Genomic compartment")+ylab("Predictor variable")+
  scale_fill_gradientn(colors=c("#0A1569","Steel Blue","White","Maroon"),breaks=c(-1,0,1),labels=c(-1,0,1),limits=c(-1.75,1))+
  theme(text=element_text(size=20),axis.text.y=element_text(color="Black"),axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,color="Black"),legend.position = "top")+
  xlab("")+ylab("")+labs(fill=expression(paste("Log"[2],"(Odds ratio)",sep="")))




################################################
### Fig 3B:
################################################
overlap<-read.table("./results/all_sites_intersected_with_chromhmm.bed")
table(overlap$V8)

names<-read.delim("./Data/chrom_hmm/chrom_hmm_info.txt")
mat<-as.data.frame(matrix(NA,nrow=15,ncol=2))
colnames(mat)<-c("name","number_segments_overlap")
mat$name<-names$MNEMONIC
overlap$V8<-as.character(overlap$V8)
for(f in 1:15){
  mat[f,2]<-length(which(overlap$V8==paste("E",f,sep = "")))
  
}
mat$lname<-names$DESCRIPTION

OR<-matrix(NA,nrow=15,ncol=4)
p<-matrix(NA,nrow=15,ncol=4)

colnames(OR)<-colnames(p)<-c("Male rank","Habitat quality","Drought (low-habitat quality)","Age")
rownames(OR)<-rownames(p)<-c(as.character(mat$name))

overlap$site<-paste(overlap$V1,overlap$V2,sep="_")

for(x in c(1:4)){
  sigs<-cbind.data.frame(rownames(results_model3)[results2_temp[,x]<.6],as.numeric(qvalue(results_temp[results2_temp[,x]<.6,x])$qvalues<.1))
  names(sigs)<-c("site","sig")
  overlap2<-merge(overlap,sigs,by="site")
  
  #Ratio of sig to non-sig in each bin
  
  for(f in 1:15){
    
    OR[f,x]<-fisher.test(matrix(c(length(which(sigs$sig==0 & (!sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")]))),
                                  length(which(sigs$sig==0 & (sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")]))),
                                  length(which(sigs$sig==1 & (!sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")]))),
                                  length(which(sigs$sig==1 & (sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")])))),
                                nrow=2))$estimate
    
    p[f,x]<-fisher.test(matrix(c(length(which(sigs$sig==0 & (!sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")]))),
                                 length(which(sigs$sig==0 & (sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")]))),
                                 length(which(sigs$sig==1 & (!sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")]))),
                                 length(which(sigs$sig==1 & (sigs$site %in% overlap2$site[overlap2$V8== paste("E",f,sep="")])))),
                               nrow=2))$p.value
  }
}


#Number of significant overlapping sites
colnames(OR)<-c("Male rank","Habitat quality","Drought (low-habitat quality)","Age")
colnames(p)<-c("Male rank","Habitat quality","Drought (low-habitat quality)","Age")
OR_melt<-melt(OR)

p_melt<-melt(p)
OR_melt3<-merge(OR_melt,p_melt,by=c("Var1","Var2"))

mat
mat$reg_activity<-NA
mat$reg_activity<-c(1,2,3,4,5,7,6,8,15,9,10,11,14,13,12)
OR_melt2<-merge(OR_melt,mat,by.x="Var1",by.y="name")


#With shading based on significance?
OR_melt3<-merge(OR_melt2,p_melt,by=c("Var1","Var2"))
OR_melt3$alpha<-as.numeric(OR_melt3$value.y<.01)


a<-OR_melt3[OR_melt3$Var2!="Age",] %>% mutate(lname=fct_reorder(lname,reg_activity)) %>%
  ggplot()+
  geom_line(aes(lname, log2(value.x),group=Var2),alpha=0.4,colour="Black")+
  geom_point(aes(lname, log2(value.x),fill=Var2,alpha=alpha),pch=21,size=5)+
  theme_classic()+
  xlab("")+
  ylab(bquote(Log[2]~'(Odds Ratio)'))+
  scale_fill_manual(values = viridis(20)[c(1,8,14,18)])+ 
  scale_colour_manual(values = viridis(20)[c(1,8,14,18)])+ 
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept=7.5,lty=2)+
  theme(axis.text.x = element_blank(),legend.title = element_blank(),legend.text = element_text(size=22,color="Black"),
        legend.background = element_rect(fill = NA),legend.position = "top",
        axis.text.y = element_text(size=16,color="Black"),plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.title.y=element_text(size=20,color="Black"))+
  scale_y_continuous(breaks = c(-2,-1,0,1,2),limits=c(-3,3))+guides(alpha=FALSE)

b<-OR_melt3[OR_melt3$Var2=="Age",] %>% mutate(lname=fct_reorder(lname,reg_activity)) %>%
  ggplot()+
  geom_line(aes(lname, log2(value.x),group=Var2),alpha=0.4,colour="Black")+
  geom_point(aes(lname, log2(value.x),fill=Var2,alpha=alpha),pch=21,size=5)+
  theme_classic()+
  xlab("")+
  ylab(bquote(Log[2]~'(Odds Ratio)'))+
  scale_fill_manual(values = viridis(5)[5])+ 
  scale_colour_manual(values = viridis(5)[5])+ 
  geom_hline(yintercept = 0,lty=2)+
  geom_vline(xintercept=7.5,lty=2)+
  theme(legend.title = element_blank(),legend.text = element_text(size=22),
        legend.background = element_rect(fill = NA),legend.position = c(.3,.85),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1,size=16,color="Black"),axis.text.y = element_text(size=16,color="Black"),plot.margin = unit(c(.1,.1,.1,.1), "cm"),
        axis.title.y=element_text(size=20,color="Black"))+
  guides(alpha=FALSE)+
  scale_y_continuous(breaks = c(-2,-1,0,1,2),limits = c(-2,2.25))


grid.arrange(a,b,nrow=2,layout_matrix=matrix(c(1,1,2,2,2),ncol=1))


rm(list = grep(ls(),invert=T,pattern="results",value = T))
rm(results_temp,results,results2_temp,results2)


####################################################################################
#Fig 4
####################################################################################


#predicted_binomial<-parallel_predict_no_norm(age_vector = info$habitat_quality,counts_file = all_imputed,no_cores = 7,nf=50,alphas=1)
predicted_binomial<-read.table("./results/habitat_quality_prediction_alpha1_50fold_binomial.txt")

info<-read.table("./Data/meta_info_n295.txt",header=T)
info$predicted_binomial_hab<-predicted_binomial$V1

info$habitat_qual<-info$habitat_quality
info$habitat_qual[info$habitat_qual==0]<-"post_shift"
info$habitat_qual[which(info$habitat_qual==1 & info$matgrp.x==1)]<-"Altos"
info$habitat_qual[which(info$habitat_qual==1 & info$matgrp.x==2)]<-"Hooks"

colors=c("#5F4B8BFF","#E69A8DFF")
p4a<-ggplot(data=info,aes(habitat_quality,predicted_binomial_hab,group=as.factor(habitat_quality),fill=as.factor(habitat_quality)))+
  geom_violin(alpha=0.5)+geom_jitter(width = .2,pch=21)+theme_bw()+
  scale_fill_manual(values = colors)+
  xlab("Habitat quality")+ylab("Predicted habitat quality")+
  guides(fill="none",alpha="none")+
  scale_x_continuous(labels= c("Low","High"),breaks = c(1,0))+
  scale_alpha_manual(values=c(.5,.9,1))+
  theme(text=element_text(size=20,color="Black"),axis.text.x=element_text(size=20,color="Black"),axis.text.y=element_text(size=20,color="Black"))


##############
# 4B
##############
# ROC curve
TPR<-1:295
FPR<-1:295

for(f in 1:295){
  thresh<-sort(info$predicted_binomial_hab,decreasing = T)[f]
  TPR[f]<-length(which(info$predicted_binomial_hab>=thresh & info$habitat_quality==1))/length(which(info$habitat_quality==1))
  FPR[f]<-length(which(info$predicted_binomial_hab>=thresh & info$habitat_quality==0))/length(which(info$habitat_quality==0))
  
}


#Area under curve?
ROC<-function(x){
  return(sapply(x,function(Z){return(TPR[min(which((FPR-Z)>=0))])}))
}

integrate(ROC,lower=0,upper=1,subdivisions = 10000)

p4b<-ggplot()+geom_line(aes(FPR,TPR),lwd=0.5)+geom_point(aes(FPR,TPR),pch=21,fill=colors[2],alpha=0.8,lwd=2)+theme_bw()+
  xlab("False positive rate (1- specificity)")+ylab("True positive rate (sensitivity)")+
  geom_abline(intercept = 0,slope=1,lty=2)+theme(text=element_text(size=20,color="Black"),axis.text=element_text(color="Black"))




##############
# 4C
##############

#Duration of time in low habitat quality?
#Time since habitat quality?
info$year_birth<-substr(info$birth.x,1,4)

max(as.numeric(info$year_birth[info$matgrp.x==1 & info$habitat_quality==1]))
max(as.numeric(info$year_birth[info$matgrp.x==2 & info$habitat_quality==1]))

info$time_since<-NA
info$time_since[info$habitat_quality==1 & info$matgrp.x==1]<-as.numeric(as.Date(info$dart_date[info$habitat_quality==1 & info$matgrp.x==1])-as.Date("1988-01-01"))
info$time_since[info$habitat_quality==1 & info$matgrp.x==2]<-as.numeric(as.Date(info$dart_date[info$habitat_quality==1 & info$matgrp.x==2])-as.Date("1992-01-01"))

info$time_since[which(info$time_since<0)]<-0

summary(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$time_since[info$habitat_quality==1]))


p4c<-ggplot()+  geom_smooth(aes(info$time_since[info$habitat_quality==1],info$predicted_binomial_hab[info$habitat_quality==1]),se=FALSE,method = "lm",lty=2,col="Black",lwd=0.5)+
  geom_point(aes(info$time_since[info$habitat_quality==1],info$predicted_binomial_hab[info$habitat_quality==1]),pch=21,fill=colors[2],alpha=0.8,lwd=2)+
  theme_bw()+
  xlab("Time since habitat shift (days)")+ylab("Predicted habitat quality")+
  theme(text=element_text(size=20),axis.text=element_text(color="Black"))


grid.arrange(p4a,p4b,p4c,ncol=3)


rm(list = grep(ls(),invert=T,pattern="results",value = T))





############################
#Figure 5A
############################

#Visualizing mSTARR pileups near a drought-associated CpG site.

#Annotations we need to include:
#Drought associated CpG sites
#Regulatory Windows
#Annotated genes
#ChromHMM annotations

#So what is a drought associated CpG site that falls within a regulatory region and is within a gene body?


library(tidyverse)
library(valr)
library(scales)
library(MetBrewer)


##############################
# Pileup
##############################
#Okay now we need to focus on a region within a gene body, 
#with MD regulatory activity, and with drought associated CpG sites within it

#Annotations we need include:
#Drought associated CpG sites
#Regulatory Windows
#Annotated genes
#ChromHMM annotations

#So what is a drought associated CpG site that falls within a regulatory region and is within a gene body?
#Preferably one that is methylation dependent?

#Alright lets focus on
#chr17: 35210300-35210400

#Are these drought sites high SE_beta?
a<-results_model3$drought_low_quality_se_beta
p<-results_model3$drought_low_quality_pvalue
ids<-rownames(results_model3)
df<-cbind.data.frame(ids,a,p)
df$qvalue<-NA
df$qvalue[which(df$a<.6)]<-qvalue(df$p[which(df$a<.6)])$qvalues
df<-df[grep(df$ids,pattern="chr17"),]

colnames(df)<-c("site","se_beta","pvalue","qvalue")
df$chr<-gsub(df$site,pattern="_.*",replacement = "")
df$pos<-gsub(df$site,pattern=".*_",replacement = "")

df<-df[as.numeric(df$pos)>35150000,]
df<-df[as.numeric(df$pos)<35250000,]
df2<-df[which(df$qvalue<.2),]



genes<-import.bed("./annotations_for_pileup_figure/genes.bed")
chmm<-import.bed("./annotations_for_pileup_figure/chmm_chr17.bed")
ms<-import.bed("./annotations_for_pileup_figure/ms.bed")
dt<-import.bed("./annotations_for_pileup_figure/drought_associated_sites.bed")


indexBam("./annotations_for_pileup_figure/rna_meth_chr17.bam")
indexBam("./annotations_for_pileup_figure/rna_sham_chr17.bam")


dTrack <- AlignmentsTrack(range="./annotations_for_pileup_figure/rna_sham_chr17.bam",
                          name = "Reads",window=-1,chromosome = "chr17",type="coverage",fill=alpha("Steel Blue",alpha = 0.3))

dTrack2 <- AlignmentsTrack(range="./annotations_for_pileup_figure/rna_meth_chr17.bam",
                           name = "Reads",window=-1,chromosome = "chr17",type="coverage",fill=alpha("Dark Green",alpha = 0.3))

aTrack <- AnnotationTrack(range = ms,
                          name = "mSTARR regulatory window", chromosome = "chr17", 
                          fill=alpha("Maroon",alpha = 0.5))

aTrack2 <- AnnotationTrack(range = chmm,
                           name = "ChromHMM", chromosome = "chr17",
                           fill=alpha("Maroon",alpha = 0.5),stacking = "full")

aTrack3 <- AnnotationTrack(range = genes,
                           name = "Gene bodies", chromosome = "chr17",
                           fill=alpha("Maroon",alpha = 0.5),shape="ellipse")

aTrack4 <- AnnotationTrack(range = dt,
                           name = "Drought-associated sites", chromosome = "chr17",
                           fill=alpha("Maroon",alpha = 0.5),stacking = "dense",shape="ellipse")

#Overlaying methylation and sham reads

ot <- OverlayTrack(trackList = list(dTrack2, dTrack))


#What are the chromHMM annotations at this location?
tt<-read.table("./annotations_for_pileup_figure/chrom_hmm_peaks_sorted.bed")
tt[tt$V1=="chr17"& 35210000>tt$V2 & 35210000 <tt$V3, ]

ttt<-read.delim("./Data/chrom_hmm/chrom_hmm_info.txt")

#And the gene name?
gg<-read.table("./annotations_for_pileup_figure/panubis1_genes_sorted.bed")
gg[which(gg$V1=="chr17"& gg$V2<35210300 & gg$V3> 35210400), ]

gg[which(gg$V1=="chr17"& gg$V2<35210300 & gg$V3> 35210400)+1, ]

#Adding DNA too
indexBam("./annotations_for_pileup_figure/dna_meth_chr17.bam")
indexBam("./annotations_for_pileup_figure/dna_sham_chr17.bam")


dTrack <- AlignmentsTrack(range="./annotations_for_pileup_figure/rna_sham_chr17.bam",
                          name = "Reads",window=-1,chromosome = "chr17",type="coverage",fill=alpha("Orange",alpha = .6),ylim=c(0,1200))

dTrack2 <- AlignmentsTrack(range="./annotations_for_pileup_figure/rna_meth_chr17.bam",
                           name = "Reads",window=-1,chromosome = "chr17",type="coverage",fill=alpha("Steel Blue",alpha = .6),ylim=c(0,1200))

dTrack3 <- AlignmentsTrack(range="./annotations_for_pileup_figure/dna_sham_chr17.bam",
                           name = "Reads",window=-1,chromosome = "chr17",type="coverage",fill=alpha("Orange",alpha = .6),ylim=c(0,350))

dTrack4 <- AlignmentsTrack(range="./annotations_for_pileup_figure/dna_meth_chr17.bam",
                           name = "Reads",window=-1,chromosome = "chr17",type="coverage",fill=alpha("Steel Blue",alpha = .6),ylim=c(0,350))


aTrack <- AnnotationTrack(range = ms,
                          name = "mSTARR regulatory window", chromosome = "chr17", 
                          fill=alpha("#2B6668",alpha = .8))

aTrack2 <- AnnotationTrack(range = chmm,
                           name = "ChromHMM", chromosome = "chr17",
                           fill=alpha("#2B6668",alpha = .8),stacking = "full")

aTrack3 <- AnnotationTrack(range = genes,
                           name = "Gene bodies", chromosome = "chr17",
                           fill=alpha("#2B6668",alpha = .8),shape="ellipse")

aTrack4 <- AnnotationTrack(range = dt,
                           name = "Drought-associated sites", chromosome = "chr17",
                           fill=alpha("#2B6668",alpha = .8),stacking = "dense",shape="ellipse")

#Overlaying methylation and sham reads
ot <- OverlayTrack(trackList = list(dTrack2, dTrack),name = "RNA counts")
ot2 <- OverlayTrack(trackList = list(dTrack3, dTrack4),name = "DNA counts")


pdf("~/Desktop/test.pdf")
plotTracks(list(ot,ot2,aTrack,aTrack2,aTrack3,aTrack4), from = 35150000,
           to = 35250000, background.title = "White",
           col.axis="#000000",col="Black",col.title="Black",collapse=TRUE,ylim=c(0,1200),lwd=0.05)

dev.off()

#Which genes are in this region?
ggg<-read.table("./annotations_for_pileup_figure/panubis1_genes_sorted.bed")
ggg[which(ggg$V1=="chr17" &ggg$V2 <35150000 & ggg$V3 >35250000),]


#Figure 5b and 5c 
m<-read.table("./results/misc_results/mstarr_results_to_plot",header=T)
m$Compartment<-gsub(m$Compartment,pattern="_",replacement = " ",fixed=T)
m$Compartment[5] <-"Repetitive CNV 2"
m$Compartment[7]<-"Heterchromatin / Low" 
m$Compartment[8] <- "Unannotated"
m$Compartment[11] <-"Repetitive CNV 1"
m$Compartment[14] <-"Strong Enhancer 2"
m$Compartment[15] <-"Strong Enhancer 1"


p5b <-m %>% mutate(Compartment = fct_reorder(Compartment, OR)) %>%
  ggplot()+ 
  geom_errorbar(aes(Compartment,ymin=log2(conf1),ymax=log2(conf2)),width=.2,lwd=0.3)+
  geom_point(aes(Compartment,log2(OR)),fill=met.brewer(name="Morgenstern",15)[13],pch=21,col="Black",size=5) +
  theme_bw()+guides(fill=FALSE)+
  theme(axis.text.x = element_text(angle=75,hjust = 1),text=element_text(size=20),plot.margin =margin(t = 1,r=20,1,1))+xlab("ChromHMM compartment")+
  ylab((bquote("Log"[2]~"(Odds ratio)")))+geom_hline(yintercept = 0,lty=2)+ylim(c(-4,4))




#Plot 5C enrichment of age/drought/rank effects in mSTARR regions
##########################
#mSTARR results
##########################
ms_full<-read.table("./results/misc_results/baboon_mstarr_summary_table_31Dec21.txt",header=T)
head(ms_full)
ms_full$reg_activity<-as.character(ms_full$reg_activity)
ms_full$meth_depend<-as.character(ms_full$meth_depend)
ms_full$reg_activity[is.na(ms_full$reg_activity)]<-"no_reg"
ms_full$meth_depend[is.na(ms_full$meth_depend)]<-"no_md"

table(ms_full$reg_activity)

tmp<-ms_full[ms_full$chr=="chrY",]
table(tmp$reg_activity)
table(tmp$reg_activity,tmp$meth_depend)


ms<-read.table("./results/misc_results/mstarr_overlap_with_tested_sites.bed",header=F)
ms$site<-paste(ms$V1,ms$V2,sep="_")
ms2<-ms[,c(7,22:28)][,c(1,5,6,7,8)]
colnames(ms2)<-c("window","reg_activity","meth_depend","source","site")
ms2$reg_activity[is.na(ms2$reg_activity)]<-"none"
ms2$meth_depend<-as.character(ms2$meth_depend)
ms2$meth_depend[is.na(ms2$meth_depend)]<-"none"

w1<-unique(ms2$window)
w3<-w2<-rep(NA,length(w1))

for(f in w1){
  w2[w1==f]<-length(which(ms2$source[ms2$window==f]=="msp1"))
  w3[w1==f]<-length(which(ms2$source[ms2$window==f]=="sheared"))
}
w4<-as.numeric(w2>0 & w3>0)

w5<-w1[w4==1]
ms3<-ms2[-which(ms2$window %in% w5),]
ms<-ms3;rm(ms2,ms3,w1,w2,w3,w4,w5)

results_temp<-results
rownames(results_temp)<-out$id

ms$duplicated<-duplicated(ms$site)
ms<-ms[!ms$duplicated,]
table(ms$reg_activity)
table(ms$meth_depend)

OR<-matrix(NA,nrow=14,ncol=3)
rownames(OR)<-colnames(results)


#Range of thresholds for rank + drought
thresh<-seq(.1,.3,by=.01)
or3<-or2<-or1<-rep(NA,length(thresh))
p3<-p2<-p1<-rep(NA,length(thresh))

for(f in thresh){
  x=1
  sigs<-cbind.data.frame(out$id[results2[,x]<.5],as.numeric(qvalue(results[results2[,x]<.5,x])$qvalues<f))
  names(sigs)<-c("site","sig")
  ms_temp<-ms
  ms2<-merge(ms_temp,sigs,by="site")
  ms2$reg<-NA
  ms2$reg[ms2$reg_activity=="none"]<-0
  ms2$reg[ms2$reg_activity!="none"]<-1
  ms2$sig<-factor(ms2$sig,levels=c(0,1))
  
  #Ratio of sig to non-sig in each bin
  or1[thresh==f]<-fisher.test(table(ms2$reg,ms2$sig))$estimate
  p1[thresh==f]<-fisher.test(table(ms2$reg,ms2$sig))$p.value
  
  x=7
  sigs<-cbind.data.frame(out$id[results2[,x]<.5],as.numeric(qvalue(results[results2[,x]<.5,x])$qvalues<f))
  names(sigs)<-c("site","sig")
  ms_temp<-ms
  ms2<-merge(ms_temp,sigs,by="site")
  ms2$reg<-NA
  ms2$reg[ms2$reg_activity=="none"]<-0
  ms2$reg[ms2$reg_activity!="none"]<-1
  ms2$sig<-factor(ms2$sig,levels=c(0,1))
  
  #Ratio of sig to non-sig in each bin
  or2[thresh==f]<-fisher.test(table(ms2$reg,ms2$sig))$estimate
  p2[thresh==f]<-fisher.test(table(ms2$reg,ms2$sig))$p.value
  
  
  x=14
  sigs<-cbind.data.frame(out$id[results2[,x]<.5],as.numeric(qvalue(results[results2[,x]<.5,x])$qvalues<f))
  names(sigs)<-c("site","sig")
  ms_temp<-ms
  ms2<-merge(ms_temp,sigs,by="site")
  ms2$reg<-NA
  ms2$reg[ms2$reg_activity=="none"]<-0
  ms2$reg[ms2$reg_activity!="none"]<-1
  ms2$sig<-factor(ms2$sig,levels=c(0,1))
  
  #Ratio of sig to non-sig in each bin
  or3[thresh==f]<-fisher.test(table(ms2$reg,ms2$sig))$estimate
  p3[thresh==f]<-fisher.test(table(ms2$reg,ms2$sig))$p.value
  
}

p5c<-ggplot()+
  geom_line(aes(thresh,log2(or1)),alpha=0.5)+
  geom_line(aes(thresh,log2(or2)),alpha=0.5)+
  geom_line(aes(thresh,log2(or3)),alpha=0.5)+
  geom_point(aes(thresh,log2(or2),fill="Drought",alpha=as.numeric(p2<.05)),pch=21,size=5)+theme_bw()+
  geom_point(aes(thresh,log2(or1),fill="Rank",alpha=as.numeric(p1<.05)),pch=21,size=5)+
  geom_point(aes(thresh,log2(or3),fill="Age",alpha=as.numeric(p3<.05)),pch=21,size=5)+
  scale_alpha_continuous(range = c(0.3,1))+
  scale_fill_manual(values = c("Grey","Maroon","Steel Blue"))+
  xlab("FDR threshold")+ylab(expression("Log"[2]~"(Odds ratio)"))+
  xlim(c(0.09,.3))+ylim(c(-0.25,1))+
  geom_hline(yintercept = 0,lty=2)+
  theme(legend.title = element_blank(),text=element_text(size=20),plot.margin =margin(t = 1,1,b=98,1),legend.text = element_text(size=22),legend.position = "bottom")+
  guides(alpha="none")

  grid.arrange(p5b,p5c,nrow=1)









####################################################################################
#Supplementary Figures
####################################################################################

################################
# Figure S1
################################
  
#When were these samples collected relative to habitat quality. 
info<-read.table("./Data/meta_info_n295.txt",header=T)
length(table(info$sname)[table(info$sname)>1])

info2<-unique(cbind.data.frame(info$sname,info$matgrp,info$birth))

colnames(info2)<-c("sname","matgrp.x","birth.x")

#SI figure of habitat quality individuals at sampling
habitat_qual<-rep(0,length(unique(info$sname)))

#Those born in or before 1987 in Alto's group
length(habitat_qual[substr(info2$birth,1,4)<=1987 & info2$matgrp==1])

#Those born in or before 1991 in Hook's group
length(habitat_qual[substr(info2$birth,1,4)<=1991 & info2$matgrp==2])

habitat_qual[substr(info2$birth,1,4)<=1987 & info2$matgrp==1]<-1
habitat_qual[substr(info2$birth,1,4)<=1991 & info2$matgrp==2]<-1

info2$habitat_qual<-habitat_qual
info2$habitat_qual[which(info2$habitat_qual==0 )]<-"Grey"
info2$habitat_qual[which(info2$habitat_qual==1 & info2$matgrp==1)]<-"#5F4B8BFF"
info2$habitat_qual[which(info2$habitat_qual==1 & info2$matgrp==2)]<-"#E69A8DFF"
table(info2$habitat_qual)

info2_plot<-info2
info2_plot[order(as.Date(info2_plot$birth)),]->info2_plot
info2_plot$order<-1:256

tmp3<-unique(cbind.data.frame(info$sname,info$dart_date))
colnames(tmp3)<-c("sname","dart_date")
info3_plot<-merge(info2_plot,tmp3,by="sname")

plot(as.Date(info3_plot$birth),info3_plot$order,
     xlim=c(min(as.Date(info3_plot$birth)),
            max(as.Date(info3_plot$dart_date))),ylim=c(0,256),
     xlab="Date",ylab="Individual",pch=20,
     col=alpha(info3_plot$habitat_qual,alpha=0.8),cex=1.5)
par(new=T)
plot(as.Date(info3_plot$dart_date),info3_plot$order,
     xlim=c(min(as.Date(info3_plot$birth)),
            max(as.Date(info3_plot$dart_date))),ylim=c(0,256),pch=4,
     xlab="",ylab="",col=alpha(info3_plot$habitat_qual,alpha=0.8),cex=1.5)
segments(x0=as.Date(info3_plot$birth),x1=as.Date(info3_plot$dart_date),
         y0=info3_plot$order,y1=info3_plot$order,
         col=alpha(info3_plot$habitat_qual,alpha=0.8),lty=2)

abline(v=as.Date("1987-12-31"),col="#5F4B8BFF",lty=2)
abline(v=as.Date("1991-12-31"),col="#E69A8DFF",lty=2)

rm(list = grep(ls(),invert=T,pattern="results",value = T))

  

################################
# Figure S2
################################

########
# A
########
#How often do these effects co-occur?
info<-read.table("./Data/meta_info_n295.txt",header=T)
early_life<-cbind.data.frame(info$sname,info$habitat_quality,info$mat_loss,info$drought,info$big_grp,info$competing_sib,info$low_rank_mom)
colnames(early_life)<-c("sname","Habitat Quality","Maternal Loss","Drought","Large Group","Close-in-age Younger Sibling","Low Rank Mom")

early_life<-unique(early_life)
library(corrplot)
p<-matrix(NA,nrow=6,ncol=6)
for(x in 1:6){
  for(y in 1:6){
    if(x!=y){
      p[x,y]<-cor.test(early_life[,-1][,x],early_life[,-1][,y])$p.value
    }
  }
}

diag(p)<-NA
colnames(p)<-rownames(p)<-colnames(early_life)[-1]

par(mfrow=c(1,1),pty='s')

corrplot(cor(early_life[,-c(1)]),
         method="ellipse",type = "lower",addCoef.col="Black",tl.col = "Black",
         col = COL2('BrBG'),order="AOE",diag=T,tl.pos = 'lt')

corrplot(cor(early_life[,-c(1)]),
         insig = "p-value",sig.level = -1,number.digits = 3,p.mat = p,
         method = "ellipse",type="upper",
         order="AOE",pch.cex=0.9,col = "White",add=T,tl.pos = FALSE,cl.pos = 'n')


########
# B
########

#What about the lower quality environment individuals themselves?
p<-matrix(NA,nrow=5,ncol=5)
for(x in 1:5){
  for(y in 1:5){
    if(x!=y){
      p[x,y]<-cor.test(early_life[which(early_life$`Habitat Quality`==1),-c(1:2)][,x],early_life[which(early_life$`Habitat Quality`==1),-c(1:2)][,y])$p.value
    }
  }
}

colnames(p)<-rownames(p)<-colnames(early_life)[-c(1:2)]

diag(p)<-NA
corrplot(cor(early_life[which(early_life$`Habitat Quality`==1),-c(1,2)]),
         method="ellipse",type = "lower",addCoef.col="Black",tl.col = "Black",
         col = COL2('BrBG'),order="AOE",diag=T,tl.pos = 'lt')


corrplot(cor(early_life[which(early_life$`Habitat Quality`==1),-c(1,2)]),
         insig = "p-value",sig.level = -1,number.digits = 3,p.mat = p,
         method = "ellipse",type="upper",
         order="AOE",pch.cex=0.9,col = "White",add=T,tl.pos = FALSE,cl.pos = 'n')

rm(list = grep(ls(),invert=T,pattern="results",value = T))



################################
# Figure S3
################################

#Current versus early life rainfall?
info<-read.table("./Data/meta_info_n295.txt",header=T)

par(mfrow=c(1,2),pty="s")
plot(info$adult_rain~info$early_rain,pch=20,
     xlab="Rainfall in the first year of life",
     ylab="Rainfall in the year leading up to darting",main="All individuals (N=295 samples)")
abline(lm(info$adult_rain~info$early_rain),lty=2)
summary(lm(info$adult_rain~info$early_rain))

#What about for individuals from the poor habitat quality?
plot(info$adult_rain[info$habitat_quality==1]~info$early_rain[info$habitat_quality==1],pch=20,
     xlab="Rainfall in the first year of life",
     ylab="Rainfall in the year leading up to darting",main="Low Habitat Quality individuals (N=64 samples)")
abline(lm(info$adult_rain[info$habitat_quality==1]~info$early_rain[info$habitat_quality==1]),lty=2)
summary(lm(info$adult_rain[info$habitat_quality==1]~info$early_rain[info$habitat_quality==1]))

rm(list = grep(ls(),invert=T,pattern="results",value = T))






###########################
#Overlap of effects between age and Rank and each of the sources of early life adversity?
###########################
overlap<-matrix(NA,nrow=2,ncol=12)
rownames(overlap)<-c('rank','age')
colnames(overlap)<-gsub(colnames(results_model3)[25:36],pattern="_pvalue",replacement = "")

OR<-overlap

for(f in 1:12){
  
  keep<-which(results_model3$age_se_beta<0.6 & results_model3[,f+12] < 0.6)
  if(sum(as.numeric(qvalue(results_model3[keep,f+24])$qvalues<.1),na.rm = T)>0){
  
  overlap[1,f]<-fisher.test(table(
    as.numeric(qvalue(results_model1$male_rank_pvalue[keep])$qvalues<.1),
    as.numeric(qvalue(results_model3[keep,f+24])$qvalues<.1)
  ))$p.value
  
  OR[1,f]<-fisher.test(table(
    as.numeric(qvalue(results_model1$male_rank_pvalue[keep])$qvalues<.1),
    as.numeric(qvalue(results_model3[keep,f+24])$qvalues<.1)
  ))$estimate
  
  
  
  overlap[2,f]<-fisher.test(table(
    as.numeric(qvalue(results_model3$age_pvalue[keep])$qvalues<.1),
    as.numeric(qvalue(results_model3[keep,f+24])$qvalues<.1)
  ))$p.value 
  
  
  OR[2,f]<-fisher.test(table(
    as.numeric(qvalue(results_model3$age_pvalue[keep])$qvalues<.1),
    as.numeric(qvalue(results_model3[keep,f+24])$qvalues<.1)
  ))$estimate
  }
}
overlap<-overlap[,c(1,3,5,7,9,11,12)]
OR<-OR[,c(1,3,5,7,9,11,12)]

library(reshape2)

o1<-melt(overlap)
o2<-melt(OR)


o3<-merge(o1,o2,by=c('Var1','Var2'))
o3<-o3[-c(1,8),]
library(stringr)
o3$Var1<-str_to_sentence(o3$Var1)  
o3$Var2<-gsub(o3$Var2,pattern="_",replacement = " ")
o3$Var2<-str_to_sentence(o3$Var2)  
o3$Var2<-gsub(o3$Var2,pattern=" low quality",replacement = "")
o3$Var2[o3$Var2=="Big grp"]<-"Group size"
o3$Var2[o3$Var2=="Competing sib"]<-"Close-in-age younger sibling"
o3$Var2[o3$Var2=="drought"]<-"Drought"
o3$Var2[o3$Var2=="Low rank mom"]<-"Maternal rank"
o3$Var2[o3$Var2=="Mat loss"]<-"Maternal loss"

o3$sig<-.2
o3$sig[o3$value.x<.05]<-.8
o3$value.x<-scientific(o3$value.x,3)
o3$value.x[o3$value.x<1*10^(-10)]<-"<1.00e-10"

ggplot(data=o3)+
  geom_tile(aes(as.factor(Var1),as.factor(Var2),fill=as.factor(sign(value.y-1))),alpha=o3$sig,col="Black")+theme_bw()+
  scale_fill_manual(values = c("Steel Blue","Maroon"))+guides(fill="none")+
  theme_minimal()+theme(axis.text.x = element_text(angle=90))+
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.3,color="Black"),axis.text=element_text(size=20,color="Black"),legend.position = "top")+
  xlab("")+ylab("")+
  geom_text(aes(as.factor(Var1),as.factor(Var2),label=paste(round(log2(value.y),3),"\n",value.x)),size=6)



################################
# Figure S5
################################

#Probability of being rank associated GE | rank associated DNAm?

rank_meth<-cbind.data.frame(beta[,12],se_beta[,12],bhat[,12],pvalue[,12])
colnames(rank_meth)<-c("beta_meth","se_beta_meth","bhat_meth","pvalue_meth")
rownames(rank_meth)<-out$id

rank_meth$qvalue_meth<-qvalue(rank_meth$pvalue_meth)$qvalues
rm(beta,bhat,bse,early_mat_n295,out,pvalue,se_beta,temp,f,files)

closest_genes<-read.table("./results/closest_genes_to_all_sites.bed")
closest_genes$site<-paste(closest_genes$V1,closest_genes$V2,sep="_")

rank_genes<-read.delim("./results/ge_results/GE_rank_effects.txt",header=T)

rank_genes2<-cbind.data.frame(rank_genes$Gene.ID,rank_genes$Male.rank.beta,rank_genes$Male.rank.var.beta.,rank_genes$Male.rank.q.value)
colnames(rank_genes2)<-c("gene","beta_ge","se_beta_ge","qvalue_ge")

closest_genes<-cbind.data.frame(closest_genes$site,closest_genes$V8,closest_genes$V9)
colnames(closest_genes)<-c("site","gene","distance")
rank_meth<-rank_meth[,c(1,5)]
rank_genes<-rank_genes2[,c(1,2,4)]

rank_combined<-merge(closest_genes,rank_meth,by.x="site",by.y="row.names")
rank_combined<-merge(rank_combined,rank_genes,by="gene")

rank_combined$sig_rank_meth<-as.numeric(rank_combined$qvalue_meth<.1)
rank_combined$sig_rank_ge<-as.numeric(rank_combined$qvalue_ge<.1)

table(rank_combined$sig_rank_meth,rank_combined$sig_rank_ge)
fisher.test(table(rank_combined$sig_rank_meth,rank_combined$sig_rank_ge))

#Focus on CpG sites within genes
rank_combined2<-rank_combined[rank_combined$distance==0,]

#A gene is more likely to be rank associated if a CpG site inside of it is rank-associated. 
table(rank_combined2$sig_rank_meth,rank_combined2$sig_rank_ge)
fisher.test(table(rank_combined2$sig_rank_meth,rank_combined2$sig_rank_ge))

#Probability of being rank associated across thresholds of significance for each
prob<-seq(to=.5,from=.01,by=.01)
or<-matrix(NA,nrow=length(prob),ncol=length(prob))

for(x in prob){
  for(y in prob){
    rank_combined2$sig_rank_meth<-as.numeric(rank_combined2$qvalue_meth<x)
    rank_combined2$sig_rank_ge<-as.numeric(rank_combined2$qvalue_ge<y)
    or[prob==x,prob==y]<-fisher.test(rank_combined2$sig_rank_meth,rank_combined2$sig_rank_ge)$estimate
  }
}

library(reshape2)
rownames(or)<-colnames(or)<-prob
or_melt<-melt(or)
or_melt$value<-log2(or_melt$value)
ggplot(data=or_melt[or_melt$Var1>=.02 & or_melt$Var1<=.3 & or_melt$Var2>=.02 & or_melt$Var2<=.3,]) +
  geom_hex(aes(Var1,Var2,fill=value),stat="identity") +
  theme_minimal()+scale_fill_viridis()+xlab("Q-value threshold: rank effects on DNAM")+
  ylab("Q-value threshold: rank effects on GE")+labs(fill = "FET Odds Ratio")+
  theme(text=element_text(size=16))

ggplot(data=or_melt[or_melt$Var1>=.02 & or_melt$Var1<=.3 & or_melt$Var2>=.02 & or_melt$Var2<=.3,]) +
  geom_raster(aes(Var1,Var2,fill=value),stat="identity",interpolate = T) +
  theme_minimal()+scale_fill_viridis(breaks=c(0.1,0.3,0.5))+xlab("FDR threshold (rank effects on DNA methylation)")+
  ylab("FDR threshold (rank effects on gene expression)")+labs(fill = expression("Log"[2]~"(OR)",sep=""))+
  theme(text=element_text(size=16),legend.position = "top")


########
# B
########
#Again focus on CpG sites within genes.

closest_genes<-closest_genes[closest_genes$distance==0,]
closest_genes2<-merge(closest_genes,rank_meth,by.x="site",by.y="row.names")
closest_genes3<-merge(closest_genes2,rank_genes,by.x="gene",by.y="gene")

dev.off()
smoothScatter(closest_genes3$beta_meth~closest_genes3$beta_ge,xlab="Rank effect gene expression",ylab="Rank effect DNAm")
abline(h=0,lty=2);abline(v=0,lty=2);abline(lm(closest_genes3$beta_meth~closest_genes3$beta_ge),lty=2)

table(sign(closest_genes3$beta_meth),sign(closest_genes3$beta_ge))
fisher.test(table(sign(closest_genes3$beta_meth),sign(closest_genes3$beta_ge)))

#Significant genes
closest_genes4<-closest_genes3[which(closest_genes3$qvalue_ge<.2 & closest_genes3$qvalue_meth<.2), ]


ggplot()+
  geom_density(aes(closest_genes4$beta_ge[which(sign(closest_genes4$beta_meth)==(-1))]),fill="Dark Blue",alpha=0.6,col="Black",bw=0.05)+
  geom_density(aes(closest_genes4$beta_ge[which(sign(closest_genes4$beta_meth)==(1))]),fill="Orange",alpha=0.6,col="Black",bw=0.05)+theme_bw()+
  xlab("Rank effect on gene expression")+ylab("Density (arbitrary units)")+theme(text=element_text(size=16))+ylim(c(0,3.5))+
  geom_vline(xintercept = 0,lty=2)

