################################################################
#This script will generate the general results referenced in text. 
################################################################

#Starting with main model results, this will produce the results as they appear within the manuscript.
rm(list = grep(ls(),invert=T,pattern="results",value = T))

#Start with results matrices and meta info table
info<-read.table("./Data/meta_info_n295.txt",header=T)

#How many unique individuals, and unique males and females?
length(unique(info$sname))

colSums(table(unique(cbind.data.frame(info$sname,info$sex))))

#How many repeated individuals?
length(table(info$sname)[table(info$sname)>1])

#How many low-quality habitat samples and individuals?
colSums(table(unique(cbind.data.frame(info$sname,info$habitat_quality))))
table(info$habitat_quality)



####################
#Results for Model 1
####################
library(qvalue)

#How many sites significant at 10% FDR excluding high SE loci?
length(which(qvalue(results_model1$habitat_quality_pvalue[which(results_model1$habitat_quality_se_beta<0.6)])$qvalues<.1))

length(which(qvalue(results_model1$male_rank_pvalue[which(results_model1$male_rank_se_beta<0.6)])$qvalues<.1))
length(which(qvalue(results_model1$female_rank_pvalue[which(results_model1$female_rank_se_beta<0.6)])$qvalues<.1))

length(which(qvalue(results_model1$age_pvalue[which(results_model1$age_se_beta<0.6)])$qvalues<.1))


#Overlap between tested loci and functional compartments
#Direction of age effects
overlap<-read.table("./results/sites_in_each_basic_compartment.bed")
overlap$V7<-as.character(overlap$V7)
overlap$V7[overlap$V7=="."]<-"Unannotated"
table(overlap$V7)
overlap$site<-paste(overlap$V1,overlap$V2,sep="_")

length(unique(overlap$site))

tmp<-cbind.data.frame(rownames(results_model1),results_model1$age_bhat,results_model1$age_pvalue)
tmp<-tmp[which(results_model1$age_se_beta<0.6),]
colnames(tmp)<-c("site","bhat","pvalue")
tmp$qvalue<-qvalue(tmp$pvalue)$qvalues
tmp2<-tmp[which(tmp$qvalue<.1),]
table(sign(tmp2$bhat))

table(sign(tmp2$bhat[tmp2$site %in% overlap$site[overlap$V7=="cpg_island"]]))/sum(table(sign(tmp2$bhat[tmp2$site %in% overlap$site[overlap$V7=="cpg_island"]])))
table(sign(tmp2$bhat[tmp2$site %in% overlap$site[overlap$V7!="cpg_island"]]))/sum(table(sign(tmp2$bhat[tmp2$site %in% overlap$site[overlap$V7!="cpg_island"]])))


#No effect of cumulative early adversity.
length(which(qvalue(results_model1$cumulative_pvalue[which(results_model1$cumulative_se_beta<0.6)])$qvalues<.1))

#Clean up
rm(list = grep(ls(),invert=T,pattern="results",value = T))




####################
#Results for model 2
####################

#Main effect of habitat quality?
length(which(qvalue(results_model2$habitat_quality_pvalue[which(results_model2$habitat_quality_se_beta<0.6)])$qvalues<.1))

#Cumulative EA effect
length(which(qvalue(results_model2$cumulative_low_quality_pvalue[which(results_model2$cumulative_low_quality_se_beta<0.6)])$qvalues<.1))
length(which(qvalue(results_model2$cumulative_high_quality_pvalue[which(results_model2$cumulative_high_quality_se_beta<0.6)])$qvalues<.1))

length(which(qvalue(results_model2$cumulative_low_quality_pvalue[which(results_model2$cumulative_low_quality_se_beta<0.6)])$qvalues<.2 | qvalue(results_model2$cumulative_high_quality_pvalue[which(results_model2$cumulative_high_quality_se_beta<0.6)])$qvalues<.2))


keep<-which(results_model2$cumulative_low_quality_se_beta<0.6 & results_model2$habitat_quality_se_beta<0.6)
tmp<-cbind.data.frame(results_model2$cumulative_low_quality_pvalue[keep],results_model2$cumulative_low_quality_bhat[keep],results_model2$habitat_quality_bhat[keep])

tmp$qvalue<-qvalue(tmp[,1])$qvalues
tmp2<-tmp[which(tmp$qvalue< 0.1),]


cor.test(tmp2[,2],tmp2[,3])

#Cumulative EA low-quality vs high-quality
keep<-which(results_model2$cumulative_low_quality_se_beta<0.6 & results_model2$cumulative_high_quality_se_beta<0.6)
tmp<-cbind.data.frame(results_model2$cumulative_low_quality_pvalue[keep],results_model2$cumulative_low_quality_bhat[keep],results_model2$cumulative_high_quality_bhat[keep])

tmp$qvalue<-qvalue(tmp[,1])$qvalues
tmp2<-tmp[which(tmp$qvalue< 0.1),]

cor.test(tmp2[,2],tmp2[,3])

rm(list = grep(ls(),invert=T,pattern="results",value = T))


####################################
info<-read.table("./Data/meta_info_n295.txt",header=T)
info$cumulative_ea<-apply(info[,7:11],1,sum)


wilcox.test(info$cumulative_ea[info$habitat_quality==0],info$cumulative_ea[info$habitat_quality==1])



#What are the effects of each individual source of adversity in the low quality environment?
library(tidyverse)

for(f in seq(5)){
  n<-str_to_title(gsub(gsub(gsub(colnames(results_model3)[f*2+25],pattern="_pvalue",replacement = ""),pattern="_",replacement=" "),pattern=" low quality",replacement = ""))
  tmp<-results_model3[,c(f*2+13,f*2+25)]
  tmp<-tmp[which(tmp[,1]<0.6),]
  print(paste(n,":",length(which(qvalue(tmp[,2])$qvalues<.1)),"loci"))
}

#No real detectable early life effects in the high quality environment.
for(f in seq(5)){
  n<-str_to_title(gsub(gsub(gsub(colnames(results_model3)[f*2+24],pattern="_pvalue",replacement = ""),pattern="_",replacement=" "),pattern="high quality",replacement = ""))
  tmp<-results_model3[,c(f*2+12,f*2+24)]
  tmp<-tmp[which(tmp[,1]<0.6),]
  print(paste(n,":",length(which(qvalue(tmp[,2])$qvalues<.1)),"loci"))
}


rm(list = grep(ls(),invert=T,pattern="results",value = T))




####################################################################
#Genomic distribution of environmental predictors of DNA methylation
####################################################################

#Overlap of early life effects-- focusing on habitat quality, drought, maternal loss, group size
tmp<-results_model3[,c(1,3,5,7)+12]
tmp2<-results_model3[,c(1,3,5,7)+24]

tmp2[tmp>=0.6]<-NA
apply(tmp2,2,function(x){return(length(which(qvalue(x)$qvalues<.1)))})

tmp3<-tmp2

for(f in 1:4){
  tmp3[,f]<-as.numeric(qvalue(tmp2[,f])$qvalues<.1)
}


or<-matrix(NA,nrow=4,ncol=4)
colnames(or)<-rownames(or)<-gsub(colnames(tmp2),pattern="_pvalue",replacement = "")
p<-or

for(x in 1:4){
  for(y in 1:4){
    or[x,y]<-or[y,x]<-fisher.test(table(tmp3[,x],tmp3[,y]))$estimate
    p[x,y]<-p[y,x]<-fisher.test(table(tmp3[,x],tmp3[,y]))$p.value
  }
}


#Overlap of drought and habitat quality?
table(tmp3[,1],tmp3[,3])
fisher.test(table(tmp3[,1],tmp3[,3]))
log2(fisher.test(table(tmp3[,1],tmp3[,3]))$estimate)

tmp4<-results_model3[,c(1,3,5,7)]

tmp4[tmp>=0.6]<-NA

keep<-which(tmp3[,1]+tmp3[,3]==2)
table(sign(tmp4[keep,1]),sign(tmp4[keep,3]))

1-(7/4000)

smoothScatter(tmp4[keep,1],tmp4[keep,3])


#Overlap between rank effects and habitat quality or drought effects
tmp<-cbind.data.frame(results_model1[,8],results_model3[,c(1,5)+12])

tmp2<-cbind.data.frame(results_model1[,13],results_model3[,c(1,5)+24])

tmp2[tmp>=0.6]<-NA

tmp3<-tmp2
for(f in 1:3){
  tmp3[,f]<-as.numeric(qvalue(tmp2[,f])$qvalues<.1)
}


or<-matrix(NA,nrow=3,ncol=3)
colnames(or)<-rownames(or)<-gsub(colnames(tmp2),pattern="_pvalue",replacement = "")
p<-or

for(x in 1:3){
  for(y in 1:3){
    or[x,y]<-or[y,x]<-fisher.test(table(tmp3[,x],tmp3[,y]))$estimate
    p[x,y]<-p[y,x]<-fisher.test(table(tmp3[,x],tmp3[,y]))$p.value
  }
}


tmp4<-cbind.data.frame(results_model1[,3],results_model3[,c(1,5)])

tmp4[tmp>=0.6]<-NA

keep<-which(tmp3[,1]+tmp3[,2]==2)
log2(fisher.test(table(sign(tmp4[keep,1]),sign(tmp4[keep,2])))$estimate)
log2(fisher.test(table(sign(tmp4[keep,1]),sign(tmp4[keep,3])))$estimate)



#Including age
tmp<-cbind.data.frame(results_model1[,8],results_model3[,c(1,5,12)+12])

tmp2<-cbind.data.frame(results_model1[,13],results_model3[,c(1,5,12)+24])

tmp2[tmp>=0.6]<-NA

tmp3<-tmp2
for(f in 1:4){
  tmp3[,f]<-as.numeric(qvalue(tmp2[,f])$qvalues<.1)
}


or<-matrix(NA,nrow=4,ncol=4)
colnames(or)<-rownames(or)<-gsub(colnames(tmp2),pattern="_pvalue",replacement = "")
p<-or

for(x in 1:4){
  for(y in 1:4){
    or[x,y]<-or[y,x]<-fisher.test(table(tmp3[,x],tmp3[,y]))$estimate
    p[x,y]<-p[y,x]<-fisher.test(table(tmp3[,x],tmp3[,y]))$p.value
  }
}

log2(or[,4])






rm(list = grep(ls(),invert=T,pattern="results",value = T))




###################################################
# Enrichment in functional compartments
# Age, early life habitat quality, drought, and male rank
###################################################

overlap<-read.table("./results/sites_in_each_basic_compartment.bed")
overlap$V7<-as.character(overlap$V7)
overlap$V7[overlap$V7=="."]<-"Unannotated"
table(overlap$V7)

OR<-matrix(NA,nrow=6,ncol=4)
p<-matrix(NA,nrow=6,ncol=4)

results_temp<-cbind.data.frame(results_model1[,13],results_model3[,c(25,29,36)])
results2_temp<-cbind.data.frame(results_model1[,8],results_model3[,c(13,17,24)])

colnames(OR)<-colnames(p)<-c("Male rank","Habitat quality","Drought (low-habitat quality)","Age")
rownames(OR)<-rownames(p)<-unique(overlap$V7)

overlap$site<-paste(overlap$V1,overlap$V2,sep="_")
u<-unique(overlap$V7)

for(x in 1:4){
  sigs<-cbind.data.frame(rownames(results_model1)[which(results2_temp[,x]<.6)],as.numeric(qvalue(results_temp[which(results2_temp[,x]<.6),x])$qvalues<.1))
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


colnames(OR)
rownames(OR)

#Enrichment of effects in gene bodies?
log2(OR[5,])
p[5,]

#Enhancers?
log2(OR[1,])
p[1,]

#Unannotated regions?
log2(OR[6,])
p[6,]

#What about age effects?
OR[,4]
log2(OR[,4])
p[4,]





###############################################
#Enrichment of effects in chromHMM annotations.
###############################################
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


colnames(OR)<-c("Male rank","Habitat quality","Drought (low-habitat quality)","Age")
colnames(p)<-c("Male rank","Habitat quality","Drought (low-habitat quality)","Age")

#Drought and dominance rank in enhancers
log2(OR[rownames(OR)=="Enh",])

#Transcription
log2(OR[rownames(OR)=="Tx",])


#Heterochromatin
log2(OR[rownames(OR)=="Het",])

#Weakly repressed, polycomb-marked
log2(OR[rownames(OR)=="ReprPCWk",])




rm(list = grep(ls(),invert=T,pattern="results",value = T))
rm(results_temp,results,results2_temp,results2)











##########################################################
# Attenuation of early life habitat quality over time. 
##########################################################

info<-read.table("./Data/meta_info_n295.txt",header=T)
predicted_binomial<-read.table("./results/habitat_quality_prediction_alpha1_50fold_binomial.txt")
info$predicted_binomial_hab<-predicted_binomial$V1

info$habitat_qual<-info$habitat_quality
info$habitat_qual[info$habitat_qual==0]<-"post_shift"
info$habitat_qual[which(info$habitat_qual==1 & info$matgrp.x==1)]<-"Altos"
info$habitat_qual[which(info$habitat_qual==1 & info$matgrp.x==2)]<-"Hooks"

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

integrate(ROC,lower=0,upper=1,subdivisions = 100000)


#Time since habitat shift? 
info$time_since<-NA
info$time_since[info$habitat_quality==1 & info$matgrp.x==1]<-as.numeric(as.Date(info$dart_date[info$habitat_quality==1 & info$matgrp.x==1])-as.Date("1988-01-01"))
info$time_since[info$habitat_quality==1 & info$matgrp.x==2]<-as.numeric(as.Date(info$dart_date[info$habitat_quality==1 & info$matgrp.x==2])-as.Date("1992-01-01"))

info$time_since[which(info$time_since<0)]<-0

summary(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$time_since[info$habitat_quality==1]))

summary(lm(info$predicted_binomial_hab~info$age))

plot(info$predicted_binomial_hab~info$age)

plot(info$predicted_binomial_hab[info$habitat_quality==1]~info$time_since[info$habitat_quality==1])

#What about sample age?
info$sample_age<-as.numeric(as.Date('2023-05-10')-as.Date(info$dart_date))

plot(info$predicted_binomial_hab[info$habitat_quality==1]~info$sample_age[info$habitat_quality==1])
abline(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$sample_age[info$habitat_quality==1]),lty=2)

plot(info$predicted_binomial_hab~info$sample_age)
plot(info$predicted_binomial_hab~info$habitat_quality)



#Duration of time in low habitat quality?
info$duration<-NA
#Alto is 1, in or before 1987
sort(info$birth.x[info$habitat_quality==1 & round(info$grp,0)==1])


info$duration[info$habitat_quality==1 & round(info$grp,0)==1]<- as.Date("1987-12-31")-as.Date(info$birth.x[info$habitat_quality==1 & round(info$grp,0)==1])

#Hooks is 2, in or before 1991
sort(info$birth.x[info$habitat_quality==1 & round(info$grp,0)==2])

info$duration[info$habitat_quality==1 & round(info$grp,0)==2]<- as.Date("1991-12-31")-as.Date(info$birth.x[info$habitat_quality==1 & round(info$grp,0)==2])


plot(info$predicted_binomial_hab[info$habitat_quality==1]~info$duration[info$habitat_quality==1])
summary(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$duration[info$habitat_quality==1]))

summary(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$duration[info$habitat_quality==1]+info$time_since[info$habitat_quality==1]+info$age[info$habitat_quality==1]))

summary(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$age[info$habitat_quality==1]))
summary(lm(info$predicted_binomial_hab[info$habitat_quality==1]~info$time_since[info$habitat_quality==1]+info$age[info$habitat_quality==1]))



#sample age?
info$sample_age<-as.Date('2023-02-20')-as.Date(info$dart_date)

summary(lm(info$predicted_binomial_hab~info$sample_age))
summary(lm(info$predicted_binomial_hab~info$habitat_quality))


plot(info$sample_age,info$predicted_binomial_hab)

par(pty="s")
plot(info$habitat_quality,info$sample_age,xlab="Habitat quality",ylab="~Sample age",col="Steel Blue",pch=20)



rm(list = grep(ls(),invert=T,pattern="results",value = T))










############################
#MSTARR results
############################
ms<-read.table("./results/misc_results/baboon_mstarr_summary_table_31Dec21.txt",header=T)

#How many of these unique windows overlap test loci?
loci<-str_split(rownames(results_model3),pattern = "_",simplify = T)
loci<-as.data.frame(loci)
loci[,2]<-as.numeric(loci[,2])


library(parallel)

ms$keep<-FALSE
no_cores<-detectCores()
cl <- makeCluster(no_cores-2)
clusterExport(cl = cl, varlist = c("loci","ms"), envir = environment())
library(parallel)

tmp<-parSapply(cl=cl,X = 1:nrow(loci),function(f){
  keep<-which(ms$chr==loci[f,1] & ms$start<= loci[f,2] & ms$end >= loci[f,2])
  tmp<-FALSE
  if(length(keep)>0){
    tmp<-TRUE
  }
  return(tmp)
})

#write.table(tmp,"~/Desktop/tmp.txt")

#We keep about 94000 loci.
tmp<-read.table("~/Desktop/tmp.txt")
table(tmp)


#mstarr results
ms<-read.table("./results/misc_results/mstarr_overlap_with_tested_sites.bed")
ms2<-read.table("./results/misc_results/baboon_mstarr_summary_table_31Dec21.txt",header=T)

colnames(ms)[4:27]<-colnames(ms2)


ms$reg_activity[is.na(ms$reg_activity)]<-"none"
table(ms$reg_activity)/sum(length(ms$reg_activity))

table(ms$meth_depend[ms$reg_activity!="none"])


#How many unique fragments?
ms3<-unique(ms[ms$site %in% ms2$site,-c(1:3)])
table(ms3$source)

table(ms3$reg_activity=="none")
table(ms3$reg_activity=="none")/length(ms3$reg_activity)


#Enrichment in chromHMM annotations?
ms_enrichment<-read.table("results/misc_results/mstarr_results_to_plot",header=T)

ms_enrichment









#Enrichment of rank, drought, and age effects across FDR thresholds.
ms_full<-read.table("./results/misc_results/baboon_mstarr_summary_table_31Dec21.txt",header=T)
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

ms$duplicated<-duplicated(ms$site)
ms<-ms[!ms$duplicated,]

#Range of thresholds for rank + drought + age
thresh<-seq(.1,.3,by=.01)
or3<-or2<-or1<-rep(NA,length(thresh))
p3<-p2<-p1<-rep(NA,length(thresh))

for(f in thresh){
  sigs<-cbind.data.frame(rownames(results_model1)[results_model1$male_rank_se_beta<.5],as.numeric(qvalue(results_model1$male_rank_pvalue[results_model1$male_rank_se_beta<.5])$qvalues<f))
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
  
  sigs<-cbind.data.frame(rownames(results_model1)[results_model3$drought_low_quality_se_beta<.5],as.numeric(qvalue(results_model3$drought_low_quality_pvalue[results_model3$drought_low_quality_se_beta<.5])$qvalues<f))
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
  
  sigs<-cbind.data.frame(rownames(results_model1)[results_model3$age_se_beta<.5],as.numeric(qvalue(results_model3$age_pvalue[results_model3$age_se_beta<.5])$qvalues<f))
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

ggplot()+
  geom_line(aes(thresh,log2(or1)),alpha=0.5)+
  geom_line(aes(thresh,log2(or2)),alpha=0.5)+
  geom_line(aes(thresh,log2(or3)),alpha=0.5)+
  geom_point(aes(thresh,log2(or2),fill="Drought",alpha=as.numeric(p2<.05)),pch=21,size=5)+theme_bw()+
  geom_point(aes(thresh,log2(or1),fill="Rank",alpha=as.numeric(p1<.05)),pch=21,size=5)+
  geom_point(aes(thresh,log2(or3),fill="Age",alpha=as.numeric(p3<.05)),pch=21,size=5)+
  scale_alpha_continuous(range = c(0.3,1))+
  scale_fill_manual(values = c("Grey","Maroon","Steel Blue"))+
  xlab("FDR threshold")+ylab(expression("Log"[2]~"(Odds ratio)"))+
  xlim(c(0.09,.3))+
  geom_hline(yintercept = 0,lty=2)+
  theme(legend.title = element_blank(),text=element_text(size=20),plot.margin =margin(t = 1,1,b=98,1),legend.text = element_text(size=22),legend.position = "bottom")+
  guides(alpha="none")





#########
#15% FDR?
f=0.15
sigs<-cbind.data.frame(rownames(results_model1)[results_model1$male_rank_se_beta<.5],as.numeric(qvalue(results_model1$male_rank_pvalue[results_model1$male_rank_se_beta<.5])$qvalues<f))
names(sigs)<-c("site","sig")
ms_temp<-ms
ms2<-merge(ms_temp,sigs,by="site")
ms2$reg<-NA
ms2$reg[ms2$reg_activity=="none"]<-0
ms2$reg[ms2$reg_activity!="none"]<-1
ms2$sig<-factor(ms2$sig,levels=c(0,1))

#Ratio of sig to non-sig in each bin
table(ms2$reg,ms2$sig)
table(ms2$reg,ms2$sig)[2,2]/sum(table(ms2$reg,ms2$sig)[,2])

fisher.test(table(ms2$reg,ms2$sig))

#Enrichment of MD regulatory in these loci?
ms2$meth_depend2<-ms2$meth_depend
ms2$meth_depend2[ms2$meth_depend!="none"]<-"yes"

table(ms2$sig,ms2$meth_depend2)
length(which(ms2$sig==1 & ms2$meth_depend2=="yes" & ms2$reg==1))/length(which(ms2$sig==1& ms2$reg==1))


fisher.test(table(ms2$sig,ms2$meth_depend!="none"))




#Drought
sigs<-cbind.data.frame(rownames(results_model1)[results_model3$drought_low_quality_se_beta<.5],as.numeric(qvalue(results_model3$drought_low_quality_pvalue[results_model3$drought_low_quality_se_beta<.5])$qvalues<f))
names(sigs)<-c("site","sig")
ms_temp<-ms
ms2<-merge(ms_temp,sigs,by="site")
ms2$reg<-NA
ms2$reg[ms2$reg_activity=="none"]<-0
ms2$reg[ms2$reg_activity!="none"]<-1
ms2$sig<-factor(ms2$sig,levels=c(0,1))

#Ratio of sig to non-sig in each bin
table(ms2$reg,ms2$sig)
table(ms2$reg,ms2$sig)[2,2]/sum(table(ms2$reg,ms2$sig)[,2])

fisher.test(table(ms2$reg,ms2$sig))

#Enrichment of MD regulatory in these loci?
ms2$meth_depend2<-ms2$meth_depend
ms2$meth_depend2[ms2$meth_depend!="none"]<-"yes"

table(ms2$sig,ms2$meth_depend2)
length(which(ms2$sig==1 & ms2$meth_depend2=="yes" & ms2$reg==1))/length(which(ms2$sig==1& ms2$reg==1))

fisher.test(table(ms2$sig,ms2$meth_depend))
fisher.test(table(ms2$sig,ms2$meth_depend!="none"))



sigs<-cbind.data.frame(rownames(results_model3)[results_model3$age_se_beta<.5],as.numeric(qvalue(results_model3$age_pvalue[results_model3$age_se_beta<.5])$qvalues<f))
names(sigs)<-c("site","sig")
ms_temp<-ms
ms2<-merge(ms_temp,sigs,by="site")
ms2$reg<-NA
ms2$reg[ms2$reg_activity=="none"]<-0
ms2$reg[ms2$reg_activity!="none"]<-1
ms2$sig<-factor(ms2$sig,levels=c(0,1))

#Ratio of sig to non-sig in each bin
table(ms2$reg,ms2$sig)
table(ms2$reg,ms2$sig)[2,2]/sum(table(ms2$reg,ms2$sig)[,2])

fisher.test(table(ms2$reg,ms2$sig))
log2(fisher.test(table(ms2$reg,ms2$sig))$estimate)






#How many mSTARR fragments are methylation dependent?
ms<-read.table("results/misc_results/baboon_mstarr_summary_table_31Dec21.txt")
tmp<-read.table("~/Desktop/tmp.txt")

ms<-ms[ms]




#Cell type models
#out<-read.table("~/Desktop/model1_chr1_w_cell_type_neutrophil_lymphocyte.assoc.txt",header=T)
out<-read.table("~/Desktop/model1_chr1_w_cell_type.assoc.txt",header=T)

hist(out$pvalue,breaks=100)

head(out)

#Calculating p-values from MACAU
beta<-cbind.data.frame(out[,seq(11,dim(out)[2],2)],out[,4])
se_beta<-cbind.data.frame(out[,seq(12,dim(out)[2],2)],out[,5])
bhat=as.matrix(beta/(1-se_beta^2))
bse=as.matrix(se_beta/sqrt(1-se_beta^2))
pvalue=1-pchisq((bhat/bse)^2,1)

dim(pvalue)
par(mfrow=c(4,4),pty="s")

cols<-colnames(pvalue)<-c("intercept","mol_ecol_dark","mol_ecol_non_dark","new_batch",
                          "november","drrbs","r21","mapped_reads",
                          "bscr","habitat_quality","cumulative","male_rank","female_rank","cell1","cell2","age")


  
  
for(f in 1:16){
  hist(pvalue[,f],xlab="pvalue",breaks=100,main=cols[f])
  
}

#how many males have blood smears
a<-read.table("~/Desktop/covariates_model1_with_smear_cell_type.txt")
table(info$sex.x[is.na(a[,15])])

par(mfrow=c(1,1),pty="s")
hist(info$rank[is.na(a[,15]) & info$sex.x=="M"])

c<-read.table("./results/macau_output/model1_chr1.assoc.txt",header=T)
smoothScatter(c$alpha12,out$alpha12)
smoothScatter(c$alpha11,out$alpha11)



b<-read.table("results/macau_output/model2_chr1.assoc.txt",header=T)
par(mfrow=c(1,1))
plot(b$pvalue,pvalue[,18])

#Age
smoothScatter(out$beta,b$beta)
summary(lm(out$beta~b$beta))

#Cumulative EA low environment







abline(0,1,lty=2)


par(mfrow=c(1,2),pty="s")
hist(pvalue[se_beta[,16]<.6,16],breaks=1000,xlab="pvalue",main="Lymphocytes")
hist(pvalue[se_beta[,17]<.6,17],breaks=1000,xlab="pvalue",main="Neutrophils")

library(qvalue)
length(which(qvalue(pvalue[se_beta[,16]<.6,16])$qvalues<.1))
length(which(qvalue(pvalue[,17])$qvalues<.1))

par(mfrow=c(1,1))
hist(pvalue[,18],breaks=1000)

fisher.test(table(qvalue(pvalue[,17])$qvalues<.1,qvalue(out$pvalue)$qvalues<.1))
fisher.test(table(qvalue(pvalue[,16])$qvalues<.1,qvalue(pvalue[,17])$qvalues<.1))

#
a<-results_model3[grep(rownames(results_model3),pattern="chr1_"),]
keep<-a$

fisher.test(table(qvalue(pvalue[,17])$qvalues<.1,qvalue(a$habitat_quality_pvalue)$qvalues<.1))


wbc<-read.csv("./Data/wbc_counts.csv",header=T)
head(wbc)




###################################################
#Male rank on gene expression versus DNA methylation
###################################################

rank_meth<-cbind.data.frame(results_model1$male_rank_bhat,results_model1$male_rank_se_beta,results_model1$male_rank_bhat,results_model1$male_rank_pvalue)
colnames(rank_meth)<-c("beta_meth","se_beta_meth","bhat_meth","pvalue_meth")
rownames(rank_meth)<-rownames(results_model1)

rank_meth$qvalue_meth<-qvalue(rank_meth$pvalue_meth)$qvalues

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


#Is a male-rank associated site closer to rank associated genes than non rank associated genes?
ks.test(rank_combined$distance[rank_combined$sig_rank_meth==0 &rank_combined$sig_rank_ge==1],
        rank_combined$distance[rank_combined$sig_rank_meth==1 &rank_combined$sig_rank_ge==1])

par(pty="s")
qqplot(rank_combined$distance[rank_combined$sig_rank_meth==0 &rank_combined$sig_rank_ge==1],
       rank_combined$distance[rank_combined$sig_rank_meth==1 &rank_combined$sig_rank_ge==1],
       xlab="Distance of non-rank sites from rank-associated genes",
       ylab="Distance of rank-associated sites from rank-associated genes",pch=20)
abline(0,1,lty=2)

mean(rank_combined$distance[rank_combined$sig_rank_meth==0 &rank_combined$sig_rank_ge==1],na.rm=T)
mean(rank_combined$distance[rank_combined$sig_rank_meth==1 &rank_combined$sig_rank_ge==1],na.rm=T)

#And much more likely to fall within a gene
fisher.test(table(rank_combined$distance[rank_combined$sig_rank_ge==1]==0,rank_combined$sig_rank_meth[rank_combined$sig_rank_ge==1]))

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
library(viridis)
rownames(or)<-colnames(or)<-prob
or_melt<-melt(or)
or_melt$value<-log2(or_melt$value)
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
fisher.test(table(sign(closest_genes4$beta_meth),sign(closest_genes4$beta_ge)))


ggplot()+
  geom_density(aes(closest_genes4$beta_ge[which(sign(closest_genes4$beta_meth)==(-1))]),fill="Dark Blue",alpha=0.6,col="Black",bw=0.05)+
  geom_density(aes(closest_genes4$beta_ge[which(sign(closest_genes4$beta_meth)==(1))]),fill="Orange",alpha=0.6,col="Black",bw=0.05)+theme_bw()+
  xlab("Rank effect on gene expression")+ylab("Density (arbitrary units)")+theme(text=element_text(size=16))+ylim(c(0,3.5))+
  geom_vline(xintercept = 0,lty=2)





#Gene set enrichment analyses of rank effects on DNA methylation
################################################################################################
############################# Gene Set Enrichment Analyses #####################################
################################################################################################
library(doParallel)
library(parallel)
library(qusage)
library(qvalue)
library(viridis)
library(ggplot2)
library(tidyverse)

########################################################################
#Source GSEA code from the BROAD institute. 
#https://github.com/GSEA-MSigDB/GSEA_R
#I do not claim to have made or maintain the "GSEA.EnrichmentScore code. 
########################################################################

GSEA.EnrichmentScore <- function(gene.list, gene.set, weighted.score.type = 1, correl.vector = NULL) {
  
  tag.indicator <- sign(match(gene.list, gene.set, nomatch = 0))  # notice that the sign is 0 (no tag) or 1 (tag)
  no.tag.indicator <- 1 - tag.indicator
  N <- length(gene.list)
  Nh <- length(gene.set)
  Nm <- N - Nh
  if (weighted.score.type == 0) {
    correl.vector <- rep(1, N)
  }
  alpha <- weighted.score.type
  correl.vector <- abs(correl.vector^alpha)
  sum.correl.tag <- sum(correl.vector[tag.indicator == 1])
  norm.tag <- 1/sum.correl.tag
  norm.no.tag <- 1/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  max.ES <- max(RES)
  min.ES <- min(RES)
  if (max.ES > -min.ES) {
    # ES <- max.ES
    ES <- signif(max.ES, digits = 5)
    arg.ES <- which.max(RES)
  } else {
    # ES <- min.ES
    ES <- signif(min.ES, digits = 5)
    arg.ES <- which.min(RES)
  }
  return(list(ES = ES, arg.ES = arg.ES, RES = RES, indicator = tag.indicator))
}

############################################################
#Read in the MSIG database Hallmark pathways
#List of 50 different hallmark pathways with associated genes
############################################################
hallmark<-read.gmt("Data/h.all.v7.2.symbols.gmt")


#For each gene set we need all of the tested CpG sites, and we need a single gene associated with each CpG site
#We'll do the latter using "bedtools closest"
#And exclude genes greater than X KB away

#Read in closest genes
closest_genes<-read.table("./results/closest_genes_to_all_sites.bed")
closest_genes$site<-paste(closest_genes$V1,closest_genes$V2,sep="_")

#If they're in two genes, make the gene NA
#If they're > 20kb, make the gene NA
hist(closest_genes$V9)

plot(ecdf(closest_genes$V9/1000),xlab="Distance of tested site from nearest gene (kb)",
     ylab="Cumulative proportion of sites",main="",xlim=c(0,1000))
ecdf(closest_genes$V9/1000)(10)
ecdf(closest_genes$V9/1000)(100)

rank_genes<-read.delim("./results/ge_results/rank_effects_ge.txt",header=T)
rank_genes<-rank_genes[rank_genes$Male.rank.q.value<.1,]

closest_genes2<-closest_genes[closest_genes$V8 %in% rank_genes$Gene.ID,]

r<-rownames(results_model1)[which(qvalue(results_model1$male_rank_pvalue)$qvalues<.1)]
r<-as.character(r)
closest_genes3<-closest_genes2[closest_genes2$site %in% r,]

qqplot(closest_genes2$V9,closest_genes3$V9,pch=20,xlab="Distance of tested sites from male rank associated genes",
       ylab="Distance of male rank associated sites from male rank associated genes")
abline(0,1,lty=2)

ks.test(closest_genes2$V9,closest_genes3$V9)


#Where are rank associated sites?
closest_genes$duplicate<-duplicated(closest_genes$site,)
closest_genes$duplicate2<-duplicated(closest_genes$site,fromLast = T)

table(closest_genes$duplicate)
closest_genes$V8[closest_genes$duplicate]<-NA
closest_genes$V8[closest_genes$duplicate2]<-NA

closest_genes<-closest_genes[-which(closest_genes$duplicate==TRUE),]
length(which(is.na(closest_genes$V8)))

#closest_genes$V8[closest_genes$V9>0]<-NA

#Enrichment of rank associated sites?
#Rank
bhats<-results_model1$male_rank_bhat
names(bhats)<-rownames(results_model1)

closest_genes<-closest_genes[order(closest_genes$site),]
names(bhats)<-closest_genes$V8
bhats<-abs(bhats)
bhats<-bhats[order(bhats)]

hall_enrich<-as.data.frame(matrix(NA,nrow=50,ncol=2))
colnames(hall_enrich)<-c("rank_ES","rank_p")

for(f in 1:50){
  bhats<-bhats[(!is.na(bhats)) & (!is.infinite(bhats))]
  hall_enrich[f,1]<-GSEA.EnrichmentScore(gene.list = names(bhats),
                                         gene.set = unlist(hallmark[[f]]),
                                         correl.vector = bhats,weighted.score.type = 1)$ES
  
  temp<-1:1000
  clus <- makeCluster(7)
  registerDoParallel(cores=7) 
  clusterExport(clus,varlist =c("hallmark",'bhats',"GSEA.EnrichmentScore","f","temp"),envir=environment())
  temp2<-parSapply(cl = clus,X = 1:1000,FUN = function(i){
    #Randomly shuffle standardized betas across genes
    rand<-sample(1:length(bhats),size= length(bhats),replace=FALSE)
    bhats_perm<-bhats[rand]
    names(bhats_perm)<-names(bhats)
    
    bhats_perm<-bhats_perm[order(bhats_perm)]
    temp[i]<-GSEA.EnrichmentScore(gene.list = names(bhats_perm), gene.set = unlist(hallmark[[f]]),correl.vector = bhats_perm,weighted.score.type = 1)$ES
    return(temp[i])
  })
  stopCluster(clus)
  
  #The pvalue is the % of abs(observed ES)< abs(permuted ES)  
  hall_enrich[f,2]<-length(which(abs(hall_enrich[f,1]) < abs(temp2)))/1000
  
  print(f)
}

rownames(hall_enrich)<-names(hallmark)
hall_enrich_rank<-hall_enrich

hall_enrich_rank$name<-rownames(hall_enrich_rank)
hall_enrich_rank %>% mutate(name= fct_reorder(name,rank_p)) %>%
  ggplot()+geom_bar(aes(x=name,y=rank_p),stat="identity",col="Black",fill=viridis(50))+coord_flip()+theme_bw()








#Correlations between drought and habitat quality is really high?
keep<-which(results_model3$habitat_quality_se_beta<.6 & results_model3$drought_low_quality_se_beta<.6)
p1<-results_model3$habitat_quality_pvalue[keep]
p2<-results_model3$drought_low_quality_pvalue[keep]
b1<-results_model3$habitat_quality_bhat[keep]
b2<-results_model3$drought_low_quality_bhat[keep]

keep2<-which(qvalue(p1)$qvalues<.1 | qvalue(p2)$qvalues<.1)


keep2<-which(qvalue(p1)$qvalues<.1 & qvalue(p2)$qvalues<.1)
table(sign(b1[keep2]),sign(b2[keep2]))

#Is this always so insane? Lets look at non-sig sites
keep2<-which(qvalue(p1)$qvalues>.6 & qvalue(p2)$qvalues>.6)

table(sign(b1[keep2]),sign(b2[keep2]))

fisher.test(table(sign(b1[keep2]),sign(b2[keep2])))








