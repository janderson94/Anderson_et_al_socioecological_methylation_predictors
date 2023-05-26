#Main model results output from MACAU (all controlling for relatedness)
#Model 1 -- Batch/technical effects + Cumulative Early Adversity + Habitat Quality
#Model 2 -- Batch/technical effects + Habitat Quality + Habitat Quality:Cumulative Early Adversity
#Model 3 -- Batch/technical effects + Habitat Quality + Habitat Quality:Individual sources of adversity

#Single dependency
library(qvalue)

########################
####### Model 1 ########
########################
files<-list.files("./results/macau_output/",pattern = "model1")[grep(x=list.files(path = "./results/macau_output/",pattern = "model1"),pattern = "w_cell_type",invert=T)]

out<-read.table(paste("./results/macau_output/",files[1],sep=""),header=T)
for(f in 2:length(files)){
  temp<-read.table(paste("./results/macau_output/",files[f],sep=""),header=T)
  out2<-rbind.data.frame(out,temp)
  out<-out2
  rm(out2)
}

#Calculating p-values from MACAU
beta<-cbind.data.frame(out[,seq(11,dim(out)[2],2)],out[,4])
se_beta<-cbind.data.frame(out[,seq(12,dim(out)[2],2)],out[,5])
bhat=as.matrix(beta/(1-se_beta^2))
bse=as.matrix(se_beta/sqrt(1-se_beta^2))
pvalue=1-pchisq((bhat/bse)^2,1)

colnames(pvalue)<-c("intercept","mol_ecol_dark","mol_ecol_non_dark","new_batch",
                    "november","drrbs","r21","mapped_reads",
                    "bscr","habitat_quality","cumulative","male_rank","female_rank","age")


#Store the necessary outputs -- Habitat quality, male rank, and female rank effects
results_model1<-cbind.data.frame(bhat[,c(10,11,12,13,14)],se_beta[,c(10,11,12,13,14)],pvalue[,c(10,11,12,13,14)])
colnames(results_model1)<-apply(expand.grid(colnames(pvalue)[c(10,11,12,13,14)], c("bhat","se_beta","pvalue")), 1, paste, collapse="_")
rownames(results_model1)<-out$id

########################################################################################################################

########################
####### Model 2 ########
########################
files<-list.files("./results/macau_output/",pattern = "model2")

out<-read.table(paste("./results/macau_output/",files[1],sep=""),header=T)
for(f in 2:length(files)){
  temp<-read.table(paste("./results/macau_output/",files[f],sep=""),header=T)
  out2<-rbind.data.frame(out,temp)
  out<-out2
  rm(out2)
}

#Calculating p-values from MACAU
beta<-cbind.data.frame(out[,seq(11,dim(out)[2],2)],out[,4])
se_beta<-cbind.data.frame(out[,seq(12,dim(out)[2],2)],out[,5])
bhat=as.matrix(beta/(1-se_beta^2))
bse=as.matrix(se_beta/sqrt(1-se_beta^2))
pvalue=1-pchisq((bhat/bse)^2,1)


colnames(pvalue)<-c("intercept","mol_ecol_dark","mol_ecol_non_dark","new_batch",
                    "november","drrbs","r21","mapped_reads",
                    "bscr","habitat_quality","cumulative_high_quality","cumulative_low_quality","age")


#Store the necessary outputs -- Habitat quality, Cumulative EA (high quality env), Cumulative EA (low quality env)
results_model2<-cbind.data.frame(bhat[,10:12],se_beta[,10:12],pvalue[,10:12])
colnames(results_model2)<-apply(expand.grid(colnames(pvalue)[10:12], c("bhat","se_beta","pvalue")), 1, paste, collapse="_")
rownames(results_model2)<-out$id


########################################################################################################################

########################
####### Model 3 ########
########################

files<-list.files("./results/macau_output/",pattern = "model3")

out<-read.table(paste("./results/macau_output/",files[1],sep=""),header=T)
for(f in 2:length(files)){
  temp<-read.table(paste("./results/macau_output/",files[f],sep=""),header=T)
  out2<-rbind.data.frame(out,temp)
  out<-out2
  rm(out2)
}

beta<-cbind.data.frame(out[,seq(11,dim(out)[2],2)],out[,4])
se_beta<-cbind.data.frame(out[,seq(12,dim(out)[2],2)],out[,5])
bhat=as.matrix(beta/(1-se_beta^2))
bse=as.matrix(se_beta/sqrt(1-se_beta^2))
pvalue=1-pchisq((bhat/bse)^2,1)

colnames(pvalue)<-c("intercept","mol_ecol_dark","mol_ecol_non_dark","new_batch",
                    "november","drrbs","r21","mapped_reads",
                    "bscr","habitat_quality","mat_loss_high_quality","mat_loss_low_quality","drought_high_quality",
                    "drought_low_quality","big_grp_high_quality","big_grp_low_quality","low_rank_mom_high_quality",
                    "low_rank_mom_low_quality","close_in_age_sib_high_quality","close_in_age_sib_low_quality","age")


results_model3<-cbind.data.frame(bhat[,10:21],se_beta[,10:21],pvalue[,10:21])
colnames(results_model3)<-apply(expand.grid(colnames(pvalue)[10:21], c("bhat","se_beta","pvalue")), 1, paste, collapse="_")
rownames(results_model3)<-out$id


#Remove everything but the results
rm(list = grep(ls(),invert=T,pattern="results",value = T))



