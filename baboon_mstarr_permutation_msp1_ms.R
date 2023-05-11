##### function to calculate fdr from permutation p values
library(cobs)
perm.fdr=function(input_df,perm_df,Pvals_col_name,name){
    
    pvals_index=which(colnames(input_df)==Pvals_col_name)
    ro<-input_df[order(input_df[,pvals_index]),]
    p_obs <- data.frame(pvalue=ro[,pvals_index])
    p_vector<-matrix(as.matrix(perm_df),ncol=1)
    p_vector=data.frame(p_vector[order(p_vector)])
    
    F<-p_obs[,1]
    F_o<-p_obs[,1]
    pi_hat<-p_obs[,1]
    
    j=1
    observed<-length(p_obs[,1])
    randoms<-length(p_vector[,1])
    
    for(i in 1:observed)
    {
        repeat
        {
            if((p_vector[j,1]<p_obs[i,1])&j<randoms){j<-j+1}else{break}
        }
        F[i]=i/observed
        F_o[i]=(j-1)/randoms
        if(F_o[i]<1){pi_hat[i]=(1-F[i])/(1-F_o[i])}else{pi_hat[i]=1}
    }
    tabla <-data.frame(pi_hat,pval=p_obs[,1])
    
    tabla[1,]=c(1,0)
    last_percentile_average=mean(tabla$pi_hat[as.integer(min((length(tabla[,1])*0.99),(nrow(tabla)-1)):length(tabla[,1]))])
    tabla[nrow(tabla),]=c(last_percentile_average,1)
    constraint_matrix=as.matrix(data.frame(c(0,2),c(0,1),c(1,0)))
    f_hat<-suppressWarnings(cobs(tabla$pval,tabla$pi_hat,constraint="convex",pointwise=constraint_matrix,maxiter=1000,print.warn=FALSE,print.mesg=FALSE))
    
    f_hat_serie=f_hat$fitted
    pi_o=f_hat_serie[length(f_hat_serie)]
    pi_o=min(pi_o,1)
    pi_o=max(pi_o,0)
    
    Fdr_ST_perm=pi_o*F_o/F
    
    for(i in 1:length(p_obs[,1]))
    {
        Fdr_ST_perm[i]=pi_o*F_o[i]/F[i]
        if(i>1)
        {
            for(j in 1:(i-1))
            {
                if(Fdr_ST_perm[i-j]>Fdr_ST_perm[i]){Fdr_ST_perm[i-j]=Fdr_ST_perm[i]}else{break}
            }
        }
        if(Fdr_ST_perm[i]>1)  Fdr_ST_perm[i]=1
    }
    
    fdrs_df <-data.frame(ro,q_ST_perm=Fdr_ST_perm)
    rownames(fdrs_df)=rownames(ro)
    colnames(fdrs_df)[ncol(fdrs_df)]=paste0("fdr_",name)
    
    return(fdrs_df)
}

#####################################
# Create PERMUTATIONS of RNA/DNA label, keeping samples paired. 
#####################################
library(tidyverse)
library(data.table)

# Permute within treatment (meth/sham)
meta <- read.csv("/data/tunglab/dlin/baboon_mSTARR/window_counts/count_mat/mSTARR_info.csv", header=TRUE)
perm.mat.m0<-matrix("a", nrow=nrow(meta)/4,ncol=100) #ncol designates number of permutations 
perm.mat.s0<-matrix("a", nrow=nrow(meta)/4,ncol=100) #ncol designates number of permutations 

perm.mat.m1 <- data.frame(apply(perm.mat.m0, 2, function(x) x=sample(c("DNA","RNA"), nrow(meta)/4, replace=TRUE)))
perm.mat.s1 <- data.frame(apply(perm.mat.s0, 2, function(x) x=sample(c("DNA","RNA"), nrow(meta)/4, replace=TRUE)))

perm.mat.m2 <- apply(perm.mat.m1,2,function(c) ifelse(c=="DNA", "RNA", "DNA")) 
perm.mat.s2 <- apply(perm.mat.s1,2,function(c) ifelse(c=="DNA", "RNA", "DNA")) 
perm.mat <- rbind(perm.mat.m1, perm.mat.s1, perm.mat.m2, perm.mat.s2)
colnames(perm.mat) <- paste0("perm_",1:ncol(perm.mat))
# write.table(perm.mat,file='mstarr_msp1_type_perm_100_13Jul21.txt', quote=F, col.names=TRUE, row.names=F)

#####################################
# Run nested model to get permuted nested betas & the respective p values.
#####################################
### == Permute 100x by RNA (maintaining paired stucture) to identify sign (positive=more RNA than DNA) regulatory activity.
# Prepare voom normalized counts file
library(tidyverse)
library(data.table)
require(limma); require(edgeR); require(statmod)

meta <- read.csv("mSTARR_info.csv", header=TRUE)
rownames(meta) <- as.vector(meta$replicate) 
meta$Treatment <- relevel(meta$Treatment, ref="sham")
meta$Type <- relevel(meta$Type, ref="DNA")
v <- readRDS("msp1_nested_mod_voom_normalized_count_mat_18Jul21.RDS")
analyzable.wins <- readRDS("msp1_500_dnacov_dna_counts_filtered_13Jul21.RDS")  
#Run permutations:
perm.mat <-read.table('mstarr_msp1_type_perm_100_13Jul21.txt', header=TRUE)
for (i in 1:100) {
  Type <-perm.mat[,i] # permuted type (DNA/RNA) 
  design2<-model.matrix(~meta$Treatment + meta$Treatment:Type)
  design2
  vfit <-lmFit(v,design2)
  vfit <- eBayes(vfit)
  se.coef <- sqrt(vfit$s2.post) * vfit$stdev.unscaled
  vfit.df <- as.data.frame(cbind(vfit$p.value[,1:4],vfit$coefficient[,2:4],se.coef[,2:4]))
  rownames(v$E)<-analyzable.wins$V1
  vfit.df$site<-analyzable.wins$V1
  colnames(vfit.df) <- c("p_intercept","p_treatment", "p_RNA_sham", "p_RNA_meth", "b_treatment", "b_RNA_sham", "b_RNA_meth", "bse_treatment", "bse_RNA_sham", "bse_RNA_meth", "site")    
  saveRDS(vfit.df,paste0('msp1_nested_mod_rna_perm_',i,'_18Jul21.RDS'))
  print(i)
}

############################################
# Subsampling p-value matrix from each permfile to calculate q values. 
############################################
## Retain windows has positive permuted RNA beta, and take subset to have consistent number of regulatory windows in nested models from empirical data.
real<-readRDS("msp1_counts/msp1_nested_mod_out_df_18Jul21.RDS") ## vfit.df
pos.rna<-subset(real,b_RNA_meth>0 | b_RNA_sham>0); dim(pos.rna) # 5878

# Read in p-values from perms:
for (p in 1:100) {
  perm<-paste0('../msp1_perm/msp1_nested_mod_rna_perm_',p,'_18Jul21.RDS')
  permfile<-readRDS(perm)
  permfile<-subset(permfile, b_RNA_meth>0 | b_RNA_sham>0)   #If reducing to windows that had pos RNA beta in permuted data
  print(dim(permfile))
  permfile<-permfile[sample(nrow(permfile), nrow(pos.rna),replace=F), ]
  saveRDS(permfile,paste0('msp1_nested_mod_rna_perm_pval',p,'_18Jul21.RDS'))
  permfile<-permfile[,c("p_RNA_meth", "p_RNA_sham")]
  
  if(p=="1")
  {
    shuffled_pvals <-data.frame(x=permfile[,c("p_RNA_meth", "p_RNA_sham")])
    rownames(shuffled_pvals)=permfile$site
  } else {
    shuffled_pvals <- cbind(shuffled_pvals,permfile)
  }
  print(p)
}  

dim(shuffled_pvals) 
shuffled_pvals_RNA_meth<-shuffled_pvals[,seq(1,199,by=2)]
shuffled_pvals_RNA_sham<-shuffled_pvals[,seq(2,200,by=2)]
length(which(pos.rna$b_RNA_meth>0)) # 2430
length(which(pos.rna$b_RNA_sham>0)) # 4544

## calculate fdr (Joaquin's code)
res_full=perm.fdr(data.frame(pos.rna),shuffled_pvals_RNA_meth,"p_RNA_meth","RNA_meth")
res_full=res_full[order(res_full$site),]
head(res_full)
res_full=perm.fdr(data.frame(res_full),shuffled_pvals_RNA_sham,"p_RNA_sham","RNA_sham")
res_full=res_full[order(order(res_full$site)),]
head(res_full)
str(res_full)
save(res_full, shuffled_pvals_RNA_meth, shuffled_pvals_RNA_sham, file="msp1_nested_mod_fdr_shuffled_pvals_26Jul21.RData")

reg.win <- res_full %>% filter((b_RNA_meth>0 & fdr_RNA_meth < 0.1)| (b_RNA_sham>0 & fdr_RNA_sham < 0.1)) %>% pull(site)
write.table(reg.win, file="msp1_500_reg_wins_fdr10_26Jul21.tmp", quote=FALSE, col.names=FALSE, row.names=FALSE)


####################################################
# Create dataframe of permuted meth for calculating interaction fdr
####################################################
# Create permutation file, permuting meth/unmeth within dna samples and matching to rna samples
# use msp1 dataset as examples
meta <- read.csv("/mSTARR_info.csv", header=TRUE)
rownames(meta) <- as.vector(meta$replicate) 
meta$Treatment <- relevel(meta$Treatment, ref="sham")
meta$Type <- relevel(meta$Type, ref="DNA")
meta.dna <-subset(meta,Type=='DNA')
cols<-as.data.frame(matrix(nrow=nrow(meta),ncol=100)) ## designate perms

for (i in 1:100) {
  random<-as.vector(sample(meta.dna$Treatment,nrow(meta.dna),replace=F))
  cols[,i]<-c(random, random)
}
colnames(cols) <- paste0("perm_",1:ncol(cols))
write.table(cols,file='mstarr_msp1_treatment_perm_100_15Jul21.txt', quote=F, col.names=TRUE, row.names=F)

####################################################
# Run nested model to get permuted interaction (meth) betas & the respective p values.
####################################################
require(limma); require(edgeR); require(statmod) 

v <- readRDS("msp1_nested_mod_voom_normalized_count_mat_18Jul21.RDS")
analyzable.wins <- readRDS("msp1_counts/msp1_500_dnacov_dna_counts_filtered_13Jul21.RDS")  

# Run permutations:
perm.mat <-read.table('mstarr_msp1_treatment_perm_100_15Jul21.txt', header=TRUE)
for (i in 1:100) {
  Treatment <-perm.mat[,i] # permuted meth/unmeth
  Type <- meta$Type
  design2<-model.matrix(~Treatment + Treatment:Type)
  design2
  vfit <-lmFit(v,design2)
  vfit <- eBayes(vfit)
  se.coef <- sqrt(vfit$s2.post) * vfit$stdev.unscaled
  vfit.df <- as.data.frame(cbind(vfit$p.value[,1:4],vfit$coefficient[,2:4],se.coef[,2:4]))
  rownames(v$E)<-analyzable.wins$V1
  vfit.df$site<-analyzable.wins$V1
  colnames(vfit.df) <- c("p_intercept","p_treatment", "p_RNA_sham", "p_RNA_meth", "b_treatment", "b_RNA_sham", "b_RNA_meth", "bse_treatment", "bse_RNA_sham", "bse_RNA_meth", "site")    
  saveRDS(vfit.df,paste0('msp1_nested_mod_meth_perm_',i,'_18Jul21.RDS'))
  print(i)
}

####################################################
# Subsampling p-value matrix from each permfile to calculate fdr for interaction (meth)
####################################################
load("msp1_nested_mod_fdr_shuffled_pvals_18Jul21.RData")
msp1.nested.fdr <- res_full
msp1.nested.fdr.b <-subset(msp1.nested.fdr,b_RNA_meth>0 | b_RNA_sham>0)

for (p in 1:100) {
  permfile<-readRDS(paste0('msp1_nested_mod_meth_perm_',p,'_18Jul21.RDS'))
  permfile$b_interaction<-permfile$b_RNA_meth-permfile$b_RNA_sham
  permfile$bse_interaction<-sqrt(permfile$bse_RNA_meth^2 + permfile$bse_RNA_sham^2)
  permfile$t_interaction=permfile$b_interaction/permfile$bse_interaction
  permfile$p_interaction=2*pt(-abs(permfile$t_interaction),df=20) 

  # subsetting to the window subsamples that have positive regulatory windows from nested models.
  permfile<-permfile %>% filter(site %in% msp1.nested.fdr.b$site) 

  if(p=="1")
  {
    shuffled_pvals_rnaxmeth <-data.frame(x=permfile[,"p_interaction"])
    rownames(shuffled_pvals_rnaxmeth)=permfile$site
  } else {
    shuffled_pvals_rnaxmeth <- cbind(shuffled_pvals_rnaxmeth,x=permfile[,"p_interaction"])
  }
  print(p)
} 
colnames(shuffled_pvals_rnaxmeth) <- paste0("perm",1:100)


# calculate qvalues
res_full=perm.fdr(data.frame(msp1.nested.fdr),shuffled_pvals_rnaxmeth,"p_interaction","interaction")
res_full=res_full[order(rownames(res_full)),]

res_full.md <- readRDS("msp1_nested_mod_md_fdr_out_18Jul21.RDS") %>% select(site, fdr_interaction)
res_full <- merge(res_full, res_full.md, by="site", all.x=TRUE)
shuffled_pvals_rnaxmeth <- readRDS('mstarr_msp1_shuffled_pvals_rnaxmeth_18Jul21.RDS')
res_full.10<-subset(res_full, (b_RNA_meth>0 & fdr_RNA_meth <0.1 | (b_RNA_sham >0 & fdr_RNA_sham <0.1)))
res_full.10.int<-subset(res_full.10, fdr_interaction <0.1)
