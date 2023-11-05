#!/usr/bin/Rscript
#rarity_demo.r
#Nazia Pathan, 2023
rm(list = ls())

ptm<-proc.time()

#dependencies
if (!require("data.table", character.only = TRUE)) {
  install.packages("data.table")
}

if (!require("dplyr", character.only = TRUE)) {
  install.packages("dplyr")
}


# Load the packages
library(data.table)
library(dplyr)

#Input Arguments
args = (commandArgs(TRUE))
geno_file_pattern<-as.character(args[1])
pheno_file<-as.character(args[2])


#required files and directories
dir.create(file.path("RESULTS_H2"), showWarnings = FALSE)
outDir<-"RESULTS_H2"
output_file<-"BLOCK_HERITABILITY.txt"
outfile<-paste0("RESULTS_H2/",output_file)

load(pheno_file)
traits<-colnames(norm_pheno)
pheno_data<-norm_pheno

geno_files<-list.files(pattern=geno_file_pattern)


#main script to calculate block-wise heritability (adj_r2)

write.table(cbind("phenotype","block","N","N_RV","r2","adj_r2","block_r2_variance","block_adj_r2_variance"),
    outfile,append=F,quote=F,row.names=F,col.names=F,sep="\t")


for (i in 1:length(geno_files))
        {
            block=gsub(".RData","",geno_files[i])
            load(geno_files[i]) 
            geno_data_block = as.matrix(norm_geno)
            rm(norm_geno)

            A = t(geno_data_block) %*% geno_data_block 
            A<-round(A,digits=5) 
            C = t(geno_data_block) %*% pheno_data 
            inv_matrix<-solve(A)
            betas= inv_matrix %*% C 
            predicted = as.numeric()
            predicted= matrix(0, nrow = dim(geno_data_block)[1], ncol = dim(pheno_data)[2]) 
            colnames(predicted)<-colnames(pheno_data)

            for(trait_num in 1:length(traits))
            {
                predicted[,traits[trait_num]] = geno_data_block%*%betas[,traits[trait_num]]    
                R2= var(predicted[,traits[trait_num]])
                #n: nb_individuals, m: nb_SNVs
                n = dim(geno_data_block)[1]
                m = dim(geno_data_block)[2]
                adj_R2= 1 - ( 1 - R2 ) * ( n - 1 ) / ( n - m - 1 )
                block_VarR2<-((4*R2)*(1-R2)^2*(n-m-1)^2)/(((n^2)-1)*(n+3))
                block_Var_adj_R2<-( (n-1) / (n-m-1) )^2 * block_VarR2
                write.table(cbind(traits[[trait_num]],block,n,m,R2,adj_R2,block_VarR2,block_Var_adj_R2), 
                outfile,append=T,quote=F,row.names=F,col.names=F,sep="\t")
            }
        }


#compute total trait heritability of all blocks

all_h2<-fread(outfile,header=T, stringsAsFactors=F)
write.table(cbind("phenotype","N","N_RV","Heritability","LCL", "UCL"), 
    paste0(outDir,"/TOTAL_HERITABILITY.txt"),append=F,quote=F,row.names=F,col.names=F,sep="\t")


for (i in 1:length(traits)) 
  {
  trait_h2<-filter(all_h2,phenotype==traits[[i]])
  N_RV<-sum(trait_h2$N_RV)
  N<-mean(trait_h2$N)
  ADJ_R2<-as.numeric(sum(trait_h2$adj_r2))
  STD2 = sqrt(sum(trait_h2$block_adj_r2_variance)) 
  LCL_adj = ADJ_R2 - 1.96 * STD2 
  UCL_adj = ADJ_R2 + 1.96 * STD2 
  write.table((cbind(traits[i],N,N_RV,ADJ_R2,LCL_adj, UCL_adj)), 
    paste0(outDir,"/TOTAL_HERITABILITY.txt"),append=T,quote=F,row.names=F,col.names=F,sep="\t")
  }

proc.time()


