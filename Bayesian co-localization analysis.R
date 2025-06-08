library(tidyverse)
library(data.table)
library(dplyr)
library(coloc)
library(gwasglue)
library(locuscomparer)

setwd()

QTLdata <- fread("2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz",sep = "\t",header = TRUE)
QTLdata <- QTLdata %>% filter(GeneSymbol == "THOC5")
ref <- fread("2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz", sep = "\t", 
             header = TRUE)
ref<-select(ref,SNP,AlleleB_all)
QTLdata <- left_join(QTLdata, ref, by = "SNP")
QTLdata$MAF <- ifelse(QTLdata$AlleleB_all<0.5,QTLdata$AlleleB_all,1-QTLdata$AlleleB_all)
QTLdata$end <- QTLdata$GenePos
QTLdata <- select(QTLdata,Pvalue,SNP,NrSamples,MAF,GeneChr,GenePos,end)
colnames(QTLdata)<-c("pvalues","snp","N","MAF","chrom","start","end")

lead <- QTLdata %>% dplyr::arrange(pvalues)
leadchr <- lead$chrom[1]
leadstart <- as.numeric(lead$start[1])
leadend <- as.numeric(lead$end[1])
data1 <- QTLdata[QTLdata$chrom==leadchr,]
data1 <- data1[data1$start>leadstart-1000000 & data1$end<leadend+1000000,]

GWAS <- fread("finngen_R10_C_DIFFICILE_ENTEROCOLITIS.txt",sep = "\t",header = TRUE)
GWAS$N <- "409432"
GWAS$MAF <- ifelse(GWAS$af_alt<0.5,GWAS$af_alt,1-GWAS$af_alt)
GWAS$varbeta <- GWAS$sebeta^2
GWAS$s <- "0.00826511"

data2 <- select(GWAS,pval,rsids,N,MAF,beta,varbeta,s)
colnames(data2)<-c("pvalues","snp","N","MAF","beta","varbeta","s")

data1 <- subset(data1, !duplicated(snp))
data2 <- subset(data2, !duplicated(snp))
data1 <- data1 %>% na.omit() %>% filter(pvalues>0)
data2 <- data2 %>% na.omit() %>% filter(pvalues>0)

sameSNP <- intersect(data1$snp,data2$snp)
data1 <- data1[data1$snp %in% sameSNP, ] %>% dplyr::arrange(snp) %>% na.omit()
data2 <- data2[data2$snp %in% sameSNP, ] %>% dplyr::arrange(snp) %>% na.omit()

result <- coloc.abf(dataset1=list(pvalues=data1$pvalues, snp=data1$snp, 
                                  type="quant", N=max(data1$N), MAF=data1$MAF), 
                    dataset2=list(pvalues=data2$pvalues, snp=data2$snp,
                                  type="quant", N=max(data2$N), MAF=data2$MAF))

result <- coloc.abf(dataset1=list(pvalues=data1$pvalues, snp=data1$snp, 
                                  type="quant", N=max(data1$N), MAF=data1$MAF), 
                    dataset2=list(pvalues=data2$pvalues, snp=data2$snp, 
                                  type="cc", s=0.00826511, N=409432, MAF=data2$MAF, 
                                  beta=data2$beta,
                                  varbeta=data2$varbeta))

need_result <- result$results %>% dplyr::arrange(desc(SNP.PP.H4))

eqtl_fn <- data1[,c('snp','pvalues')] %>%
  dplyr::rename(rsid = snp, pval = pvalues)
gwas_fn <- data2[,c('snp','pvalues')] %>%
  dplyr::rename(rsid = snp, pval = pvalues)

print(locuscompare(in_fn2 = eqtl_fn, 
                   in_fn1 = gwas_fn, 
                   title2 = 'eQTL', 
                   title1 = 'GWAS'))
