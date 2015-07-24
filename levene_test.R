#/usr/bin/env Rscript

# what metabolite to analyze?
args <- commandArgs(TRUE)
metabolite <- args[1]
  
# returns p value for differences in variance between groups
# zero genotype samples not used in levene.test
# only SNPs with >100 people in each group are tested
# also there need to be more than one group
process_levene_a <- function(x){
  data <- x
  data <- gsub("00",NA,data)
  data <- data[!is.na(data)]
  temp <- merge(select_phenotype,data.frame(data),by="row.names")
  if(length(table(temp[,ncol(temp)]))>1 & all(table(temp[,ncol(temp)])>50)){
    stats <- tidy(levene.test(temp$MMA,temp[,ncol(temp)]))
    return(stats)
  } else{
  return()
  }
}

# packages
library(broom)
library(reshape2)
library(data.table)
library(lawstat)
library(qvalue)
library(parallel)

# load, clean, and prepare data 
trinity_tped <- data.frame(fread('~/Desktop/TRINITY_2232_757577_w_sex.collapsed.tped'))
pheno<-read.csv("~/Desktop/TSS_GWAS_Nov10.csv",header=T)
row.names(pheno)<-pheno$UniqueID
samples<-read.table("~/Desktop/TRINITY_2232_757577_w_sex.tfam")
pheno<-pheno[intersect(row.names(pheno),samples$V1),]
select_phenotype<-pheno[,c("UniqueID",metabolite)]
select_phenotype[,2] <- as.numeric(as.character(select_phenotype[,2]))
select_phenotype<-select_phenotype[complete.cases(select_phenotype[,2]),]
first<-c("Chr","SNP","0","Position")
titles<-c(first,as.character(samples$V1))
colnames(trinity_tped)<-titles
row.names(trinity_tped)<-trinity_tped$SNP

# initialize multi-core and analyze
c1<-makeCluster(6)
system.time(raw_results<-unlist(parApply(c1,trinity_tped[1:10000,row.names(select_phenotype)],1,function(x) process_levene_a(x))))
head(raw_results)

# pull out p values, display, and save
pvals <- raw_results[grep("value",names(raw_results))]
head(pvals)
print(summary(pvals))
print(qsummary(qvalue(pvals)))
name0<-paste(metabolite,"_pvals.Rdata")
save(pvals,file=name0)

# ggplot the top 9 hits by raw value
top<-sapply(strsplit(names(head(sort(pvals),n=9)),"\\."),'[',1)
hits<-merge(select_phenotype,data.frame(t(trinity_tped[top,5:ncol(trinity_tped)])),by="row.names")
hits<-(melt(hits,id.vars=c("Row.names","UniqueID",metabolite)))
hits<-subset(hits,value!='00')
name1<-paste(metabolite,"_raw.pdf",sep="")
pdf(name1)
print(ggplot(hits,aes(x=value,y=as.numeric(MMA))) + facet_wrap(~variable,scale="free") + geom_violin() + geom_boxplot(width=0.1))
dev.off()

# ggplot the top 9 hits by abs(distance) from median
median_dist <- function(x) {
  data <- x
  data <- gsub("00",NA,data)
  data <- data[!is.na(data)]
  temp <- merge(select_phenotype,data.frame(data),by="row.names")
  temp <- sapply(split(temp[,metabolite],temp[,ncol(temp)]),function(x) abs(x-median(x)))
  return(temp)
}
med_dist<-apply(trinity_tped[top,5:ncol(trinity_tped)],1,function(x) median_dist(x))
med_dist<-melt(med_dist)
name2<-paste(metabolite,"_distance.pdf",sep="")
pdf(name2)
print(ggplot(med_dist,aes(x=abs(value),colour=L2)) + geom_density() + theme_bw() + facet_wrap(~L1))
dev.off()

