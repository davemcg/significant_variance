#!/usr/bin/env Rscript

library(reshape2)
library(data.table)
library(ggplot2)
samples<-read.table("TRINITY_2232_757577_w_sex.tfam")
first<-c("Chr","SNP","0","Position")
titles<-c(first,as.character(samples$V1))
rdata <- list.files(pattern='*.Rdata')
setwd("/data/mcgaugheyd/projects/brody/trinity_variance")
pheno<-read.csv("TSS_GWAS_Nov10.csv",header=T)
row.names(pheno)<-pheno$UniqueID
samples<-read.table("TRINITY_2232_757577_w_sex.tfam")
pheno<-pheno[intersect(row.names(pheno),samples$V1),]
trinity_tped <- data.frame(fread('TRINITY_2232_757577_w_sex.collapsed.tped'))
colnames(trinity_tped)<-titles
row.names(trinity_tped)<-trinity_tped$SNP

for (i in rdata) {
    metabolite <- strsplit(i,"_")[[1]][1]
    select_phenotype<-pheno[,c("UniqueID",metabolite)]
    select_phenotype[,2] <- as.numeric(as.character(select_phenotype[,2]))
    select_phenotype<-select_phenotype[complete.cases(select_phenotype[,2]),]
    row.names(select_phenotype) <- select_phenotype$UniqueID
    load(i)
    top<-sapply(strsplit(names(head(sort(pvals),n=9)),"\\."),'[',1)
    hits<-merge(select_phenotype,data.frame(t(trinity_tped[top,5:ncol(trinity_tped)])),by="row.names")
    hits<-(melt(hits,id.vars=c("Row.names","UniqueID",metabolite)))
    hits<-subset(hits,value!='00')
    counts<-aggregate(hits[,3]~variable+value,data=hits,FUN=length)
    rs_levels<-factor(hits$variable)
	counts<-with(counts,counts[order(rs_levels),])
	hit_pvals<-data.frame(cbind(head(sort(pvals),n=9),top))
    colnames(hit_pvals) <- c("value","variable")
    # set y dimension for p values
    temp<-aggregate(hits[,3]~variable+value,data=hits,FUN=max)
    yPmax<-min(c(max(temp[,3])*1.2,max(temp[,3])+10))
    #yPmin<-min(c(-min(temp[,3])*1.2,min(temp[,3])-10))
    export<-paste(metabolite,"_rawB.pdf",sep="")
    #setwd("~/Desktop")
    pdf(export)
    print(ggplot(hits,aes(x=value,y=hits[,metabolite])) + facet_wrap(~variable,scale="free") + geom_violin() + geom_boxplot(width=0.1) + geom_text(data=hit_pvals,aes(x=0,y=yPmax,label=as.character(value),colour='red'),parse=TRUE,hjust=0) + geom_text(data=counts,aes(x=value,y=-0.1,label=as.character(counts[,3]),parse=TRUE,size=4,colour='red')) + ylab(metabolite) + theme(legend.position="none"))
    dev.off()
}
