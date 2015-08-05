library(NCBI2R)
library(qvalue)
rdata <- list.files(pattern='*.Rdata')

for (i in rdata) { 
  load(i)
  metabolite <- strsplit(i,"_")[[1]][1]
  #all_sig <- length(pvals[pvals<1.939962e-09]) bonferonni cutoff for all 50 tests for 517061 snps
  q <- qvalue(pvals)
  q <- q$qvalues
  all_sig <- length(q[q<0.01])
  if (all_sig==0){
    out<-paste("No sig hits for",metabolite)
    print(out)
    print("___________________________")
    next
  }
  rs <- names(head(sort(pvals),n=all_sig))
  rs <- sapply(rs,function(x) strsplit(x,"\\.")[[1]][1])
  rs <- tolower(rs)
  info <- GetSNPInfo(rs)
  genes<-info$genesymbol
  temp<-''
  for (i in 1:length(genes)){temp<-c(temp,strsplit(genes[i],",")[[1]])}
  u_genes<-unique(temp)
  print(metabolite)
  print(sort(u_genes))
  entrezID <- info$locusID
  temp<-''
  for (i in 1:length(entrezID)){temp<-c(temp,strsplit(entrezID[i],",")[[1]])}
  u_entrezID<-unique(temp)
  print(u_entrezID)
  print("____________________________________")
}
