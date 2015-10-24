#/usr/bin/env Rscript

# what metabolite to analyze?
args <- commandArgs(TRUE)
met <- args[1]


# packages
library(car)
library(dglm)  
library(geoR)

# from Ronnegard and Valdar, Genetics 2011
dglm.Pvalues <- function(dglm.fit){
   P.disp = anova.dglm(dglm.fit)$Adj.P[2]
   P.mean = summary(dglm.fit)$coef[2,4]
   list(P.mean=P.mean, P.disp=P.disp)
}


load("cov_AllMetablites_AllGenotypes012_TSS.Rdata")
# Calculate boxcox lambda
bc_lambda <- boxcoxfit(cov_met_geno[,met][complete.cases(cov_met_geno[,met])])$lambda


snps <- names(cov_met_geno[grepl("RS",names(cov_met_geno))])
snps <- snps[1:100]

cov <- c("Ht","Wt","Sex","Sibcode")


run_dglm <- function(geno) {

  if (length(table(geno) > 0) & all(data.frame(table(geno))$Freq>1)) {  

    fitting_data <- cbind(cov_met_geno[,met],geno,cov_met_geno[,cov])
    fitting_data <- fitting_data[complete.cases(fitting_data),]
    fitting_data <<- fitting_data # dglm needs this object to be global
    d.fit <- dglm( formula =
                 bcPower(fitting_data[,1], bc_lambda ) ~ fitting_data[,2] +
                      (fitting_data[,"Ht"] +
                      fitting_data[,"Wt"] +
                      fitting_data[,"Sex"] +
                      fitting_data[,"Sibcode"]),
                      dformula = ~ fitting_data[,2] )
    return(dglm.Pvalues(d.fit))
	}
}


pvalues <- unlist(apply(cov_met_geno[,snps],2, function(x) run_dglm(x)))





output_name<- paste(met, "dglm_Pvalues.Rdata",sep="_")
save(pvalues,file=output_name)
