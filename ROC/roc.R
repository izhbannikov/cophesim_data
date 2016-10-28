# cophesim tests #

library(GenABEL)

calculate_roc <- function(df, cost_of_fp=NULL, cost_of_fn=NULL, n=100, column="P") {
  tpr <- function(df, threshold) {
    sum(df[[column]] <= threshold & df$true == 1, na.rm=T) / sum(df$true == 1, na.rm = T)
  }
  
  fpr <- function(df, threshold) {
    sum(df[[column]] <= threshold & df$true == 0, na.rm = T) / sum(df$true == 0, na.rm=T)
  }
  
  roc <- data.frame(threshold = seq(0,1,length.out=n), tpr=NA, fpr=NA)
  roc$tpr <- sapply(roc$threshold, function(th) tpr(df, th))
  roc$fpr <- sapply(roc$threshold, function(th) fpr(df, th))
  
  return(roc)
}

## ROC curves ##
N.SNP <- 1000 # Total # of snps
N.SNP.C <- 100 # Causal snps
### Simulate some genotypic data ###
system("plink --simulate-ncases 5000 --simulate-ncontrols 5000 --simulate /Users/ilya/Projects/cophesim_test/ROC/wgas.sim --out /Users/ilya/Projects/cophesim_test/ROC/sim.plink --make-bed")
system("plink --bfile /Users/ilya/Projects/cophesim_test/ROC/sim.plink --recode --out /Users/ilya/Projects/cophesim_test/ROC/sim.plink")

### Simulate phenotypes for this data ###
#### Make file effects.txt ####
effects <- runif(n=N.SNP, min=-2, max=2)
sink("/Users/ilya/Projects/cophesim_test/ROC/effects.txt")
for(i in seq(0, (N.SNP.C-1))) {
  cat(paste(i, ":", effects[i+1], "\n", sep=""))
}
sink()

#### Simulate phenotypes with cophesim ####
system("python /Users/ilya/Projects/cophesim/cophesim.py -i /Users/ilya/Projects/cophesim_test/ROC/sim.plink -o /Users/ilya/Projects/cophesim_test/ROC/testout -itype plink -otype plink -c -ce /Users/ilya/Projects/cophesim_test/ROC/effects.txt -s -gomp")

### Association tests with Plink ###
#### Logistic ####
system("plink --file /Users/ilya/Projects/cophesim_test/ROC/testout_pheno_bin.txt --1 --logistic --noweb --allow-no-sex --out /Users/ilya/Projects/cophesim_test/ROC/out")
#### Linear ####
system("plink --file /Users/ilya/Projects/cophesim_test/ROC/testout_pheno_cont.txt --linear --noweb --allow-no-sex --out /Users/ilya/Projects/cophesim_test/ROC/out")
#### Survival with GenABEL####
convert.snp.ped(pedfile = "/Users/ilya/Projects/cophesim_test/ROC/testout_pheno_surv.txt.ped",
                mapfile = "/Users/ilya/Projects/cophesim_test/ROC/testout_pheno_surv.txt.map",
                outfile = "/Users/ilya/Projects/cophesim_test/ROC/testout_pheno_surv.txt.dat")

d <- load.gwaa.data(phenofile = "/Users/ilya/Projects/cophesim_test/ROC/testout_pheno_surv.txt",
                    genofile = "/Users/ilya/Projects/cophesim_test/ROC/testout_pheno_surv.txt.dat")

assoc.res.survival <- mlreg(GASurv(age,case)~1,d)



### Handling results ###
assoc.res.logistic <- read.table("/Users/ilya/Projects/cophesim_test/ROC/out.assoc.logistic", header = T)
assoc.res.linear <- read.table("/Users/ilya/Projects/cophesim_test/ROC/out.assoc.linear", header = T)
assoc.res.survival <- assoc.res.survival@results

### Adding classification column ###
#### Classified ###
assoc.res.logistic <- cbind(assoc.res.logistic, classified=ifelse(assoc.res.logistic$P <= 0.05, 1, 0))
assoc.res.linear <- cbind(assoc.res.linear, classified=ifelse(assoc.res.linear$P <= 0.05, 1, 0))
assoc.res.survival <- cbind(assoc.res.survival, classified=ifelse(assoc.res.survival$P1df <= 0.05, 1, 0))
#### Reality ####
cc <- replicate(0, n=dim(assoc.res.logistic)[1])
cc[1:N.SNP.C] <- 1
assoc.res.logistic <- cbind(assoc.res.logistic, true=cc)
assoc.res.linear <- cbind(assoc.res.linear, true=cc)
assoc.res.survival <- cbind(assoc.res.survival, true=cc[1:(N.SNP-1)])

### plotting ROC curves ###
roc.logistic <- data.frame(calculate_roc(assoc.res.logistic, n=N.SNP, column="P"))
roc.linear <- calculate_roc(assoc.res.linear, n=N.SNP)
roc.survival <- calculate_roc(assoc.res.survival, n=N.SNP-1, column="P1df")

plot(x=roc.logistic$fpr, y=roc.logistic$tpr, type = "l", main="ROC curve", ylab="TPR", xlab="FPR", col="red")
lines(x=roc.linear$fpr, y=roc.linear$tpr, type = "l", col="blue")
lines(x=roc.survival$fpr, y=roc.survival$tpr, type = "l", col="green")

#### Plotting with ggplot2 ####
library(ggplot2)

idx_threshold = which.min(abs(roc.logistic$threshold-0.05))
p_roc <- ggplot(roc.logistic,aes(x=fpr)) + 
  geom_line(aes(y=tpr, color="Dichotomous"), alpha=0.3) +
  geom_point(aes(y=tpr, color="Dichotomous"), size=1, alpha=0.5) +
  coord_fixed() +
  geom_line(aes(threshold,threshold), linetype="dashed") +
  geom_line(aes(x=fpr,y=tpr, color="Continous"), data = roc.linear, alpha=0.3) +
  geom_point(aes(y=tpr, color="Continous"), size=1, alpha=0.5, data = roc.linear) +
  geom_line(aes(x=fpr,y=tpr, color="Survival"), data = roc.survival, alpha=0.3) +
  geom_point(aes(y=tpr, color="Survival"), size=1, alpha=0.5, data = roc.survival) +
  labs(title = sprintf("ROC")) + xlab("FPR") + ylab("TPR") +
  scale_colour_manual(values=c("#D55E00", "#0072B2", "#009E73"), 
                      name="Trait type")
p_roc


