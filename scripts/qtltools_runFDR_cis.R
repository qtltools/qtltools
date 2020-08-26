#!/usr/bin/env Rscript
#Load qvalue package
suppressMessages(library(qvalue))

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) != 3)  stop("Incorrect number of arguments\n usage\t qtltools_runFDR_cis.R INPUT FDR OUTPUT")
opt_input  <- args[1]
opt_fdr    <- as.numeric(args[2])
opt_output <- args[3]
opt_col    <- 20
opt_dof1   <- 12
opt_dof2   <- 13
opt_bml1   <- 14
opt_bml2   <- 15
temp <- read.table(opt_input,nrows = 1)
header <- FALSE
if (temp$V1 == "phe_id" || temp$V1 == "grp_id"){
  header  <- TRUE
  opt_col <- which(temp == "adj_beta_pval")
  if(length(opt_col) != 1)  stop("Header present but cannot find column adj_beta_pval")
  opt_dof1 <- which(temp == "dof1")
  if(length(opt_dof1) != 1)  stop("Header present but cannot find column dof1")
  opt_dof2 <- which(temp == "dof2")
  if(length(opt_dof2) != 1)  stop("Header present but cannot find column dof2")
  opt_bml1 <- which(temp == "bml1")
  if(length(opt_bml1) != 1)  stop("Header present but cannot find column bml1")
  opt_bml2 <- which(temp == "bml2")
  if(length(opt_bml2) != 1)  stop("Header present but cannot find column bml2")
}

#Verbose
cat("\nProcessing QTLtools output\n");
cat("  * Input  = [", opt_input, "]\n");
cat("  * FDR    = ", opt_fdr, "\n");
cat("  * Output = [", opt_output, "]\n");

#Read data
cat("\nRead Input data\n");
D = read.table(opt_input,hea=header, stringsAsFactors=FALSE)
if (!header){
  if(ncol(D) == 21){
    opt_col <- 21
    if (!all(D[,16] >= 0 & D[,16] <= 1,na.rm = T)){
      cat("Assuming --grp-mean\n")
      opt_dof1 <- 13
      opt_dof2 <- 14
      opt_bml1 <- 15
      opt_bml2 <- 16
    }else{
      cat("Assuming --std-err\n")
    }
  }else if (ncol(D) == 22){
    opt_col <- 22
    if(!all(D[,17] >= 0 & D[,17] <= 1,na.rm = T)){
      cat("Assuming --grp-{pca,best}\n")
      opt_dof1 <- 14
      opt_dof2 <- 15
      opt_bml1 <- 16
      opt_bml2 <- 17
    }else{
      cat("Assuming --grp-mean and --std-err\n")
      opt_dof1 <- 13
      opt_dof2 <- 14
      opt_bml1 <- 15
      opt_bml2 <- 16
    }
  }else if (ncol(temp) == 23){
    cat("Assuming --grp-{pca,best} and --std-err\n")
    opt_col    <- 23
    opt_dof1   <- 14
    opt_dof2   <- 15
    opt_bml1   <- 16
    opt_bml2   <- 17
  }else if (ncol(temp) != 20){
    stop("Unknown file type")
  }
}

cat("  * Selected columns for dof1 dof2 bml1 bml2 adj_beta_pval\n")
print(head(D[c(opt_dof1,opt_dof2,opt_bml1,opt_bml2,opt_col)]))

cat("  * Number of molecular phenotypes =" , nrow(D), "\n")
cat("  * Number of NA lines =" , sum(is.na(D[,opt_col - 1])), "\n")
cat("  * Correlation between Beta approx. and Empirical p-values =", round(cor(D[, opt_col-1], D[, opt_col],use = "pairwise.complete.obs"), 6), "\n")

#Run qvalue on pvalues for best signals
cat("\nProcess Input data with Qvalue\n")
Q <- qvalue(D[,opt_col])
D$qval <- Q$qvalues
cat("  * Proportion of significant phenotypes =" , round((1 - Q$pi0) * 100, 2), "%\n")

#Determine significance threshold
cat("\nDetermine significance thresholds\n");
set0 = D[which(D$qval <= opt_fdr),]
set1 = D[which(D$qval > opt_fdr),]
pthreshold = (sort(set1[,opt_col])[1] - sort(-1.0 * set0[,opt_col])[1]) / 2
cat("  * Corrected p-value threshold = ", pthreshold, "\n")
pval0 = qbeta(pthreshold, D[,opt_bml1], D[,opt_bml2], ncp = 0, lower.tail = TRUE, log.p = FALSE)
test0 = qf(pval0, 1, D[,opt_dof2], ncp = 0, lower.tail = FALSE, log.p = FALSE)
corr0 = sqrt(test0 / (D[,opt_dof2] + test0))
test1 = D[,opt_dof1] * corr0 * corr0 / (1 - corr0 * corr0)
pval1 = pf(test1, 1, D[,opt_dof1], ncp = 0, lower.tail = FALSE, log.p = FALSE)
cat("  * pval0 = ", mean(pval0,na.rm = T), " +/- ", sd(pval0,na.rm = T), "\n")
cat("  * test0 = ", mean(test0,na.rm = T), " +/- ", sd(test0,na.rm = T), "\n")
cat("  * corr0 = ", mean(corr0,na.rm = T), " +/- ", sd(corr0,na.rm = T), "\n")
cat("  * test1 = ", mean(test1,na.rm = T), " +/- ", sd(test1,na.rm = T), "\n")
cat("  * pval1 = ", mean(pval1,na.rm = T), " +/- ", sd(pval1,na.rm = T), "\n")
D$nthresholds = pval1

#Write significant hits
fout1=paste(opt_output, "significant.txt", sep=".")
cat("\nWrite significant hits in [", fout1, "]\n");
write.table(D[which(D$qval <= opt_fdr),], fout1, quote=FALSE, row.names=FALSE, col.names=FALSE)

#Write thresholds
fout2=paste(opt_output, "thresholds.txt", sep=".")
cat("\nWrite nominal thresholds in [", fout2, "]\n");
write.table(D[,c(1,ncol(D))], fout2, quote=FALSE, row.names=FALSE, col.names=FALSE)
