#Load qvalue package
suppressMessages(library(qvalue))

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript runFDR.R INPUT FDR OUTPUT"))
opt_input  = args[1];
opt_fdr    = as.numeric(args[2]);
opt_output = args[3];

#Verbose
cat("\nProcessing fastQTL output\n");
cat("  * Input  = [", opt_input, "]\n");
cat("  * FDR    = ", opt_fdr, "\n");
cat("  * Output = [", opt_output, "]\n");

#Read data
cat("\nRead Input data\n");
D = read.table(opt_input,hea=FALSE, stringsAsFactors=FALSE)
exon_offset = ifelse(ncol(D) == 19, 0, 2)
if (exon_offset == 2) cat("  * Gene level correction detected\n")
MASK=!is.na(D[,18+exon_offset])
Dnas=D[!MASK,]
D = D[MASK,]
cat("  * Number of molecular phenotypes =" , nrow(D), "\n")
cat("  * Number of NA lines =" , nrow(Dnas), "\n")
cat("  * Correlation between Beta approx. and Empirical p-values =", round(cor(D[, 18+exon_offset], D[, 19+exon_offset]), 4), "\n")

#Run qvalue on pvalues for best signals
cat("\nProcess Input data with Qvalue\n");
MASK=!is.na(D[,18+exon_offset])
Q = qvalue(D[MASK,19+exon_offset]);
D$qval = NA;
D$qval[MASK] = Q$qvalue;
cat("  * Proportion of significant phenotypes =" , round((1 - Q$pi0) * 100, 2), "%\n")

#Determine significance threshold
cat("\nDetermine significance thresholds\n");
set0 = D[which(D$qval <= opt_fdr),]
set1 = D[which(D$qval > opt_fdr),]
pthreshold = (sort(set1[,19+exon_offset])[1] - sort(-1.0 * set0[,19+exon_offset])[1]) / 2
cat("  * Corrected p-value threshold = ", pthreshold, "\n")
pval0 = qbeta(pthreshold, D[,14+exon_offset], D[,15+exon_offset], ncp = 0, lower.tail = TRUE, log.p = FALSE)
test0 = qf(pval0, 1, D[,13+exon_offset], ncp = 0, lower.tail = FALSE, log.p = FALSE)
corr0 = sqrt(test0 / (D[,13+exon_offset] + test0))
test1 = D[,12+exon_offset] * corr0 * corr0 / (1 - corr0 * corr0)
pval1 = pf(test1, 1, D[,12+exon_offset], ncp = 0, lower.tail = FALSE, log.p = FALSE)
cat("  * pval0 = ", mean(pval0), " +/- ", sd(pval0), "\n")
cat("  * test0 = ", mean(test0), " +/- ", sd(test0), "\n")
cat("  * corr0 = ", mean(corr0), " +/- ", sd(corr0), "\n")
cat("  * test1 = ", mean(test1), " +/- ", sd(test1), "\n")
cat("  * pval1 = ", mean(pval1), " +/- ", sd(pval1), "\n")
D$nthresholds = pval1

#Write significant hits
fout1=paste(opt_output, "significant.txt", sep=".")
cat("\nWrite significant hits in [", fout1, "]\n");
write.table(D[D$qval <= opt_fdr,], fout1, quote=FALSE, row.names=FALSE, col.names=FALSE)

#Write thresholds
fout2=paste(opt_output, "thresholds.txt", sep=".")
cat("\nWrite nominal thresholds in [", fout2, "]\n");
D1=D[, c(1, 21+exon_offset)]
D2=Dnas[, c(1,18)]
names(D2)=names(D1)
D3=rbind(D1, D2)  
write.table(D3, fout2, quote=FALSE, row.names=FALSE, col.names=FALSE)
