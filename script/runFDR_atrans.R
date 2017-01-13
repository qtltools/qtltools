#Load qvalue package
suppressMessages(library(qvalue))

#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 4) stop("Incorrect number of arguments, usage> Rscript runFDR_atrans.R adjusted.best.txt adjusted.hits.txt FDR output.txt"))

cat("\nProcessing QTLtools approximate trans output\n");
cat("  * File best  = [", args[1], "]\n");
cat("  * File hits  = [", args[2], "]\n");
cat("  * FDR        = [", args[3], "]\n");
cat("  * Output     = [", args[4], "]\n");


B = read.table(args[1], head=FALSE, stringsAsFactors=FALSE)
H = read.table(args[2], head=FALSE, stringsAsFactors=FALSE)
FDR = as.numeric(args[3])
B$qval = qvalue(B$V2)$qval
threshold = min(B$V2[which(B$qval > FDR)])
cat("  * Threshold of significance for adjusted P-values =" , threshold, "\n")

cat("\nFiltering hits and output results\n");
S = H[which(H$V8 <= threshold), ]
cat("  * " , nrow(S) , " are significante out of ", nrow(H), "\n")
write.table(S, args[4], quote=FALSE, row.names=FALSE, col.names=FALSE)
