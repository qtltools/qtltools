#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 3) stop("Incorrect number of arguments, usage> Rscript runFDR_ftrans.R nominal.hits.txt.gz permutation.hits.txt.gz output.txt"))

cat("\nProcessing QTLtools full trans output\n");
cat("  * File hits nominal = [", args[1], "]\n");
cat("  * File hits permute = [", args[2], "]\n");
cat("  * Output            = [", args[3], "]\n");

#Read data
Nh=read.table(args[1], head=FALSE, stringsAsFactors=FALSE)
Ph=read.table(args[2], head=FALSE, stringsAsFactors=FALSE)

#Sort best hits
Nh=Nh[order(Nh$V7), ]
Ph=Ph[order(Ph$V7), ]

#Estimate FDR
Nh$fdr=1.0
for (h in 1:nrow(Nh)) {
	Nh$fdr[h] = sum(Ph$V7 <= Nh$V7[h]) / h
}

#OUTPUT
write.table(Nh, args[3], quote=FALSE, row.names=FALSE, col.names=FALSE)
