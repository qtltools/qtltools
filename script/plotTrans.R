#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) < 3) stop("Incorrect number of arguments, usage> Rscript plotTrans.R output.pdf nominals.txt permutation1.txt ... permutationN.txt"))

N=read.table(args[2], head=FALSE, stringsAsFactors=FALSE)
N=N[order(N$V7), ]
N$logpva = -log10(N$V7)

P=read.table(args[3], head=FALSE, stringsAsFactors=FALSE)
P=P[order(P$V7), ]
P$logpva = -log10(P$V7)

pdf(args[1], 6, 6)
plot(0, 0, type="n", main="QQplot", xlab="Expected -log10(P-value)", ylab="Observed -log10(P-value)", xlim=c(min(P$logpva), max(P$logpva)+1), ylim=c(min(N$logpva), max(N$logpva)+1))
for (p in 3:length(args)) {
	P=read.table(args[p], head=FALSE, stringsAsFactors=FALSE)
	P=P[order(P$V7), ]
	P$logpva = -log10(P$V7)
	MIN = min(nrow(P), nrow(N))
	points(P$logpva[1:MIN], N$logpva[1:MIN], type="l", col=rgb(0,0,1, alpha=0.1))
}
abline(0, 1, col="red")
dev.off()



