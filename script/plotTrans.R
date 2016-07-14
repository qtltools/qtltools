#Read command line arguments
args <- commandArgs(trailingOnly = TRUE)
try(if(length(args) != 5) stop("Incorrect number of arguments, usage> Rscript transQQplot.R qqplot.pdf nominal.hits.txt.gz nominal.bins.txt.gz permutation.hits.txt.gz permutation.bins.txt.gz"))

#Read data
Nh=read.table(args[2], head=FALSE, stringsAsFactors=FALSE)
Nb=read.table(args[3], head=FALSE, stringsAsFactors=FALSE)
Ph=read.table(args[4], head=FALSE, stringsAsFactors=FALSE)
Pb=read.table(args[5], head=FALSE, stringsAsFactors=FALSE)

#Sort best hits
Nh=Nh[order(-Nh$V7), ]
Ph=Ph[order(-Ph$V7), ]

#Counts
n_bins = nrow(Nb)
n_tests = nrow(Nh) + sum(Nb$V6)

#Cumulative sums
Nb$cs0=c(0, cumsum(Nb$V6[1:(n_bins-1)]))
Nb$cs1=cumsum(Nb$V6)
Pb$cs0=c(0, cumsum(Pb$V6[1:(n_bins-1)]))
Pb$cs1=cumsum(Pb$V6)
 
#Get the null p-value of a rank
pvalue_estimate <- function(idx, Ph, Pb) {
	pvalue = -1.0;
	if (idx > Pb$cs1[n_bins]) {
		pvalue = Ph$V7[idx - Pb$cs1[n_bins]]
	} else {
		nidx = which(Pb$cs0 < idx & idx <= Pb$cs1)

		cs0 = Pb$cs0[nidx]
		cs1 = Pb$cs1[nidx]
		pv0 = Pb$V4[nidx]
		pv1 = Pb$V5[nidx]

		pvalue = pv0 - (pv0 - pv1) * (idx - cs0) / (cs1 - cs0)
	}
	return (pvalue);
}

#Build the QQplot data
MP = matrix(0, 2 * rep (n_bins + nrow(Nh)), ncol=2)
for (b in 1:n_bins) {
	MP[b, 1] = (Nb$V4[b] + Nb$V5[b]) / 2
	MP[b, 2] = pvalue_estimate(round((Nb$cs0[b] + Nb$cs1[b]) / 2), Ph, Pb)
}
for (h in 1:nrow(Nh)) {
	MP[h+n_bins, 1] = Nh$V7[h]
	MP[h+n_bins, 2] = pvalue_estimate(sum(Nb$V6) + h, Ph, Pb)
}

#Plot the QQplot
pdf(args[1], 5, 5)
plot(-log10(MP[, 2]), -log10(MP[, 1]), xlab="-log10(permutation P-values)", ylab="-log10(nominal P-values)", main="QQplot")
abline(0, 1, col="red")
dev.off()

