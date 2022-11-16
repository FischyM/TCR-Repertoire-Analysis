# https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html

# BiocManager::install()

library("Biostrings")
library("stringr")
library("ape")
library("phangorn")

# read multiple aligned sequence
seqs_msa <- Biostrings::readDNAMultipleAlignment("testing/full-seqs.clustalo.aln.fa", format="fasta")

# convert fasta MSA file to Biostrings file type
alignTCRs <- Biostrings::DNAStringSet(seqs_msa)

# get unaligned strings for distance/identity matrix
unalignTCRs <- stringr::str_replace_all(as.character(alignTCRs), "-", "")
unalignTCRs <- Biostrings::DNAStringSet(unalignTCRs)

# calculate levenshtein distance matrix
# aln.distMat <- ape::dist.dna(as.DNAbin(alignTCRs), pairwise.deletion=T)
unaln.distMat <- Biostrings::stringDist(unalignTCRs)

# create unrooted nj tree
# initialTree1 <- ape::bionj(unaln.distMat)
initialTree1 <- ape::fastme.ols(unaln.distMat, nni=T)

# set negative branch length to 0
initialTree1$edge.length[initialTree1$edge.length < 0] <- 0
ape::write.tree(initialTree1, file='testing/full-seqs.clustalo.aln-fastme.ols.nwk', tree.names=T)

# plot tree
ape::plot.phylo(initialTree1, type="cladogram")
ape::plot.phylo(initialTree1, type="tidy")



initialTree2 <- ape::nj(unaln.distMat)
ape::write.tree(initialTree2, file='testing/full-seqs.unaln-nj.nwk', tree.names=T)
ape::plot.phylo(initialTree2, type="cladogram")
ape::plot.phylo(initialTree2, type="tidy")





save.image("C:/Users/Mathew/My Drive/LauraRogers/TCR_Repertoire/tree_construction.RData")

load("C:/Users/Mathew/My Drive/LauraRogers/TCR_Repertoire/tree_construction.RData")





# maximum likelihood
dna2 <- phangorn::as.phyDat(alignTCRs)
phangorn::pml(initialTree1, dna2, k=4)

# if loglikelihood is NaN, likely there are missing seq data
table(as.character(dna2))
na.posi <- which(apply(as.character(seqs_msa),2, function(e) any(!e %in% c("a","t","g","c"))))

# plot freq of non ATCG bases over position in seqs
temp <- apply(as.character(seqs_msa),2, function(e) sum(!e %in% c("a","t","g","c")))
plot(temp, type="l", col="blue", xlab="Position in HA segment", ylab="Number of NAs")

# exclude missing data
dna3 <- seqs_msa[,-na.posi]
table(as.character(dna3))

# convert back to phyDat
dna4 <- phangorn::as.phyDat(dna3)

# remake NJ tree and use pml to calc likelihood
tre.ini <- ape::bionj(dist.dna(dna3, model="TN93"))
fit.ini <- phangorn::pml(tre.ini, dna4, k=4)
fit.ini

# optimize tree
fit <- phangorn::optim.pml(fit.ini, optNni=TRUE, optBf=TRUE, optQ=TRUE, optGamma=TRUE)
fit

# compare optimized tree to the nj tree using an anova
anova(fit.ini, fit)

# lower AIC is better
AIC(fit.ini)
AIC(fit)

# plot tree again
tre4 <- root(fit$tree,1)
tre4 <- ladderize(tre4)
plot(tre4, show.tip=FALSE, edge.width=2)
title("Maximum-likelihood tree")
tiplabels(annot$year, bg=transp(num2col(annot$year, col.pal=myPal),.7), cex=.5, fg="transparent")
axisPhylo()
temp <- pretty(1993:2008, 5)
legend("topright", fill=transp(num2col(temp, col.pal=myPal),.7), leg=temp, ncol=2)



# analysis using Jalview
# http://www.jalview.org/sites/jalview.org/files/TheJalviewTutorial-screen.pdf
