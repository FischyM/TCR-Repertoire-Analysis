# Run multiple sequence alignment on commandline after downloading clustalo
#
# cd /mnt/c/Users/Mathew/My\ Drive/LauraRogers/TCR_Repertoire
# ~/./clustalo-1.2.4 -i testing/full-seqs.fa -o testing/full-seqs.clustalo.aln.fa --threads 8 -v --distmat-out testing/full-seqs.clustalo.distmat --percent-id --full --force
#
# Use Jalview to view the alignments


library("Biostrings")
library("msa")

fasta_file <- "testing/full-seqs.fa"

# load fasta file
unaln_seqs <- Biostrings::readDNAStringSet(fasta_file)

# multiple alignment sequences
aln_seqs <- msa::msa(unaln_seqs, method="ClustalOmega")

# save msa file
Biostrings::writeXStringSet(Biostrings::DNAStringSet(aln_seqs), filepath="testing/full-seqs.clustalo.aln.fa")


# save.image("C:/Users/Mathew/My Drive/LauraRogers/TCR_Repertoire/run_msa.RData")


