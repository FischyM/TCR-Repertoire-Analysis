# https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0243927
# https://a-little-book-of-r-for-bioinformatics.readthedocs.io/en/latest/src/chapter5.html
# https://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf
# http://ape-package.ird.fr/

# https://fuzzyatelin.github.io/bioanth-stats/module-24/module-24.html

# BiocManager::install()

library("dplyr")
library("tidyr")
library("Biostrings")

# prepare input data
file_name <- "scRepertoire_outs/scaledClonotypeAbundance.csv"

seqs <- read.csv(file_name)
seqs <- seqs %>% dplyr::select(- "X")
seqs <- seqs %>% tidyr::separate("CTstrict", c('alpha regions', 'alpha chain', 'beta regions', 'beta chain'), sep='_')
seqs <- seqs %>% tidyr::unite(seq, c("alpha chain", "beta chain"), remove=F, sep='')
seqs <- seqs %>% tidyr::unite(names, c("values", "alpha regions", "beta regions"), remove=F, sep='-')

# make the fasta file
fasta_seqs <- Biostrings::DNAStringSet(seqs$seq)
fasta_seqs@ranges@NAMES <- make.unique(seqs$names, sep="|")

# save fasta file
Biostrings::writeXStringSet(fasta_seqs, filepath="testing/full-seqs.fa")


# save.image("C:/Users/Mathew/My Drive/LauraRogers/TCR_Repertoire/process_data_to_fasta.RData")

