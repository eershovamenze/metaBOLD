library(ape)
BiocManager::install("DECIPHER")
library(DECIPHER)
setwd("~/Documents/R projects/metaBOLD/src/data_preprocessing")

motu_init <- read.csv('nrws_motu.csv')

library(Biostrings)
seqs <- DNAStringSet(motu_init$sequence)
names(seqs) <- motu_init$id
alignment <- AlignSeqs(seqs)
summary(alignment)
alignment_dnabin <- as.DNAbin(alignment)
dist_matrix <- dist.dna(alignment_dnabin, model = "K80")
dist_mat <- as.matrix(dist_matrix)
close_pairs <- which(dist_mat <= 0.03, arr.ind = TRUE)
close_pairs <- close_pairs[close_pairs[,1] != close_pairs[,2], ]
motu_ids <- rownames(dist_mat)
result <- data.frame(
  MOTU1 = motu_ids[close_pairs[,1]],
  MOTU2 = motu_ids[close_pairs[,2]],
  Distance = dist_mat[close_pairs]
)
bad_indices <- which(!is.finite(dist_mat), arr.ind = TRUE)

# Get unique problematic rows/columns
bad_rows <- unique(bad_indices[,1])

# Remove these from dist_mat
dist_mat_clean <- dist_mat[-bad_rows, -bad_rows]

heatmap(dist_mat_clean, 
        main = "Pairwise Distance Heatmap", 
        xlab = "MOTUs", 
        ylab = "MOTUs",
        col = hcl.colors(100, "YlOrRd", rev = TRUE))

seq_list <- setNames(as.list(motu_init$sequence), motu_init$id)
seq_list_split <- lapply(seq_list, function(x) unlist(strsplit(x, "")))
# Convert to DNAbin format:
# Here we assume all sequences are the same length and aligned.
# We'll first bind them by rows into a matrix.
dna_matrix <- as.matrix(as.DNAbin(seq_list_split))

# Now dna_matrix is a DNAbin object containing all sequences.
# You can calculate the distance matrix:
dist_matrix <- dist.dna(dna_matrix, model = "K80")