library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)

args = commandArgs(trailingOnly=TRUE)

id <- args[1]

file_fasta <- args[2]

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

# --- MASKING PATCH START ---
# Load metadata (MANDATORY)
boundaries_path <- args[3]
ortholog_table_path <- args[4]

source("./PHACE_Codes/load_metadata.R")
protein_boundaries <- load_protein_boundaries(boundaries_path)
ortholog_matrix <- load_ortholog_table(ortholog_table_path, rownames(msa))
cat("Metadata loaded. Artificial gaps will be encoded as missing data.\n")
# --- MASKING PATCH END ---

for (i in 1:length(msa[1,])){
  # --- MASKING PATCH START ---
  # Determine which leaves are missing this protein
  mask_leaves <- get_masked_leaves_for_position(i, protein_boundaries, ortholog_matrix)
  # --- MASKING PATCH END ---
  
  gap <- which(msa[,i]=="-")
  oth <- which(msa[,i]!="-")
  
  msa[gap,i] <- "G"
  msa[oth,i] <- "C"
  
  # --- MASKING PATCH START ---
  # Overwrite artificial gaps with missing data
  if (length(mask_leaves) > 0) {
    msa[mask_leaves, i] <- "?"
  }
  # --- MASKING PATCH END ---
}

num_leaves <- length(msa[,1])
vect <- matrix(0, (num_leaves*2), 1)
k <- 0
for (i in seq(1, (num_leaves*2), 2)){
  k <- k + 1
  vect[i, 1] <- paste(">", row.names(msa)[k], sep = "")
}
k <- 0
for (i in seq(2, (num_leaves*2), 2)){
  k <- k + 1
  vect[i, 1] <- paste(msa[k,], collapse = "")
}

write.table(vect, quote = F, col.names = F, row.names = F, sprintf("MSA2/%s_MSA2.fasta", id))
