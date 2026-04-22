library(ape)
library(Biostrings)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)

args = commandArgs(trailingOnly=TRUE)
id <- args[1]

res <- read.csv(sprintf("ToleranceScores/%s.csv", id))
eps <- 10^(-15)
res[,2:21] <- 1-log(res[,2:21] +eps)/log(eps)

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
  sc <- res[i, 2:21]
  
  # --- MASKING PATCH START ---
  # Exclude both "-" AND "?" from dominant state computation
  valid_rows <- which(msa[,i] != "-" & msa[,i] != "?")
  obs <- table(msa[valid_rows, i])
  # --- MASKING PATCH END ---
  
  dom <- names(sort(obs, decreasing = T))[1]
  
  lc <- which(colnames(sc)==dom)
  dom_p2 <- names(sc)[which(as.numeric(sc)>=as.numeric(sc[lc]))]
  
  dom_loc <- c()
  for (aa in c(dom, dom_p2)){
    dom_loc <- c(dom_loc, which(msa[,i]==aa))
  }
 
  msa[dom_loc,i] <- "XX"
  
  upd <- intersect(which(msa[,i]!="XX"), which(msa[,i]!="-"))
  msa[upd,i] <- "A"
  msa[dom_loc, i] <- "C"
  
  # --- MASKING PATCH START ---
  # Overwrite artificial gaps with missing data AFTER C/A/- encoding
  mask_leaves <- get_masked_leaves_for_position(i, protein_boundaries, ortholog_matrix)
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

write.table(vect, quote = F, col.names = F, row.names = F, sprintf("MSA1/%s_MSA1.fasta", id))
