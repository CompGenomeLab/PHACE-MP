 #!/usr/bin/env Rscript
library(ape)
library(tidytree)
library(stringr)
library(dplyr)
library(bio3d)
library(Peptides)
source("./PHACE_Codes/position_score.R")

args = commandArgs(trailingOnly=TRUE)

aa_to_num <- function(aa) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  num <- sapply(aa, function(a){ifelse(sum(amino_acids %in% a) == 1, as.numeric(which(amino_acids %in% a)), 21)})
  # num <- ifelse(sum(amino_acids %in% aa) == 1, as.numeric(which(amino_acids %in% aa)), 21)
  return(num)
}

num_to_aa <- function(num) {
  amino_acids <- c("G", "A", "L", "M", "F", "W", "K", "Q", "E", "S", "P", "V", "I", "C", "Y", "H", "R", "N", "D", "T")
  aa <- ifelse(num == 21, 21, amino_acids[num])
  return(aa)
}

id <- args[1]
file_fasta <- args[2]
file_nwk <- args[3]
file_rst <- args[4]

output_name <- id

pos_chosen <- "all"

tr_org <- read.tree(file_nwk)

# --- MASKING PATCH START ---
# Load metadata (MANDATORY for leaf-level masking)
boundaries_path <- args[5]
ortholog_table_path <- args[6]

source("./PHACE_Codes/load_metadata.R")
protein_boundaries <- load_protein_boundaries(boundaries_path)
ortholog_matrix <- load_ortholog_table(ortholog_table_path, tr_org$tip.label)
cat("Metadata loaded. Leaf-level masking enabled.\n")
# --- MASKING PATCH END ---

# Read ancestral state file
x <- read.table(file = file_rst, sep = '\t', header = TRUE, fill = TRUE)
x[,1] <- str_remove(x[,1], "Node")
colnames(x)[4:ncol(x)] <- gsub("p_", replacement = "", x = colnames(x)[4:ncol(x)], fixed = TRUE)

# --- PARTITIONED STATE FIX START ---
if ("Part" %in% colnames(x)) {
  # optional safety check
  if (max(x$Part, na.rm = TRUE) > nrow(protein_boundaries)) {
    stop("State file Part indices exceed number of entries in protein_boundaries")
  }

  x$LocalSite <- x$Site

  # map partition-local sites to global concatenated coordinates
  x$Site <- protein_boundaries$start[x$Part] + x$LocalSite - 1

  # reorder by node and global site
  x <- x[order(as.numeric(x[,1]), x$Site), ]

  cat("Partitioned .state file detected. Site column remapped to global concatenated positions.\n")
}
# --- PARTITIONED STATE FIX END ---

tree_info <- as.data.frame(as_tibble(tr_org))

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali
# #region agent log
msa_nrow <- if (is.null(dim(msa))) length(msa) else nrow(msa)
msa_ncol <- if (is.null(dim(msa))) NA else ncol(msa)


# connections_1: Parent node, connections_2: connected node/leaf
connections_1 <- tree_info$parent
connections_2 <- tree_info$node

# Names of leaves
names_all <- tr_org[["tip.label"]]
# --- TREE/MSA ORDER PATCH START ---
msa <- reorder_msa_to_tree_tips(msa, names_all, msa_label = "AA MSA (ToleranceScore)")
# --- TREE/MSA ORDER PATCH END ---

# Number of total leaves&nodes
num_leaves <- length(tr_org[["tip.label"]])
num_nodes <- tr_org[["Nnode"]]

nodes_raxml <- as.numeric(gsub(pattern = "Node", replacement = "", x = tree_info[num_leaves+1:num_nodes, "label"])) #Node or Branch
names(nodes_raxml) <- tree_info[num_leaves+1:num_nodes, "node"]
nodes_raxml_prev <- nodes_raxml
num_nodes_prev <- num_nodes
num_leaves_prev <- num_leaves

# Total number of positions from ancestralProbs file
total_pos <- max(x$Site)
# #region agent log
site_min <- min(x$Site, na.rm = TRUE)
site_max <- max(x$Site, na.rm = TRUE)
site_unique <- length(unique(x$Site))
positions <- 1:total_pos
score_all <- matrix(0, total_pos, 21)

####################################################
####################################################
chosen_leaves <- tree_info[1:num_leaves,c("parent", "node")]
chosen_nodes <- tree_info[(num_leaves+2):(num_leaves +num_nodes),c("parent", "node")]
leaf_names <- tree_info$label

chosen_nodes2 <- chosen_nodes

scores <- t(mapply(function(ps){
  # --- MASKING PATCH START ---
  mask_leaves <- get_masked_leaves_for_position(ps, protein_boundaries, ortholog_matrix)
  # --- MASKING PATCH END ---
  
  position_score(ps, x, msa, trim_final, names_all, tr_org, num_nodes, num_leaves, tree_info, num_nodes, nodes_raxml, num_leaves, total_pos, nodes_raxml, mask_leaves = mask_leaves)
}, rep(positions)))

tolerance_scores <- cbind(positions, scores)
colnames(tolerance_scores) <- c("Pos/AA", num_to_aa(1:20))

write.csv(tolerance_scores, quote = F, row.names = F, paste("ToleranceScores/", output_name, ".csv", sep = ""))
