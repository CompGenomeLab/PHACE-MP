library(ape)
library(Biostrings)
library(tidytree) 
library(stringr)
library(dplyr) 
library(bio3d)
library(mltools)
library("irr")

# Removed aa_to_num and num_to_aa functions as they were not used.

wccc <- function(Y, X, W) {
  mean_Y_w <- weighted.mean(Y, W)
  mean_X_w <- weighted.mean(X, W)
  cov_w <- sum(W * (Y - mean_Y_w) * (X - mean_X_w))
  sd_Y_w_sq <- sum(W * (Y - mean_Y_w)^2)
  sd_X_w_sq <- sum(W * (X - mean_X_w)^2)
  
  wcc <- 2 * cov_w / (sd_Y_w_sq + sd_X_w_sq + (mean_Y_w - mean_X_w)^2)
  return(wcc)
}

args <- commandArgs(trailingOnly = TRUE)
id <- args[1]
slurm_array_task_id <- as.integer(args[2])

file_fasta <- sprintf("MSA1/%s_MSA1.fasta", id)
file_fasta_org <- args[6]
file_nwk <- sprintf("MSA_ASRs/%s_MSA1.treefile", id)

tr_org <- read.tree(file_nwk)
tree_info <- as.data.frame(as_tibble(tr_org))

# --- TREE/MSA ORDER PATCH START ---
# Load helper(s) to enforce MSA row order = tree tip order (strict check).
source("./PHACE_Codes/load_metadata.R")
# --- TREE/MSA ORDER PATCH END ---

# Read fasta file, MSA
fasta <- read.fasta(file = file_fasta)
msa <- fasta$ali

fasta_org <- read.fasta(file = file_fasta_org)
msa_org <- fasta_org$ali

# Names of leaves
names_all <- tr_org[["tip.label"]]
# --- TREE/MSA ORDER PATCH START ---
msa <- reorder_msa_to_tree_tips(msa, names_all, msa_label = "MSA1 (masked)")
msa_org <- reorder_msa_to_tree_tips(msa_org, names_all, msa_label = "Original AA MSA")
# --- TREE/MSA ORDER PATCH END ---
# Number of total leaves&nodes
num_leaves <- length(tr_org[["tip.label"]])
num_nodes <- tr_org[["Nnode"]]

# Total number of positions from ancestralProbs file
total_pos <- ncol(msa)

# --- MASKING PATCH START ---
# Load metadata (MANDATORY for leaf-level masking)
ortholog_table_path <- args[4]
boundaries_path <- args[5]

protein_boundaries <- load_protein_boundaries(boundaries_path)
ortholog_matrix <- load_ortholog_table(ortholog_table_path, tr_org$tip.label)
cat(sprintf("Job %d: Metadata loaded, leaf-level masking enabled\n", slurm_array_task_id))
# --- MASKING PATCH END ---

# --- Parallelization setup ---

# OPTIMIZATION: Avoid generating all pairs in every job.
# Instead, calculate the indices for the pairs this job should process.
total_pairs <- choose(total_pos, 2) # Total number of pairs without generating them
num_windows <- as.integer(args[3]) # Total number of jobs (should match SLURM_ARRAY_TASK_COUNT)

pairs_per_job <- ceiling(total_pairs / num_windows)
start_idx <- (slurm_array_task_id - 1) * pairs_per_job + 1
end_idx <- min(slurm_array_task_id * pairs_per_job, total_pairs)

# This function correctly maps a 1-based index 'k'
# to its (i, j) pair for 'n' items.
get_pair_from_index <- function(n, k) {
  # Convert k to a 0-based index
  k_idx <- k - 1 
  
  # This formula finds the first element 'i' (0-based)
  i <- floor(((2 * n - 1) - sqrt((2 * n - 1)^2 - 8 * k_idx)) / 2)
  
  # This formula finds the second element 'j' (0-based)
  j <- k_idx - (i * (2 * n - i - 1) / 2)
  
  # Return the pair as a 1-based index
  return(c(i + 1, i + j + 2))
}

# Generate only the pairs needed for this specific job.
job_pairs <- lapply(start_idx:end_idx, function(k) get_pair_from_index(total_pos, k))
# --- End of Parallelization setup ---

num_branch <- length(tree_info$branch.length)-1

load(sprintf("totalChanges/%s_TotalChange.RData", id))
mat_total_change <- mat_of_changes$total_change
mat_aa_change <- mat_of_changes$aa_change
mat_aff_branch_num <- mat_of_changes$aff_branch_num
mat_aff_branches <- mat_of_changes$aff_branches
mat_info <- mat_of_changes$info

# --- MASKING PATCH START ---
# Internal-branch masking support:
# Precompute descendant tip membership for every branch row once per run,
# then cache masked-row vectors per protein-pair to avoid per-pair tree traversals.
branch_labels <- mat_info[,2]

desc_tips_by_branch <- compute_descendant_tip_membership(tr_org, branch_labels)
stopifnot(identical(rownames(ortholog_matrix), tr_org$tip.label))

# Cache: position -> protein_id (boundaries lookup is expensive in tight loops)
pos_to_protein <- vapply(seq_len(total_pos), function(p) {
  get_protein_for_position(p, protein_boundaries)
}, character(1))
pos_to_protein[is.na(pos_to_protein)] <- NA_character_

# Cache: protein-pair -> internal masked row indices (in mat_info row space)
mask_rows_internal_cache <- new.env(parent = emptyenv())
get_internal_mask_rows_for_protein_pair <- function(prot1, prot2) {
  if (is.na(prot1) || is.na(prot2)) {
    return(integer(0))
  }
  if (!(prot1 %in% colnames(ortholog_matrix)) || !(prot2 %in% colnames(ortholog_matrix))) {
    return(integer(0))
  }
  key <- paste(sort(c(prot1, prot2)), collapse = "||")
  if (exists(key, envir = mask_rows_internal_cache, inherits = FALSE)) {
    return(get(key, envir = mask_rows_internal_cache, inherits = FALSE))
  }

  present_tips <- ortholog_matrix[, prot1] & ortholog_matrix[, prot2]
  present_tips[is.na(present_tips)] <- FALSE

  present_counts <- as.numeric(desc_tips_by_branch %*% as.integer(present_tips))
  mask_rows <- which(present_counts == 0)

  assign(key, mask_rows, envir = mask_rows_internal_cache)
  return(mask_rows)
}
# --- MASKING PATCH END ---

mat_of_changes_combined_num <- mat_total_change
mat_of_changes_combined_num <-  matrix(as.numeric(mat_of_changes_combined_num), nrow(mat_of_changes_combined_num), ncol(mat_of_changes_combined_num))
row_means <- rowMeans(mat_of_changes_combined_num)
general_mean <- mean(mat_of_changes_combined_num)
row_means[which(row_means<general_mean)] <- general_mean
row_means <- row_means/general_mean
base_weight_branch <- 1/row_means

# Pre-computation outside the loop
base_data <- as.data.frame(mat_info)
base_data$node <- mat_info[,2]
base_data$branch_len <-  mat_info[,1]
base_data <- base_data[,-c(1:2)]

# --- OPTIMIZATION: Pre-computation for loops ---
# 1. Create a "hash map" (named vector) for fast row name lookups
#    This replaces slow `which(rownames(msa_org) == aff_br)`
msa_org_row_indices <- setNames(seq_along(rownames(msa_org)), rownames(msa_org))

# 2. Pre-calculate all amino acid frequencies for all columns
#    This replaces slow `length(which(msa_org[,i1] == alt))`
msa_org_freqs <- apply(msa_org, 2, table)
# --- End of Optimization ---


# Pre-allocate a list to store results (This is already good)
results_list <- vector("list", length(job_pairs))

# Loop over the pairs assigned to this job
for (k in seq_along(job_pairs)) {
  
    # The "%%" (modulo) operator gives the remainder of a division.
    # k %% 1000 == 0 will be TRUE for k = 1000, 2000, 3000, etc.
    if (k %% 10000 == 0) {
    # Print a USEFUL message, not just the numbers
    cat(sprintf("Job %d: Processed %d / %d pairs (%s)\n",
                slurm_array_task_id, k, length(job_pairs), Sys.time()))
    }
  
    pair <- job_pairs[[k]]
    i1 <- pair[1]
    i2 <- pair[2]
    
    # --- MASKING PATCH START ---
    # Determine masked leaves for this pair (union of missing species)
    mask_leaves <- get_masked_leaves_for_pair(i1, i2, protein_boundaries, ortholog_matrix)
    # --- MASKING PATCH END ---

    
    # print(c(i1, i2)) # REMOVED: Printing in a loop is very slow
    
    upgr <- 0
    weight_branch <- base_weight_branch
    
    com_gap <- length(unique(c(names(msa[which(msa[,i1]=="-"),i1]), names(msa[which(msa[,i2]=="-"),i2]))))
    cons1 <- length(unique(msa_org[,i1]))
    cons2 <- length(unique(msa_org[,i2]))
    
    data <- base_data
    data$dif1 <- mat_total_change[,i1]
    data$dif2 <- mat_total_change[,i2]
    data$num_eff1 <- mat_aff_branch_num[,i1]
    data$num_eff2 <- mat_aff_branch_num[,i2]
    data$pos1_aff_br <- mat_aff_branches[,i1]
    data$pos2_aff_br <- mat_aff_branches[,i2]
    data$dif_aa1 <- mat_aa_change[,i1]
    data$dif_aa2 <- mat_aa_change[,i2]
    
    # --- MASKING PATCH START ---
    # Leaf masking (already present) + NEW internal-branch masking:
    # mask any branch whose descendant subtree contains zero "present" tips for this protein pair.
    mask_matrix_rows_leaf <- integer(0)
    if (length(mask_leaves) > 0) {
      mask_matrix_rows_leaf <- map_tips_to_branch_rows(mask_leaves, branch_labels, names_all)
    }

    prot1 <- pos_to_protein[i1]
    prot2 <- pos_to_protein[i2]
    mask_matrix_rows_internal <- get_internal_mask_rows_for_protein_pair(prot1, prot2)

    mask_matrix_rows_all <- sort(unique(c(mask_matrix_rows_leaf, mask_matrix_rows_internal)))
    if (length(mask_matrix_rows_all) > 0) {
      data$dif1[mask_matrix_rows_all] <- 0
      data$dif2[mask_matrix_rows_all] <- 0
      data$dif_aa1[mask_matrix_rows_all] <- "--"
      data$dif_aa2[mask_matrix_rows_all] <- "--"
      # CRITICAL: Also zero weights so masked branches contribute zero to wccc()
      weight_branch[mask_matrix_rows_all] <- 0
    }
    # --- MASKING PATCH END ---

    
    w1 <- which(data$dif_aa1=="X-")
    w2 <- which(data$dif_aa2=="X-")
    
    br1 <- setdiff(unlist(strsplit(data$pos1_aff_br[w1],"-")), "0")
    br2 <- setdiff(unlist(strsplit(data$pos2_aff_br[w2],"-")), "0")
    ww <- length(unique(c(br1, br2)))
    com_ww <- length(intersect(br1, br2))
    
    weight_inc1 <- max(1 - (ww-com_ww)/num_branch,0)
    
    pos1_totalscore <- as.numeric(data$dif1)
    pos2_totalscore <- as.numeric(data$dif2)
    
    score3 <- -1
    
    if (min(cons1, cons2) > 1){
        
        dif_tot <- pos1_totalscore - pos2_totalscore
        loc1 <- which(dif_tot>=0.5)
        loc2 <- which(dif_tot<=-0.5)
        
        # --- OPTIMIZATION: Vectorized loop 1 ---
        # This replaces the slow `for (lx_idx in seq_along(loc1))` loop
        if (length(loc1) > 0) {
        aff_lfs <- as.numeric(data$num_eff1[loc1])
        indices_to_check <- loc1[aff_lfs == 1]
        
        if (length(indices_to_check) > 0) {
            brs_to_check <- data$pos1_aff_br[indices_to_check]
            row_idxs <- msa_org_row_indices[brs_to_check]
            
            valid_rows <- !is.na(row_idxs)
            indices_to_check <- indices_to_check[valid_rows]
            row_idxs <- row_idxs[valid_rows]
            
            if (length(indices_to_check) > 0) {
   
            alts_to_check <- msa_org[row_idxs, i1]

            col_freqs <- msa_org_freqs[[i1]]
            
            obs_orgs <- sapply(alts_to_check, function(alt) {
                count <- col_freqs[alt]
                ifelse(is.na(count), 0, count)
            })
            

            indices_to_zero_out <- indices_to_check[obs_orgs == 1]
            
            if(length(indices_to_zero_out) > 0) {
                upgr <- upgr + length(indices_to_zero_out)
                pos1_totalscore[indices_to_zero_out] <- 0
                pos2_totalscore[indices_to_zero_out] <- 0
            }
            }
        }
        }
        # --- End of optimization ---
        
        # --- OPTIMIZATION: Vectorized loop 2 ---
        # This replaces the slow `for (lx_idx in seq_along(loc2))` loop
        if (length(loc2) > 0) {
        aff_lfs <- as.numeric(data$num_eff2[loc2])
        indices_to_check <- loc2[aff_lfs == 1]
        
        if (length(indices_to_check) > 0) {
            brs_to_check <- data$pos2_aff_br[indices_to_check]
            row_idxs <- msa_org_row_indices[brs_to_check]
            
            valid_rows <- !is.na(row_idxs)
            indices_to_check <- indices_to_check[valid_rows]
            row_idxs <- row_idxs[valid_rows]
            
            if (length(indices_to_check) > 0) {
            alts_to_check <- msa_org[row_idxs, i2]
            col_freqs <- msa_org_freqs[[i2]]
            
            obs_orgs <- sapply(alts_to_check, function(alt) {
                count <- col_freqs[alt]
                ifelse(is.na(count), 0, count)
            })
            
            indices_to_zero_out <- indices_to_check[obs_orgs == 1]
            
            if(length(indices_to_zero_out) > 0) {
                upgr <- upgr + length(indices_to_zero_out)
                pos1_totalscore[indices_to_zero_out] <- 0
                pos2_totalscore[indices_to_zero_out] <- 0
            }
            }
        }
        }
        # --- End of optimization ---

        data$dif1 <- pos1_totalscore
        data$dif2 <- pos2_totalscore
        
        inc_gap1 <- intersect(intersect(which(data$dif_aa1=="X-"), which(data$dif_aa2!="--")),
                                which(data$dif_aa2!="X-"))
        inc_gap2 <- intersect(intersect(which(data$dif_aa2=="X-"), which(data$dif_aa1!="--")),
                                which(data$dif_aa1!="X-"))
        
        if (length(inc_gap1)>0){
            data$dif1[inc_gap1] <- 0
        }
        if (length(inc_gap2)>0){
            data$dif2[inc_gap2] <- 0
        }
        
        pos1_totalscore <- as.numeric(data$dif1)
        pos2_totalscore <- as.numeric(data$dif2)
        
        if (sum(pos1_totalscore)>=1 && sum(pos2_totalscore)>=1 && max(pos1_totalscore)>=0.5 && max(pos2_totalscore)>=0.5){
        df <- abs(pos1_totalscore-pos2_totalscore)
        plc <- which(df>=0.5)
        weight_branch[plc] <- 1
        
        cc <- pmax(pos1_totalscore, pos2_totalscore) # pmax is fast
        cc[cc==0] <- 1
        weight_branch2 <- sqrt(weight_branch*cc)
        score3 <- wccc(pos1_totalscore, pos2_totalscore, weight_branch2)*weight_inc1
          
                
        }
    }
    
    results_list[[k]] <- c(i1, i2, score3)
}

mat <- as.data.frame(do.call(rbind, results_list))
colnames(mat) <- c("Pos1", "Pos2", "PHACE_Score")
save(mat, file= sprintf("PHACE_scores/%s/%s_PHACE_part%d.RData", id, id , slurm_array_task_id))
cat(sprintf("Job %d completed at %s\n", slurm_array_task_id, Sys.time()))
