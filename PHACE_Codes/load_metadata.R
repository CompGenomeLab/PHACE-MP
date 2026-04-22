# load_metadata.R
# Helper functions for PHACE gap masking
# Loads metadata and provides position → protein → masked species mapping

# Load and validate protein boundaries
load_protein_boundaries <- function(filepath) {
  if (!file.exists(filepath)) {
    stop(sprintf("Protein boundaries file not found: %s", filepath))
  }
  
  boundaries <- read.csv(filepath, stringsAsFactors = FALSE)
  
  # Validate required columns
  required_cols <- c("protein_name", "start", "end")
  missing_cols <- setdiff(required_cols, colnames(boundaries))
  if (length(missing_cols) > 0) {
    stop(sprintf("Protein boundaries file missing required columns: %s", 
                 paste(missing_cols, collapse = ", ")))
  }
  
  # Validate no gaps or overlaps
  boundaries <- boundaries[order(boundaries$start), ]
  for (i in 1:(nrow(boundaries) - 1)) {
    if (boundaries$end[i] >= boundaries$start[i + 1]) {
      stop(sprintf("Overlapping protein boundaries: %s (end=%d) and %s (start=%d)",
                   boundaries$protein_name[i], boundaries$end[i],
                   boundaries$protein_name[i + 1], boundaries$start[i + 1]))
    }
  }
  
  cat(sprintf("Loaded %d protein boundaries\n", nrow(boundaries)))
  return(boundaries)
}

# Load and validate ortholog table
load_ortholog_table <- function(filepath, tip_labels) {
  if (!file.exists(filepath)) {
    stop(sprintf("Ortholog table file not found: %s", filepath))
  }
  
  ortholog_table <- read.table(filepath, header = TRUE, sep = "\t", 
                                stringsAsFactors = FALSE, row.names = 1,
                                check.names = FALSE)
  
  # Validate species names match tip labels
  table_species <- rownames(ortholog_table)
  missing_in_table <- setdiff(tip_labels, table_species)
  missing_in_tree <- setdiff(table_species, tip_labels)
  
  if (length(missing_in_table) > 0) {
    cat(sprintf("WARNING: %d species in tree but not in ortholog table. They will be treated as missing all proteins.\n",
                length(missing_in_table)))
    cat("Missing species:", paste(head(missing_in_table, 10), collapse = ", "), "\n")
  }
  
  if (length(missing_in_tree) > 0) {
    cat(sprintf("WARNING: %d species in ortholog table but not in tree (will be ignored)\n",
                length(missing_in_tree)))
  }
  
  # Reorder rows to match tip_labels and add missing species as all-NA
  ortholog_matrix <- matrix(NA, nrow = length(tip_labels), ncol = ncol(ortholog_table),
                            dimnames = list(tip_labels, colnames(ortholog_table)))
  
  for (sp in tip_labels) {
    if (sp %in% table_species) {
      ortholog_matrix[sp, ] <- as.matrix(ortholog_table[sp, ])
    }
    # else: remains NA (treated as missing all proteins)
  }
  
  # Convert to logical: TRUE = has ortholog, FALSE/NA = missing
  # Assuming non-empty/non-NA values mean "has ortholog"
  ortholog_logical <- !is.na(ortholog_matrix) & ortholog_matrix != ""
  
  cat(sprintf("Ortholog table: %d species × %d proteins\n",
              nrow(ortholog_logical), ncol(ortholog_logical)))
  
  return(ortholog_logical)
}

# Map position to protein
get_protein_for_position <- function(pos, boundaries) {
  idx <- which(boundaries$start <= pos & boundaries$end >= pos)
  
  if (length(idx) == 0) {
    # Warn once per unmapped position (use a global flag)
    if (!exists(".warned_unmapped_positions", envir = .GlobalEnv)) {
      assign(".warned_unmapped_positions", TRUE, envir = .GlobalEnv)
      cat(sprintf("WARNING: Position %d does not map to any protein block. No masking for this position.\n", pos))
    }
    return(NA)
  }
  
  if (length(idx) > 1) {
    stop(sprintf("Position %d maps to multiple protein blocks: %s",
                 pos, paste(boundaries$protein_name[idx], collapse = ", ")))
  }
  
  # Return protein_id (not protein_name) because ortholog table columns use IDs
  return(boundaries$protein_id[idx])
}

# Get masked leaf indices for a position
get_masked_leaves_for_position <- function(pos, boundaries, ortholog_matrix) {
  protein <- get_protein_for_position(pos, boundaries)
  
  if (is.na(protein)) {
    return(integer(0))  # No masking for unmapped positions
  }
  
  # Check if protein exists in ortholog matrix columns
  if (!protein %in% colnames(ortholog_matrix)) {
    cat(sprintf("WARNING: Protein %s not found in ortholog table columns. No masking for position %d.\n",
                protein, pos))
    return(integer(0))
  }
  
  # Find species (row indices) that are missing this protein
  # ortholog_matrix[species, protein] = TRUE means has ortholog
  # We want indices where it's FALSE or NA
  has_protein <- ortholog_matrix[, protein]
  mask_indices <- which(!has_protein | is.na(has_protein))
  
  return(mask_indices)
}

# Get masked leaf indices for a position pair (union)
get_masked_leaves_for_pair <- function(pos1, pos2, boundaries, ortholog_matrix) {
  mask1 <- get_masked_leaves_for_position(pos1, boundaries, ortholog_matrix)
  mask2 <- get_masked_leaves_for_position(pos2, boundaries, ortholog_matrix)
  
  # Union: species missing at least one of the two proteins
  return(unique(c(mask1, mask2)))
}

# Map tip label indices to branch_info row indices for safe masking
# tip_indices: 1-based indices into tip_labels (e.g., from get_masked_leaves_for_*)
# branch_labels: vector of branch labels from mat_info[,2] or tree_info$label
# tip_labels: tr_org$tip.label (tip labels in tree order)
# Returns: row indices into the change matrices where branch_labels match tip_labels[tip_indices]
# Uses match() to guarantee 1:1 mapping and detect missing/duplicate labels
map_tips_to_branch_rows <- function(tip_indices, branch_labels, tip_labels) {
  if (length(tip_indices) == 0) {
    return(integer(0))
  }
  
  # Get the tip labels to mask
  tips_to_mask <- tip_labels[tip_indices]
  
  # Use match() for safe 1:1 mapping: returns one index per tip (or NA if not found)
  matched_rows <- match(tips_to_mask, branch_labels)
  
  # Check for unmapped tips (NAs)
  unmapped_mask <- is.na(matched_rows)
  if (any(unmapped_mask)) {
    unmapped_tips <- tips_to_mask[unmapped_mask]
    warning(sprintf("Could not map %d tip labels to branch rows: %s", 
                    length(unmapped_tips), paste(head(unmapped_tips, 5), collapse=", ")))
    # Remove NAs from result
    matched_rows <- matched_rows[!unmapped_mask]
  }
  
  # Check for potential duplicates in branch_labels (would cause silent mis-masking)
  if (any(duplicated(branch_labels[matched_rows]))) {
    warning("Duplicate branch labels detected in matched rows - masking may be incorrect")
  }
  
  return(matched_rows)
}

# Enforce: MSA rownames must be exactly the tree tip labels (set equality),
# then reorder deterministically to tip order.
#
# This prevents brittle "silent assumptions" where missing/extra taxa are
# silently dropped/introduced as NA rows by matrix subsetting.
reorder_msa_to_tree_tips <- function(msa, tip_labels, msa_label = "MSA") {
  if (is.null(msa)) {
    stop(sprintf("%s is NULL; cannot reorder to tree tips.", msa_label))
  }
  if (is.null(rownames(msa))) {
    stop(sprintf("%s has no rownames; expected species names as rownames.", msa_label))
  }
  if (any(duplicated(rownames(msa)))) {
    dups <- unique(rownames(msa)[duplicated(rownames(msa))])
    stop(sprintf("%s has duplicate rownames (species labels). Example(s): %s",
                 msa_label, paste(head(dups, 10), collapse = ", ")))
  }
  missing_in_msa <- setdiff(tip_labels, rownames(msa))
  extra_in_msa <- setdiff(rownames(msa), tip_labels)
  if (length(missing_in_msa) > 0 || length(extra_in_msa) > 0) {
    stop(sprintf(
      paste0(
        "%s species labels do not match tree tip labels (strict check failed).\n",
        "- Missing in %s (present in tree): %d. Example(s): %s\n",
        "- Extra in %s (absent from tree): %d. Example(s): %s\n",
        "Fix input names so the sets are identical."
      ),
      msa_label, msa_label, length(missing_in_msa), paste(head(missing_in_msa, 10), collapse = ", "),
      msa_label, length(extra_in_msa), paste(head(extra_in_msa, 10), collapse = ", ")
    ))
  }
  msa <- msa[tip_labels, , drop = FALSE]
  return(msa)
}


# Precompute descendant-tip membership for each branch label (matrix row label).
#
# In PHACE change matrices, each row corresponds to the branch leading into the
# node whose label is stored in `branch_labels` (tips + internal nodes, excluding root).
#
# Returns a logical matrix with:
# - rows in the same order as `branch_labels`
# - columns in `tr_org$tip.label` order
#
# This representation is compact (O((2N-2)*N) booleans) and allows fast masking via
# `%*%` against a present-tips vector.
compute_descendant_tip_membership <- function(tr_org, branch_labels) {
  if (is.null(tr_org) || class(tr_org)[1] != "phylo") {
    stop("compute_descendant_tip_membership: tr_org must be a phylo object.")
  }
  if (is.null(branch_labels) || length(branch_labels) == 0) {
    stop("compute_descendant_tip_membership: branch_labels is empty.")
  }

  num_tips <- length(tr_org$tip.label)
  num_nodes <- tr_org$Nnode
  num_all <- num_tips + num_nodes

  # Map branch labels to child node indices in the phylo numbering.
  tip_idx <- match(branch_labels, tr_org$tip.label)
  node_idx <- match(branch_labels, tr_org$node.label)

  child_nodes <- rep(NA_integer_, length(branch_labels))
  child_nodes[!is.na(tip_idx)] <- tip_idx[!is.na(tip_idx)]
  child_nodes[is.na(tip_idx) & !is.na(node_idx)] <- num_tips + node_idx[is.na(tip_idx) & !is.na(node_idx)]

  if (any(is.na(child_nodes))) {
    bad <- branch_labels[is.na(child_nodes)]
    stop(sprintf(
      paste0(
        "Could not map %d branch labels to tree nodes for internal masking.\n",
        "Example(s): %s\n",
        "Expected each label to be either a tip label or an internal node label in the Newick tree."
      ),
      length(bad), paste(head(bad, 10), collapse = ", ")
    ))
  }

  # Compute descendant tips for every node using a postorder dynamic program.
  tr_post <- ape::reorder.phylo(tr_org, order = "postorder")
  edge <- tr_post$edge

  desc_tips_by_node <- matrix(FALSE, nrow = num_all, ncol = num_tips)
  desc_tips_by_node[seq_len(num_tips), seq_len(num_tips)] <- diag(TRUE, num_tips)

  for (e in seq_len(nrow(edge))) {
    parent <- edge[e, 1]
    child <- edge[e, 2]
    desc_tips_by_node[parent, ] <- desc_tips_by_node[parent, ] | desc_tips_by_node[child, ]
  }

  desc_tips_by_branch <- desc_tips_by_node[child_nodes, , drop = FALSE]
  colnames(desc_tips_by_branch) <- tr_org$tip.label
  return(desc_tips_by_branch)
}
