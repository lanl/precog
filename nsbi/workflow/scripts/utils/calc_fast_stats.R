#!/usr/bin/env Rscript

#' Calculate Fast Tree Statistics
#' 
#' Optimized version that computes fast statistics (< 0.1 seconds for 10k tips)
#' with optional tree size normalization for coalescent/epidemic trees.
#' 
#' Based on comprehensive benchmarking:
#' - Only includes statistics with mean_time < 0.1s on 10,000 tip trees
#' - Excludes slow stats: vpd (0.54s), phylogenetic_div (0.97s), imbalance_steps (0.54s)
#' - Excludes ultrametric-only stats: gamma, nltt_base (return NA for non-ultrametric trees)
#' - Includes proper function name mappings (e.g., beta_statistic, entropy_j)
#' 
#' Size Agnostic Mode (NEW - STRICT INTERPRETATION):
#' - When size_agnostic=TRUE: Returns ONLY 14 inherently size-agnostic statistics
#' - These are statistics that are independent of tree size by their nature (ratios, branch lengths, entropy)
#' - Does NOT use Yule/PDA normalization (those methods are not truly tree size-agnostic)
#' - When size_agnostic=FALSE: Returns all 51 statistics (original behavior)
#' 
#' Performance:
#' - size_agnostic=TRUE: ~0.02s per 10k tip tree (14 stats)
#' - size_agnostic=FALSE: ~0.055s per 10k tip tree (51 stats)
#' - For 10,000 trees: 5-15 minutes total
#' 
#' @param tree A phylo object (from ape package)
#' @param size_agnostic Logical. If TRUE, returns only inherently size-agnostic stats (14 stats). Default: TRUE
#' @param normalization Character. Ignored when size_agnostic=TRUE. For size_agnostic=FALSE, not currently used. Default: "both"
#' @return A named vector of tree statistics
#' @export

library(treestats)

calc_fast_stats <- function(tree, size_agnostic = TRUE, normalization = "both") {
  
  # When size_agnostic=TRUE, we compute ONLY inherently size-agnostic statistics
  # (no Yule/PDA normalization, as those methods are not truly size-agnostic)
  # This gives us 14 statistics that are independent of tree size by their nature.
  
  # When size_agnostic=FALSE, we compute all 51 statistics (original behavior, no normalization)
  # The normalization parameter is ignored in both modes
  
  if (size_agnostic) {
    # Strict size-agnostic mode: only inherently agnostic stats
    compute_yule <- FALSE
    compute_pda <- FALSE
    compute_none <- FALSE
    compute_agnostic_only <- TRUE
  } else {
    # Size-dependent mode: compute all stats WITHOUT normalization (original behavior)
    compute_yule <- FALSE
    compute_pda <- FALSE
    compute_none <- TRUE
    compute_agnostic_only <- FALSE
  }
  
  stats <- list()
  
  # Helper function to safely compute a statistic
  try_stat <- function(stat_name, func_call) {
    tryCatch({
      func_call
    }, error = function(e) {
      NA
    })
  }
  
  # Helper function to compute a stat with multiple normalizations
  compute_stat_normalized <- function(base_name, stat_func, supports_norm = TRUE) {
    # In strict agnostic mode, skip all normalized stats
    if (compute_agnostic_only) {
      return()  # Don't compute this stat
    }
    
    if (!supports_norm) {
      # Stat doesn't support normalization, compute once
      stats[[base_name]] <<- try_stat(base_name, stat_func("none"))
    } else if (compute_none) {
      # Size-dependent mode: use no normalization
      stats[[base_name]] <<- try_stat(base_name, stat_func("none"))
    } else {
      # Size-agnostic mode with normalization: compute with requested normalizations
      if (compute_yule) {
        stats[[paste0(base_name, "_yule")]] <<- try_stat(paste0(base_name, "_yule"), stat_func("yule"))
      }
      if (compute_pda) {
        stats[[paste0(base_name, "_pda")]] <<- try_stat(paste0(base_name, "_pda"), stat_func("pda"))
      }
    }
  }
  
  # ===== BALANCE & IMBALANCE METRICS =====
  # Most require Yule/PDA normalization - skip in strict agnostic mode
  compute_stat_normalized("colless", function(norm) treestats::colless(tree, normalization = norm))
  compute_stat_normalized("colless_corr", function(norm) treestats::colless_corr(tree, normalization = norm))
  compute_stat_normalized("colless_quad", function(norm) treestats::colless_quad(tree, normalization = norm))
  compute_stat_normalized("sackin", function(norm) treestats::sackin(tree, normalization = norm))
  compute_stat_normalized("b1", function(norm) treestats::b1(tree, normalization = norm))
  compute_stat_normalized("b2", function(norm) treestats::b2(tree, normalization = norm))
  compute_stat_normalized("rogers", function(norm) treestats::rogers(tree, normalization = norm))
  compute_stat_normalized("rquartet", function(norm) treestats::rquartet(tree, normalization = norm))
  
  # These are inherently size-agnostic (computed in all modes)
  stats$ew_colless <- try_stat("ew_colless", treestats::ew_colless(tree))
  stats$beta <- try_stat("beta", treestats::beta_statistic(tree))
  stats$root_imbalance <- try_stat("root_imbalance", treestats::root_imbalance(tree))
  
  # blum: uses TRUE/FALSE for normalization parameter - skip in strict agnostic mode
  if (compute_agnostic_only) {
    # Skip blum in strict agnostic mode
  } else if (compute_none) {
    stats$blum <- try_stat("blum", treestats::blum(tree, normalization = FALSE))
  } else {
    if (compute_yule) {
      stats$blum_yule <- try_stat("blum_yule", treestats::blum(tree, normalization = TRUE))
    }
    if (compute_pda) {
      stats$blum_pda <- try_stat("blum_pda", treestats::blum(tree, normalization = TRUE))
    }
  }
  
  # ===== TOPOLOGY METRICS =====
  # cherries and pitchforks require normalization - skip in strict agnostic mode
  compute_stat_normalized("cherries", function(norm) treestats::cherries(tree, normalization = norm))
  compute_stat_normalized("pitchforks", function(norm) treestats::pitchforks(tree, normalization = norm))
  
  # These topology metrics do NOT support normalization - exclude in all agnostic modes
  if (!compute_agnostic_only && compute_none) {
    stats$double_cherries <- try_stat("double_cherries", treestats::double_cherries(tree))
    stats$four_prong <- try_stat("four_prong", treestats::four_prong(tree))
    stats$stairs <- try_stat("stairs", treestats::stairs(tree))
    stats$stairs2 <- try_stat("stairs2", treestats::stairs2(tree))
    stats$avg_ladder <- try_stat("avg_ladder", treestats::avg_ladder(tree))
    stats$max_ladder <- try_stat("max_ladder", treestats::max_ladder(tree))
  }
  
  # il_number and symmetry_nodes require normalization - skip in strict agnostic mode
  compute_stat_normalized("il_number", function(norm) treestats::ILnumber(tree, normalization = norm))
  compute_stat_normalized("symmetry_nodes", function(norm) treestats::sym_nodes(tree, normalization = norm))
  
  # tot_path doesn't have normalization parameter - exclude in all agnostic modes
  if (!compute_agnostic_only && compute_none) {
    stats$tot_path <- try_stat("tot_path", treestats::tot_path_length(tree))
  }
  
  # ===== BRANCH LENGTH STATISTICS =====
  # These are inherently size-agnostic (temporal measures, not count-based)
  stats$mean_branch_length <- try_stat("mean_branch_length", treestats::mean_branch_length(tree))
  stats$mean_branch_length_int <- try_stat("mean_branch_length_int", treestats::mean_branch_length_int(tree))
  stats$mean_branch_length_ext <- try_stat("mean_branch_length_ext", treestats::mean_branch_length_ext(tree))
  stats$var_branch_length <- try_stat("var_branch_length", treestats::var_branch_length(tree))
  stats$var_branch_length_int <- try_stat("var_branch_length_int", treestats::var_branch_length_int(tree))
  stats$var_branch_length_ext <- try_stat("var_branch_length_ext", treestats::var_branch_length_ext(tree))
  stats$treeness <- try_stat("treeness", treestats::treeness(tree))
  
  # ===== TREE SHAPE & SIZE METRICS =====
  # Most require normalization - skip in strict agnostic mode
  compute_stat_normalized("max_depth", function(norm) treestats::max_depth(tree, normalization = norm))
  compute_stat_normalized("max_width", function(norm) treestats::max_width(tree, normalization = norm))
  compute_stat_normalized("average_leaf_depth", function(norm) treestats::average_leaf_depth(tree, normalization = norm))
  compute_stat_normalized("max_del_width", function(norm) treestats::max_del_width(tree, normalization = norm))
  compute_stat_normalized("var_depth", function(norm) treestats::var_leaf_depth(tree, normalization = norm))
  
  # These don't support normalization - exclude in all agnostic modes
  if (!compute_agnostic_only && compute_none) {
    stats$diameter <- try_stat("diameter", treestats::diameter(tree))
    stats$avg_vert_depth <- try_stat("avg_vert_depth", treestats::avg_vert_depth(tree))
    stats$number_of_lineages <- try_stat("number_of_lineages", treestats::number_of_lineages(tree))
  }
  
  # mw_over_md is a ratio - inherently size-agnostic (computed in all modes)
  stats$mw_over_md <- try_stat("mw_over_md", treestats::mw_over_md(tree))
  
  # ===== DISTANCE METRICS =====
  # Most require normalization - skip in strict agnostic mode
  compute_stat_normalized("tot_coph", function(norm) treestats::tot_coph(tree, normalization = norm))
  compute_stat_normalized("area_per_pair", function(norm) treestats::area_per_pair(tree, normalization = norm))
  compute_stat_normalized("mpd", function(norm) treestats::mean_pair_dist(tree, normalization = norm))
  
  # tot_internal_path doesn't support normalization - exclude in all agnostic modes
  if (!compute_agnostic_only && compute_none) {
    stats$tot_internal_path <- try_stat("tot_internal_path", treestats::tot_internal_path(tree))
  }
  
  # wiener: uses TRUE/FALSE for normalization - skip in strict agnostic mode
  if (compute_agnostic_only) {
    # Skip wiener in strict agnostic mode
  } else if (compute_none) {
    stats$wiener <- try_stat("wiener", treestats::wiener(tree, normalization = FALSE))
  } else {
    if (compute_yule) {
      stats$wiener_yule <- try_stat("wiener_yule", treestats::wiener(tree, normalization = TRUE))
    }
    if (compute_pda) {
      stats$wiener_pda <- try_stat("wiener_pda", treestats::wiener(tree, normalization = TRUE))
    }
  }
  
  # ===== CENTRALITY METRICS =====
  # All require normalization - skip in strict agnostic mode
  compute_stat_normalized("max_betweenness", function(norm) treestats::max_betweenness(tree, normalization = norm))
  compute_stat_normalized("max_closeness", function(norm) treestats::max_closeness(tree, weight = FALSE, normalization = norm))
  compute_stat_normalized("max_closenessW", function(norm) treestats::max_closeness(tree, weight = TRUE, normalization = norm))
  
  # ===== OTHER METRICS =====
  # These are information-theoretic/correlation-based (inherently size-agnostic)
  stats$i_stat <- try_stat("i_stat", treestats::mean_i(tree))
  stats$j_one <- try_stat("j_one", treestats::j_one(tree))
  stats$j_stat <- try_stat("j_stat", treestats::entropy_j(tree))
  
  # psv requires normalization - skip in strict agnostic mode
  compute_stat_normalized("psv", function(norm) treestats::psv(tree, normalization = norm))
  
  # Convert to named vector (matching calc_all_stats output format)
  stats <- unlist(stats)
  stats <- stats[order(names(stats))]
  
  return(stats)
}


# If running as script
if (!interactive()) {
  cat("========================================\n")
  cat("calc_fast_stats() function loaded\n")
  cat("========================================\n\n")
  cat("UPDATED: Strict size-agnostic mode!\n\n")
  cat("Usage:\n")
  cat("  library(ape)\n")
  cat("  source('calc_fast_stats.R')\n")
  cat("  tree <- read.tree('tree.nwk')\n\n")
  cat("  # Size-agnostic mode (default) - ONLY inherently agnostic stats:\n")
  cat("  stats <- calc_fast_stats(tree)\n")
  cat("  # Returns 14 stats (no Yule/PDA normalization)\n\n")
  cat("  # Size-dependent mode - all 51 original stats:\n")
  cat("  stats <- calc_fast_stats(tree, size_agnostic=FALSE)\n")
  cat("  # Returns 51 stats (original behavior)\n\n")
  cat("Performance:\n")
  cat("  size_agnostic=TRUE: ~0.02s per 10k tip tree (14 stats)\n")
  cat("  size_agnostic=FALSE: ~0.055s per 10k tip tree (51 stats)\n")
  cat("  For 10,000 trees: 5-15 minutes total\n")
  cat("========================================\n")
}

