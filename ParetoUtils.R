# pareto front function

is_pareto_efficient <- function(fitness_matrix) {
  n <- nrow(fitness_matrix)
  is_efficient <- rep(TRUE, n)
  for (i in seq_len(n)) {
    if (is_efficient[i]) {
      is_dominated <- apply(fitness_matrix, 1, function(x) {
        all(x >= fitness_matrix[i, ]) && any(x > fitness_matrix[i, ])
      })
      is_dominated[i] <- FALSE
      if (any(is_dominated)) {
        is_efficient[i] <- FALSE
      }
    }
  }
  return(is_efficient)
}

# specialist/generalist function

compute_epsilon <- function(x, y, pareto_df) {
  distances <- sqrt((x - pareto_df$mean_Glycolysis)^2 + (y - pareto_df$mean_OxPhos)^2)
  return(min(distances))
}

# interpatient function

run_pairwise_regressions <- function(listA, listB, group_label) {
  results <- data.frame(Pair = character(), Slope = numeric(), R2 = numeric(), stringsAsFactors = FALSE)
  
  for (a_name in names(listA)) {
    for (b_name in names(listB)) {
      a_df <- listA[[a_name]]
      b_df <- listB[[b_name]]
      
      shared_k <- intersect(a_df$k, b_df$k)
      if (length(shared_k) == 0) next  
      
      a_shared <- a_df[a_df$k %in% shared_k, ]
      b_shared <- b_df[b_df$k %in% shared_k, ]
      
      a_shared <- a_shared[order(a_shared$k), ]
      b_shared <- b_shared[order(b_shared$k), ]
      
      comp_df <- data.frame(mean_A = a_shared$mean, mean_B = b_shared$mean)
      fit_mat <- as.matrix(comp_df)
      pareto_mask <- is_pareto_efficient(fit_mat)
      comp_df <- comp_df[pareto_mask, ]
      
      if (nrow(comp_df) < 2) next
      
      lm_model <- lm(mean_B ~ mean_A, data = comp_df)
      slope <- coef(lm_model)[2]
      r2 <- summary(lm_model)$r.squared
      
      results <- rbind(results, data.frame(
        Pair = paste0(a_name, "_vs_", b_name),
        Slope = round(slope, 4),
        R2 = round(r2, 4),
        Group = group_label,
        stringsAsFactors = FALSE
      ))
    }
  }
  return(results)
}
