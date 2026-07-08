# Multi-seed stability check for the Powell/Brent coordinate-wise optimizer.
#
# Context: JSS editor/referee comments ask how frequently near-flat log-likelihood
# regions / boundary solutions occur in practice, and what role random
# initialization plays in the reported optimum (see REVISION_NOTES.md sections 7,
# 7bis and 8, and Rev_2_3.md). A single seed cannot answer "how frequently" -
# this script refits the paper's own worked example (same formula/time_axis/
# parameter ranges as ReplicationCode.R) from several different random seeds and
# reports, for each one: the optimized log-likelihood, convergence status/run
# count, the estimated regression coefficients, mu1/nu, the first-order
# optimality (gradient) diagnostic and its boundary-adjacent parameters, and
# (where the interior Hessian is positive definite) the full-Hessian standard
# error correction. Comparing these across seeds is direct evidence for whether
# the reported optimum - and the boundary/ridge behaviour discussed in the
# revision notes - is stable, rather than an artifact of one particular seed.
#
# Usage: Rscript Examples/MultiSeedStabilityCheck.R <seed>
# Each fit takes several minutes; intended to be run once per seed, in parallel.

args <- commandArgs(trailingOnly = TRUE)
seed <- suppressWarnings(as.integer(args[1]))
if (is.na(seed)) stop("Usage: Rscript MultiSeedStabilityCheck.R <seed>")

devtools::load_all("/Users/alessandragni/Documents/DATA/POLITECNICO/PHD/CODE_REPO/TimeDepFrail", quiet = TRUE)

data(data_dropout)
formula <- time_to_event ~ Gender + CFUP + cluster(group)
time_axis <- c(1.0, 1.4, 1.8, 2.3, 3.1, 3.8, 4.3, 5.0, 5.5, 5.8, 6.0)
eps <- 1e-10
categories_range_min <- c(-8, -2, eps, eps, eps)
categories_range_max <- c(-eps, 0.5, 1 - eps, 1, 10)

cat("=== seed =", seed, "===\n")
set.seed(seed)
t0 <- Sys.time()
result <- AdPaikModel(formula, data_dropout, time_axis,
                      categories_range_min, categories_range_max,
                      full_hessian_se = TRUE)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("Fit time:", elapsed, "sec\n\n")

cat("Loglikelihood:", result$Loglikelihood, "\n")
cat("AIC:", result$AIC, "\n")
cat("Status:", result$Status, " NRun:", result$NRun, "\n")
cat("beta (GenderMale, CFUP):", result$OptimalParameters[11:12], "\n")
cat("mu1, nu:", result$OptimalParameters[13:14], "\n\n")

cat("GradientCheck: GradientNorm =", result$GradientCheck$GradientNorm,
   " GradientMaxAbs =", result$GradientCheck$GradientMaxAbs, "\n")
cat("BoundaryAdjacent indices:", result$GradientCheck$BoundaryAdjacent, "\n\n")

cat("HessianCheck IsPositiveDefinite:", result$HessianCheck$IsPositiveDefinite, "\n")
cat("HessianCheck Mu1NuBlockPositiveDefinite:", result$HessianCheck$Mu1NuBlockPositiveDefinite, "\n")
cat("SE_ratio (GenderMale, CFUP):", result$HessianCheck$SE_ratio[11:12], "\n")

cat("\n--- one-line summary (for cross-seed comparison table) ---\n")
cat(sprintf("SEED=%d LL=%.6f NRun=%d Status=%s beta_GenderMale=%.6f beta_CFUP=%.6f mu1=%.6g nu=%.6f GradNorm=%.6f BoundaryAdjacent=[%s] SEratio_GenderMale=%.4f SEratio_CFUP=%.4f\n",
           seed, result$Loglikelihood, result$NRun, result$Status,
           result$OptimalParameters[11], result$OptimalParameters[12],
           result$OptimalParameters[13], result$OptimalParameters[14],
           result$GradientCheck$GradientNorm,
           paste(result$GradientCheck$BoundaryAdjacent, collapse = ","),
           result$HessianCheck$SE_ratio[11], result$HessianCheck$SE_ratio[12]))

cat("\nDONE seed", seed, "\n")
