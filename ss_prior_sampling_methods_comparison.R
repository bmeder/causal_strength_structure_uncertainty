# ============================================================
# Sanity check for Lu et al. (2008) / Meder et al. (2014)
# SS prior implementation: direct truncated-exponential mixture
# versus uniform-sampling + SS prior weights.
# ============================================================
#
# This file compares three implementations of the
# generative Lu et al. (2008) SS strength model.
#
# Lu, H., Yuille, A. L., Liljeholm, M., Cheng, P. W., & Holyoak, K. J. (2008).
# Bayesian generic priors for causal learning. Psychological Review, 115, 955-984.
# Matlab code: https://cvl.psych.ucla.edu/publications/
#
# Generative SS prior:
#   p(wc, wa) proportional to
#     exp[-alpha*wa - alpha*(1-wc)] +
#     exp[-alpha*(1-wa) - alpha*wc]
#
# Mapping:
#   Lu et al. notation:   w0 = background strength, w1 = candidate strength
#   SI/Meder notation:    wa = background strength, wc = candidate strength
#
# 1. direct_SS_prior
#    Draws (wc, wa) directly from the joint SS prior using a two-component
#    truncated-exponential mixture. Posterior predictions are then computed
#    by weighting samples by the likelihood.
#
# 2. uniform_plus_SS_weights
#    Draws wc and wa from independent Uniform(0,1) priors, evaluates the
#    Lu et al. SS prior density at each sampled point, and computes posterior
#    predictions using likelihood * SS-prior weights (importance sampling).
#
# 3. grid 
#.   This is one of two implementations used by Lu et al. (2008). 
#    Matlab code available at: 
#    Evaluates the SS prior and likelihood on a deterministic grid over
#    wc, wa and computes posterior predictions as weighted
#    averages over the grid. This is deterministic numerical integration
#    and serves as a reference check for the Monte Carlo methods.
#
# The three methods estimate the same posterior quantities. 
# Differences should be small and reflect Monte Carlo error or grid resolution.
# ============================================================

# -----------------------------
# 1. Settings
# -----------------------------
library(tidyverse)
set.seed(0815)

alpha <- 5
m <- 1e6
n_grid <- 501

# Example data rows in SI convention:
#   N_00 = N(C=0,E=0)
#   N_01 = N(C=0,E=1)
#   N_10 = N(C=1,E=0)
#   N_11 = N(C=1,E=1)
# These include zero-contingency and positive generative cases.

data_examples <- rbind(
  zero_low_density  = c(N_00 = 15, N_01 = 5,  N_10 = 15, N_11 = 5),
  zero_high_density = c(N_00 = 5,  N_01 = 15, N_10 = 5,  N_11 = 15),
  gen_moderate      = c(N_00 = 15, N_01 = 5,  N_10 = 8,  N_11 = 12),
  gen_strong        = c(N_00 = 15, N_01 = 5,  N_10 = 2,  N_11 = 18),
  ceiling_zero      = c(N_00 = 0,  N_01 = 20, N_10 = 0,  N_11 = 20)
)


# -----------------------------
# 2. Truncated exponential sampler
# -----------------------------
# Draws x in [0,1] from density proportional to exp(-rate*x).
#   rate > 0: mass near 0
#   rate < 0: mass near 1
#   rate = 0: uniform

rtrunc_exp01 <- function(n, rate) {
  u <- runif(n)

  if (abs(rate) < 1e-12) {
    return(u)
  }

  -log(1 - u * (1 - exp(-rate))) / rate
}


# -----------------------------
# 3. Direct sampler from generative SS prior
# -----------------------------
# Eq. 10 / Appendix B Eq. B1 as a two-component mixture.

rSS_gen_direct <- function(n, alpha = 5) {
  out <- data.frame(
    wc = numeric(n),
    wa = numeric(n),
    component = integer(n)
  )

  comp1 <- runif(n) < 0.5
  out$component[comp1] <- 1L
  out$component[!comp1] <- 2L

  # Component 1: exp[-alpha*wa - alpha*(1-wc)]
  # Peak at wa = 0, wc = 1
  out$wa[comp1] <- rtrunc_exp01(sum(comp1),  alpha)
  out$wc[comp1] <- rtrunc_exp01(sum(comp1), -alpha)

  # Component 2: exp[-alpha*(1-wa) - alpha*wc]
  # Peak at wa = 1, wc = 0
  out$wa[!comp1] <- rtrunc_exp01(sum(!comp1), -alpha)
  out$wc[!comp1] <- rtrunc_exp01(sum(!comp1),  alpha)

  out
}


# -----------------------------
# 4. Uniform proposal + SS prior weight
# -----------------------------
# This implements the SS equation literally as an importance weight.

SS_gen_weight <- function(wc, wa, alpha = 5) {
  exp(-alpha * wa - alpha * (1 - wc)) +
    exp(-alpha * (1 - wa) - alpha * wc)
}

rSS_gen_uniform_weighted <- function(n, alpha = 5) {
  out <- data.frame(
    wc = runif(n),
    wa = runif(n)
  )

  out$prior_weight <- SS_gen_weight(out$wc, out$wa, alpha = alpha)
  out
}

weighted_mean <- function(x, w) {
  sum(x * w) / sum(w)
}

weighted_var <- function(x, w) {
  mu <- weighted_mean(x, w)
  sum(w * (x - mu)^2) / sum(w)
}

weighted_cor <- function(x, y, w) {
  mux <- weighted_mean(x, w)
  muy <- weighted_mean(y, w)
  cov_xy <- sum(w * (x - mux) * (y - muy)) / sum(w)
  cov_xy / sqrt(weighted_var(x, w) * weighted_var(y, w))
}


# -----------------------------
# 5. Likelihood and posterior summaries
# -----------------------------
# Conditional likelihood of E given C under S1/noisy-OR.
# bc is omitted because it cancels for posterior summaries of wc, wa
# when P(bc) factorizes from P(wc,wa).

likelihood_E_given_C <- function(data, wc, wa) {
  N_00 <- data[1]
  N_01 <- data[2]
  N_10 <- data[3]
  N_11 <- data[4]

  p_E_noC <- wa
  p_E_C <- wc + wa - wc * wa

  # Clamp to avoid numerical issues at exact 0/1.
  eps <- .Machine$double.xmin
  p_E_noC <- pmin(pmax(p_E_noC, eps), 1 - eps)
  p_E_C <- pmin(pmax(p_E_C, eps), 1 - eps)

  (1 - p_E_noC)^N_00 *
    p_E_noC^N_01 *
    (1 - p_E_C)^N_10 *
    p_E_C^N_11
}

posterior_from_direct_prior <- function(data, prior_draws) {
  lik <- likelihood_E_given_C(data, prior_draws$wc, prior_draws$wa)

  c(
    wc_pp = weighted_mean(prior_draws$wc, lik),
    wa_pp = weighted_mean(prior_draws$wa, lik),
    pEC_pp = weighted_mean(prior_draws$wc + prior_draws$wa - prior_draws$wc * prior_draws$wa, lik),
    pEnoC_pp = weighted_mean(prior_draws$wa, lik)
  )
}

posterior_from_uniform_weighted <- function(data, uniform_draws) {
  lik <- likelihood_E_given_C(data, uniform_draws$wc, uniform_draws$wa)
  w <- lik * uniform_draws$prior_weight

  c(
    wc_pp = weighted_mean(uniform_draws$wc, w),
    wa_pp = weighted_mean(uniform_draws$wa, w),
    pEC_pp = weighted_mean(uniform_draws$wc + uniform_draws$wa - uniform_draws$wc * uniform_draws$wa, w),
    pEnoC_pp = weighted_mean(uniform_draws$wa, w)
  )
}

posterior_from_grid <- function(data, alpha = 5, n_grid = 501) {
  wc_values <- seq(0, 1, length.out = n_grid)
  wa_values <- seq(0, 1, length.out = n_grid)

  grid <- expand.grid(wc = wc_values, wa = wa_values)

  prior_w <- SS_gen_weight(grid$wc, grid$wa, alpha = alpha)
  lik <- likelihood_E_given_C(data, grid$wc, grid$wa)
  w <- prior_w * lik

  c(
    wc_pp = weighted_mean(grid$wc, w),
    wa_pp = weighted_mean(grid$wa, w),
    pEC_pp = weighted_mean(grid$wc + grid$wa - grid$wc * grid$wa, w),
    pEnoC_pp = weighted_mean(grid$wa, w)
  )
}


# -----------------------------
# 6. Prior sanity checks
# -----------------------------

cat("\n=== Prior sanity check ===\n")
cat("alpha =", alpha, "\n")
cat("m =", format(m, scientific = FALSE), "\n\n")

direct_prior <- rSS_gen_direct(m, alpha = alpha)
uniform_weighted <- rSS_gen_uniform_weighted(m, alpha = alpha)

prior_check <- rbind(
  direct_SS_prior = c(
    mean_wc = mean(direct_prior$wc),
    mean_wa = mean(direct_prior$wa),
    cor_wc_wa = cor(direct_prior$wc, direct_prior$wa),
    q10_wc = unname(quantile(direct_prior$wc, .10)),
    q50_wc = unname(quantile(direct_prior$wc, .50)),
    q90_wc = unname(quantile(direct_prior$wc, .90)),
    q10_wa = unname(quantile(direct_prior$wa, .10)),
    q50_wa = unname(quantile(direct_prior$wa, .50)),
    q90_wa = unname(quantile(direct_prior$wa, .90))
  ),
  uniform_plus_SS_weights = c(
    mean_wc = weighted_mean(uniform_weighted$wc, uniform_weighted$prior_weight),
    mean_wa = weighted_mean(uniform_weighted$wa, uniform_weighted$prior_weight),
    cor_wc_wa = weighted_cor(uniform_weighted$wc, uniform_weighted$wa, uniform_weighted$prior_weight),
    q10_wc = NA_real_,
    q50_wc = NA_real_,
    q90_wc = NA_real_,
    q10_wa = NA_real_,
    q50_wa = NA_real_,
    q90_wa = NA_real_
  )
)

print(round(prior_check, 4))

# -----------------------------
# 7. Posterior sanity checks
# -----------------------------
# Comparing direct SS-prior sampler, uniform+SS weights, and deterministic grid.

results <- list()

for (i in seq_len(nrow(data_examples))) {
  row_name <- rownames(data_examples)[i]
  data_i <- data_examples[i, ]

  direct_i <- posterior_from_direct_prior(data_i, direct_prior)
  weighted_i <- posterior_from_uniform_weighted(data_i, uniform_weighted)
  grid_i <- posterior_from_grid(data_i, alpha = alpha, n_grid = n_grid)

  results[[row_name]] <- rbind(
    direct_SS_prior = direct_i,
    uniform_plus_SS_weights = weighted_i,
    grid = grid_i,
    direct_minus_weighted = direct_i - weighted_i,
    direct_minus_grid = direct_i - grid_i
  )
}

for (nm in names(results)) {
  cat("\n---", nm, "---\n")
  print(round(results[[nm]], 4))
}


# -----------------------------
# 8. visual sanity check for truncated sampler
# -----------------------------
# plot samples from truncated mixtures; should show two prior clouds: 
# one near (wa=0,wc=1) and 
# one near (wa=1,wc=0).

idx <- round(seq(1, m, length.out = 10000))

prior_cloud <- data.frame(
  wa = direct_prior$wa[idx],
  wc = direct_prior$wc[idx]
)

ggplot(prior_cloud, aes(x = wa, y = wc)) +
  geom_point(size = 0.6, alpha = 0.5) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  labs(
    x = "wa = background strength",
    y = "wc = candidate strength",
    title = "Direct samples from generative SS prior"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = "plots/ss_truncated_sanity_check_prior_cloud.png",
  dpi=300,
  width = 6,
  height = 6
)

