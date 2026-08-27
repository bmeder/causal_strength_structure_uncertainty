# Lu et al. (2008) / Meder et al. (2014) SS strength model ------------------
# Generative, strength-only implementation structured to parallel the SI code.
#
# Scope:
#   - Structure S1 only: C and background A can generate E by a noisy-OR.
#   - No S0 branch, no structure posterior, no Bayesian model averaging.
#   - Main prediction is the posterior mean of wc, E[wc | D, S1, SS prior].
#
# Paper mapping:
#   - Lu et al. (2008), Eq. 10, p. 959:
#       Generative sparse-and-strong (SS) prior over w0,w1.
#       In SI/Meder notation: w0 = wa, w1 = wc.
#
#   - Meder, Mayrhofer, & Waldmann (2014), Appendix B, Eq. B1, p. 300:
#       Same generative SS prior over wc, wa.
#
#   - Meder et al. (2014), Appendix B, Eq. B3, p. 301:
#       Monte Carlo approximation for the Bayesian Power PC / SS model.
#       Their diagnostic prediction averages P(c|e; theta); the strength
#       version below averages wc instead.
#
# Input data convention, identical to generate_SI_preds():
#   data[,1] = N_00 = N(C=0, E=0)
#   data[,2] = N_01 = N(C=0, E=1)
#   data[,3] = N_10 = N(C=1, E=0)
#   data[,4] = N_11 = N(C=1, E=1)
#
# Parameters:
#   bc = base rate of candidate cause C
#   wc = candidate-cause strength
#   wa = background-cause strength
#
# Three numerical implementations are available:
#
#   method = "importance"  [default]
#     Draws bc, wc, wa from independent Uniform(0,1) priors, evaluates the
#     SS prior density at each sampled (wc, wa), and computes posterior
#     predictions using likelihood * SS-prior weights. Technically, this is
#     self-normalized importance sampling with a uniform proposal.
#
#   method = "direct"
#     Draws (wc, wa) directly from the joint SS prior using a two-component
#     truncated-exponential mixture, and draws bc from Uniform(0,1). Posterior
#     predictions are computed by weighting samples by likelihood only. This is
#     closest to the wording of Meder et al. Appendix B, Eq. B3.
#
#   method = "grid"
#     Evaluates the prior and likelihood on a deterministic grid over
#     bc, wc, wa in [0,1]^3, then computes posterior predictions as weighted
#     averages over the grid. This is deterministic numerical integration and
#     is mainly intended as a validation/reference method. It can be slow.
#
# Recommended use:
#   - Use method = "importance" for the main SI-style model-comparison table.
#   - Use method = "direct" or method = "grid" as validation checks.


# ---------------------------------------------------------------------------
# Basic checks and helpers
# ---------------------------------------------------------------------------

check_SS_data <- function(data) {
  if (!is.matrix(data) && !is.data.frame(data)) {
    stop("data must be a matrix or data.frame with columns N_00, N_01, N_10, N_11.")
  }
  data <- as.matrix(data)
  if (ncol(data) != 4) {
    stop("data must have exactly 4 columns: N_00, N_01, N_10, N_11.")
  }
  if (any(data < 0, na.rm = TRUE)) {
    stop("data contains negative counts.")
  }
  data
}

weighted_mean_safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w >= 0
  if (!any(ok)) return(NA_real_)
  sw <- sum(w[ok])
  if (!is.finite(sw) || sw <= 0) return(NA_real_)
  sum(x[ok] * w[ok]) / sw
}

importance_ess <- function(w) {
  ok <- is.finite(w) & w >= 0
  if (!any(ok)) return(NA_real_)
  sw <- sum(w[ok])
  sw2 <- sum(w[ok]^2)
  if (!is.finite(sw) || !is.finite(sw2) || sw2 <= 0) return(NA_real_)
  sw^2 / sw2
}


# ---------------------------------------------------------------------------
# MLE summaries, copied conceptually from the SI code
# ---------------------------------------------------------------------------

generate_MLE_est_SS <- function(data) {
  data <- check_SS_data(data)

  est <- matrix(
    data = 0,
    nrow = nrow(data),
    ncol = 10,
    dimnames = list(
      seq_len(nrow(data)),
      c("N_00", "N_01", "N_10", "N_11",
        "deltaP_MLE", "wc_MLE", "bc_MLE", "wa_MLE",
        "pEC_MLE", "pCE_MLE")
    )
  )
  est <- as.data.frame(est)

  est$N_00 <- data[, 1]
  est$N_01 <- data[, 2]
  est$N_10 <- data[, 3]
  est$N_11 <- data[, 4]

  est$bc_MLE  <- rowSums(data[, 3:4, drop = FALSE]) / rowSums(data)
  est$wa_MLE  <- data[, 2] / (data[, 1] + data[, 2])
  est$pEC_MLE <- data[, 4] / (data[, 3] + data[, 4])
  est$pCE_MLE <- data[, 4] / (data[, 2] + data[, 4])

  est$deltaP_MLE <- est$pEC_MLE - est$wa_MLE

  # Generative causal power. For negative contingencies this function returns
  # the preventive power estimate, as the SI helper does. The SS model below is
  # still generative-only; this MLE column is only a descriptive summary.
  est$wc_MLE <- ifelse(
    is.nan(est$deltaP_MLE), NaN,
    ifelse(
      est$deltaP_MLE == 0, 0,
      ifelse(
        est$deltaP_MLE > 0,
        est$deltaP_MLE / (1 - est$wa_MLE),
        -est$deltaP_MLE / est$wa_MLE
      )
    )
  )

  est
}


# ---------------------------------------------------------------------------
# Likelihood and derived quantities under S1
# ---------------------------------------------------------------------------

likelihood_data_SS_S1 <- function(data, theta) {
  # Same S1 noisy-OR likelihood as in the SI code.
  # theta[,1] = bc, theta[,2] = wc, theta[,3] = wa.
  ((1 - theta[, 1]) * (1 - theta[, 3]))^data[1] *
    ((1 - theta[, 1]) * theta[, 3])^data[2] *
    (theta[, 1] * (1 - theta[, 2]) * (1 - theta[, 3]))^data[3] *
    ((theta[, 2] + theta[, 3] - theta[, 2] * theta[, 3]) * theta[, 1])^data[4]
}

SS_p_EC_theta <- function(theta) {
  theta[, 2] + theta[, 3] - theta[, 2] * theta[, 3]
}

SS_p_EnoC_theta <- function(theta) {
  theta[, 3]
}

SS_p_CE_theta <- function(theta) {
  pEC <- SS_p_EC_theta(theta)
  pEnoC <- SS_p_EnoC_theta(theta)
  num <- pEC * theta[, 1]
  den <- num + pEnoC * (1 - theta[, 1])
  out <- num / den
  out[!is.finite(out)] <- NA_real_
  out
}


# ---------------------------------------------------------------------------
# SS prior density / weight
# ---------------------------------------------------------------------------

SS_prior_weight_gen <- function(theta, alpha = 5) {
  # Lu et al. Eq. 10 / Meder et al. Appendix B Eq. B1.
  # In this code: wc = theta[,2], wa = theta[,3].
  wc <- theta[, 2]
  wa <- theta[, 3]

  exp(-alpha * wa - alpha * (1 - wc)) +
    exp(-alpha * (1 - wa) - alpha * wc)
}


# ---------------------------------------------------------------------------
# Method 1: self-normalized importance sampling with uniform proposal
# ---------------------------------------------------------------------------

draw_theta_SS_importance <- function(m) {
  theta <- matrix(data = 0, nrow = m, ncol = 3)
  colnames(theta) <- c("bc", "wc", "wa")

  theta[, 1] <- rbeta(m, 1, 1)  # bc
  theta[, 2] <- rbeta(m, 1, 1)  # wc
  theta[, 3] <- rbeta(m, 1, 1)  # wa

  theta
}


# ---------------------------------------------------------------------------
# Method 2: direct sampling from the joint SS prior
# ---------------------------------------------------------------------------

rtrunc_exp01 <- function(n, rate) {
  # Draw from density proportional to exp(-rate*x), truncated to x in [0,1].
  # rate > 0 favors x near 0; rate < 0 favors x near 1.
  u <- runif(n)
  if (abs(rate) < 1e-12) return(u)
  -log(1 - u * (1 - exp(-rate))) / rate
}

draw_theta_SS_direct <- function(m, alpha = 5) {
  theta <- matrix(data = 0, nrow = m, ncol = 3)
  colnames(theta) <- c("bc", "wc", "wa")

  theta[, 1] <- rbeta(m, 1, 1)  # bc; flat prior as in Meder Appendix B

  comp1 <- runif(m) < 0.5

  # Component 1: exp[-alpha*wa - alpha*(1-wc)]
  # peak at wa = 0, wc = 1.
  theta[comp1, 3] <- rtrunc_exp01(sum(comp1),  alpha)  # wa near 0
  theta[comp1, 2] <- rtrunc_exp01(sum(comp1), -alpha)  # wc near 1

  # Component 2: exp[-alpha*(1-wa) - alpha*wc]
  # peak at wa = 1, wc = 0.
  theta[!comp1, 3] <- rtrunc_exp01(sum(!comp1), -alpha) # wa near 1
  theta[!comp1, 2] <- rtrunc_exp01(sum(!comp1),  alpha) # wc near 0

  theta
}


# ---------------------------------------------------------------------------
# Method 3: deterministic grid integration
# ---------------------------------------------------------------------------

draw_theta_SS_grid <- function(n_grid = 101) {
  if (n_grid > 151) {
    warning("Grid over bc,wc,wa has n_grid^3 points; large n_grid values can be slow and memory-intensive.")
  }

  grid_values <- seq(0, 1, length.out = n_grid)

  theta <- expand.grid(
    bc = grid_values,
    wc = grid_values,
    wa = grid_values
  )

  as.matrix(theta)
}


# ---------------------------------------------------------------------------
# Posterior summaries for one data row
# ---------------------------------------------------------------------------

posterior_SS_row <- function(data_row, theta, alpha = 5,
                             method = c("importance", "direct", "grid")) {
  method <- match.arg(method)

  lik <- likelihood_data_SS_S1(data_row, theta)

  if (method == "importance" || method == "grid") {
    prior_w <- SS_prior_weight_gen(theta, alpha = alpha)
    post_w <- lik * prior_w
  } else {
    # theta was already drawn from the SS prior.
    prior_w <- rep(1, nrow(theta))
    post_w <- lik
  }

  pEC <- SS_p_EC_theta(theta)
  pEnoC <- SS_p_EnoC_theta(theta)
  pCE <- SS_p_CE_theta(theta)

  c(
    SS_pCE = weighted_mean_safe(pCE, post_w),
    SS_pEC = weighted_mean_safe(pEC, post_w),
    SS_pEnoC = weighted_mean_safe(pEnoC, post_w),
    SS_wc_pp = weighted_mean_safe(theta[, 2], post_w),
    SS_bc_pp = weighted_mean_safe(theta[, 1], post_w),
    SS_wa_pp = weighted_mean_safe(theta[, 3], post_w),
    SS_mean_pCE = weighted_mean_safe(pCE, post_w),
    SS_mean_pEC = weighted_mean_safe(pEC, post_w),
    SS_mean_pEnoC = weighted_mean_safe(pEnoC, post_w),
    SS_mean_wc = weighted_mean_safe(theta[, 2], post_w),
    SS_mean_bc = weighted_mean_safe(theta[, 1], post_w),
    SS_mean_wa = weighted_mean_safe(theta[, 3], post_w),
    SS_weight_sum = sum(post_w, na.rm = TRUE),
    SS_weight_ess = importance_ess(post_w)
  )
}


# ---------------------------------------------------------------------------
# Main function
# ---------------------------------------------------------------------------

generate_SS_strength_preds <- function(data, m = 1e6, alpha = 5,
                                       method = c("importance", "direct", "grid"),
                                       n_grid = 101) {
  method <- match.arg(method)
  data <- check_SS_data(data)

  if (method %in% c("importance", "direct")) {
    if (missing(m) || is.null(m) || length(m) != 1L || !is.finite(m) || m <= 0) {
      stop("For method = 'importance' or 'direct', m must be a positive finite number.")
    }
    m <- as.integer(m)
  }

  pred_cols <- c(
    "N", "N_00", "N_01", "N_10", "N_11",
    "deltaP_MLE", "wc_MLE", "bc_MLE", "wa_MLE", "pEC_MLE", "pCE_MLE",
    "SS_method", "SS_alpha", "SS_m", "SS_n_grid",
    "SS_pCE", "SS_pEC", "SS_pEnoC",
    "SS_wc_pp", "SS_bc_pp", "SS_wa_pp",
    "SS_mean_pCE", "SS_mean_pEC", "SS_mean_pEnoC",
    "SS_mean_wc", "SS_mean_bc", "SS_mean_wa",
    "SS_weight_sum", "SS_weight_ess"
  )

  pred <- as.data.frame(matrix(NA, nrow = nrow(data), ncol = length(pred_cols)))
  names(pred) <- pred_cols

  # Observed frequencies
  pred$N <- rowSums(data)
  pred$N_00 <- data[, 1]
  pred$N_01 <- data[, 2]
  pred$N_10 <- data[, 3]
  pred$N_11 <- data[, 4]

  # MLE summaries
  mle <- generate_MLE_est_SS(data)
  pred$deltaP_MLE <- mle$deltaP_MLE
  pred$wc_MLE <- mle$wc_MLE
  pred$bc_MLE <- mle$bc_MLE
  pred$wa_MLE <- mle$wa_MLE
  pred$pEC_MLE <- mle$pEC_MLE
  pred$pCE_MLE <- mle$pCE_MLE

  pred$SS_method <- method
  pred$SS_alpha <- alpha
  pred$SS_m <- if (method == "grid") NA_integer_ else m
  pred$SS_n_grid <- if (method == "grid") n_grid else NA_integer_

  # Draw/evaluate theta once, as in the SI implementation.
  if (method == "importance") {
    theta <- draw_theta_SS_importance(m)
  } else if (method == "direct") {
    theta <- draw_theta_SS_direct(m, alpha = alpha)
  } else {
    theta <- draw_theta_SS_grid(n_grid = n_grid)
  }

  # Row-wise predictions
  for (i in seq_len(nrow(data))) {
    row_summary <- posterior_SS_row(
      data_row = data[i, ],
      theta = theta,
      alpha = alpha,
      method = method
    )

    for (nm in names(row_summary)) {
      pred[i, nm] <- row_summary[[nm]]
    }
  }

  pred
}


# ---------------------------------------------------------------------------
# Example use
# ---------------------------------------------------------------------------
# data <- matrix(data = 0, nrow = 5, ncol = 4)
# data[1, ] <- c(5, 15, 5, 15)    # zero contingency, high effect density
# data[2, ] <- c(15, 5, 15, 5)    # zero contingency, low effect density
# data[3, ] <- c(10, 10, 0, 20)   # generative moderate/strong
# data[4, ] <- c(20, 0, 10, 10)   # generative with low background
# data[5, ] <- c(0, 20, 0, 20)    # ceiling zero-contingency case
#
# set.seed(1)
# ss_imp <- generate_SS_strength_preds(data, m = 1e6, alpha = 5, method = "importance")
# ss_dir <- generate_SS_strength_preds(data, m = 1e6, alpha = 5, method = "direct")
# ss_grid <- generate_SS_strength_preds(data, alpha = 5, method = "grid", n_grid = 101)
