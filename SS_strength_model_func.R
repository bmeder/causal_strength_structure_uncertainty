# Lu et al. (2008) / Meder et al. (2014) SS strength model ------------------
# Direct-sampling implementation with generative/preventive model_type argument.
#
# Scope:
#   - Structure S1 only.
#   - No S0 branch, no structure posterior, no Bayesian model averaging.
#   - Main prediction is posterior mean strength:
#       E[wc | D, S1, SS prior]
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
#        - if model_type = "noisy_or_generative", wc is generative strength
#        - if model_type = "noisy_and_not_preventive", wc is preventive strength
#   wa = background-cause generative strength
#
# Paper mapping:
#   Lu et al. (2008), Eq. 10, p. 959:
#     Generative SS prior over w0,w1.
#     In this code: w0 = wa, w1 = wc.
#
#   Lu et al. (2008), Eq. 11, p. 959:
#     Preventive SS prior over w0,w1.
#     In this code: w0 = wa, w1 = wc.
#
#   Meder et al. (2014), Appendix B, Eq. B1/B3:
#     Direct Monte Carlo with (wc, wa) drawn from the joint SS prior and
#     bc drawn from Beta(1,1). Meder et al. average P(c|e;theta); the
#     strength version here averages wc.
#
# Numerical method:
#   Direct prior-sampling Monte Carlo.
#   - Draw bc ~ Beta(1,1).
#   - Draw (wc, wa) directly from the model_type-specific joint SS prior.
#   - Compute posterior summaries by weighting samples by the S1 likelihood.



# ---------------------------------------------------------------------------
# S1 likelihood functions
# ---------------------------------------------------------------------------

SS_likelihood_data_S1_gen <- function(data, theta) {
  # Generative S1: noisy-OR.
  #
  # P(E=1 | C=0) = wa
  # P(E=1 | C=1) = wc + wa - wc*wa
  
  bc <- theta[, "bc"]
  wc <- theta[, "wc"]
  wa <- theta[, "wa"]
  
  pEC <- wc + wa - wc * wa
  
  ((1 - bc) * (1 - wa))^data[1] *
    ((1 - bc) * wa)^data[2] *
    (bc * (1 - pEC))^data[3] *
    (bc * pEC)^data[4]
}


SS_likelihood_data_S1_prev <- function(data, theta) {
  # Preventive S1: noisy-AND-NOT.
  #
  # P(E=1 | C=0) = wa
  # P(E=1 | C=1) = wa * (1 - wc)
  
  bc <- theta[, "bc"]
  wc <- theta[, "wc"]
  wa <- theta[, "wa"]
  
  pEC <- wa * (1 - wc)
  
  ((1 - bc) * (1 - wa))^data[1] *
    ((1 - bc) * wa)^data[2] *
    (bc * (1 - pEC))^data[3] *
    (bc * pEC)^data[4]
}


# ---------------------------------------------------------------------------
# Conditional probabilities
# ---------------------------------------------------------------------------

SS_p_EC_gen <- function(theta) {
  # Generative: P(E=1 | C=1, theta, S1)
  theta[, "wc"] + theta[, "wa"] - theta[, "wc"] * theta[, "wa"]
}


SS_p_EC_prev <- function(theta) {
  # Preventive: P(E=1 | C=1, theta, S1)
  theta[, "wa"] * (1 - theta[, "wc"])
}


SS_p_EnoC <- function(theta) {
  # Both directions: P(E=1 | C=0, theta, S1)
  theta[, "wa"]
}


SS_p_CE_gen <- function(theta) {
  # Generative diagnostic probability:
  # P(C=1 | E=1, theta, S1)
  
  pEC <- SS_p_EC_gen(theta)
  pEnoC <- SS_p_EnoC(theta)
  
  (pEC * theta[, "bc"]) /
    (pEC * theta[, "bc"] + pEnoC * (1 - theta[, "bc"]))
}


SS_p_CE_prev <- function(theta) {
  # Preventive diagnostic probability:
  # P(C=1 | E=1, theta, S1)
  
  pEC <- SS_p_EC_prev(theta)
  pEnoC <- SS_p_EnoC(theta)
  
  (pEC * theta[, "bc"]) /
    (pEC * theta[, "bc"] + pEnoC * (1 - theta[, "bc"]))
}


# ---------------------------------------------------------------------------
# Posterior means under S1
# ---------------------------------------------------------------------------

SS_S1_wc_pp <- function(S1_lik, theta) {
  sum(S1_lik * theta[, "wc"]) / sum(S1_lik)
}


SS_S1_wa_pp <- function(S1_lik, theta) {
  sum(S1_lik * theta[, "wa"]) / sum(S1_lik)
}


SS_S1_bc_pp <- function(S1_lik, theta) {
  sum(S1_lik * theta[, "bc"]) / sum(S1_lik)
}


SS_S1_pEC_pp <- function(S1_lik, pEC) {
  sum(S1_lik * pEC) / sum(S1_lik)
}


SS_S1_pEnoC_pp <- function(S1_lik, pEnoC) {
  sum(S1_lik * pEnoC) / sum(S1_lik)
}


SS_S1_pCE_pp <- function(S1_lik, pCE) {
  sum(S1_lik * pCE) / sum(S1_lik)
}


# ---------------------------------------------------------------------------
# Direct sampling from model_type-specific SS priors
# ---------------------------------------------------------------------------

SS_rtrunc_exp01 <- function(n, alpha) {
  # Draw x in [0,1] from density proportional to exp(-alpha * x).
  #
  # alpha > 0: values near 0.
  # alpha < 0: values near 1.
  # alpha = 0: uniform on [0,1].
  #
  # the function is called with either:
  #   alpha  =  alpha_SS   # draw near 0
  #   alpha  = -alpha_SS   # draw near 1
  
  u <- runif(n)
  
  if (abs(alpha) < 1e-12) {
    return(u)
  }
  
  -log(1 - u * (1 - exp(-alpha))) / alpha
}


SS_draw_direct <- function(m, 
                           alpha = 5,
                           model_type = c("noisy_or_generative", "noisy_and_not_preventive"),
                           mixture_weight = 0.5) # default, both components weighted equally
  {
  # Draw theta = (bc, wc, wa).
  #
  # bc:
  #   Flat Beta(1,1) prior.
  #
  # wc, wa:
  #   Drawn from the model_type-specific sparse-and-strong prior.
  #
  # Generative SS prior:
  #   p(wc, wa) ∝ exp[-alpha*wa - alpha*(1-wc)] +
  #               exp[-alpha*(1-wa) - alpha*wc]
  #   Peaks: (wa=0,wc=1) and (wa=1,wc=0).
  #
  # Preventive SS prior:
  #   p(wc, wa) ∝ exp[-alpha*(1-wa) - alpha*(1-wc)] +
  #               exp[-alpha*(1-wa) - alpha*wc]
  #   Peaks: (wa=1,wc=1) and (wa=1,wc=0).
  
  model_type <- match.arg(model_type)
  
  theta <- matrix(NA_real_, nrow = m, ncol = 3)
  colnames(theta) <- c("bc", "wc", "wa")
  
  theta[, "bc"] <- rbeta(m, 1, 1)
  
  # deault in Lu et al's 2008 SS priors are equal-mass two-component mixtures
  # adjust to include arbitrary mixture weights
  comp1 <- runif(m) < mixture_weight
  
  if (model_type == "noisy_or_generative") {
    
    # Component 1: wa near 0, wc near 1.
    theta[comp1, "wa"] <- SS_rtrunc_exp01(sum(comp1),  alpha)
    theta[comp1, "wc"] <- SS_rtrunc_exp01(sum(comp1), -alpha)
    
    # Component 2: wa near 1, wc near 0.
    theta[!comp1, "wa"] <- SS_rtrunc_exp01(sum(!comp1), -alpha)
    theta[!comp1, "wc"] <- SS_rtrunc_exp01(sum(!comp1),  alpha)
    
  } else {
    
    # Component 1: wa near 1, wc near 1.
    theta[comp1, "wa"] <- SS_rtrunc_exp01(sum(comp1), -alpha)
    theta[comp1, "wc"] <- SS_rtrunc_exp01(sum(comp1), -alpha)
    
    # Component 2: wa near 1, wc near 0.
    theta[!comp1, "wa"] <- SS_rtrunc_exp01(sum(!comp1), -alpha)
    theta[!comp1, "wc"] <- SS_rtrunc_exp01(sum(!comp1),  alpha)
  }
  
  theta
}


# ---------------------------------------------------------------------------
# Main function: SS strength predictions
# ---------------------------------------------------------------------------

generate_SS_strength_preds <- function(data, 
                                       m = 1e6, 
                                       alpha = 5,
                                       model_type = c("noisy_or_generative", "noisy_and_not_preventive"),
                                       mixture_weight = 0.5) {
  model_type <- match.arg(model_type)
  
  data <- as.matrix(data)
  m <- as.integer(m)
  
  pred_cols <- c(
    "N", "N_00", "N_01", "N_10", "N_11",
    "deltaP_MLE", "wc_MLE", "bc_MLE", "wa_MLE", "pEC_MLE", "pCE_MLE",
    "SS_direction", "SS_method", 
    "SS_alpha", "SS_mixture_weight",
    "SS_marglik", "SS_logmarglik",
    "SS_wc_pp", "SS_wa_pp", "SS_bc_pp",
    "SS_pEC", "SS_pEnoC", "SS_pCE"
  )
  
  pred <- matrix(
    data = NA,
    nrow = nrow(data),
    ncol = length(pred_cols),
    dimnames = list(seq_len(nrow(data)), pred_cols)
  )
  
  pred <- as.data.frame(pred)
  
  
  # -------------------------------------------------------------------------
  # Observed frequencies
  # -------------------------------------------------------------------------
  
  pred$N <- rowSums(data)
  
  pred$N_00 <- data[, 1]
  pred$N_01 <- data[, 2]
  pred$N_10 <- data[, 3]
  pred$N_11 <- data[, 4]
  
  # -------------------------------------------------------------------------
  # MLE summaries
  # -------------------------------------------------------------------------
  
  pred$bc_MLE   <- (data[,3] + data[,4]) /(data[,1] + data[,2] + data[,3] + data[,4]) # base rate of cause C, P(c)
  pred$wa_MLE   <- data[,2]/(data[,1]+data[,2]) # strength of background cause, equals the likelihood of the effect if the cause is absent, P(e|no c)
  pred$pEC_MLE  <- data[,4]/(data[,3]+data[,4]) # likelihood of the effect if the cause is present
  pred$pCE_MLE  <- data[,4]/(data[,2]+data[,4]) # probability of the cause if the effect is present
  
  # maximum likelihood estimate of contingency delta P
  pred$deltaP_MLE <- pred$pEC_MLE - pred$wa_MLE
  
  # Descriptive causal-power MLE.
  # Positive deltaP: generative power.
  # Negative deltaP: preventive power.
  # Zero deltaP: zero strength.
  pred$wc_MLE <- ifelse(
    is.nan(pred$deltaP_MLE), NaN,
    ifelse(
      pred$deltaP_MLE == 0, 0,
      ifelse(
        pred$deltaP_MLE > 0,
        pred$deltaP_MLE / (1 - pred$wa_MLE),
        -pred$deltaP_MLE / pred$wa_MLE
      )
    )
  )
  
  # -------------------------------------------------------------------------
  # Model metadata
  # -------------------------------------------------------------------------
  
  pred$SS_direction <- model_type
  pred$SS_alpha <- alpha
  pred$SS_mixture_weight <- mixture_weight
  pred$SS_prior_mean_wc <-
    mixture_weight * SS_component_mean_high(alpha) +
    (1 - mixture_weight) * SS_component_mean_low(alpha)
  pred$SS_method <- "direct sampling"
  
  # -------------------------------------------------------------------------
  # Draw theta once and reuse for all data sets
  # -------------------------------------------------------------------------
  
  theta_S1 <- SS_draw_direct(
    m = m,
    alpha = alpha,
    model_type = model_type,
    mixture_weight = mixture_weight
  )
  # Theta-level quantities that do not depend on the data 
  if (model_type == "noisy_or_generative") {
    pEC <- SS_p_EC_gen(theta_S1)
    pCE <- SS_p_CE_gen(theta_S1)
  } else {
    pEC <- SS_p_EC_prev(theta_S1)
    pCE <- SS_p_CE_prev(theta_S1)
  }
  
  pEnoC <- SS_p_EnoC(theta_S1)
  
  # -------------------------------------------------------------------------
  # posterior summaries for each data set
  # -------------------------------------------------------------------------
  
  for (i in seq_len(nrow(data))) {
    
    if (model_type == "noisy_or_generative") {
      S1_lik <- SS_likelihood_data_S1_gen(data[i, ], theta_S1)
    } else {
      S1_lik <- SS_likelihood_data_S1_prev(data[i, ], theta_S1)
    }
    
    # Monte Carlo marginal likelihood estimate:
    # E_prior[P(D | theta, S1)] ≈ mean(P(D | theta_k, S1)).
    pred$SS_marglik[i] <- mean(S1_lik)
    pred$SS_logmarglik[i] <- log(pred$SS_marglik[i])
    
    # Posterior parameter means
    pred$SS_wc_pp[i] <- sum(S1_lik * theta_S1[, "wc"]) / sum(S1_lik)
    pred$SS_wa_pp[i] <- sum(S1_lik * theta_S1[, "wa"]) / sum(S1_lik)
    pred$SS_bc_pp[i] <- sum(S1_lik * theta_S1[, "bc"]) / sum(S1_lik)
    
    # Posterior conditional-probability means
    pred$SS_pEC[i] <- sum(S1_lik * pEC) / sum(S1_lik)
    pred$SS_pEnoC[i] <- sum(S1_lik * pEnoC) / sum(S1_lik)
    pred$SS_pCE[i] <- sum(S1_lik * pCE) / sum(S1_lik)
  }
  
  pred
}

# -------------------------------------------------------------------------
# SS prior component means
# -------------------------------------------------------------------------

# Mean of the truncated exponential component concentrated near zero:
# density proportional to exp(-alpha * x), 0 <= x <= 1
SS_component_mean_low <- function(alpha = 5) {
  
  if (!is.numeric(alpha) || length(alpha) != 1L ||
      !is.finite(alpha) || alpha < 0) {
    stop("alpha must be a single finite non-negative number.")
  }
  
  if (alpha < 1e-12) {
    return(0.5)
  }
  
  1 / alpha - 1 / expm1(alpha)
}


# Mean of the mirrored component concentrated near one
SS_component_mean_high <- function(alpha = 5) {
  1 - SS_component_mean_low(alpha)
}

SS_component_mean_low(alpha = 5)
SS_component_mean_high(alpha = 5)


# ---------------------------------------------------------------------------
# Example use
# ---------------------------------------------------------------------------

# data <- matrix(data = 0, nrow = 5, ncol = 4)
# data[1, ] <- c(5, 15, 5, 15)    # zero contingency, high effect density
# data[2, ] <- c(15, 5, 15, 5)    # zero contingency, low effect density
# data[3, ] <- c(10, 10, 0, 20)   # generative example
# data[4, ] <- c(20, 0, 10, 10)   # generative with low background
# data[5, ] <- c(5, 15, 15, 5)    # preventive example
# 
# set.seed(1)
# ss_gen <- generate_SS_strength_preds(
#   data = data,
#   m = 10^5,
#   alpha = 5,
#   model_type = "noisy_or_generative"
# )
# 
# set.seed(1)
# ss_prev <- generate_SS_strength_preds(
#   data = data,
#   m = 10^5,
#   alpha = 5,
#   model_type = "noisy_and_not_preventive"
# )