# Structure Induction Model -----------------------------------------------
# Meder, B., Mayrhofer, R., & Waldmann, M. R. (2014).
# Structure induction in diagnostic causal reasoning. Psychological Review, 121, 277–301.
#
# Revised 2025/26 by Björn Meder.
#
# This implementation generates predictions for a two-structure SI model:
#
#   S0 = no causal link from C to E
#   S1 = causal link from C to E
#
# The argument `model_type` specifies the parameterization of S1:
#
#   model_type = "noisy_or_generative"  -> S1 is a generative cause, using noisy-OR
#   model_type = "noisy_and_not_preventive"  -> S1 is a preventive cause, using noisy-AND-NOT

# For ΔP=0, the generative and preventive SI variants typically converge toward w_c = 0, producing qualitatively similar predictions. 
# They are not guaranteed to be numerically identical, however, because noisy-OR and noisy-AND-NOT induce different likelihoods of S1.

# The same structure prior is used in both directional versions:
#
#   P(S1) = S1_prior
#   P(S0) = 1 - S1_prior
#
# By default:
#
#   S1_prior = 0.5
#
# so that:
#
#   P(S1) = P(S0) = 0.5
#
# For model_type = "noisy_or_generative", wc is generative strength.
# For model_type = "noisy_and_not_preventive", wc is preventive strength.


###########################################
# Generating predictions
###########################################

# Number of Monte Carlo samples:
#
#   m <- 10^6

# Structure prior:
#
#   S1_prior <- 0.5
#
# This gives uniform priors over the two structures:
#
#   P(S1) = P(S0) = 0.5

###########################################
# Input data format
###########################################

# Input: `data`, a matrix or data frame of contingency tables.
# Each row is one contingency table with columns:
#
#   data[,1] = N(C = 0, E = 0)
#   data[,2] = N(C = 0, E = 1)
#   data[,3] = N(C = 1, E = 0)
#   data[,4] = N(C = 1, E = 1)

###########################################
# Example data
###########################################

# data <- matrix(data = 0, nrow = 5, ncol = 4)
# data[1,] <- c(5, 15, 5, 15)    # zero contingency, high effect density
# data[2,] <- c(15, 5, 15, 5)    # zero contingency, low effect density
# data[3,] <- c(15, 5, 5, 15)    # positive/generative contingency
# data[4,] <- c(5, 15, 15, 5)    # negative/preventive contingency
# data[5,] <- c(10, 10, 0, 20)   # deterministic positive example

###########################################
# Example predictions
###########################################

# Generative S1 vs S0:
#
# pred_gen <- generate_SI_preds(
#   data = data,
#   m = m,
#   S1_prior = 0.5,
#   model_type = "noisy_or_generative"
# )

# Preventive S1 vs S0:
#
# pred_prev <- generate_SI_preds(
#   data = data,
#   m = m,
#   S1_prior = 0.5,
#   model_type = "noisy_and_not_preventive"
# )


################################################################################
## functions implmenting the SI model
################################################################################


## Likelihood for Structure S0 given the data
## (i.e., given a specific contingency table)
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##
likelihood_data_S0 <- function(data, theta)
{
  ((1-theta[,1])*(1-theta[,3]))^data[1]*
    ((1-theta[,1])*theta[,3])^data[2]*
    (theta[,1]*(1-theta[,3]))^data[3]*
    (theta[,1]*theta[,3])^data[4]
}

## Likelihood for Structure S1 given the data
## (i.e., given a specific contingency table)
##
## noisy-OR parameterization for *generative* causes
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##
likelihood_data_S1 <- function(data, theta)
{
  ((1-theta[,1])*(1-theta[,3]))^data[1]*
    ((1-theta[,1])*theta[,3])^data[2]*
    (theta[,1]*(1-theta[,2])*(1-theta[,3]))^data[3]*
    ((theta[,2]+theta[,3]-theta[,2]*theta[,3])*theta[,1])^data[4]
  
  # bc <- theta[,1]  # base rate of candidate cause, P(C)
  # wc <- theta[,2]  # preventive strength of C
  # wa <- theta[,3]  # background generative strength
  # 
  # ((1-bc)*(1-wa))^data[1]*
  #   ((1-bc)*wa)^data[2]*
  #   (bc*(1-wc)*(1-wa))^data[3]*
  #   ((wc+wa-wc*wa)*bc)^data[4]
}

## Likelihood for Structure S1 given the data
## (i.e., given a specific contingency table)
##
## noisy-AND-NOT parameterization for *preventive* causes
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##
likelihood_data_S1_prev <- function(data, theta)
{
  bc <- theta[,1]  # base rate of candidate cause, P(C)
  wc <- theta[,2]  # preventive strength of C
  wa <- theta[,3]  # background generative strength
  
  ((1 - bc) * (1 - wa))^data[1] *
    ((1 - bc) * wa)^data[2] *
    (bc * (1 - wa * (1 - wc)))^data[3] *
    (bc * wa * (1 - wc))^data[4]
}

# -------------------------------------------------------------------------
# Log likelihoods for stable model evidence
# -------------------------------------------------------------------------
# Optional log-likelihood versions for later numerical checks.
# They are not currently used by generate_SI_preds().

likelihood_data_S0_log <- function(data, theta)
{
  ifelse(
    data[1] == 0,
    0,
    data[1] * log((1 - theta[,1]) * (1 - theta[,3]))
  ) +
    ifelse(
      data[2] == 0,
      0,
      data[2] * log((1 - theta[,1]) * theta[,3])
    ) +
    ifelse(
      data[3] == 0,
      0,
      data[3] * log(theta[,1] * (1 - theta[,3]))
    ) +
    ifelse(
      data[4] == 0,
      0,
      data[4] * log(theta[,1] * theta[,3])
    )
}


likelihood_data_S1_log <- function(data, theta)
{
  ifelse(
    data[1] == 0,
    0,
    data[1] * log((1 - theta[,1]) * (1 - theta[,3]))
  ) +
    ifelse(
      data[2] == 0,
      0,
      data[2] * log((1 - theta[,1]) * theta[,3])
    ) +
    ifelse(
      data[3] == 0,
      0,
      data[3] * log(theta[,1] * (1 - theta[,2]) * (1 - theta[,3]))
    ) +
    ifelse(
      data[4] == 0,
      0,
      data[4] * log((theta[,2] + theta[,3] - theta[,2] * theta[,3]) * theta[,1])
    )
}


likelihood_data_S1_prev_log <- function(data, theta)
{
  bc <- theta[,1]
  wc <- theta[,2]
  wa <- theta[,3]
  
  ifelse(
    data[1] == 0,
    0,
    data[1] * log((1 - bc) * (1 - wa))
  ) +
    ifelse(
      data[2] == 0,
      0,
      data[2] * log((1 - bc) * wa)
    ) +
    ifelse(
      data[3] == 0,
      0,
      data[3] * log(bc * (1 - wa * (1 - wc)))
    ) +
    ifelse(
      data[4] == 0,
      0,
      data[4] * log(bc * wa * (1 - wc))
    )
}

## Diagnostic Probability P(C|E) given Structure S0 and the data
## (i.e., given a specific contingency table)
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##
p_CE_S0 <- function(data, S0_lik, theta)
{
  pCE <- theta[,1]
  return(mean(pCE*S0_lik)/mean(S0_lik))
  # return(mean(pCE*likelihood_data_S0(data, theta))/mean(likelihood_data_S0(data, theta)))
}

## Diagnostic Probability P(C|E) given Structure S1 and the data
## (i.e., given a specific contingency table)
##
## noisy-OR parameterization for *generative* causes
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##
p_CE_S1 <- function(data, S1_lik, theta)
{
  pCE <- (theta[,2]+theta[,3]-theta[,2]*theta[,3])*theta[,1]/( (theta[,2]+theta[,3]-theta[,2]*theta[,3])*theta[,1] + theta[,3]*(1-theta[,1]))
  return(mean(pCE*S1_lik)/mean(S1_lik))
}


## Diagnostic Probability P(C|E) given Structure S1 and the data
## (i.e., given a specific contingency table)
##
## noisy-AND-NOT parameterization for *preventive* causes
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (preventive causal strength)
##           theta[,3] = w_a (strength of background) 
##
p_CE_S1_prev <- function(data, S1_lik, theta)
{
  pCE <- (theta[,3] * (1 - theta[,2])) * theta[,1] /
    ((theta[,3] * (1 - theta[,2])) * theta[,1] + theta[,3] * (1 - theta[,1]))
  
  return(mean(pCE * S1_lik) / mean(S1_lik))
}


## Predictive Probability P(E|C) given Structure S0 and the data
## (i.e., given a specific contingency table)
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##
p_EC_S0 <- function(data, S0_lik, theta)
{
  pEC <- theta[,3]
  return(mean(pEC*S0_lik)/mean(S0_lik))
}

# probability of effect given no C, P(E|no C), given S0 and data
# identical to p_EC_S0 under both preventive and generative causes but explicit name for clarity
p_EnoC_S0 <- p_EC_S0
p_EnoC_S0_prev <- p_EnoC_S0

## Predictive Probability P(E|C) given Structure S1 and the data
## (i.e., given a specific contingency table)
##
## noisy-OR parameterization for *generative* causes
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##

# probability of effect E given C under S1
p_EC_S1 <- function(data, S1_lik, theta)
{
  pEC <- theta[,2] + theta[,3] - theta[,2]*theta[,3]
  
  return(mean(pEC*S1_lik)/mean(S1_lik))
}

## Predictive Probability P(E|C) given Structure S1 and the data
## (i.e., given a specific contingency table)
##
## noisy-AND-NOT parameterization for *preventive* causes
##
## Input: "data", a contingency table with
##           data[1] = N(C=0, E=0)
##           data[2] = N(C=0, E=1)
##           data[3] = N(C=1, E=0)
##           data[4] = N(C=1, E=1)
##
##        "theta", a set of parameter vectors (samples from prior)
##           theta[,1] = b_c (base rate of cause)
##           theta[,2] = w_c (causal strength)
##           theta[,3] = w_a (strength of background) 
##

## probability of effect E given C under S1
## E occurs if wa causes the effect and C does not prevent it
p_EC_S1_prev <- function(data, S1_lik, theta)
{
  pEC <- theta[,3] * (1- theta[,2])
  
  return(mean(pEC*S1_lik)/mean(S1_lik))
}


## probability of effect E given no C, P(E|no C), given S1 and data
## under both generative and preventive S1, when C is absent, only the background cause can generate E.
p_EnoC_S1 <- function(data, S1_lik, theta)
{
  # P(E=1 | C=0, w_a, S1) = w_a  (only background case wa)
  pEnoC <- theta[,3]
  
  return(mean(pEnoC*S1_lik)/mean(S1_lik))
  # return(mean(pEC*likelihood_data_S0(data, theta))/mean(likelihood_data_S0(data, theta)))
}

## Function to compute Bayesian estimate of wc (causal strength)
## only done for structure S1; under S0 wc is fixed to 0
## generative case: posterior mean generative strength
## preventive case: posterior mean preventive strength
S1_wc_pp <- function(data, S1_lik, theta)
{
  return( sum(S1_lik*theta[,2])/sum(S1_lik))
}

## Function to compute Bayesian estimate of bc (base rate of cause event)
## identical for generative and preventive case
S0_bc_pp <- function(data, S0_lik, theta)
{
  return( sum(S0_lik*theta[,1])/sum(S0_lik))
}

S1_bc_pp <- function(data, S1_lik, theta)
{
  return( sum(S1_lik*theta[,1])/sum(S1_lik))
}

## Function to compute Bayesian estimate of wa (strength of background cause, estimated from P(e|no c))
## identical for generative and preventive case
S0_wa_pp <- function(data, S0_lik, theta)
{
  return( sum(S0_lik*theta[,3])/sum(S0_lik))
}

S1_wa_pp <- function(data, S1_lik, theta)
{
  return( sum(S1_lik*theta[,3])/sum(S1_lik))
}


## MAIN FUNCTION - generates predictions of the structure induction (SI) model
## 
## Input: "data", a set of contingency tables (rows) with
##           data[,1] = N(C=0, E=0)
##           data[,2] = N(C=0, E=1)
##           data[,3] = N(C=1, E=0)
##           data[,4] = N(C=1, E=1)
##
##        "m", number of bootstrap samples to be drawn
##
## Output: a data.frame (created from the pred matrix) with columns:
##
##  Data (observed frequencies)
##    N_00        = N(C=0, E=0)   (data[,1])
##    N_01        = N(C=0, E=1)   (data[,2])
##    N_10        = N(C=1, E=0)   (data[,3])
##    N_11        = N(C=1, E=1)   (data[,4])
##
##  MLE summaries from observed frequencies
##    deltaP_MLE  = ΔP = P(E|C) - P(E|¬C) (MLE from data)
##    wc_MLE      = causal power (MLE; generative vs. preventive formula depending on sign of ΔP)
##    bc_MLE      = base rate of cause, P(C) (MLE from data)
##    wa_MLE      = background strength, P(E|¬C) (MLE from data)
##    pEC_MLE     = P(E|C) (MLE from data)
##    pCE_MLE     = P(C|E) (MLE from data; “simple Bayes” from frequencies)
##
##  Structure-marginal likelihoods (prior-free model evidence)
##    S0_marglik     = Monte-Carlo estimate of P(D|S0)
##    S1_marglik     = Monte-Carlo estimate of P(D|S1)
##    S0_logmarglik  = log P(D|S0)  (numerically stabilized)
##    S1_logmarglik  = log P(D|S1)  (numerically stabilized)
##
##  Causal support / Bayes factor (Griffiths & Tenenbaum, 2005)
##    logBF10     = log BF_10 = log(P(D|S1) / P(D|S0))
##    BF10   = exp(logBF10)
##
##  Structure priors (user-specified; may be scalar or per-row)
##    S0_prior    = P(S0)
##    S1_prior    = P(S1)
##
##  Structure posteriors given data
##    S0_pp       = P(S0|D)
##    S1_pp       = P(S1|D)
##
##  Structure-conditional posterior predictions (Bayesian parameter averaging within structure)
##    S0_pCE      = P(C|E, D, S0)
##    S1_pCE      = P(C|E, D, S1)
##    S0_pEC      = P(E|C, D, S0)
##    S1_pEC      = P(E|C, D, S1)
##    S0_pEnoC    = P(E|¬C, D, S0)
##    S1_pEnoC    = P(E|¬C, D, S1)
##
##  Structure-conditional posterior parameter estimates
##    S0_wc_pp    = E[w_c | D, S0]  (fixed to 0 by design)
##    S1_wc_pp    = E[w_c | D, S1]
##    S0_bc_pp    = E[b_c | D, S0]
##    S1_bc_pp    = E[b_c | D, S1]
##    S0_wa_pp    = E[w_a | D, S0]
##    S1_wa_pp    = E[w_a | D, S1]
##
##  SI model predictions (Bayesian model averaging across structures)
##    SI_mean_pEC    = P(E|C, D) averaged over S0/S1 using structure posterior
##    SI_mean_pEnoC  = P(E|¬C, D) averaged over S0/S1 using structure posterior
##    SI_mean_pCE    = P(C|E, D) averaged over S0/S1 using structure posterior
##    SI_mean_wc     = E[w_c | D] averaged over S0/S1 (note: S0 contributes 0)
##    SI_mean_bc     = E[b_c | D] averaged over S0/S1
##    SI_mean_wa     = E[w_a | D] averaged over S0/S1

generate_SI_preds <- function(data, # contingency data (observed frequencies)
                  m, # number of Monte-Carlo samples
                  S1_prior = 0.5, # prior of S1, with S0 = 1-S1_prior
                  model_type = c("noisy_or_generative", "noisy_and_not_preventive"), # S1 is generative (noisy-OR) or preventive (noisy-AND-NOT)
                  wc_prior_alpha = 1, # default shape parameter 1 of beta distributions (alpha)
                  wc_prior_beta = 1 # # default shape parameter 2 of beta distributions (beta)
                  ) 
{
  
  model_type <- match.arg(model_type)
  
  if (model_type == "noisy_or_generative") {
    likelihood_S1_fun     <- likelihood_data_S1
    likelihood_S1_log_fun <- likelihood_data_S1_log
    p_CE_S1_fun           <- p_CE_S1
    p_EC_S1_fun           <- p_EC_S1
  } else if (model_type == "noisy_and_not_preventive") {
    likelihood_S1_fun     <- likelihood_data_S1_prev
    likelihood_S1_log_fun <- likelihood_data_S1_prev_log
    p_CE_S1_fun           <- p_CE_S1_prev
    p_EC_S1_fun           <- p_EC_S1_prev
  }
  
  pred <- matrix(data=0, nrow = dim(data)[1], ncol=40,
                 dimnames = list(1:dim(data)[1],
                                 c("N", "N_00", "N_01", "N_10", "N_11", # data (observed frequencies)
                                   "deltaP_MLE", "wc_MLE", "bc_MLE", "wa_MLE", # maximum likelihood estimates of parameters
                                   "pEC_MLE", "pCE_MLE", # # maximum likelihood estimates of the conditional probabilities P(e|c) and P(c|e)
                                   "model_type", # 
                                   "S0_marglik", "S1_marglik", "S0_logmarglik", "S1_logmarglik", # (log) marginal likelihood of data under S1 and S0
                                   "logBF10", "BF10", # Bayes Factor, aka "causal support" (Griffiths & Tenenbaum, 2005)
                                   "S0_prior", "S1_prior", # prior probabilities of causal structures S1 and S0
                                   "S0_pp", "S1_pp", # posterior probabilities of causal structures S1 and S0
                                   "S0_pCE", "S1_pCE",   # posterior conditional probability of P(c|e) under causal structures S1 and S0 and mean weighted by posterior of structures
                                   "S0_pEC", "S1_pEC",  # posterior conditional probability of P(e|c) under causal structures S1 and S0 and mean weighted by posterior of structures
                                   "S0_pEnoC", "S1_pEnoC",  # posterior conditional probability of P(e|c) under causal structures S1 and S0 and mean weighted by posterior of structures
                                   "S0_wc_pp", "S1_wc_pp","S0_bc_pp", "S1_bc_pp", "S0_wa_pp", "S1_wa_pp", # posterior estimates of the graphs' parameters, for each structure
                                   "SI_mean_pEC", "SI_mean_pEnoC",  "SI_mean_pCE", # mean posterior conditional probabilities P(c|e) and P(e|c) weighted by posterior probability of causal structures S1 and S0 (i.e. Bayesian model averaging)
                                   "SI_mean_wc", "SI_mean_bc", "SI_mean_wa") # posterior mean parameter estimates
                 )
  ) 
  
  
  # make df
  pred <- as.data.frame(pred)
  
  # precision
  eps <- .Machine$double.xmin
  
  # generative vs preventive model
  pred$model_type <- model_type
  
  # observed frequencies
  pred$N <- data[,1] + data[,2] + data[,3] + data[,4]
  pred$N_00 <- data[,1]
  pred$N_01 <- data[,2]
  pred$N_10 <- data[,3]
  pred$N_11 <- data[,4]
  
  # maximum likelihood estimates of conditional and unconditional probabilities
  pred$bc_MLE   <- (data[,3] + data[,4]) /(data[,1] + data[,2] + data[,3] + data[,4]) # base rate of cause C, P(c)
  pred$wa_MLE   <- data[,2]/(data[,1]+data[,2]) # strength of background cause, equals the likelihood of the effect if the cause is absent, P(e|no c)
  pred$pEC_MLE  <- data[,4]/(data[,3]+data[,4]) # likelihood of the effect if the cause is present
  pred$pCE_MLE  <- data[,4]/(data[,2]+data[,4]) # probability of the cause if the effect is present
  
  # maximum likelihood estimate of contingency delta P
  pred$deltaP_MLE <- (data[,4]/(data[,3]+data[,4]) - data[,2]/(data[,1]+data[,2])) 
  
  # # maximum likelihood estimate of causal power (with different calculation depending on whether the contingency is positive (=generative cause) or negative (inhibitory cause))
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
  
  #draw flat priors over parameters of structure S1
  # theta_S1 <- matrix(data=0, nrow=m, ncol=3)
  # theta_S1[,1] <- rbeta(m, 1, 1)
  # theta_S1[,2] <- rbeta(m, 1, 1)
  # theta_S1[,3] <- rbeta(m, 1, 1)
  
  # Draw parameter priors under S1
  theta_S1 <- matrix(
    data = 0,
    nrow = m,
    ncol = 3
  )
  
  # Base rate of the candidate cause:
  # bc ~ Beta(1,1)
  theta_S1[,1] <- rbeta(
    n = m,
    shape1 = 1,
    shape2 = 1
  )
  
  # Causal strength:
  # wc ~ Beta(wc_prior_alpha, wc_prior_beta)
  theta_S1[,2] <- rbeta(
    n = m,
    shape1 = wc_prior_alpha,
    shape2 = wc_prior_beta
  )
  
  # Strength of the background cause:
  # wa ~ Beta(1,1)
  theta_S1[,3] <- rbeta(
    n = m,
    shape1 = 1,
    shape2 = 1
  )
  
  #draw flat priors for S0; ignore wc by setting it to zero
  theta_S0 <- matrix(data=0, nrow=m, ncol=3)
  theta_S0[,1] <- rbeta(n=m, shape1 = 1, shape2 = 1)
  theta_S0[,2] <- rep(0, m)
  theta_S0[,3] <- rbeta(n=m, shape1 = 1, shape2 = 1)
  
  # causal structure priors
  # default: P(S1)=P(S0)=0.5 i.e. uniform prior
  # S0_prior <- 1 - S1_prior
  
  # allow scalar or per-row priors
  if (length(S1_prior) == 1L) {
    S1_prior <- rep(S1_prior, nrow(data))
  } else if (length(S1_prior) != nrow(data)) {
    stop("S1_prior must be a single value or a vector of length nrow(data).")
  }
  
  #for each contingency table:
  for (i in 1:dim(data)[1])
  {
    S1p <- S1_prior[i]
    S0p <- 1 - S1p
    
    pred$S1_prior[i] <- S1p
    pred$S0_prior[i] <- S0p
    
    # likelihood of data under S0 and S1
    S0_lik <- likelihood_data_S0(data[i, ], theta_S0)
    S1_lik <- likelihood_S1_fun(data[i, ], theta_S1)
    # S1_lik <- likelihood_data_S1(data[i, ], theta_S1)

    # marginal likelihood of data under S0 and S1
    pred$S0_marglik[i] <- mean(S0_lik)
    pred$S1_marglik[i] <- mean(S1_lik)

    # log marginal likelihood of data under S0 and S1
    pred$S0_logmarglik[i] <- log(pmax(pred$S0_marglik[i], eps))
    pred$S1_logmarglik[i] <- log(pmax(pred$S1_marglik[i], eps))

    # causal support (log Bayes factor)
    pred$logBF10[i]   <- pred$S1_logmarglik[i] - pred$S0_logmarglik[i]
    pred$BF10[i] <- exp(pred$logBF10[i])

    # unnormalized posterior of S0
    # pred$S0_pp[i] <- mean(likelihood_data_S0(data[i,], theta_S0))* S0p
    pred$S0_pp[i] <- pred$S0_marglik[i] * S0p

    # unnormalized posterior of S1
    # pred$S1_pp[i] <- mean(likelihood_data_S1(data[i,], theta_S1))* S1p
    pred$S1_pp[i] <-  pred$S1_marglik[i] * S1p
    
   
    
    # diagnostic probs given each structure
    pred$S0_pCE[i] <- p_CE_S0(data[i,], S0_lik, theta_S0)
    pred$S1_pCE[i] <- p_CE_S1_fun(data[i,], S1_lik, theta_S1)
    # pred$S1_pCE[i] <- p_CE_S1(data[i,], S1_lik, theta_S1)
    
    # predictive probs given each structure
    pred$S0_pEC[i] <- p_EC_S0(data[i,], S0_lik, theta_S0)
    pred$S1_pEC[i] <- p_EC_S1_fun(data[i,], S1_lik, theta_S1)
    # pred$S1_pEC[i] <- p_EC_S1(data[i,], S1_lik, theta_S1)
    
    pred$S0_pEnoC[i] <- p_EnoC_S0(data[i,], S0_lik, theta_S0)
    pred$S1_pEnoC[i] <- p_EnoC_S1(data[i,], S1_lik, theta_S1)
    
    # Compute Bayesian posterior estimates of graph's parameters for each data set under S1 and S0
    pred$S0_wc_pp[i] <- 0
    pred$S1_wc_pp[i] <- S1_wc_pp(data[i,], S1_lik, theta_S1)     
    
    pred$S0_bc_pp[i] <- S0_bc_pp(data[i,], S0_lik, theta_S0)
    pred$S1_bc_pp[i] <- S1_bc_pp(data[i,], S1_lik, theta_S1)  
    
    pred$S0_wa_pp[i] <- S0_wa_pp(data[i,], S0_lik, theta_S0)
    pred$S1_wa_pp[i] <- S1_wa_pp(data[i,], S1_lik, theta_S1)  
    
  }
  
  # normalize posterior over structures 
  pred[,c("S0_pp", "S1_pp")] <- pred[,c("S0_pp", "S1_pp")]/rowSums(pred[,c("S0_pp", "S1_pp")])
  
  # calculate Bayesian model average over diagnostic probs (i.e., prediction of SI model)
  pred$SI_mean_pCE <- pred$S0_pp*pred$S0_pCE + pred$S1_pp*pred$S1_pCE
  
  # calculate Bayesian model average over predictive probs (i.e., prediction of SI model)
  pred$SI_mean_pEC <- pred$S0_pp*pred$S0_pEC + pred$S1_pp*pred$S1_pEC
  pred$SI_mean_pEnoC <- pred$S0_pp*pred$S0_pEnoC + pred$S1_pp*pred$S1_pEnoC
  
  
  # calculate Bayesian model average over graph parameters (i.e., prediction of SI model)
  pred$SI_mean_wc <- pred$S0_pp * pred$S0_wc_pp + pred$S1_pp*pred$S1_wc_pp
  pred$SI_mean_bc <- pred$S0_pp * pred$S0_bc_pp + pred$S1_pp*pred$S1_bc_pp
  pred$SI_mean_wa <- pred$S0_pp * pred$S0_wa_pp + pred$S1_pp*pred$S1_wa_pp
  
  return(pred)
}




