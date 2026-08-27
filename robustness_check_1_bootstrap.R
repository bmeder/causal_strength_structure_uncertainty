# ==============================================================================
# Robustness check 1: Bootstrap model comparison
# ==============================================================================
# Assess whether the relative performance of the models is stable across studies.
# In each of 10,000 bootstrap iterations, sample 12 studies with replacement and
# compute the correlation and RMSE for each model. Track how often each model
# yields the highest correlation and lowest RMSE.


source("structure_induction_func.R") # functions for generating predictions of the structure induction (SI) model (Meder et al., 2014, 2026)
source("SS_strength_model_func.R") # functions for generating causal strength estimates under S1 given sparse and strong (SS) priors (Lu et al, 2008, Psych Review)
source("meta_analysis_preprocessing.R") # reads data for meta analysis from various data files stored in directory "data"



set.seed(0815)

B <- 10000 # number of siomulations


# --- get data for meta analysis --------------
df_data <- meta_analysis_preprocessing()

# Model comparison: get model predictions

# --- make df with all experimental conditions (cell frequencies) --------------
# cell input order for models: [d,c,b,a] = [N_00,N_01,N_10,N_11]
df_freq <- df_data %>% 
  select(N_00,  N_01, N_10, N_11) %>% 
  as.matrix()

# --- Generate SI predictions in the same row order as df_data --------------
# generate SI predictions under uniform priors for structure, P(S1)=P(S0)=0.5

m <- 10^6

# generative (noisy-OR)
set.seed(0815)
SI_gen <- generate_SI_preds(
  df_freq,
  m = m,
  S1_prior = 0.5,
  model_type = "noisy_or_generative"
) %>%
  mutate(id = factor(row_number())) %>%
  select(
    id,
    N, N_00, N_01, N_10, N_11,
    BF10, logBF10,
    S0_marglik, S1_marglik,
    S0_logmarglik, S1_logmarglik,
    S0_pp, S1_pp,
    S1_wc_pp,
    SI_mean_wc
  ) %>%
  rename(
    BF10_gen = BF10,
    logBF10_gen = logBF10,
    S0_marglik_gen = S0_marglik,
    S1_marglik_gen = S1_marglik,
    S0_logmarglik_gen = S0_logmarglik,
    S1_logmarglik_gen = S1_logmarglik,
    S0_pp_gen = S0_pp,
    S1_pp_gen = S1_pp,
    S1_wc_pp_gen = S1_wc_pp,
    SI_mean_wc_gen = SI_mean_wc
  )

# preventive (noisy-AND-NOT)
set.seed(0815)
SI_prev <- generate_SI_preds(
  df_freq,
  m = m,
  S1_prior = 0.5,
  model_type = "noisy_and_not_preventive"
) %>%
  mutate(id = factor(row_number())) %>%
  select(
    id,
    N, N_00, N_01, N_10, N_11,
    BF10, logBF10,
    S0_marglik, S1_marglik,
    S0_logmarglik, S1_logmarglik,
    S0_pp, S1_pp,
    S1_wc_pp,
    SI_mean_wc
  ) %>%
  rename(
    BF10_prev = BF10,
    logBF10_prev = logBF10,
    S0_marglik_prev = S0_marglik,
    S1_marglik_prev = S1_marglik,
    S0_logmarglik_prev = S0_logmarglik,
    S1_logmarglik_prev = S1_logmarglik,
    S0_pp_prev = S0_pp,
    S1_pp_prev = S1_pp,
    S1_wc_pp_prev = S1_wc_pp,
    SI_mean_wc_prev = SI_mean_wc
  )

# --- Generate SS priors prediction: posterior causal strength under S1 given sparse and strong priors (Lu et al., 2008) --------------
# only structure S1 is considered
# generative or preventive SS prior over wc and wa:
# base rate bc sampled from uniform Beta(1,1) prior as in SI model

# generative (noisy-OR)
set.seed(0815)
SS_gen <- generate_SS_strength_preds(
  df_freq,
  m = m,
  alpha = 5,
  model_type = "noisy_or_generative"
) %>%
  mutate(id = factor(row_number())) %>%
  select(id, 
         N, N_00, N_01, N_10, N_11,
         SS_wc_pp, SS_marglik, SS_logmarglik) %>%
  rename(SS_wc_pp_gen = SS_wc_pp,
         SS_marglik_gen = SS_marglik,
         SS_logmarglik_gen = SS_logmarglik)

# preventive (noisy-AND-NOT)
set.seed(0815)
SS_prev <- generate_SS_strength_preds(
  df_freq,
  m = m,
  alpha = 5,
  model_type = "noisy_and_not_preventive"
) %>%
  mutate(id = factor(row_number())) %>%
  select(id, 
         N, N_00, N_01, N_10, N_11,
         SS_wc_pp, SS_marglik, SS_logmarglik) %>%
  rename(SS_wc_pp_prev = SS_wc_pp,
         SS_marglik_prev = SS_marglik,
         SS_logmarglik_prev = SS_logmarglik)

# --- reference table with all relevant info --------------------------
model_predictions <- df_data %>%
  #select(study, experiment, id, a, b, c, d, deltaP, type_deltaP, rating) %>%
  select(-reference) %>%
  left_join(SI_gen, by = c("id", "N_00",  "N_01", "N_10", "N_11")) %>%
  left_join(SI_prev, by = c("id", "N_00",  "N_01", "N_10", "N_11", "N")) %>%
  left_join(SS_gen, by = c("id", "N_00",  "N_01", "N_10", "N_11", "N")) %>%
  left_join(SS_prev, by = c("id", "N_00",  "N_01", "N_10", "N_11", "N")) %>%
  rename(human = rating)


# --- generate model predictions -----------------------------------------
# specify which functional form is used (e.g. generative via noisy-OR or preventive via noisy-AND-NOT)
# using the framing of the causal inference context
# generative: participants' task was to infer a generative relation (e.g. on a scale from 0-100)
# preventive: participants' task was to infer a preventive relation (e.g. on a scale from -100 to 0)
# bidirectional: participants' were informed that the relation could be positive or negative (e.g. on a scale from -100  to +100)
# In bidirectional conditions, generative and preventive predictions are combined using their posterior or marginal-likelihood weights

model_predictions <- model_predictions %>%
  mutate(
    
    # ------------------------------------------------------------------
    # Posterior link probability and signed response prediction
    # ------------------------------------------------------------------
    # For bidirectional framing, following Griffiths and Tenenbaum
    # (2009, pp. 681–682), jointly compare three structures using a
    # uniform prior:
    #   P(S0)            = 1/3
    #   P(S1_generative) = 1/3
    #   P(S1_preventive) = 1/3
    
    # Unnormalized structure posteriors:
    # marginal likelihood multiplied by structure prior.
    SI_bidir_S0_unnorm = (1/3) * S0_marglik_gen,
    SI_bidir_S1_gen_unnorm = (1/3) * S1_marglik_gen,
    SI_bidir_S1_prev_unnorm = (1/3) * S1_marglik_prev,
    
    # Common normalizing constant for the three structures.
    SI_bidir_normalizer = SI_bidir_S0_unnorm + SI_bidir_S1_gen_unnorm +  SI_bidir_S1_prev_unnorm,
    
    # Posterior probabilities of the three structures.
    SI_bidir_S0_pp = SI_bidir_S0_unnorm / SI_bidir_normalizer, # structure S0
    SI_bidir_S1_gen_pp = SI_bidir_S1_gen_unnorm / SI_bidir_normalizer, # S1 generative
    SI_bidir_S1_prev_pp = SI_bidir_S1_prev_unnorm / SI_bidir_normalizer, # S1 preventive
    
    # check
    # model_predictions$SI_bidir_S1_gen_pp + model_predictions$SI_bidir_S1_prev_pp + model_predictions$SI_bidir_S0_pp
    
    # Posterior probability that some causal link exists,
    # irrespective of whether it is generative or preventive.
    SI_bidir_link_pp = SI_bidir_S1_gen_pp + SI_bidir_S1_prev_pp,
    
    # Equivalent:SI_bidir_link_pp = 1 - SI_bidir_S0_pp
    
    # Prediction on the signed -100 to +100 response scale.
    # For bidirectional framing, Griffiths and Tenenbaum map the
    # three-structure posterior to:
    # -100 * P(S1_preventive | D) + 100 * P(S1_generative | D).
    S1_pp_prediction = case_when(
      framing == "generative" ~ S1_pp_gen * 100,
      framing == "preventive" ~ -S1_pp_prev * 100,
      framing == "bidirectional" ~ -100 * SI_bidir_S1_prev_pp + 100 * SI_bidir_S1_gen_pp
    ),
    
    # ------------------------------------------------------------------
    # Mean posterior strength under S1 with sparse-and-strong (SS) priors over w_a and w_c (Lu et al, 2008)
    # ------------------------------------------------------------------
    
    # Version 1: sparse-and-strong (SS) priors over w_a and w_c (Lu et al, 2008)
    # for generative and preventive framing the corresponding S1 model (noisy-or vs. noisy-and-not) are used
    # For bidirectional judgments, following Lu et al. (2008, Eq. 15, p.968),
    # the positive generative and negative preventive posterior mean
    # strengths are averaged using P(g | D), with g={S1_generative, S1_preventive}
    # obtained by normalizing the generative and preventive S1 marginal likelihoods under SS priors.
    SS_gen_direction_pp = SS_marglik_gen / (SS_marglik_gen + SS_marglik_prev),
    SS_prev_direction_pp = SS_marglik_prev / (SS_marglik_gen + SS_marglik_prev),
    
    SS_wc_prediction = case_when(
      framing == "generative" ~ SS_wc_pp_gen * 100,
      framing == "preventive" ~ SS_wc_pp_prev * -100,
      framing == "bidirectional" ~
        SS_gen_direction_pp  * (SS_wc_pp_gen  * 100) +
        SS_prev_direction_pp * (SS_wc_pp_prev *-100)
    ),
    
    # ------------------------------------------------------------------
    # SI mode: Structure-averaged strength
    # ------------------------------------------------------------------
    # this approach considers S0 and S1 (generative and preventive parameterization, respectively)
    # for generative and preventive framing, with uniform prior over structures, P(S0)=P(S1)=0.5
    # For bidirectional framing, following Griffiths & Tenenbaum (2009) we use a uniform prior over all three structures:
    #   P(S0)            = 1/3
    #   P(S1_generative) = 1/3
    #   P(S1_preventive) = 1/3
    
    # S0 is the same model in the generative and preventive SI runs,
    # so only one S0 marginal likelihood is used. sanity check:
    # all.equal(model_predictions$S0_marglik_gen, model_predictions$S0_marglik_prev)
    
    # these quantitities are computed above 
    # Compute the unnormalized posterior probability of each structure
    # as marginal likelihood × structure prior.
    # SI_bidir_S0_unnorm =  1/3 * S0_marglik_gen,
    # SI_bidir_S1_gen_unnorm =  1/3 * S1_marglik_gen,
    # SI_bidir_S1_prev_unnorm =  1/3 * S1_marglik_prev,
    # 
    # Normalize across all three structures so that their posterior
    # probabilities sum to 1.
    # SI_bidir_normalizer = SI_bidir_S0_unnorm + SI_bidir_S1_gen_unnorm + SI_bidir_S1_prev_unnorm,
    
    # SI_bidir_S0_pp = SI_bidir_S0_unnorm / SI_bidir_normalizer, 
    # SI_bidir_S1_gen_pp = SI_bidir_S1_gen_unnorm / SI_bidir_normalizer,
    # SI_bidir_S1_prev_pp =  SI_bidir_S1_prev_unnorm / SI_bidir_normalizer,
    # 
    # Bayesian model average of causal strength across the three structures:
    # generative strength enters positively, preventive strength negatively,
    # and S0 contributes zero.
    
    SI_mean_wc_prediction = case_when(
      framing == "generative" ~ SI_mean_wc_gen * 100,
      framing == "preventive" ~ SI_mean_wc_prev * -100,
      framing == "bidirectional" ~
        SI_bidir_S1_gen_pp  * ( S1_wc_pp_gen  * 100) +
        SI_bidir_S1_prev_pp * (S1_wc_pp_prev * -100) +
        SI_bidir_S0_pp      * 0
    )
  )

models <- c(
  "SI_mean_wc_prediction",
  "S1_pp_prediction",
  "SS_wc_prediction"
)

model_labels <- c(
  SI_mean_wc_prediction = "Structure-averaged\nstrength",
  S1_pp_prediction      = "Posterior link\nprobability",
  SS_wc_prediction      = "Bayesian strength\n(SS priors)"
)

# -------------------------------------------------------------------------
# Bootstrap dataset
# -------------------------------------------------------------------------

df_boot <- model_predictions %>%
  select(
    study,
    human,
    all_of(models)
  ) %>%
  filter(
    !is.na(human),
    if_all(all_of(models), ~ !is.na(.x))
  )

studies <- unique(df_boot$study)
n_studies <- length(studies)


# -------------------------------------------------------------------------
# Compute study-level fit values
# -------------------------------------------------------------------------
# Each study contributes one correlation and one RMSE per model.
# Pearson correlations are computed directly with cor(), including studies
# with only three conditions. RMSE is computed directly from squared errors.

study_fit <- df_boot %>%
  group_by(study) %>%
  summarise(
    across(
      all_of(models),
      list(
        r = ~ cor(
          human,
          .x,
          use = "pairwise.complete.obs"
        ),
        rmse = ~ sqrt(
          mean(
            (human - .x)^2,
            na.rm = TRUE
          )
        )
      ),
      .names = "{.col}_{.fn}"
    ),
    .groups = "drop"
  )

r_cols <- paste0(models, "_r")
rmse_cols <- paste0(models, "_rmse")

study_r_matrix <- as.matrix(study_fit[, r_cols])
study_rmse_matrix <- as.matrix(study_fit[, rmse_cols])

colnames(study_r_matrix) <- models
colnames(study_rmse_matrix) <- models

# -------------------------------------------------------------------------
# Bootstrap
# -------------------------------------------------------------------------

boot_r_values <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(models),
  dimnames = list(NULL, models)
)

boot_rmse_values <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(models),
  dimnames = list(NULL, models)
)

boot_r_winner <- character(B)
boot_rmse_winner <- character(B)

for (b in seq_len(B)) {
  
  sampled_studies <- sample(
    seq_len(n_studies),
    size = n_studies,
    replace = TRUE
  )
  
  r_b <- colMeans(
    study_r_matrix[sampled_studies, , drop = FALSE],
    na.rm = TRUE
  )
  
  rmse_b <- colMeans(
    study_rmse_matrix[sampled_studies, , drop = FALSE],
    na.rm = TRUE
  )
  
  boot_r_values[b, ] <- r_b
  boot_rmse_values[b, ] <- rmse_b
  
  boot_r_winner[b] <- models[which.max(r_b)]
  boot_rmse_winner[b] <- models[which.min(rmse_b)]
}

# -------------------------------------------------------------------------
# Winner proportions
# -------------------------------------------------------------------------

bootstrap_winner_proportions <- bind_rows(
  tibble(
    metric = "Correlation",
    winner = boot_r_winner
  ),
  tibble(
    metric = "RMSE",
    winner = boot_rmse_winner
  )
) %>%
  count(metric, winner, name = "n_best") %>%
  group_by(metric) %>%
  complete(
    winner = models,
    fill = list(n_best = 0)
  ) %>%
  mutate(
    B = B,
    prop_best = n_best / B,
    model = recode(winner, !!!model_labels),
    model = factor(model, levels = unname(model_labels)),
    metric = factor(metric, levels = c("Correlation", "RMSE"))
  ) %>%
  ungroup()


# results of bootstrap analysis -------------------------------------------

bootstrap_winner_proportions

















