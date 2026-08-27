# =============================================================================
# Robustness check 2: Original Perales & Shanks (2007) data only
# =============================================================================
# Repeat the main model comparison using only the original contingency-table
# compilation reported by Perales and Shanks (2007).
#
# Models:
#   1. Posterior link probability
#   2. Bayesian strength under S1 with sparse-and-strong (SS) priors
#   3. Structure-averaged strength
#
# Fit is evaluated as in the main analysis:
#   - Pearson correlation calculated separately within each study
#   - RMSE calculated separately within each study
#   - Plain average across studies, giving each study equal weight
# =============================================================================


# -----------------------------------------------------------------------------
# Setup
# -----------------------------------------------------------------------------

library(tidyverse)
library(readr)

source("structure_induction_func.R")
source("SS_strength_model_func.R")


# -----------------------------------------------------------------------------
# Load and preprocess original Perales & Shanks (2007) data
# -----------------------------------------------------------------------------

df_data <- read_csv(
  "data/dat_Perales_Shanks_2007.csv",
  show_col_types = FALSE
) %>%
  filter(
    !if_all(everything(), is.na)
  ) %>%
  select(
    where(~ !all(is.na(.)))
  ) %>%
  rename(
    N_11 = a,
    N_10 = b,
    N_01 = c,
    N_00 = d
  ) %>%
  mutate(
    experiment = as.character(experiment),
    
    id = factor(row_number()),
    
    N =
      N_11 +
      N_10 +
      N_01 +
      N_00,
    
    pC =
      (N_11 + N_10) / N,
    
    pEC =
      N_11 /
      (N_11 + N_10),
    
    pEnoC =
      N_01 /
      (N_01 + N_00),
    
    deltaP =
      pEC - pEnoC,
    
    power = case_when(
      is.nan(deltaP) ~ NaN,
      deltaP == 0 ~ 0,
      deltaP > 0 ~
        deltaP / (1 - pEnoC),
      deltaP < 0 ~
        -deltaP / pEnoC
    ),
    
    sample_direction = case_when(
      deltaP > 0 ~ "generative",
      deltaP == 0 ~ "zero",
      deltaP < 0 ~ "preventive"
    )
  ) %>%
  select(
    id,
    study,
    experiment,
    n,
    framing,
    N,
    N_11,
    N_10,
    N_01,
    N_00,
    pEC,
    pEnoC,
    deltaP,
    power,
    sample_direction,
    rating
  )


# -----------------------------------------------------------------------------
# Contingency-table input for model predictions
# -----------------------------------------------------------------------------
# Cell order:
# [d, c, b, a] = [N_00, N_01, N_10, N_11]

df_freq <- df_data %>%
  select(
    N_00,
    N_01,
    N_10,
    N_11
  ) %>%
  as.matrix()


# -----------------------------------------------------------------------------
# Structure induction: noisy-OR (generative)
# -----------------------------------------------------------------------------

set.seed(0815)

SI_gen <- generate_SI_preds(
  df_freq,
  m = 10^6,
  S1_prior = 0.5,
  model_type = "noisy_or_generative"
) %>%
  mutate(
    id = factor(row_number())
  ) %>%
  select(
    id,
    N,
    N_00,
    N_01,
    N_10,
    N_11,
    BF10,
    logBF10,
    S0_marglik,
    S1_marglik,
    S0_logmarglik,
    S1_logmarglik,
    S0_pp,
    S1_pp,
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


# -----------------------------------------------------------------------------
# Structure induction: noisy-AND-NOT (preventive)
# -----------------------------------------------------------------------------

set.seed(0815)

SI_prev <- generate_SI_preds(
  df_freq,
  m = 10^6,
  S1_prior = 0.5,
  model_type = "noisy_and_not_preventive"
) %>%
  mutate(
    id = factor(row_number())
  ) %>%
  select(
    id,
    N,
    N_00,
    N_01,
    N_10,
    N_11,
    BF10,
    logBF10,
    S0_marglik,
    S1_marglik,
    S0_logmarglik,
    S1_logmarglik,
    S0_pp,
    S1_pp,
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


# -----------------------------------------------------------------------------
# Bayesian strength under S1 with SS priors: noisy-OR
# -----------------------------------------------------------------------------

set.seed(0815)

SS_gen <- generate_SS_strength_preds(
  df_freq,
  m = 10^6,
  alpha = 5,
  model_type = "noisy_or_generative"
) %>%
  mutate(
    id = factor(row_number())
  ) %>%
  select(
    id,
    N,
    N_00,
    N_01,
    N_10,
    N_11,
    SS_wc_pp,
    SS_marglik,
    SS_logmarglik
  ) %>%
  rename(
    SS_wc_pp_gen = SS_wc_pp,
    SS_marglik_gen = SS_marglik,
    SS_logmarglik_gen = SS_logmarglik
  )


# -----------------------------------------------------------------------------
# Bayesian strength under S1 with SS priors: noisy-AND-NOT
# -----------------------------------------------------------------------------

set.seed(0815)

SS_prev <- generate_SS_strength_preds(
  df_freq,
  m = 10^6,
  alpha = 5,
  model_type = "noisy_and_not_preventive"
) %>%
  mutate(
    id = factor(row_number())
  ) %>%
  select(
    id,
    N,
    N_00,
    N_01,
    N_10,
    N_11,
    SS_wc_pp,
    SS_marglik,
    SS_logmarglik
  ) %>%
  rename(
    SS_wc_pp_prev = SS_wc_pp,
    SS_marglik_prev = SS_marglik,
    SS_logmarglik_prev = SS_logmarglik
  )


# -----------------------------------------------------------------------------
# Merge model outputs with observed judgments
# -----------------------------------------------------------------------------

model_predictions <- df_data %>%
  left_join(
    SI_gen,
    by = c(
      "id",
      "N_00",
      "N_01",
      "N_10",
      "N_11",
      "N"
    )
  ) %>%
  left_join(
    SI_prev,
    by = c(
      "id",
      "N_00",
      "N_01",
      "N_10",
      "N_11",
      "N"
    )
  ) %>%
  left_join(
    SS_gen,
    by = c(
      "id",
      "N_00",
      "N_01",
      "N_10",
      "N_11",
      "N"
    )
  ) %>%
  left_join(
    SS_prev,
    by = c(
      "id",
      "N_00",
      "N_01",
      "N_10",
      "N_11",
      "N"
    )
  ) %>%
  rename(
    human = rating
  )


# -----------------------------------------------------------------------------
# Map model quantities to the judgment scale
# -----------------------------------------------------------------------------

model_predictions <- model_predictions %>%
  mutate(
    
    # Three-structure posterior for bidirectional framing:
    # S0, S1-generative, S1-preventive
    
    SI_bidir_S0_unnorm =
      (1 / 3) * S0_marglik_gen,
    
    SI_bidir_S1_gen_unnorm =
      (1 / 3) * S1_marglik_gen,
    
    SI_bidir_S1_prev_unnorm =
      (1 / 3) * S1_marglik_prev,
    
    SI_bidir_normalizer =
      SI_bidir_S0_unnorm +
      SI_bidir_S1_gen_unnorm +
      SI_bidir_S1_prev_unnorm,
    
    SI_bidir_S0_pp =
      SI_bidir_S0_unnorm /
      SI_bidir_normalizer,
    
    SI_bidir_S1_gen_pp =
      SI_bidir_S1_gen_unnorm /
      SI_bidir_normalizer,
    
    SI_bidir_S1_prev_pp =
      SI_bidir_S1_prev_unnorm /
      SI_bidir_normalizer,
    
    
    # Posterior link probability
    
    S1_pp_prediction = case_when(
      
      framing == "generative" ~
        S1_pp_gen * 100,
      
      framing == "preventive" ~
        -S1_pp_prev * 100,
      
      framing == "bidirectional" ~
        100 * SI_bidir_S1_gen_pp -
        100 * SI_bidir_S1_prev_pp
    ),
    
    
    # Bayesian strength under S1 with SS priors
    
    SS_gen_direction_pp =
      SS_marglik_gen /
      (SS_marglik_gen + SS_marglik_prev),
    
    SS_prev_direction_pp =
      SS_marglik_prev /
      (SS_marglik_gen + SS_marglik_prev),
    
    SS_wc_prediction = case_when(
      
      framing == "generative" ~
        SS_wc_pp_gen * 100,
      
      framing == "preventive" ~
        SS_wc_pp_prev * -100,
      
      framing == "bidirectional" ~
        SS_gen_direction_pp *
        (SS_wc_pp_gen * 100) +
        SS_prev_direction_pp *
        (SS_wc_pp_prev * -100)
    ),
    
    
    # Structure-averaged strength
    
    SI_mean_wc_prediction = case_when(
      
      framing == "generative" ~
        SI_mean_wc_gen * 100,
      
      framing == "preventive" ~
        SI_mean_wc_prev * -100,
      
      framing == "bidirectional" ~
        SI_bidir_S1_gen_pp *
        (S1_wc_pp_gen * 100) +
        SI_bidir_S1_prev_pp *
        (S1_wc_pp_prev * -100)
    )
  )


# =============================================================================
# Correlations
# =============================================================================


# -----------------------------------------------------------------------------
# Correlations by study
# -----------------------------------------------------------------------------

model_correlations_by_study <- model_predictions %>%
  select(
    study,
    human,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(
      study,
      human
    ),
    names_to = "model",
    values_to = "prediction"
  ) %>%
  filter(
    !is.na(human),
    !is.na(prediction)
  ) %>%
  mutate(
    model = recode(
      model,
      S1_pp_prediction =
        "Posterior link probability",
      SS_wc_prediction =
        "Bayesian strength (SS priors)",
      SI_mean_wc_prediction =
        "Structure-averaged strength"
    )
  ) %>%
  group_by(
    study,
    model
  ) %>%
  summarise(
    n = n(),
    
    r = cor(
      human,
      prediction,
      use = "pairwise.complete.obs"
    ),
    
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# Mean correlation across studies
# -----------------------------------------------------------------------------

model_correlations_mean_study <- model_correlations_by_study %>%
  filter(
    !is.na(r),
    is.finite(r)
  ) %>%
  group_by(
    model
  ) %>%
  summarise(
    mean_r = mean(r),
    n_studies = n(),
    .groups = "drop"
  )


# =============================================================================
# RMSE
# =============================================================================


# -----------------------------------------------------------------------------
# RMSE by study
# -----------------------------------------------------------------------------

model_rmse_by_study <- model_predictions %>%
  select(
    study,
    human,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(
      study,
      human
    ),
    names_to = "model",
    values_to = "prediction"
  ) %>%
  filter(
    !is.na(human),
    !is.na(prediction)
  ) %>%
  mutate(
    sq_error =
      (human - prediction)^2,
    
    model = recode(
      model,
      S1_pp_prediction =
        "Posterior link probability",
      SS_wc_prediction =
        "Bayesian strength (SS priors)",
      SI_mean_wc_prediction =
        "Structure-averaged strength"
    )
  ) %>%
  group_by(
    study,
    model
  ) %>%
  summarise(
    n = n(),
    
    mse = mean(
      sq_error,
      na.rm = TRUE
    ),
    
    rmse = sqrt(mse),
    
    .groups = "drop"
  )


# -----------------------------------------------------------------------------
# Mean RMSE across studies
# -----------------------------------------------------------------------------

model_rmse_mean_study <- model_rmse_by_study %>%
  filter(
    !is.na(rmse),
    is.finite(rmse)
  ) %>%
  group_by(
    model
  ) %>%
  summarise(
    mean_rmse = mean(rmse),
    n_studies = n(),
    .groups = "drop"
  )


# =============================================================================
# Final comparison table
# =============================================================================

model_comparison_perales <- model_correlations_mean_study %>%
  select(
    model,
    mean_r
  ) %>%
  left_join(
    model_rmse_mean_study %>%
      select(
        model,
        mean_rmse
      ),
    by = "model"
  ) %>%
  mutate(
    model = factor(
      model,
      levels = c(
        "Structure-averaged strength",
        "Posterior link probability",
        "Bayesian strength (SS priors)"
      )
    )
  ) %>%
  arrange(
    model
  )


model_comparison_perales