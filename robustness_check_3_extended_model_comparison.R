# =============================================================================
# Robustness check 3: Extended model comparison
# =============================================================================
# Compare the three main models with two additional variants:
#
#   1. Posterior link probability
#   2. Bayesian strength under S1 with uniform priors
#   3. Bayesian strength under S1 with sparse-and-strong (SS) priors
#   4. Structure-averaged strength with uniform priors
#   5. Structure-averaged strength with sparse-and-strong (SS) priors
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
library(stringr)

source("structure_induction_func.R")
source("SS_strength_model_func.R")


# =============================================================================
# Load and preprocess full data set
# =============================================================================


# -----------------------------------------------------------------------------
# Perales & Shanks (2007) compilation
# -----------------------------------------------------------------------------
# Buehner et al. (2003) and White (2003) are removed here because their
# original data are loaded separately below.

dat_Perales_Shanks_2007 <- read_csv(
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
    experiment = as.character(experiment)
  ) %>%
  filter(
    !str_detect(study, "Buehner"),
    !str_detect(study, "White")
  ) %>%
  select(
    study,
    experiment,
    n,
    framing,
    N_11,
    N_10,
    N_01,
    N_00,
    rating
  )


# -----------------------------------------------------------------------------
# Buehner, Cheng, & Clifford (2003)
# -----------------------------------------------------------------------------

dat_Buehner_Cheng_Clifford_2003 <- read_csv(
  "data/dat_buehner_et_al_2003.csv",
  show_col_types = FALSE
) %>%
  mutate(
    experiment = as.character(experiment),
    
    # Put preventive judgments on the negative side of the common scale
    rating = case_when(
      framing == "preventive" ~ -rating,
      TRUE ~ rating
    )
  ) %>%
  filter(
    experiment != 4
  ) %>%
  select(
    study,
    experiment,
    n,
    framing,
    N_11,
    N_10,
    N_01,
    N_00,
    rating
  )


# -----------------------------------------------------------------------------
# Liljeholm & Cheng (2009)
# -----------------------------------------------------------------------------

dat_Liljeholm_Cheng_2009 <- read_csv(
  "data/dat_liljeholm_cheng_2009.csv",
  show_col_types = FALSE
) %>%
  mutate(
    experiment = as.character(experiment),
    
    # Put preventive judgments on the negative side of the common scale
    rating = case_when(
      framing == "preventive" ~ -rating,
      TRUE ~ rating
    )
  ) %>%
  select(
    study,
    experiment,
    n,
    framing,
    N_11,
    N_10,
    N_01,
    N_00,
    rating
  )


# -----------------------------------------------------------------------------
# White (2004)
# -----------------------------------------------------------------------------

dat_White_2004 <- read_csv(
  "data/dat_white_2004.csv",
  show_col_types = FALSE
) %>%
  filter(
    !if_all(everything(), is.na)
  ) %>%
  select(
    where(~ !all(is.na(.)))
  ) %>%
  mutate(
    experiment = as.character(experiment)
  ) %>%
  rename(
    N_11 = a,
    N_10 = b,
    N_01 = c,
    N_00 = d
  ) %>%
  select(
    study,
    experiment,
    n,
    framing,
    N_11,
    N_10,
    N_01,
    N_00,
    rating
  )


# -----------------------------------------------------------------------------
# White (2003a): JEP:LMC
# -----------------------------------------------------------------------------

dat_White_2003_jep <- read_csv(
  "data/dat_white_2003_making_causal_judgments_from_proportion_confirming_instances.csv",
  show_col_types = FALSE
) %>%
  filter(
    !if_all(everything(), is.na)
  ) %>%
  select(
    where(~ !all(is.na(.)))
  ) %>%
  mutate(
    study = "White (2003a)",
    experiment = as.character(experiment)
  ) %>%
  rename(
    N_11 = a,
    N_10 = b,
    N_01 = c,
    N_00 = d
  ) %>%
  select(
    study,
    experiment,
    n,
    framing,
    N_11,
    N_10,
    N_01,
    N_00,
    rating
  )


# -----------------------------------------------------------------------------
# White (2003b): QJEP
# -----------------------------------------------------------------------------

dat_White_2003_qjep <- read_csv(
  "data/dat_white_2003_causal_judgement_evaluation_evidence.csv",
  show_col_types = FALSE
) %>%
  filter(
    !if_all(everything(), is.na)
  ) %>%
  select(
    where(~ !all(is.na(.)))
  ) %>%
  mutate(
    study = "White (2003b)",
    experiment = as.character(experiment)
  ) %>%
  rename(
    N_11 = a,
    N_10 = b,
    N_01 = c,
    N_00 = d
  ) %>%
  select(
    study,
    experiment,
    n,
    framing,
    N_11,
    N_10,
    N_01,
    N_00,
    rating
  )


# -----------------------------------------------------------------------------
# Combine studies and derive contingency measures
# -----------------------------------------------------------------------------

df_data <- bind_rows(
  dat_Perales_Shanks_2007,
  dat_Buehner_Cheng_Clifford_2003,
  dat_Liljeholm_Cheng_2009,
  dat_White_2003_qjep,
  dat_White_2003_jep,
  dat_White_2004
) %>%
  mutate(
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
# Check size of full data set
# -----------------------------------------------------------------------------

message(
  "Number of studies: ",
  n_distinct(df_data$study)
)

message(
  "Number of experiments: ",
  n_distinct(
    paste(
      df_data$study,
      df_data$experiment
    )
  )
)

message(
  "Number of conditions: ",
  nrow(df_data)
)


# =============================================================================
# Generate model predictions
# =============================================================================


# -----------------------------------------------------------------------------
# Contingency-table input
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
# Structure induction with uniform priors: noisy-OR (generative)
# -----------------------------------------------------------------------------

set.seed(0815)

SI_gen <- generate_SI_preds(
  df_freq,
  m = 10^5,
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
# Structure induction with uniform priors: noisy-AND-NOT (preventive)
# -----------------------------------------------------------------------------

set.seed(0815)

SI_prev <- generate_SI_preds(
  df_freq,
  m = 10^5,
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
# Bayesian strength under S1 with SS priors: noisy-OR (generative)
# -----------------------------------------------------------------------------

set.seed(0815)

SS_gen <- generate_SS_strength_preds(
  df_freq,
  m = 10^5,
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
# Bayesian strength under S1 with SS priors: noisy-AND-NOT (preventive)
# -----------------------------------------------------------------------------

set.seed(0815)

SS_prev <- generate_SS_strength_preds(
  df_freq,
  m = 10^5,
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

model_predictions_extended <- df_data %>%
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


# =============================================================================
# Map model quantities to the judgment scale
# =============================================================================

model_predictions_extended <- model_predictions_extended %>%
  mutate(
    
    # -------------------------------------------------------------------------
    # Uniform-prior structure posterior for bidirectional judgments
    # -------------------------------------------------------------------------
    # Three structures:
    #   S0
    #   S1-generative
    #   S1-preventive
    #
    # Each structure receives prior probability 1/3.
    
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
    
    
    # -------------------------------------------------------------------------
    # 1. Posterior link probability
    # -------------------------------------------------------------------------
    
    S1_pp_prediction = case_when(
      
      framing == "generative" ~
        S1_pp_gen * 100,
      
      framing == "preventive" ~
        S1_pp_prev * -100,
      
      framing == "bidirectional" ~
        SI_bidir_S1_gen_pp * 100 +
        SI_bidir_S1_prev_pp * -100
    ),
    
    
    # -------------------------------------------------------------------------
    # 2. Bayesian strength under S1 with uniform priors
    # -------------------------------------------------------------------------
    # For bidirectional judgments, average the generative and preventive
    # strength estimates according to their relative marginal likelihoods.
    
    S1_gen_direction_pp =
      S1_marglik_gen /
      (S1_marglik_gen + S1_marglik_prev),
    
    S1_prev_direction_pp =
      S1_marglik_prev /
      (S1_marglik_gen + S1_marglik_prev),
    
    S1_wc_prediction = case_when(
      
      framing == "generative" ~
        S1_wc_pp_gen * 100,
      
      framing == "preventive" ~
        S1_wc_pp_prev * -100,
      
      framing == "bidirectional" ~
        S1_gen_direction_pp *
        (S1_wc_pp_gen * 100) +
        S1_prev_direction_pp *
        (S1_wc_pp_prev * -100)
    ),
    
    
    # -------------------------------------------------------------------------
    # 3. Bayesian strength under S1 with SS priors
    # -------------------------------------------------------------------------
    
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
    
    
    # -------------------------------------------------------------------------
    # 4. Structure-averaged strength with uniform priors
    # -------------------------------------------------------------------------
    
    SI_mean_wc_prediction = case_when(
      
      framing == "generative" ~
        SI_mean_wc_gen * 100,
      
      framing == "preventive" ~
        SI_mean_wc_prev * -100,
      
      framing == "bidirectional" ~
        SI_bidir_S1_gen_pp *
        (S1_wc_pp_gen * 100) +
        SI_bidir_S1_prev_pp *
        (S1_wc_pp_prev * -100) +
        SI_bidir_S0_pp * 0
    ),
    
    
    # -------------------------------------------------------------------------
    # SS structure posterior: unidirectional judgments
    # -------------------------------------------------------------------------
    # S0 uses the same no-link marginal likelihood.
    # S1 is evaluated under the SS prior.
    
    SS_S0_gen_unnorm =
      0.5 * S0_marglik_gen,
    
    SS_S1_gen_unnorm =
      0.5 * SS_marglik_gen,
    
    SS_gen_normalizer =
      SS_S0_gen_unnorm +
      SS_S1_gen_unnorm,
    
    SS_S1_pp_gen =
      SS_S1_gen_unnorm /
      SS_gen_normalizer,
    
    
    SS_S0_prev_unnorm =
      0.5 * S0_marglik_prev,
    
    SS_S1_prev_unnorm =
      0.5 * SS_marglik_prev,
    
    SS_prev_normalizer =
      SS_S0_prev_unnorm +
      SS_S1_prev_unnorm,
    
    SS_S1_pp_prev =
      SS_S1_prev_unnorm /
      SS_prev_normalizer,
    
    
    # -------------------------------------------------------------------------
    # SS structure posterior: bidirectional judgments
    # -------------------------------------------------------------------------
    
    SS_bidir_S0_unnorm =
      (1 / 3) * S0_marglik_gen,
    
    SS_bidir_S1_gen_unnorm =
      (1 / 3) * SS_marglik_gen,
    
    SS_bidir_S1_prev_unnorm =
      (1 / 3) * SS_marglik_prev,
    
    SS_bidir_normalizer =
      SS_bidir_S0_unnorm +
      SS_bidir_S1_gen_unnorm +
      SS_bidir_S1_prev_unnorm,
    
    SS_bidir_S0_pp =
      SS_bidir_S0_unnorm /
      SS_bidir_normalizer,
    
    SS_bidir_S1_gen_pp =
      SS_bidir_S1_gen_unnorm /
      SS_bidir_normalizer,
    
    SS_bidir_S1_prev_pp =
      SS_bidir_S1_prev_unnorm /
      SS_bidir_normalizer,
    
    
    # -------------------------------------------------------------------------
    # 5. Structure-averaged strength with SS priors
    # -------------------------------------------------------------------------
    
    SS_mean_wc_prediction = case_when(
      
      framing == "generative" ~
        SS_S1_pp_gen *
        (SS_wc_pp_gen * 100),
      
      framing == "preventive" ~
        SS_S1_pp_prev *
        (SS_wc_pp_prev * -100),
      
      framing == "bidirectional" ~
        SS_bidir_S1_gen_pp *
        (SS_wc_pp_gen * 100) +
        SS_bidir_S1_prev_pp *
        (SS_wc_pp_prev * -100) +
        SS_bidir_S0_pp * 0
    )
  )


# =============================================================================
# Correlations
# =============================================================================


# -----------------------------------------------------------------------------
# Correlations by study
# -----------------------------------------------------------------------------

model_correlations_by_study_extended <- model_predictions_extended %>%
  select(
    study,
    human,
    S1_pp_prediction,
    S1_wc_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction,
    SS_mean_wc_prediction
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
      
      S1_wc_prediction =
        "Bayesian strength (uniform priors)",
      
      SS_wc_prediction =
        "Bayesian strength (SS priors)",
      
      SI_mean_wc_prediction =
        "Structure-averaged strength (uniform priors)",
      
      SS_mean_wc_prediction =
        "Structure-averaged strength (SS priors)"
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
# Plain average of study-level correlations, giving each study equal weight.

model_correlations_mean_study_extended <- model_correlations_by_study_extended %>%
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

model_rmse_by_study_extended <- model_predictions_extended %>%
  select(
    study,
    human,
    S1_pp_prediction,
    S1_wc_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction,
    SS_mean_wc_prediction
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
      
      S1_wc_prediction =
        "Bayesian strength (uniform priors)",
      
      SS_wc_prediction =
        "Bayesian strength (SS priors)",
      
      SI_mean_wc_prediction =
        "Structure-averaged strength (uniform priors)",
      
      SS_mean_wc_prediction =
        "Structure-averaged strength (SS priors)"
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
# Plain average of study-level RMSEs, giving each study equal weight.

model_rmse_mean_study_extended <- model_rmse_by_study_extended %>%
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
# Final extended model-comparison table
# =============================================================================

model_comparison_extended <- model_correlations_mean_study_extended %>%
  select(
    model,
    mean_r
  ) %>%
  left_join(
    model_rmse_mean_study_extended %>%
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
        "Structure-averaged strength (uniform priors)",
        "Bayesian strength (uniform priors)",
        "Posterior link probability",
        "Bayesian strength (SS priors)",
        "Structure-averaged strength (SS priors)"
      )
    )
  ) %>%
  arrange(
    model
  )


model_comparison_extended