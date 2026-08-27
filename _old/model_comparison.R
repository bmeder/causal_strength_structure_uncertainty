library("tidyverse")
library("readxl")

source("structure_induction_func.R") # functions for generating predictions of the structure induction (SI) model (Meder et al., 2014, 2026)
source("SS_strength_model_func.R") # functions for generating causal strength estimates under S1 given sparse and strong (SS) priors (Lu et al, 2008, Psych Review)
source("meta_analysis_preprocessing.R") # reads data for meta analysis from various csv files stored in directory "data"


# --- get data for meta analysis --------------
df_data <- meta_analysis_preprocessing()

# Studies: distinct values of study
# Experiments: distinct combinations of study and experiment
# Conditions: rows in the final dataset
# Proportions: proportions of experimental conditions, not participants

# number of publications
df_data %>%
  summarise(n = n_distinct(reference)) %>%
  pull(n)

# number of conditions across all studies
df_data %>%
  summarise(n = n_distinct(id)) 


# number of contingency conditions
# contingency_summary <-
  df_data %>%
  count(sample_direction, name = "n_conditions") %>%
  mutate(
    proportion = n_conditions / sum(n_conditions),
    percent = round(100 * proportion, 1)
  )

# framing_summary <-
  df_data %>%
  count(framing, name = "n_conditions") %>%
  mutate(
    proportion = n_conditions / sum(n_conditions),
    percent = round(100 * proportion, 1)
  )

# total number of participants: sum over experiments
    df_data  %>%
    distinct(study, experiment, n) %>%
    group_by(study) %>%
    summarise(
      total_n = sum(n, na.rm = TRUE),
      n_experiments = n(),
      n_missing = sum(is.na(n)),
      .groups = "drop"
    ) %>% 
    summarise(total_n = sum(total_n))
  


# --- make df with all experimental conditions (cell frequencies) --------------
# cell input order for models: [d,c,b,a] = [N_00,N_01,N_10,N_11]
df_freq <- df_data %>% 
  select(N_00,  N_01, N_10, N_11) %>% 
  as.matrix()

# --- Generate SI predictions in the same row order as df_data --------------
# generate SI predictions under uniform priors for structure, P(S1)=P(S0)=0.5
# note that in this case causal support (Bayes Factor, BF10) and posterior link probability, S1_pp, yield identical prediction (because the prior odds are 1 it reduces to BF10)
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

# generative (noisy-NOT-AND)
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

# preventive (noisy-NOT-AND)
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
# in the bidirectional case, the average of the preventive and generative version of each model is taken (cf. Lu et al., 2008)



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
    
    # Equivalent identity:
    # SI_bidir_link_pp = 1 - SI_bidir_S0_pp
    
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
    # Structure-averaged strength
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


# -------------------------------------------------------------------------
# Correlations between model predictions and human ratings
# -------------------------------------------------------------------------

# --- Overall pooled correlation per model ------------------------------------
# This treats each condition as one observation and weights larger studies more.

model_correlations_overall <- model_predictions %>%
  select(
    human,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -human,
    names_to = "model",
    values_to = "prediction"
  ) %>%
  filter(!is.na(human), !is.na(prediction)) %>%
  mutate(
    model = recode(
      model,
      S1_pp_prediction      = "Posterior link\nprobability",
      SS_wc_prediction      = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction = "Structure-averaged\nstrength"
    )
  ) %>%
  group_by(model) %>%
  summarise(
    study = "All studies pooled",
    n = n(),
    test = list(cor.test(human, prediction, method = "pearson")),
    .groups = "drop"
  ) %>%
  mutate(
    r = map_dbl(test, ~ unname(.x$estimate)),
    CI_low = map_dbl(test, ~ .x$conf.int[1]),
    CI_high = map_dbl(test, ~ .x$conf.int[2]),
    row_type = "Overall pooled"
  ) %>%
  select(-test)

# --- Correlations by study ---------------------------------------------------
# Note: Fisher-z confidence intervals require n > 3.
# Studies with too few conditions, e.g. Shanks & Perales (2004), are retained but have no CI.

# helper function: returns NA if correlation/CI cannot be computed
safe_cor_test <- function(x, y) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  
  if (length(x) < 4 || sd(x) == 0 || sd(y) == 0) {
    return(list(
      r = NA_real_,
      CI_low = NA_real_,
      CI_high = NA_real_,
      p_value = NA_real_
    ))
  }
  
  test <- cor.test(x, y, method = "pearson")
  
  list(
    r = unname(test$estimate),
    CI_low = unname(test$conf.int[1]),
    CI_high = unname(test$conf.int[2]),
    p_value = unname(test$p.value)
  )
}
# -------------------------------------------------------------------------
# Mean correlation across studies: plain average --------------------------
# -------------------------------------------------------------------------
# This gives each study equal weight and averages raw study-level correlations.

model_correlations_by_study <- model_predictions %>%
  select(
    study,
    human,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(study, human),
    names_to = "model",
    values_to = "prediction"
  ) %>%
  filter(!is.na(human), !is.na(prediction)) %>%
  mutate(
    model = recode(
      model,
      S1_pp_prediction      = "Posterior link\nprobability",
      SS_wc_prediction      = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction = "Structure-averaged\nstrength"
    )
  ) %>%
  group_by(study, model) %>%
  summarise(
    n = n(),
    r = cor(
      human,
      prediction,
      use = "pairwise.complete.obs"
    ),
    test = list(
      if (
        n() >= 4 &&
        sd(human, na.rm = TRUE) > 0 &&
        sd(prediction, na.rm = TRUE) > 0
      ) {
        cor.test(
          human,
          prediction,
          method = "pearson"
        )
      } else {
        NULL
      }
    ),
    .groups = "drop"
  ) %>%
  mutate(
    CI_low = map_dbl(
      test,
      \(x) if (is.null(x)) NA_real_ else unname(x$conf.int[1])
    ),
    CI_high = map_dbl(
      test,
      \(x) if (is.null(x)) NA_real_ else unname(x$conf.int[2])
    )
  ) %>%
  select(-test) %>%
  arrange(model, desc(r))

# --- Mean correlation across studies: plain average --------------------------
# This gives each study equal weight and averages raw study-level correlations.

model_correlations_mean_study <- model_correlations_by_study %>%
  filter(!is.na(r), is.finite(r)) %>%
  group_by(model) %>%
  summarise(
    study = "Mean across studies",
    r_mean = mean(r),
    r_sd = sd(r),
    n = n(),
    r_se = r_sd / sqrt(n),
    CI_low = r_mean - 1.96 * r_se,
    CI_high = r_mean + 1.96 * r_se,
    n = sum(n, na.rm = TRUE),
    row_type = "mean across studies",
    .groups = "drop"
  ) %>%
  rename(r = r_mean)


# --- Combine study-level, pooled, and mean across studies ------------
study_order <- model_correlations_by_study %>%
  filter(model == "Structure-averaged\nstrength") %>%
  mutate(study_year = as.numeric(stringr::str_extract(study, "\\d{4}"))) %>%
  arrange(desc(study_year), desc(study)) %>%
  pull(study)


model_correlations_forest <- model_correlations_by_study %>%
  mutate(row_type = "study") %>%
  bind_rows(
    model_correlations_overall %>%
      mutate(
        study = "All studies pooled",
        row_type = "overall pooled"
      )
  ) %>%
  bind_rows(select(model_correlations_mean_study, -c(r_sd,r_se))) %>%
  mutate(
    study = factor(
      study,
      levels = c(
        "All studies pooled",
        "Mean across studies",
        study_order
      )
    ),
    overall_label = case_when(
      row_type == "overall pooled" ~ paste0("r = ", sprintf("%.2f", r)),
      row_type == "mean across studies" ~ paste0("r = ", sprintf("%.2f", r)),
      TRUE ~ NA_character_
    ),
    label_x = case_when(
      row_type %in% c("overall pooled", "mean across studies") & r > 0.75 ~ r - 0.18,
      row_type %in% c("overall pooled", "mean across studies") ~ r + 0.08,
      TRUE ~ NA_real_
    ),
    label_hjust = case_when(
      row_type %in% c("overall pooled", "mean across studies") & r > 0.75 ~ 1,
      row_type %in% c("overall pooled", "mean across studies") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>% 
  mutate(
    model = recode(
      model,
      S1_pp_prediction            = "Posterior link\nprobability",
      SS_wc_prediction            = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction       = "Structure-averaged\nstrength"
    )
  ) %>%
  mutate(
    model = factor(model, levels = c(
      "Structure-averaged\nstrength",
      "Posterior link\nprobability",
      "Bayesian strength\n(SS priors)"
    ))) 

# -------------------------------------------------------------------------
# Plot correlations by study --------------------------------------
# -------------------------------------------------------------------------
# order models

colors_models <- c(
  "Posterior link\nprobability" = "#009E73",
  "Structure-averaged\nstrength" = "#D55E00",
  "Bayesian strength\n(SS priors)" = "#56B4E9"
)


plot_r_by_study <- ggplot(model_correlations_forest, aes(x = r, y = study, color = model)) +
  facet_wrap(~ model, nrow = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.15,
    linewidth = 0.6,
    na.rm = TRUE,
    orientation = "y"
  ) +
  geom_point(
    aes(shape = row_type, size = row_type),
    na.rm = TRUE
  ) +
  geom_text(
    aes(x = label_x, label = overall_label, hjust = label_hjust),
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +
  coord_cartesian(xlim = c(-1, 1), clip = "off") +
  scale_color_manual(values = colors_models) +
  scale_shape_manual(
    values = c(
      "study" = 16,
      "overall pooled" = 18,
      "mean across studies" = 15
    )
  ) +
  scale_size_manual(
    values = c(
      "study" = 2.2,
      "overall pooled" = 3.8,
      "mean across studies" = 3.5
    )
  ) +
  labs(
    x = "Correlation with human ratings (r)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

ggsave("plots/model_comparison_cor.png", plot_r_by_study, dpi = 300, width= 13, height = 4)

# -------------------------------------------------------------------------
# Root mean squared error (RMSE) between model predictions and human ratings
# -------------------------------------------------------------------------


# --- Overall pooled RMSE per model -------------------------------------------
# This treats each condition as one observation and thus weights larger studies more.
# The RMSE CI is approximated with the delta method applied to MSE.

model_rmse_overall <- model_predictions %>%
  select(
    study,
    human,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(study, human),
    names_to = "model",
    values_to = "prediction"
  ) %>% 
  filter(!is.na(human), !is.na(prediction)) %>%
  mutate(
    sq_error = (human - prediction)^2,
    model = recode(
      model,
      S1_pp_prediction            = "Posterior link\nprobability",
      SS_wc_prediction            = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction       = "Structure-averaged\nstrength"
    )
  ) %>%
  group_by(model) %>%
  summarise(
    study = "All studies pooled",
    mse = mean(sq_error, na.rm = TRUE),
    mse_sd = sd(sq_error, na.rm = TRUE),
    n = n(),
    mse_se = mse_sd / sqrt(n),
    rmse = sqrt(mse),
    rmse_se = mse_se / (2 * rmse),
    CI_low = rmse - 1.96 * rmse_se,
    CI_high = rmse + 1.96 * rmse_se,
    row_type = "overall pooled",
    .groups = "drop"
  )


# --- RMSE by study ------------------------------------------------------------

model_rmse_by_study <- model_predictions %>%
  select(
    study,
    human,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(study, human),
    names_to = "model",
    values_to = "prediction"
  )  %>%
  filter(!is.na(human), !is.na(prediction)) %>%
  mutate(
    sq_error = (human - prediction)^2,
    model = recode(
      model,
      S1_pp_prediction            = "Posterior link\nprobability",
      SS_wc_prediction            = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction       = "Structure-averaged\nstrength"
    )
  ) %>%
  group_by(study, model) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    mse_sd = sd(sq_error, na.rm = TRUE),
    n = n(),
    mse_se = mse_sd / sqrt(n),
    rmse = sqrt(mse),
    rmse_se = mse_se / (2 * rmse),
    CI_low = rmse - 1.96 * rmse_se,
    CI_high = rmse + 1.96 * rmse_se,
    .groups = "drop"
  ) %>%
  arrange(model, rmse)


# --- Mean RMSE across studies ---------------------------------
# This gives each study equal weight and averages study-level RMSEs.
# The interval is a descriptive 95% interval around the mean across studies.

model_rmse_mean_study <- model_rmse_by_study %>%
  group_by(model) %>%
  summarise(
    study = "Mean across studies",
    rmse_mean = mean(rmse, na.rm = TRUE),
    rmse_sd = sd(rmse, na.rm = TRUE),
    n_studies = n(),
    rmse_se = rmse_sd / sqrt(n_studies),
    CI_low = rmse_mean - 1.96 * rmse_se,
    CI_high = rmse_mean + 1.96 * rmse_se,
    n = sum(n, na.rm = TRUE),
    row_type = "mean across studies",
    .groups = "drop"
  ) %>%
  rename(rmse = rmse_mean)


# --- Combine study-level, pooled, and mean-across-study estimates ------------

df_model_comparison_rmse <- 
  model_rmse_by_study %>%
  mutate(row_type = "study") %>%
  bind_rows(model_rmse_overall) %>%
  bind_rows(model_rmse_mean_study) %>%
  mutate(
    study = factor(
      study,
      levels = c(
        "All studies pooled",
        "Mean across studies",
        study_order
      )
    ),
    overall_label = case_when(
      row_type == "overall pooled" ~ paste0("RMSE = ", sprintf("%.1f", rmse)),
      row_type == "mean across studies" ~ paste0("RMSE = ", sprintf("%.1f", rmse)),
      TRUE ~ NA_character_
    ),
    label_x = case_when(
      row_type %in% c("overall pooled", "mean across studies") ~ rmse + 5,
      TRUE ~ NA_real_
    ),
    label_hjust = case_when(
      row_type %in% c("overall pooled", "mean across studies") ~ 0,
      TRUE ~ NA_real_
    )
  ) %>% 
  mutate(
    model = factor(model, levels = c(
      "Structure-averaged\nstrength",
      "Posterior link\nprobability",
      "Bayesian strength\n(SS priors)"
    ))) 


# --- Plot RMSE by study -------------------------------------------------------
plot_rmse_by_study <- ggplot(
  df_model_comparison_rmse,
  aes(x = rmse, y = study, fill = model, color = model)
) +
  facet_wrap(~ model, nrow = 1) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.18,
    linewidth = 0.3,
    orientation = "y",
    na.rm = TRUE
  ) +
  geom_col(
    width = 0.65,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  geom_text(
    aes(x = label_x, label = overall_label, hjust = label_hjust),
    size = 3,
    color = "black",
    na.rm = TRUE
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors_models) +
  scale_color_manual(values = colors_models) +
  labs(
    x = "Root mean squared error (RMSE)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

ggsave("plots/model_comparison_rmse.png", plot_rmse_by_study, dpi = 300, width = 12, height = 4)


# combine correlation and RMSE by study plots -------------------------------------

library(cowplot)

plot_model_comparison <- plot_grid(
  plot_r_by_study,
  plot_rmse_by_study + theme(axis.text.y = element_blank()),
  ncol = 2,
  labels = c("a) Correlations", "b) Prediction error (RMSE)"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0,
  label_y = 1,
  hjust = 0,
  vjust = 1
)

plot_model_comparison <- ggdraw() +
  draw_label(
    "Model comparison",
    x = 0.11,
    y = 0.98,
    fontface = "bold",
    size = 14,
    hjust = 1
  ) +
  draw_plot(
    plot_model_comparison,
    x = 0,
    y = 0,
    width = 1,
    height = 0.93
  )


ggsave("plots/model_comparison.png", plot_model_comparison, dpi = 300, width = 16, height = 6)
ggsave("plots/model_comparison.pdf", plot_model_comparison, width = 16, height = 6)



# # Model comparison: cor and RMSE by contingency type --------------------

# -------------------------------------------------------------------------
# Mean correlation across contingency types -------------------------------
# -------------------------------------------------------------------------
model_correlations_by_contingency_type <- model_predictions %>%
  select(
    study,
    human,
    sample_direction,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(study, human, sample_direction),
    names_to = "model",
    values_to = "prediction"
  ) %>% 
  filter(!is.na(human), !is.na(prediction)) %>%
  mutate(
    model = recode(
      model,
      S1_pp_prediction            = "Posterior link\nprobability",
      SS_wc_prediction            = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction       = "Structure-averaged\nstrength"
    )
  ) %>%
  mutate(
    model = factor(model, levels = c(
      "Structure-averaged\nstrength",
      "Posterior link\nprobability",
      "Bayesian strength\n(SS priors)"
    ))) %>% 
  mutate(
    sample_direction = recode(
      sample_direction,
      "zero"            = "Zero\ncontingency",
      "generative"            = "Positive\ncontingency",
      "preventive"            = "Negative\ncontingency"
    )
  ) %>%
  group_by(sample_direction, model) %>%
  summarise(
    n = n(),
    # Always compute the point estimate when possible
    r = cor(human, prediction, use = "pairwise.complete.obs"),
    # cor.test() only for CI/p-value when n is large enough and both variables vary
    test = list(
      if (
        n() >= 4 &&
        sd(human, na.rm = TRUE) > 0 &&
        sd(prediction, na.rm = TRUE) > 0
      ) {
        cor.test(human, prediction, method = "pearson")
      } else {
        NULL
      }
    ),
    .groups = "drop"
  ) %>%
  mutate(
    CI_low   = map_dbl(test, \(x) if (is.null(x)) NA_real_ else unname(x$conf.int[1])),
    CI_high  = map_dbl(test, \(x) if (is.null(x)) NA_real_ else unname(x$conf.int[2]))
  ) %>%
  select(-test) %>%
  arrange(model, desc(r))


# Plot correlation across contingency type -----------------------------------------

# colors_models2 <- c(
#   "Structure-averaged\nstrength" = "#009E73",
#   "Posterior link\nprobability"          = "#CC79A7", 
#   "Bayesian strength\n(SS priors)" = "#0072B2"
# )

plot_r_by_contingency <- ggplot(
  model_correlations_by_contingency_type,
  aes(x = r, y = sample_direction, color = model)
) +
  facet_wrap( ~ model, nrow = 1) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.1,
    linewidth = 0.6,
    na.rm = TRUE,
    orientation = "y"
  ) +
  geom_text(
    aes(label = sprintf("r = %.2f", r)),
    nudge_x = -0.43,
    #x=0.15,
    hjust = 0,
    size = 2.5,
    color = "black",
    na.rm = TRUE
  ) +
  geom_point(na.rm = TRUE) +
  coord_cartesian(xlim = c(0, 1), clip = "off") +
  scale_color_manual(values = colors_models) +
  scale_x_continuous(breaks = c(0,0.5,1), labels = c("0", "0.5", "1")) +
  labs(
    x = "Correlation with human ratings (Pearson r)",
    y = ""
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 12),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )

ggsave("plots/model_comparison_cor_by_contingency_type.png", plot_r_by_contingency, dpi = 300, width= 7, height = 2.5)
ggsave("plots/model_comparison_cor_by_contingency_type.pdf", plot_r_by_contingency,  width= 7, height = 2.5)



# --- RMSE by contingency type ------------------------------------------------------------
model_rmse_by_contingency <-
  model_predictions %>%
  select(
    study,
    human,
    sample_direction,
    S1_pp_prediction,
    SS_wc_prediction,
    SI_mean_wc_prediction
  ) %>%
  pivot_longer(
    cols = -c(study, human, sample_direction),
    names_to = "model",
    values_to = "prediction"
  )  %>%
  filter(!is.na(human), !is.na(prediction)) %>%
  mutate(
    sq_error = (human - prediction)^2,
    model = recode(
      model,
      S1_pp_prediction            = "Posterior link\nprobability",
      SS_wc_prediction            = "Bayesian strength\n(SS priors)",
      SI_mean_wc_prediction       = "Structure-averaged\nstrength"
    )
  ) %>%
  mutate(
    model = factor(model, levels = c(
      "Structure-averaged\nstrength",
      "Posterior link\nprobability",
      "Bayesian strength\n(SS priors)"
    ))) %>% 
  mutate(
    sample_direction = recode(
      sample_direction,
      "zero"            = "Zero\ncontingency",
      "generative"            = "Positive\ncontingency",
      "preventive"            = "Negative\ncontingency"
    )
  ) %>%
  group_by(sample_direction, model) %>%
  summarise(
    mse = mean(sq_error, na.rm = TRUE),
    mse_sd = sd(sq_error, na.rm = TRUE),
    n = n(),
    mse_se = mse_sd / sqrt(n),
    rmse = sqrt(mse),
    rmse_se = mse_se / (2 * rmse),
    CI_low = rmse - 1.96 * rmse_se,
    CI_high = rmse + 1.96 * rmse_se,
    .groups = "drop"
  ) %>%
  arrange(model, rmse)

# Plot RMSE across contingency types -----------------------------------------


p_rmse_by_contingency <- ggplot(
  model_rmse_by_contingency,
  aes(x = rmse, y = sample_direction, fill = model, color = model)
) +
  facet_wrap(~ model, nrow = 1) +
  geom_errorbar(
    aes(xmin = CI_low, xmax = CI_high),
    height = 0.18,
    linewidth = 0.3,
    orientation = "y",
    na.rm = TRUE
  ) +
  geom_col(
    width = 0.65,
    alpha = 0.85,
    na.rm = TRUE
  ) +
  geom_text(
    aes(label = sprintf("RMSE = %.1f", rmse)),
    nudge_x = 5,
    #x=0.15,
    hjust = 0,
    size = 2.5,
    color = "black",
    na.rm = TRUE
  ) +
  coord_cartesian(clip = "off") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = colors_models) +
  scale_color_manual(values = colors_models) +
  labs(
    x = "Root mean squared error (RMSE)",
    y = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.y = element_text(size = 8),
    plot.margin = margin(5.5, 35, 5.5, 5.5)
  )


ggsave("plots/model_comparison_rmse_by_contingency_type.png", p_rmse_by_contingency, dpi = 300, width= 8, height = 2.5)
ggsave("plots/model_comparison_rmse_by_contingency_type.pdf", p_rmse_by_contingency, width= 8, height = 2.5)



# combine plots for comparison by contingency type ------------------------


plot_model_comparison <- plot_grid(
  plot_r_by_study,
  plot_rmse_by_study + theme(axis.text.y = element_blank()),
  ncol = 2,
  labels = c("a) Correlations", "b) Prediction error (RMSE)"),
  label_size = 12,
  label_fontface = "bold",
  label_x = 0,
  label_y = 1,
  hjust = 0,
  vjust = 1
)

plot_model_comparison <- ggdraw() +
  draw_label(
    "Model comparison",
    x = 0.11,
    y = 0.98,
    fontface = "bold",
    size = 14,
    hjust = 1
  ) +
  draw_plot(
    plot_model_comparison,
    x = 0,
    y = 0,
    width = 1,
    height = 0.93
  )


ggsave("plots/model_comparison.png", plot_model_comparison, dpi = 300, width = 16, height = 6)
ggsave("plots/model_comparison.pdf", plot_model_comparison, width = 16, height = 6)



# --- Best-fitting model per study -------------------------------------------
# r: higher is better
# RMSE: lower is better

best_r_by_study <- model_correlations_by_study %>%
  filter(!is.na(r), is.finite(r)) %>%
  group_by(study) %>%
  filter(r == max(r, na.rm = TRUE)) %>%
  summarise(
    best_r_model = paste(model, collapse = "; "),
    best_r = max(r, na.rm = TRUE),
    .groups = "drop"
  )

best_rmse_by_study <- model_rmse_by_study %>%
  filter(!is.na(rmse), is.finite(rmse)) %>%
  group_by(study) %>%
  filter(rmse == min(rmse, na.rm = TRUE)) %>%
  summarise(
    best_rmse_model = paste(model, collapse = "; "),
    best_rmse = min(rmse, na.rm = TRUE),
    .groups = "drop"
  )

best_models_by_study <- best_r_by_study %>%
  full_join(best_rmse_by_study, by = "study") %>%
  arrange(study)

best_models_by_study

best_r_counts <- best_r_by_study %>%
  count(best_r_model, name = "n_best_r") %>%
  arrange(desc(n_best_r))

best_rmse_counts <- best_rmse_by_study %>%
  count(best_rmse_model, name = "n_best_rmse") %>%
  arrange(desc(n_best_rmse))

best_r_counts
best_rmse_counts



# -------------------------------------------------------------------------
# Correlations by condition type and model
# -------------------------------------------------------------------------

model_correlations_by_type <- 
  model_predictions %>%
  select(
    sample_direction,
    human,
    all_of(models_cor)
  ) %>%
  pivot_longer(
    cols = all_of(models_cor),
    names_to = "model_raw",
    values_to = "prediction"
  ) %>%
  filter(
    !is.na(sample_direction),
    !is.na(human),
    !is.na(prediction)
  ) %>%
  mutate(
    model = recode(model_raw, !!!model_labels),
    model = factor(model, levels = unname(model_labels[models_cor]))
  ) %>%
  group_by(sample_direction, model) %>%
  summarise(
    r = cor(human, prediction, use = "pairwise.complete.obs"),
    n = n(),
    .groups = "drop"
  ) 

model_correlations_by_type

# -------------------------------------------------------------------------
# RMSE by condition type and model
# -------------------------------------------------------------------------

model_rmse_by_type <- model_predictions %>%
  select(
    sample_direction,
    human,
    all_of(models_rmse)
  ) %>%
  pivot_longer(
    cols = all_of(models_rmse),
    names_to = "model_raw",
    values_to = "prediction"
  ) %>%
  filter(
    !is.na(sample_direction),
    !is.na(human),
    !is.na(prediction)
  ) %>%
  mutate(
    sq_error = (human - prediction)^2,
    model = recode(model_raw, !!!model_labels),
    model = factor(model, levels = unname(model_labels[models_rmse]))
  ) %>%
  group_by(sample_direction, model) %>%
  summarise(
    rmse = sqrt(mean(sq_error, na.rm = TRUE)),
    n = n(),
    .groups = "drop"
  )

model_rmse_by_type



# -------------------------------------------------------------------------
# Combine r and RMSE into one plotting table
# -------------------------------------------------------------------------

df_type_fit_plot <- bind_rows(
  model_correlations_by_type %>%
    transmute(
      sample_direction,
      model,
      metric = "Correlation",
      fit = r,
      n
    ),
  model_rmse_by_type %>%
    transmute(
      sample_direction,
      model,
      metric = "RMSE",
      fit = rmse,
      n
    )
) %>%
  mutate(
    metric = factor(metric, levels = c("Correlation", "RMSE"))
  )

df_type_fit_plot %>% 
  group_by(sample_direction) %>% 
  summarise(n=n())
# -------------------------------------------------------------------------
# Plot model fit by condition type
# -------------------------------------------------------------------------

ggplot(
  df_type_fit_plot,
  aes(x = model, y = fit, fill = model)
) +
  geom_col(
    width = 0.75,
    alpha = 0.9
  ) +
  geom_text(
    aes(label = sprintf("%.2f", fit)),
    vjust = -0.25,
    size = 3
  ) +
  facet_grid(metric ~ sample_direction, scales = "free_y") +
  scale_fill_manual(values = colors_models, drop = FALSE) +
  labs(
    x = NULL,
    y = "Model fit"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10)
  )

ggsave(
  "plots/model_fit_by_condition_type.png",
  dpi = 300,
  width = 11,
  height = 6
)


# -------------------------------------------------------------------------
# Bootstrap model-selection analysis
# -------------------------------------------------------------------------
# Goal:
#   Estimate how often each model is the best-fitting model across bootstrap
#   samples.
#
# Metrics:
#   r:    higher values indicate better fit.
#   RMSE: lower values indicate better fit.
#
# Bootstrap variants:
#   1. Study-balanced bootstrap:
#      Studies are sampled with replacement. Model fit is computed within
#      each sampled study and then averaged across studies.
#
#   2. Condition-level bootstrap:
#      Conditions are sampled with replacement. Model fit is computed across
#      the sampled conditions.
#
# Causal support is included for r because correlations do not require a
# shared response scale. Causal support is omitted for RMSE because it is not
# on the signed 0-100 response scale.

set.seed(0815)

B <- 10000


# -------------------------------------------------------------------------
# Models and labels
# -------------------------------------------------------------------------

models_cor <- c(
  # "support_prediction_fitted",
  "S1_pp_prediction",
  "S1_wc_prediction",
  "SS_wc_prediction",
  "SI_mean_wc_prediction"
)

models_rmse <- c(
  "S1_pp_prediction",
  "S1_wc_prediction",
  "SS_wc_prediction",
  "SI_mean_wc_prediction"
)

model_labels <- c(
  # support_prediction_fitted = "Causal support",
  S1_pp_prediction          = "Posterior link probability",
  S1_wc_prediction          = "Bayesian Power PC (uniform priors)",
  SS_wc_prediction          = "Bayesian Power PC (SS priors)",
  SI_mean_wc_prediction     = "Structure-averaged strength"
)

# -------------------------------------------------------------------------
# Bootstrap datasets
# -------------------------------------------------------------------------
# df_boot_r:
#   Complete cases for the correlation analysis.
#
# df_boot_rmse:
#   Complete cases for the RMSE analysis.

df_boot_r <- model_predictions %>%
  select(
    id,
    study,
    human,
    all_of(models_cor)
  ) %>%
  filter(
    !is.na(human),
    if_all(all_of(models_cor), ~ !is.na(.x))
  )

df_boot_rmse <- model_predictions %>%
  select(
    id,
    study,
    human,
    all_of(models_rmse)
  ) %>%
  filter(
    !is.na(human),
    if_all(all_of(models_rmse), ~ !is.na(.x))
  )

# Explicit sample sizes for transparency.
n_conditions_r <- nrow(df_boot_r)
n_conditions_rmse <- nrow(df_boot_rmse)

studies_r <- unique(df_boot_r$study)
studies_rmse <- unique(df_boot_rmse$study)

n_studies_r <- length(studies_r)
n_studies_rmse <- length(studies_rmse)

n_conditions_r
n_conditions_rmse
n_studies_r
n_studies_rmse

# -------------------------------------------------------------------------
# Helper function: compute model fit
# -------------------------------------------------------------------------

compute_model_fit <- function(dat, models, metric) {
  
  out <- numeric(length(models))
  names(out) <- models
  
  for (m in models) {
    
    if (metric == "r") {
      out[m] <- cor(dat$human, dat[[m]], use = "pairwise.complete.obs")
    }
    
    if (metric == "rmse") {
      out[m] <- sqrt(mean((dat$human - dat[[m]])^2, na.rm = TRUE))
    }
  }
  
  out
}


# -------------------------------------------------------------------------
# Study-balanced bootstrap: precompute study-level r
# -------------------------------------------------------------------------
# One row per study, one column per model.

study_r_values <- matrix(
  NA_real_,
  nrow = n_studies_r,
  ncol = length(models_cor)
)

rownames(study_r_values) <- studies_r
colnames(study_r_values) <- models_cor

for (s in studies_r) {
  
  dat_s <- df_boot_r %>%
    filter(study == s)
  
  study_r_values[s, ] <- compute_model_fit(
    dat = dat_s,
    models = models_cor,
    metric = "r"
  )
}

study_r_values

# -------------------------------------------------------------------------
# Study-balanced bootstrap: precompute study-level RMSE
# -------------------------------------------------------------------------
# One row per study, one column per model.

study_rmse_values <- matrix(
  NA_real_,
  nrow = n_studies_rmse,
  ncol = length(models_rmse)
)

rownames(study_rmse_values) <- studies_rmse
colnames(study_rmse_values) <- models_rmse

for (s in studies_rmse) {
  
  dat_s <- df_boot_rmse %>%
    filter(study == s)
  
  study_rmse_values[s, ] <- compute_model_fit(
    dat = dat_s,
    models = models_rmse,
    metric = "rmse"
  )
}

study_rmse_values

# -------------------------------------------------------------------------
# Study-balanced bootstrap: storage
# -------------------------------------------------------------------------

boot_study_r_values <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(models_cor)
)

colnames(boot_study_r_values) <- models_cor

boot_study_r_winner <- character(B)


boot_study_rmse_values <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(models_rmse)
)

colnames(boot_study_rmse_values) <- models_rmse

boot_study_rmse_winner <- character(B)


# -------------------------------------------------------------------------
# Study-balanced bootstrap: r
# -------------------------------------------------------------------------
# Each bootstrap sample contains n_studies_r sampled studies.
# Fit is averaged across sampled studies.
# Winner = model with highest mean r.

for (b in seq_len(B)) {
  
  sampled_studies <- sample(
    studies_r,
    size = n_studies_r,
    replace = TRUE
  )
  
  r_b <- colMeans(
    study_r_values[sampled_studies, , drop = FALSE],
    na.rm = TRUE
  )
  
  boot_study_r_values[b, ] <- r_b
  
  boot_study_r_winner[b] <- models_cor[which.max(r_b)]
}

# -------------------------------------------------------------------------
# Study-balanced bootstrap: RMSE
# -------------------------------------------------------------------------
# Each bootstrap sample contains n_studies_rmse sampled studies.
# Fit is averaged across sampled studies.
# Winner = model with lowest mean RMSE.

for (b in seq_len(B)) {
  
  sampled_studies <- sample(
    studies_rmse,
    size = n_studies_rmse,
    replace = TRUE
  )
  
  rmse_b <- colMeans(
    study_rmse_values[sampled_studies, , drop = FALSE],
    na.rm = TRUE
  )
  
  boot_study_rmse_values[b, ] <- rmse_b
  
  boot_study_rmse_winner[b] <- models_rmse[which.min(rmse_b)]
}

# -------------------------------------------------------------------------
# Condition-level bootstrap: storage
# -------------------------------------------------------------------------

boot_condition_r_values <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(models_cor)
)

colnames(boot_condition_r_values) <- models_cor

boot_condition_r_winner <- character(B)


boot_condition_rmse_values <- matrix(
  NA_real_,
  nrow = B,
  ncol = length(models_rmse)
)

colnames(boot_condition_rmse_values) <- models_rmse

boot_condition_rmse_winner <- character(B)


# -------------------------------------------------------------------------
# Condition-level bootstrap: r
# -------------------------------------------------------------------------
# Each bootstrap sample contains n_conditions_r sampled conditions.
# Winner = model with highest r.

for (b in seq_len(B)) {
  
  sampled_rows <- sample(
    seq_len(n_conditions_r),
    size = n_conditions_r,
    replace = TRUE
  )
  
  dat_b <- df_boot_r[sampled_rows, ]
  
  r_b <- compute_model_fit(
    dat = dat_b,
    models = models_cor,
    metric = "r"
  )
  
  boot_condition_r_values[b, ] <- r_b
  
  boot_condition_r_winner[b] <- models_cor[which.max(r_b)]
}


# -------------------------------------------------------------------------
# Condition-level bootstrap: RMSE
# -------------------------------------------------------------------------
# Each bootstrap sample contains n_conditions_rmse sampled conditions.
# Winner = model with lowest RMSE.

for (b in seq_len(B)) {
  
  sampled_rows <- sample(
    seq_len(n_conditions_rmse),
    size = n_conditions_rmse,
    replace = TRUE
  )
  
  dat_b <- df_boot_rmse[sampled_rows, ]
  
  rmse_b <- compute_model_fit(
    dat = dat_b,
    models = models_rmse,
    metric = "rmse"
  )
  
  boot_condition_rmse_values[b, ] <- rmse_b
  
  boot_condition_rmse_winner[b] <- models_rmse[which.min(rmse_b)]
}

# -------------------------------------------------------------------------
# Helper: winner proportions
# -------------------------------------------------------------------------

winner_proportions <- function(winners, models, metric_label) {
  
  out <- tibble(
    model = factor(winners, levels = models)
  ) %>%
    count(model, name = "n_best", .drop = FALSE) %>%
    mutate(
      model = recode(as.character(model), !!!model_labels),
      model = factor(model, levels = unname(model_labels[models])),
      B = B,
      prop_best = n_best / B,
      metric = metric_label
    ) %>%
    arrange(desc(prop_best))
  
  out
}

# -------------------------------------------------------------------------
# Winner proportions: study-balanced bootstrap
# -------------------------------------------------------------------------

bootstrap_study_r_winners <- winner_proportions(
  winners = boot_study_r_winner,
  models = models_cor,
  metric_label = "r"
) %>%
  mutate(bootstrap_unit = "study-balanced")

bootstrap_study_rmse_winners <- winner_proportions(
  winners = boot_study_rmse_winner,
  models = models_rmse,
  metric_label = "RMSE"
) %>%
  mutate(bootstrap_unit = "study-balanced")

bootstrap_study_r_winners
bootstrap_study_rmse_winners

# -------------------------------------------------------------------------
# Winner proportions: condition-level bootstrap
# -------------------------------------------------------------------------

bootstrap_condition_r_winners <- winner_proportions(
  winners = boot_condition_r_winner,
  models = models_cor,
  metric_label = "r"
) %>%
  mutate(bootstrap_unit = "condition-level")

bootstrap_condition_rmse_winners <- winner_proportions(
  winners = boot_condition_rmse_winner,
  models = models_rmse,
  metric_label = "RMSE"
) %>%
  mutate(bootstrap_unit = "condition-level")

bootstrap_condition_r_winners
bootstrap_condition_rmse_winners

# -------------------------------------------------------------------------
# Combined winner-proportion table
# -------------------------------------------------------------------------

bootstrap_winner_proportions <- bind_rows(
  bootstrap_study_r_winners,
  bootstrap_study_rmse_winners,
  bootstrap_condition_r_winners,
  bootstrap_condition_rmse_winners
) %>%
  select(bootstrap_unit, metric, model, n_best, B, prop_best)

bootstrap_winner_proportions

# -------------------------------------------------------------------------
# plot winner-proportion
# -------------------------------------------------------------------------

df_bootstrap_winners_plot <- bootstrap_winner_proportions %>%
  mutate(
    model = recode(
      as.character(model),
      "Structure-averaged strength" = "Structure-averaged strength\n",
      "Causal support" = "Causal support\n",
      "Posterior link probability" = "Posterior link probability\n",
      "Bayesian Power PC (uniform priors)" = "Bayesian Power PC\n(uniform priors)",
      "Bayesian Power PC (SS priors)" = "Bayesian strength\n(SS priors)"
    ),
    model = factor(model, levels = names(colors_models)),
    metric = factor(
      metric,
      levels = c("r", "RMSE"),
      labels = c("Correlation", "RMSE")
    ),
    bootstrap_unit = factor(
      bootstrap_unit,
      levels = c("study-balanced", "condition-level"),
      labels = c("Study-balanced bootstrap", "Condition-level bootstrap")
    )
  )

ggplot(
  df_bootstrap_winners_plot,
  aes(x = model, y = prop_best, fill = model)
) +
  geom_col(width = 0.75, alpha = 0.9) +
  geom_text(
    aes(label = sprintf("%.2f", prop_best)),
    vjust = -0.3,
    size = 3
  ) +
  # facet_grid(bootstrap_unit ~ metric) +
  facet_wrap(bootstrap_unit ~ metric, nrow=1) +
  scale_fill_manual(values = colors_models, drop = FALSE) +
  scale_y_continuous(
    limits = c(0, 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = NULL,
    y = "Proportion best-fitting"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 10, face = "bold"),
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    axis.title.y = element_text(size = 10)
  )




