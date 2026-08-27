# -------------------------------------------------------------------------
# Illustrative comparison:
# posterior link probability, conditional SS strength,
# and structure-averaged strength
# -------------------------------------------------------------------------

library("tidyverse")
library("cowplot")

source("structure_induction_func.R")
source("SS_strength_model_func.R")


# --- Example data --------------------------------------------------------
# Cell order: [N_00, N_01, N_10, N_11]
#
# Positive contingency: noisy-OR (generative)
# Zero contingency: noisy-OR (generative)
# Negative contingency: noisy-AND-NOT (preventive)

df_examples <- tibble(
  condition = factor(
    c(
      "Positive contingency",
      "Zero contingency",
      "Negative contingency"
    ),
    levels = c(
      "Positive contingency",
      "Zero contingency",
      "Negative contingency"
    )
  ),
  model_type = c(
    "generative",
    "generative",
    "preventive"
  ),
  N_00 = c(16, 5, 4),
  N_01 = c(4, 15, 16),
  N_10 = c(8, 5, 10),
  N_11 = c(12, 15, 10)
) %>%
  mutate(id = factor(row_number()))


# --- Cell-frequency matrix ----------------------------------------------
# Input order required by the model functions:
# [N_00, N_01, N_10, N_11]

df_freq <- df_examples %>%
  select(N_00, N_01, N_10, N_11) %>%
  as.matrix()


# --- Monte Carlo samples -------------------------------------------------

m <- 10^6


# --- Structure-induction predictions ------------------------------------
# Uniform structure prior: P(S1) = P(S0) = 0.5
# Uniform parameter priors within each structure

# Generative link structure: noisy-OR
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
    deltaP_MLE, wc_MLE,
    S1_pp,
    S1_wc_pp,
    SI_mean_wc
  ) %>%
  rename(
    S1_pp_gen = S1_pp,
    S1_wc_pp_gen = S1_wc_pp,
    SI_mean_wc_gen = SI_mean_wc
  )


# Preventive link structure: noisy-AND-NOT
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
    deltaP_MLE, wc_MLE,
    S1_pp,
    S1_wc_pp,
    SI_mean_wc
  ) %>%
  rename(
    S1_pp_prev = S1_pp,
    S1_wc_pp_prev = S1_wc_pp,
    SI_mean_wc_prev = SI_mean_wc
  )


# --- Conditional strength under S1 with SS priors ------------------------

# Generative link structure: noisy-OR
set.seed(0815)
SS_gen <- generate_SS_strength_preds(
  df_freq,
  m = m,
  alpha = 5,
  model_type = "noisy_or_generative"
) %>%
  mutate(id = factor(row_number())) %>%
  select(
    id,
    N, N_00, N_01, N_10, N_11,
    deltaP_MLE, wc_MLE,
    SS_wc_pp
  ) %>%
  rename(
    SS_wc_pp_gen = SS_wc_pp
  )


# Preventive link structure: noisy-AND-NOT
set.seed(0815)
SS_prev <- generate_SS_strength_preds(
  df_freq,
  m = m,
  alpha = 5,
  model_type = "noisy_and_not_preventive"
) %>%
  mutate(id = factor(row_number())) %>%
  select(
    id,
    N, N_00, N_01, N_10, N_11,
    deltaP_MLE, wc_MLE,
    SS_wc_pp
  ) %>%
  rename(
    SS_wc_pp_prev = SS_wc_pp
  )


# --- Reference table with all predictions --------------------------------

model_predictions_examples <- df_examples %>%
  left_join(
    SI_gen,
    by = c("id", "N_00", "N_01", "N_10", "N_11")
  ) %>%
  left_join(
    SI_prev,
    by = c("id", "N", "N_00", "N_01", "N_10", "N_11", "deltaP_MLE", "wc_MLE")
  ) %>%
  left_join(
    SS_gen,
    by = c("id", "N", "N_00", "N_01", "N_10", "N_11", "deltaP_MLE", "wc_MLE")
  ) %>%
  left_join(
    SS_prev,
    by = c("id", "N", "N_00", "N_01", "N_10", "N_11", "deltaP_MLE", "wc_MLE")
  )


# --- Select the appropriate functional form ------------------------------
# The positive and zero-contingency examples use the noisy-OR link model.
# The negative-contingency example uses the noisy-AND-NOT link model.

model_predictions_examples <- model_predictions_examples %>%
  mutate(
    posterior_link = case_when(
      model_type == "generative" ~ S1_pp_gen,
      model_type == "preventive" ~ S1_pp_prev
    ),
    
    strength_SS = case_when(
      model_type == "generative" ~ SS_wc_pp_gen,
      model_type == "preventive" ~ SS_wc_pp_prev
    ),
    
    structure_averaged = case_when(
      model_type == "generative" ~ SI_mean_wc_gen,
      model_type == "preventive" ~ SI_mean_wc_prev
    )
  )


# --- Reshape predictions for plotting ------------------------------------

plot_df <- model_predictions_examples %>%
  select(
    condition,
    posterior_link,
    structure_averaged,
    strength_SS
  ) %>%
  pivot_longer(
    cols = c(
      posterior_link,
      structure_averaged,
      strength_SS
    ),
    names_to = "model",
    values_to = "prediction"
  ) %>%
  mutate(
    model = recode(
      model,
      posterior_link =
        "Posterior link\nprobability",
      structure_averaged =
        "Structure-averaged\nstrength",
      strength_SS =
        "Bayesian strength\n(SS priors)"
    ),
    model = factor(
      model,
      levels = c(
        "Posterior link\nprobability",
        "Bayesian strength\n(SS priors)",
        "Structure-averaged\nstrength"
      )
    )
  )


# --- Plot ----------------------------------------------------------------

colors_models <- c(
  "Posterior link\nprobability" = "#009E73",
  "Structure-averaged\nstrength" = "#D55E00",
  "Bayesian strength\n(SS priors)" = "#56B4E9"
)


plot_model_predictions <- ggplot(
  plot_df,
  aes(
    x = model,
    y = prediction,
    fill = model
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = sprintf("%.2f", prediction)),
    hjust = -0.25,
    size = 5,
    colour = "black"
  ) +
  facet_wrap(
    ~ condition,
    ncol = 1
  ) +
  coord_flip(
    clip = "off"
  ) +
  scale_fill_manual(
    values = colors_models
  ) +
  scale_y_continuous(
    name = "Model prediction",
    limits = c(0, 1.08),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    name = NULL
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    # strip.text = element_text(size = 20,face = "bold" ),
    strip.text = element_blank(),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(5.5, 25, 5.5, 5.5)
  )


plot_model_predictions


ggsave(
  filename = "figures/fig2c_model_predictions.png",
  plot = plot_model_predictions,
  dpi = 300,
  width = 6,
  height = 9
)



# -------------------------------------------------------------------------
# Plot observed conditional effect probabilities
# for the three illustrative contingency conditions
# -------------------------------------------------------------------------

learning_data_plot <- model_predictions_examples %>%
  transmute(
    condition,
    pEC = N_11 / (N_10 + N_11),
    pEnoC = N_01 / (N_00 + N_01)
  ) %>%
  pivot_longer(
    cols = c(pEC, pEnoC),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    y_position = case_when(
      measure == "pEC"   ~ 2,
      measure == "pEnoC" ~ 1
    ),
    measure = factor(
      measure,
      levels = c("pEC", "pEnoC")
    )
  )


p_learning_data <- 
  ggplot(
  learning_data_plot,
  aes(
    x = value,
    y = y_position
  )
) +
  geom_col(
    fill = scales::alpha("grey90", alpha = 0.6),
    colour = "black",
    orientation = "y",
    width = 0.68,
    linewidth = 1
  ) +
  geom_text(
    aes(
      x = value + 0.03,
      label = sprintf("%.2f", value)
    ),
    hjust = 0,
    size = 6
  ) +
  facet_wrap(
    ~ condition,
    ncol = 1
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = c(2, 1),
    labels = c(
      expression(italic(P)(italic(e) * "|" * italic(c))),
      expression(italic(P)(italic(e) * "|" * "\u00AC" * italic(c)))
    ),
    limits = c(0.5, 2.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    name = "Sample probability",
    limits = c(0, 1.08),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(0, 1.08),
    clip = "off"
  ) +
  theme_classic() +
  theme(
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    # axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    # strip.text = element_text(size = 20,face = "bold" ),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(3, 22, 3, 5)
  )


p_learning_data


ggsave(
  filename = "figures/fig2b_sample_probabilities.png",
  plot = p_learning_data,
  dpi = 300,
  width = 4,
  height = 7
)



# -------------------------------------------------------------------------
# Shared formatting
# -------------------------------------------------------------------------

shared_plot_theme <- theme_classic() +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_blank(),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 20),
    axis.ticks = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(5.5, 15, 5.5, 5.5)
  )


# -------------------------------------------------------------------------
# Left column: observed probabilities
# -------------------------------------------------------------------------

p_learning_data <- ggplot(
  learning_data_plot,
  aes(
    x = value,
    y = y_position
  )
) +
  geom_col(
    fill = scales::alpha("grey90", alpha = 0.6),
    colour = "black",
    orientation = "y",
    width = 0.68,
    linewidth = 0.6
  ) +
  geom_text(
    aes(
      x = value + 0.03,
      label = sprintf("%.2f", value)
    ),
    hjust = 0,
    size = 7
  ) +
  facet_wrap(
    ~ condition,
    ncol = 1
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = c(2, 1),
    labels = c(
      expression(italic(P)(italic(e) * "|" * italic(c))),
      expression(italic(P)(italic(e) * "|" * "\u00AC" * italic(c)))
    ),
    limits = c(0.5, 2.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    name = "Observed probability",
    limits = c(0, 1.08),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  shared_plot_theme


# -------------------------------------------------------------------------
# Right column: model predictions
# -------------------------------------------------------------------------

plot_model_predictions <- ggplot(
  plot_df,
  aes(
    x = model,
    y = prediction,
    fill = model
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.6
  ) +
  geom_text(
    aes(label = sprintf("%.2f", prediction)),
    hjust = -0.25,
    size = 7,
    colour = "black"
  ) +
  facet_wrap(
    ~ condition,
    ncol = 1
  ) +
  scale_fill_manual(
    values = colors_models
  ) +
  scale_y_continuous(
    name = "Model prediction",
    limits = c(0, 1.08),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1"),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_discrete(
    name = NULL
  ) +
  coord_flip(
    clip = "off"
  ) +
  shared_plot_theme


# -------------------------------------------------------------------------
# Align plot panels and axes
# -------------------------------------------------------------------------

aligned_plots <- align_plots(
  p_learning_data,
  plot_model_predictions,
  align = "hv",
  axis = "tblr"
)

spacer <- ggplot() + theme_void()

fig_model_comparison <- plot_grid(
  aligned_plots[[1]],
  spacer,
  aligned_plots[[2]],
  ncol = 3,
  rel_widths = c(0.9, 0.2, 1.1),
  # labels = c("A", "", "B"),
  label_size = 18,
  label_fontface = "bold",
  hjust = 0,
  vjust = 1
)


fig_model_comparison

ggsave(
  filename = "figures/fig2_model_predictions.png",
  plot = fig_model_comparison,
  dpi = 300,
  width = 12,
  height = 8
)

ggsave(
  filename = "figures/fig2_model_predictions.pdf",
  plot = fig_model_comparison,
  width = 12,
  height = 8
)

