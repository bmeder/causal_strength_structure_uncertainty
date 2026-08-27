# -------------------------------------------------------------------------
# Zero-contingency simulations across outcome density
# -------------------------------------------------------------------------
# All conditions satisfy:
# P(e | c) = P(e | not c)
#
# The link hypothesis is noisy-OR (generative) S1 versus S0.
#
# Model predictions:
# 1. Posterior link probability
# 2. Conditional Bayesian strength under S1 with SS priors
# 3. Structure-averaged strength
# -------------------------------------------------------------------------

source("structure_induction_func.R") # functions for generating predictions of the structure induction (SI) model (Meder et al., 2014, 2026)
source("SS_strength_model_func.R") # functions for generating causal strength estimates under S1 given sparse and strong (SS) priors (Lu et al, 2008, Psych Review)

packages <- c('gridExtra', 'tidyverse', "cowplot", "kableExtra", "scales")
install.packages(setdiff(packages, rownames(installed.packages())))
sapply(packages, require, character.only = TRUE)

# --- Simulation parameters ----------------------------------------------

m <- 10^5

outcome_density_levels <- seq(0.1, 0.9, by = 0.1)
sample_sizes <- c(40, 80, 160)
S1_priors <- c(0.1, 0.5, 0.9)

# -------------------------------------------------------------------------
# Construct zero-contingency data
# -------------------------------------------------------------------------
# P(c) = P(not c) = .5
#
# Cell order:
# [N_00, N_01, N_10, N_11]
#
# N_00 = N(not c, not e)
# N_01 = N(not c, e)
# N_10 = N(c, not e)
# N_11 = N(c, e)

design_data <- expand_grid(
  outcome_density = outcome_density_levels,
  N = sample_sizes
) %>%
  mutate(
    N_c = N / 2,
    N_not_c = N / 2,
    
    N_11 = round(N_c * outcome_density),
    N_10 = N_c - N_11,
    
    N_01 = round(N_not_c * outcome_density),
    N_00 = N_not_c - N_01
  ) %>%
  select(
    outcome_density,
    N,
    N_00, N_01, N_10, N_11
  )


df_freq <- design_data %>%
  select(N_00, N_01, N_10, N_11) %>%
  as.matrix()


# -------------------------------------------------------------------------
# Structure-induction predictions
# -------------------------------------------------------------------------
# SI predictions are computed separately for each structure prior.
# Parameter priors within S0 and S1 are uniform.
#
# The zero-contingency hypothesis space contains:
# S0: no link
# S1: noisy-OR generative link

  SI_density <- 
  map_dfr(
    S1_priors,
    function(prior_i) {
      
      set.seed(0815)
      
      generate_SI_preds(
        data = df_freq,
        m = m,
        S1_prior = prior_i,
        model_type = "noisy_or_generative"
      ) %>%
        bind_cols(
          design_data %>%
            select(outcome_density)
        ) %>%
        mutate(
          S1_prior = prior_i
        ) %>%
        select(
          outcome_density,
          S1_prior,
          N,
          N_00, N_01, N_10, N_11,
          S0_marglik, S1_marglik,
          S0_pp, S1_pp,
          S1_wc_pp,
          SI_mean_wc
        )
    }
)

# -------------------------------------------------------------------------
# Conditional strength under S1 with sparse-and-strong priors
# -------------------------------------------------------------------------
# This model conditions on S1 and therefore does not use a structure prior.
# It is computed once for each outcome-density × sample-size condition.

  set.seed(0815)
  
  SS_density <- generate_SS_strength_preds(
    data = df_freq,
    m = m,
    alpha = 5,
    model_type = "noisy_or_generative"
  ) %>%
    bind_cols(
      design_data %>%
        select(outcome_density)
    ) %>%
    select(
      outcome_density,
      N,
      N_00, N_01, N_10, N_11,
      SS_wc_pp,
      SS_marglik,
      SS_logmarglik
    )
  


# -------------------------------------------------------------------------
# Combine model predictions
# -------------------------------------------------------------------------

df_density <- SI_density %>%
  left_join(
    SS_density,
    by = c(
      "outcome_density",
      "N",
      "N_00", "N_01", "N_10", "N_11"
    )
  ) %>%
  mutate(
    posterior_link = S1_pp,
    strength_SS = SS_wc_pp,
    structure_averaged = SI_mean_wc
  )
  
  # -------------------------------------------------------------------------
  # Prepare plotting data
  # -------------------------------------------------------------------------
  
  lab_S1_low  <- "Prior~italic(P)(italic(S)[1]) == 0.1"
  lab_S1_med  <- "Prior~italic(P)(italic(S)[1]) == 0.5"
  lab_S1_high <- "Prior~italic(P)(italic(S)[1]) == 0.9"
  
  
  df_density_plot <- df_density %>%
    filter(N == 40) %>%
    mutate(
      S1_prior_label = case_when(
        S1_prior == 0.1 ~ lab_S1_low,
        S1_prior == 0.5 ~ lab_S1_med,
        S1_prior == 0.9 ~ lab_S1_high
      ),
      S1_prior_label = factor(
        S1_prior_label,
        levels = c(
          lab_S1_low,
          lab_S1_med,
          lab_S1_high
        )
      )
    )
  
  
  df_density_plot_long <- df_density_plot %>%
    select(
      N,
      S1_prior,
      S1_prior_label,
      outcome_density,
      posterior_link,
      strength_SS,
      structure_averaged
    ) %>%
    pivot_longer(
      cols = c(
        posterior_link,
        strength_SS,
        structure_averaged
      ),
      names_to = "model",
      values_to = "prediction"
    ) %>%
    mutate(
      model = recode(
        model,
        posterior_link =
          "Posterior link probability",
        strength_SS =
          "Bayesian Power (SS priors)",
        structure_averaged =
          "Structure-averaged strength"
      ),
      model = factor(
        model,
        levels = c(
          "Structure-averaged strength",
          "Bayesian Power (SS priors)",
          "Posterior link probability"
        )
      )
    )
  
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  colors_models_density <- c(
    "Posterior link probability" = "#009E73",
    "Bayesian Power (SS priors)" = "#56B4E9",
    "Structure-averaged strength" = "#D55E00"
  )
  
  shapes_models_density <- c(
    "Posterior link probability" = 21,
    "Bayesian Power (SS priors)" = 22,
    "Structure-averaged strength" = 24
  )
  
  p_outcome_density <- 
    ggplot(
      df_density_plot_long,
      aes(
        x = outcome_density,
        y = prediction * 100,
        group = model,
        colour = model,
        fill = model,
        shape = model
      )
    ) +
    facet_wrap(
      ~ S1_prior_label,
      labeller = labeller(
        S1_prior_label = label_parsed
      )
    ) +
    geom_line(
      linewidth = 1
    ) +
    geom_point(
      size = 3,
      colour = "black",
      stroke = 0.5
    ) +
    scale_x_continuous(
      name = expression(
        "Outcome density:" ~
          italic(P)(italic(e)) * "=" *
          italic(P)(italic(e) ~ "|" ~ italic(c)) * "=" *
          italic(P)(italic(e) ~ "|" ~ "¬" ~ italic(c))
      ),
      breaks = outcome_density_levels,
      labels = outcome_density_levels
    ) +
    scale_y_continuous(
      name = "Model prediction",
      limits = c(0, 100),
      breaks = c(0, 25, 50, 75, 100),
      expand = expansion(mult = c(0, 0.02))
    ) +
    scale_colour_manual(
      values = colors_models_density
    ) +
    scale_fill_manual(
      values = colors_models_density
    ) +
    scale_shape_manual(
      values = shapes_models_density
    ) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 18),
      axis.text = element_text(color = "black", size =12),  
      axis.title = element_text(color = "black", size =14),
      axis.text.x = element_text(),  
      panel.grid = element_blank(),  
      axis.ticks.y = element_blank(), 
      axis.line.y = element_line(color = "black"),
      strip.text =  element_text(size = 14, colour="black"), 
      strip.background =element_blank(),
      legend.key.height = unit(0.9, "lines"),
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.box.background = element_rect(fill = "transparent", colour = NA),
      legend.title = element_blank(),
      legend.position="inside",
      legend.text = element_text(size = 12),
      legend.position.inside = c(0.15, 0.85),
      panel.spacing = unit(1.2, "lines")
    ) 
  
  p_outcome_density
  
  
  ggsave(
    filename = "plots/outcome_density_predictions.png",
    plot = p_outcome_density,
    dpi = 300,
    width = 10,
    height = 4
  )
