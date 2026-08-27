# -------------------------------------------------------------------------
# Zero-contingency simulations with matched parameter-prior means
# -------------------------------------------------------------------------
#
# Target prior means for causal strength wc:
#   0.2, 0.5, 0.8
#
# SI model:
#   wc | S1 ~ Beta(alpha_c, beta_c)
#   The prior mean varies while concentration is held fixed.
#
# SS model:
#   alpha is held fixed at 5.
#   The mixture weight is varied so that the marginal prior mean of wc
#   matches the SI Beta-prior means.
#
# Both models use the same zero-contingency data.
# -------------------------------------------------------------------------
source("structure_induction_func.R") # functions for generating predictions of the structure induction (SI) model (Meder et al., 2014, 2026)
source("SS_strength_model_func.R") # functions for generating causal strength estimates under S1 given sparse and strong (SS) priors (Lu et al, 2008, Psych Review)

# -------------------------------------------------------------------------
# Step 1: Define common target prior means
# -------------------------------------------------------------------------

wc_prior_means <- c(0.2, 0.5, 0.8)


# -------------------------------------------------------------------------
# Step 2: Define Beta priors for the SI model
# -------------------------------------------------------------------------

# Total concentration of the Beta prior.
# Holding this constant means that the three priors differ in mean,
# but have the same total pseudo-count.
wc_prior_kappa <- 10


prior_design_SI <- tibble(
  wc_prior_mean = wc_prior_means
) %>%
  mutate(
    wc_prior_alpha = wc_prior_mean * wc_prior_kappa,
    wc_prior_beta = (1 - wc_prior_mean) * wc_prior_kappa
  )

prior_design_SI

# The resulting parameter priors are:
#
# mean = 0.2: Beta(2, 8)
# mean = 0.5: Beta(5, 5)
# mean = 0.8: Beta(8, 2)


# -------------------------------------------------------------------------
# Step 3: Construct zero-contingency data
# -------------------------------------------------------------------------

outcome_density_levels <- seq(
  from = 0.1,
  to = 0.9,
  by = 0.1
)

sample_sizes <- c(40, 80, 160)


design_data <- expand_grid(
  outcome_density = outcome_density_levels,
  N = sample_sizes
) %>%
  mutate(
    # Cause is present on half of the trials
    N_c = N / 2,
    
    # Cause is absent on half of the trials
    N_not_c = N / 2,
    
    # Outcome-present trials when the cause is present
    N_11 = round(N_c * outcome_density),
    
    # Outcome-absent trials when the cause is present
    N_10 = N_c - N_11,
    
    # Outcome-present trials when the cause is absent
    N_01 = round(N_not_c * outcome_density),
    
    # Outcome-absent trials when the cause is absent
    N_00 = N_not_c - N_01
  ) %>%
  select(
    outcome_density,
    N,
    N_00,
    N_01,
    N_10,
    N_11
  )

design_data

# Convert the contingency tables to the matrix format expected by
# generate_SI_preds() and generate_SS_strength_preds().
df_freq <- design_data %>%
  select(
    N_00,
    N_01,
    N_10,
    N_11
  ) %>%
  as.matrix()


# -------------------------------------------------------------------------
# Step 4: Generate SI predictions under different Beta priors over wc
# -------------------------------------------------------------------------

m <- 10^6

# Hold the structure prior fixed while varying the parameter prior.
S1_prior_fixed <- 0.5


SI_density_parameter <- map(
  seq_len(nrow(prior_design_SI)),
  function(i) {
    
    prior_i <- prior_design_SI[i, ]
    
    # Resetting the seedto make Monte Carlo draws comparable
    set.seed(0815)
    
    generate_SI_preds(
      data = df_freq,
      m = m,
      S1_prior = S1_prior_fixed,
      model_type = "noisy_or_generative",
      wc_prior_alpha = prior_i$wc_prior_alpha,
      wc_prior_beta = prior_i$wc_prior_beta
    ) %>%
      bind_cols(
        design_data %>%
          select(outcome_density)
      ) %>%
      mutate(
        wc_prior_mean = prior_i$wc_prior_mean,
        wc_prior_alpha = prior_i$wc_prior_alpha,
        wc_prior_beta = prior_i$wc_prior_beta
      ) %>%
      select(
        outcome_density,
        wc_prior_mean,
        wc_prior_alpha,
        wc_prior_beta,
        S1_prior,
        N,
        N_00,
        N_01,
        N_10,
        N_11,
        S0_marglik,
        S1_marglik,
        S0_pp,
        S1_pp,
        S1_wc_pp,
        SI_mean_wc
      )
  }
) %>%
  list_rbind() 

rownames(SI_density_parameter) <- NULL 


# Inspect SI predictions for N = 40
SI_density_parameter %>%
  filter(N == 40) %>%
  select(
    outcome_density,
    wc_prior_mean,
    S1_pp,
    S1_wc_pp,
    SI_mean_wc
  )


# -------------------------------------------------------------------------
# Step 5: Define the SS mixture-weight matching function
# -------------------------------------------------------------------------
#
# These functions assume that SS_component_mean_low() and
# SS_component_mean_high() are available from the sourced SS model file.
#
# With alpha = 5:
#
#   low component mean  ≈ 0.1932
#   high component mean ≈ 0.8068
#
# The target means 0.2, 0.5, and 0.8 are therefore attainable by changing
# the mixture weight.


SS_mixture_weight_from_mean <- function(
    target_mean,
    alpha = 5
) {
  
  mu_low <- SS_component_mean_low(alpha)
  mu_high <- SS_component_mean_high(alpha)
  
  if (any(!is.finite(target_mean))) {
    stop("target_mean must contain only finite values.")
  }
  
  if (any(target_mean < mu_low | target_mean > mu_high)) {
    stop(
      sprintf(
        paste0(
          "With alpha = %.2f, target means must lie between ",
          "%.4f and %.4f."
        ),
        alpha,
        mu_low,
        mu_high
      )
    )
  }
  
  (target_mean - mu_low) /
    (mu_high - mu_low)
}


# -------------------------------------------------------------------------
# Step 6: Construct matched SS prior settings
# -------------------------------------------------------------------------

SS_alpha_fixed <- 5

prior_design_SS <- tibble(
  wc_prior_mean = wc_prior_means
) %>%
  mutate(
    SS_alpha = SS_alpha_fixed,
    SS_mixture_weight = SS_mixture_weight_from_mean(
      target_mean = wc_prior_mean,
      alpha = SS_alpha_fixed
    ),
    SS_implied_mean_wc =
      SS_mixture_weight *
      SS_component_mean_high(SS_alpha_fixed) +
      (1 - SS_mixture_weight) *
      SS_component_mean_low(SS_alpha_fixed)
  )

prior_design_SS


# Expected mixture weights are approximately:
#
# target mean 0.2: mixture weight 0.0111
# target mean 0.5: mixture weight 0.5000
# target mean 0.8: mixture weight 0.9889


# -------------------------------------------------------------------------
# Step 7: Generate SS predictions under matched prior means
# -------------------------------------------------------------------------

SS_density_parameter <- map(
  seq_len(nrow(prior_design_SS)),
  function(i) {
    
    prior_i <- prior_design_SS[i, ]
    
    set.seed(0815)
    
    generate_SS_strength_preds(
      data = df_freq,
      m = m,
      alpha = prior_i$SS_alpha,
      model_type = "noisy_or_generative",
      mixture_weight = prior_i$SS_mixture_weight
    ) %>%
      bind_cols(
        design_data %>%
          select(outcome_density)
      ) %>%
      mutate(
        wc_prior_mean = prior_i$wc_prior_mean,
        SS_mixture_weight = prior_i$SS_mixture_weight,
        SS_implied_mean_wc = prior_i$SS_implied_mean_wc
      ) %>%
      select(
        outcome_density,
        wc_prior_mean,
        SS_alpha,
        SS_mixture_weight,
        SS_implied_mean_wc,
        N,
        N_00,
        N_01,
        N_10,
        N_11,
        SS_marglik,
        SS_logmarglik,
        SS_wc_pp,
        SS_wa_pp,
        SS_bc_pp
      )
  }
) %>%
  list_rbind()


# Inspect SS predictions for N = 40
SS_density_parameter %>%
  filter(N == 40) %>%
  select(
    outcome_density,
    wc_prior_mean,
    SS_mixture_weight,
    SS_wc_pp,
    SS_wa_pp
  )


# -------------------------------------------------------------------------
# Step 8: Optional verification of the sampled SS prior means
# -------------------------------------------------------------------------
#
# This checks the prior sampler directly, without conditioning on data.

SS_prior_check <- map(
  seq_len(nrow(prior_design_SS)),
  function(i) {
    
    prior_i <- prior_design_SS[i, ]
    
    set.seed(0815)
    
    theta_i <- SS_draw_direct(
      m = 10^6,
      alpha = prior_i$SS_alpha,
      model_type = "noisy_or_generative",
      mixture_weight = prior_i$SS_mixture_weight
    )
    
    tibble(
      target_mean_wc = prior_i$wc_prior_mean,
      mixture_weight = prior_i$SS_mixture_weight,
      sampled_mean_wc = mean(theta_i[, "wc"]),
      sampled_mean_wa = mean(theta_i[, "wa"])
    )
  }
) %>%
  list_rbind()


SS_prior_check

# -------------------------------------------------------------------------
# Step 9: Combine SI and SS predictions
# -------------------------------------------------------------------------
#
# For each matched prior mean:
#
# posterior_link:
#   P(S1 | D) from the SI model
#
# strength_SS:
#   E[wc | D, S1] under the matched SS prior
#
# structure_averaged:
#   P(S1 | D) * E[wc | D, S1] from the SI model
#
# Note that the SS and SI priors are matched on their marginal prior
# mean for wc, but differ in shape and in their treatment of wa.

df_density_parameter <- SI_density_parameter %>%
  left_join(
    SS_density_parameter %>%
      select(
        outcome_density,
        wc_prior_mean,
        N,
        N_00,
        N_01,
        N_10,
        N_11,
        SS_mixture_weight,
        SS_implied_mean_wc,
        SS_wc_pp
      ),
    by = c(
      "outcome_density",
      "wc_prior_mean",
      "N",
      "N_00",
      "N_01",
      "N_10",
      "N_11"
    )
  ) %>%
  mutate(
    posterior_link = S1_pp,
    strength_SS = SS_wc_pp,
    structure_averaged = SI_mean_wc
  )


df_density_parameter %>%
  filter(N == 40) %>%
  select(
    outcome_density,
    wc_prior_mean,
    posterior_link,
    strength_SS,
    structure_averaged
  )


# -------------------------------------------------------------------------
# Step 10: Prepare plotting data
# -------------------------------------------------------------------------

lab_wc_low <-
  "Prior~italic(E)(italic(w)[c]) == 0.2"

lab_wc_med <-
  "Prior~italic(E)(italic(w)[c]) == 0.5"

lab_wc_high <-
  "Prior~italic(E)(italic(w)[c]) == 0.8"


df_density_parameter_plot <- df_density_parameter %>%
  filter(N == 40) %>%
  mutate(
    wc_prior_label = case_when(
      wc_prior_mean == 0.2 ~ lab_wc_low,
      wc_prior_mean == 0.5 ~ lab_wc_med,
      wc_prior_mean == 0.8 ~ lab_wc_high,
      TRUE ~ NA_character_
    ),
    wc_prior_label = factor(
      wc_prior_label,
      levels = c(
        lab_wc_low,
        lab_wc_med,
        lab_wc_high
      )
    )
  )


df_density_parameter_plot_long <-
  df_density_parameter_plot %>%
  select(
    N,
    wc_prior_mean,
    wc_prior_label,
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

# inspection of some key predictions
df_density_parameter_plot %>% 
  filter(
    wc_prior_mean == 0.2,
    outcome_density %in% c(0.1, 0.5, 0.9)
  ) %>%
  select(
    outcome_density,
    wc_prior_mean,
    S1_wc_pp,
    posterior_link,
    strength_SS, 
    structure_averaged
  ) 

# -------------------------------------------------------------------------
# Step 11: Plot parameter-prior effects
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


p_outcome_density_parameter <-
  ggplot(
    df_density_parameter_plot_long,
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
    ~ wc_prior_label,
    labeller = labeller(
      wc_prior_label = label_parsed
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
    expand = expansion(
      mult = c(0, 0.02)
    )
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
    plot.title = element_text(
      size = 18
    ),
    axis.text = element_text(
      colour = "black",
      size = 12
    ),
    axis.title = element_text(
      colour = "black",
      size = 14
    ),
    axis.text.x = element_text(),
    panel.grid = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(
      colour = "black"
    ),
    strip.text = element_text(
      size = 14,
      colour = "black"
    ),
    strip.background = element_blank(),
    legend.key.height = unit(
      0.9,
      "lines"
    ),
    legend.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.box.background = element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.title = element_blank(),
    legend.position = "inside",
    legend.text = element_text(
      size = 12
    ),
    legend.position.inside = c(
      0.15,
      0.85
    ),
    panel.spacing = unit(
      1.2,
      "lines"
    )
  )


p_outcome_density_parameter

# -------------------------------------------------------------------------
# Step 12: Save plot
# -------------------------------------------------------------------------

ggsave(
  filename =
    "plots/outcome_density_parameter_prior_predictions.png",
  plot = p_outcome_density_parameter,
  dpi = 300,
  width = 10,
  height = 4
)
