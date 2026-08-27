# Housekeeping -------------------------------------------------------------

rm(list = ls())

packages <- c("tidyverse", "gridExtra")
invisible(lapply(packages, require, character.only = TRUE))

set.seed(0815)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("structure_induction_func.r")
source("SS_strength_model_func.R") # functions for generating causal strength estimates under S1 given sparse and strong (SS) priors (Lu et al, 2008, Psych Review)

# Model predictions --------------------------------------------------------
m <- 10^6

# Input order:
# data[, 1] = N(C = 0, E = 0)
# data[, 2] = N(C = 0, E = 1)
# data[, 3] = N(C = 1, E = 0)
# data[, 4] = N(C = 1, E = 1)

# data <- matrix(0, nrow = 5, ncol = 4)
# data[1, ] <- c(5, 15, 5, 15)   # P(e|c) = P(e|¬c) = .75
# data[2, ] <- c(15, 5, 15, 5)   # P(e|c) = P(e|¬c) = .25
# data[3, ] <- c(15, 5, 5, 15)   # P(e|c) = 0.75, P(e|¬c) = .25
# data[4, ] <- c(20, 0, 10, 10)   # P(e|c) = 0.5, P(e|¬c) = 0
# data[5, ] <- c(16, 4, 10, 10)   # P(e|c) = 0.5, P(e|¬c) = 0.2


data <- matrix(0, nrow = 1, ncol = 4)
data[1, ] <- c(16, 4, 10, 10)   # P(e|c) = 0.5, P(e|¬c) = .2

# Structure Induction Model:
set.seed(0815)
SI_predictions <- generate_SI_preds(
  data = data,
  m = m,
  S1_prior = 0.5,
  model_type = "noisy_or_generative"
)

# =========================================================================
# Figure 1a: Sample Probabilities and strength estimates  -----------------
# =========================================================================
learning_data_plot <- SI_predictions %>%
  slice(1) %>%
  transmute(
    pEC          = pEC_MLE,
    pEnoC        = wa_MLE,
    deltaP       = deltaP_MLE,
    causal_power = wc_MLE
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    y_position = c(5, 4, 2, 1),
    group = factor(
      c("Probability", "Probability", "Derived measure", "Derived measure"),
      levels = c("Probability", "Derived measure")
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
    aes(fill = group, colour = group),
    orientation = "y",
    width = 0.68,
    linewidth = 1
  ) +
  scale_fill_manual(
    values = c(
      # "Probability" = scales::alpha("darkseagreen", 0.45),
      # "Derived measure" = scales::alpha("deeppink", 0.45)
      "Probability" = scales::alpha("grey80", alpha=0.6),
      "Derived measure" = scales::alpha("grey80", alpha=0.6)
    )
  ) +
  scale_colour_manual(
    values = c(
      "Probability" = "black", #darkseagreen4
      "Derived measure" = "black"  #deeppink4
    )
  ) + 
  geom_text(
    aes(
      x = value + 0.03,
      label = sprintf("%.2f", value)
    ),
    hjust = 0,
    size = 6
  ) +
  scale_y_continuous(
    name = NULL,
    breaks = c(5, 4, 2, 1),
    labels = c(
      expression(italic(P)(italic(e)*"|"*italic(c))),
      expression(italic(P)(italic(e)*"|"*"\u00AC"*italic(c))),
      "\u0394P",
      "power"
    ),
    limits = c(0.5, 5.5),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_x_continuous(
    name = NULL,
    breaks = c(0, 0.5, 1),
    labels = NULL,
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    xlim = c(0, 0.6),
    clip = "off"
  ) +
  theme_classic() +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_text(size = 18),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    axis.line=element_blank(),
    legend.position = "none",
    plot.margin = margin(3, 22, 3, 5)
  )


p_learning_data

ggsave(
  filename = "figures/fig1a_data_MLEs.png",
  plot = p_learning_data,
  dpi = 300,
  width = 3.5,
  height = 2
)



get_SI_parameter_posteriors <- function(
    data,
    m = 10^6,
    model_type = c(
      "noisy_or_generative",
      "noisy_and_not_preventive"
    )
) {
  
  model_type <- match.arg(model_type)
  
  # Accept a numeric vector or a one-row matrix/data frame
  if (is.matrix(data) || is.data.frame(data)) {
    if (nrow(data) != 1L) {
      stop("Supply exactly one contingency table.")
    }
    data <- as.numeric(data[1, ])
  } else {
    data <- as.numeric(data)
  }
  
  if (length(data) != 4L) {
    stop("Data must contain N_00, N_01, N_10, and N_11.")
  }
  
  # Same S1 likelihood selection as generate_SI_preds()
  if (model_type == "noisy_or_generative") {
    likelihood_S1_fun <- likelihood_data_S1
  } else {
    likelihood_S1_fun <- likelihood_data_S1_prev
  }
  
  # Same prior sampling as generate_SI_preds()
  theta_S1 <- matrix(data = 0, nrow = m, ncol = 3)
  theta_S1[, 1] <- rbeta(m, 1, 1)
  theta_S1[, 2] <- rbeta(m, 1, 1)
  theta_S1[, 3] <- rbeta(m, 1, 1)
  
  theta_S0 <- matrix(data = 0, nrow = m, ncol = 3)
  theta_S0[, 1] <- rbeta(m, 1, 1)
  theta_S0[, 2] <- rep(0, m)
  theta_S0[, 3] <- rbeta(m, 1, 1)
  
  # Same likelihood functions as generate_SI_preds()
  S0_lik <- likelihood_data_S0(data, theta_S0)
  S1_lik <- likelihood_S1_fun(data, theta_S1)
  
  if (!all(is.finite(S0_lik)) || !all(is.finite(S1_lik))) {
    stop("Non-finite likelihood values were generated.")
  }
  
  if (sum(S0_lik) == 0 || sum(S1_lik) == 0) {
    stop("All likelihood values are zero for at least one structure.")
  }
  
  # Normalize likelihoods to posterior weights
  S0_weight <- S0_lik / sum(S0_lik)
  S1_weight <- S1_lik / sum(S1_lik)
  
  # Same posterior means as generate_SI_preds()
  posterior_means <- tibble(
    structure = c("S1", "S0"),
    bc = c(
      S1_bc_pp(data, S1_lik, theta_S1),
      S0_bc_pp(data, S0_lik, theta_S0)
    ),
    wc = c(
      S1_wc_pp(data, S1_lik, theta_S1),
      0
    ),
    wa = c(
      S1_wa_pp(data, S1_lik, theta_S1),
      S0_wa_pp(data, S0_lik, theta_S0)
    )
  )
  
  # Keep all original prior draws with their posterior weights
  posterior_draws <- bind_rows(
    tibble(
      structure = "S1",
      bc = theta_S1[, 1],
      wc = theta_S1[, 2],
      wa = theta_S1[, 3],
      posterior_weight = S1_weight
    ),
    tibble(
      structure = "S0",
      bc = theta_S0[, 1],
      wc = theta_S0[, 2],
      wa = theta_S0[, 3],
      posterior_weight = S0_weight
    )
  )
  
  list(
    draws = posterior_draws,
    means = posterior_means
  )
}

set.seed(0815)
SI_parameter_post <- get_SI_parameter_posteriors(
  data = data,
  m = m,
  model_type = "noisy_or_generative"
)
# =========================================================================
# Prepare posterior draws and posterior means
# =========================================================================

posterior_long <- SI_parameter_post$draws %>%
  pivot_longer(
    cols = c(bc, wc, wa),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  mutate(
    structure = factor(
      structure,
      levels = c("S1", "S0"),
      labels = c(
        "Link structure S1",
        "No-link structure S0"
      )
    ),
    parameter = factor(
      parameter,
      levels = c("bc", "wc", "wa"),
      labels = c(
        "Base rate",
        "Causal strength",
        "Background strength"
      )
    )
  )


means_long <- SI_parameter_post$means %>%
  pivot_longer(
    cols = c(bc, wc, wa),
    names_to = "parameter",
    values_to = "posterior_mean"
  ) %>%
  mutate(
    structure = factor(
      structure,
      levels = c("S1", "S0"),
      labels = c(
        "Link structure S1",
        "No-link structure S0"
      )
    ),
    parameter = factor(
      parameter,
      levels = c("bc", "wc", "wa"),
      labels = c(
        "Base rate",
        "Causal strength",
        "Background strength"
      )
    )
  )


# =========================================================================
# Facet labels
# =========================================================================

parameter_labels <- c(
  "Base rate" =
    "Base~rate~italic(b)[c]",
  "Causal strength" =
    "Strength~italic(w)[c]",
  "Background strength" =
    "Background~italic(w)[a]"
)


structure_labels <- c(
  "Link structure S1" =
    "Link~structure~italic(S)[1]",
  "No-link structure S0" =
    "No-link~structure~italic(S)[0]"
)


# =========================================================================
# Figure 1c: All parameter posteriors under S1 and S0 ---------------------
# =========================================================================

# Exclude w_c under S0 from the density layer because it is fixed at zero.
posterior_continuous <- posterior_long %>%
  filter(
    !(structure == "No-link structure S0" &
        parameter == "Causal strength")
  )


# Row used to place the vertical line at w_c = 0 in the correct facet
wc_S0 <- means_long %>%
  filter(
    structure == "No-link structure S0",
    parameter == "Causal strength"
  )

p_parameter_posteriors <- ggplot() +
  
  # Continuous weighted posterior densities
  geom_density(
    data = posterior_continuous,
    aes(
      x = value,
      weight = posterior_weight,
      fill = structure,
      colour = structure
    ),
    alpha = 0.4,
    linewidth = 1,
    adjust = 1.2
  ) +
  
  # Degenerate posterior for w_c under S0:
  # w_c is fixed at zero
  geom_vline(
    data = wc_S0,
    aes(xintercept = posterior_mean),
    colour = "turquoise4",
    linewidth = 1,
    linetype = "solid"
  ) +
  
  facet_grid(
    rows = vars(structure),
    cols = vars(parameter),
    scales = "free_y",
    axes = "all_x",
    
    axis.labels = "margins",
    labeller = labeller(
      parameter = as_labeller(
        parameter_labels,
        default = label_parsed
      ),
      structure = as_labeller(
        structure_labels,
        default = label_parsed
      )
    )
  ) +
  
  scale_fill_manual(
    values = c(
      # "Link structure S1" = "darkorange1", #7570b3
      "Link structure S1" = "#0072B2", ##7570b3
      "No-link structure S0" = "#CC79A7"
    )
  ) +
  
  scale_colour_manual(
    values = c(
      "Link structure S1" = "#0072B2",
      "No-link structure S0" = "#CC79A7"
    )
  ) +
  scale_x_continuous(
    name = "Parameter value",
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1.0"),
    
    # Keeps the line at zero visible inside the plotting area
    expand = expansion(mult = c(0.025, 0.01))
  ) +
  scale_y_continuous(
    name = NULL,
    expand = expansion(mult = c(0, 0.05))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    axis.title.y = element_blank(),
    axis.text = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    strip.text = element_text(size = 18),
    strip.background = element_blank(),
    strip.text.y = element_blank(),
    panel.spacing.x = unit(2, "lines"),
    panel.spacing.y = unit(2, "lines"),
    plot.margin = margin(5, 12, 5, 5)
  )


p_parameter_posteriors


ggsave(
  filename = "figures/fig1c_parameter_posteriors.png",
  plot = p_parameter_posteriors,
  dpi = 300,
  width = 6,
  height = 3
)


# =========================================================================
# Figure 1c: Causal strength posterior over w_c  under S1 ------------------
# =========================================================================

wc_draws_S1 <- posterior_long %>%
  filter(
    structure == "Link structure S1",
    parameter == "Causal strength"
  )

wc_mean_S1 <- means_long %>%
  filter(
    structure == "Link structure S1",
    parameter == "Causal strength"
  ) %>%
  pull(posterior_mean)


p_wc_S1 <- ggplot(
  wc_draws_S1,
  aes(
    x = value,
    weight = posterior_weight
  )
) +
  
  geom_density(
    colour = "#0072B2",
    fill = "#0072B2",
    linewidth = 1,
    alpha = 0.2,
    adjust = 1.2
  ) +
  
  # geom_vline(
  #   xintercept = wc_mean_S1,
  #   linetype = "dashed",
  #   linewidth = 0.8,
  #   colour = "black"
  # ) +
  # 
  # annotate(
  #   geom = "text",
  #   x = wc_mean_S1,
  #   y = Inf,
  #   label = sprintf(
  #     "E(w[c]~'|'~D*','*S[1]) == %.2f",
  #     wc_mean_S1
  #   ),
  #   parse = TRUE,
  #   hjust = -0.05,
  #   vjust = 0.9,
  #   size = 5
  # ) +
  
  scale_x_continuous(
    name = expression(
      "Causal strength"~italic(w)[c]
    ),
    limits = c(0, 1),
    breaks = c(0, 0.5, 1),
    labels = c("0", "0.5", "1.0"),
    expand = expansion(mult = c(0, 0))
  ) +
  
  scale_y_continuous(
    name = "Posterior density",
    expand = expansion(mult = c(0, 0.05))
  ) +
  
  theme_classic() +
  
  theme(
    axis.title = element_text(size = 18),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks.x = element_line(),
    axis.ticks.y = element_blank(),
    plot.margin = margin(5, 12, 5, 5)
  )

p_wc_S1

ggsave(
  filename = "figures/fig1c_distribution_wc.png",
  plot = p_wc_S1,
  dpi = 300,
  width = 4,
  height = 4
)

# =========================================================================
# Figure 1d: Structure posteriors ------------------
# =========================================================================

posteriors_S1_S0 <- SI_predictions %>% 
  select(pEC_MLE, wa_MLE, deltaP_MLE, S0_pp, S1_pp) %>% 
  pivot_longer(cols = c(S0_pp, S1_pp), names_to = "structure", values_to = "posterior") %>% 
  mutate(structure = recode(structure, "S0_pp" = "S0", "S1_pp" = "S1"))


ggplot(posteriors_S1_S0, aes(x = structure, y = posterior, fill = structure)) +
  #geom_hline(yintercept = 0.5, linetype = "dashed") +
  # geom_bar(stat = "identity", position = position_dodge(), colour = c("turquoise4", "darkorange3"), fill = c("turquoise4", "darkorange1"), linewidth = 1.5,  alpha=0.5) + 
  geom_bar(stat = "identity", position = position_dodge(), colour = c("#CC79A7", "#0072B2"), fill = c("#CC79A7", "#0072B2"), linewidth = 1.5,  alpha=0.5) + 
  scale_x_discrete(name = "", labels = c(expression(italic(S)[0]), expression(italic(S)[1]))) +
  scale_y_continuous(name="Posterior Probability", expand = c(0,0), limits = c(0,1), breaks=c(0,0.5,1)) +
  coord_flip() +
  theme_classic() + 
  theme(legend.position = "none",
        axis.title = element_text(size=18),
        axis.text = element_text(size=18),
        strip.text = element_blank(),  
        strip.background = element_blank(),
        panel.spacing = unit(2, "lines"),
        plot.margin = margin(5, 12, 5, 5)
  )  

ggsave("figures/fig1d_structure_posteriors.png", dpi=300, width = 4, height = 4)

 
# =========================================================================
# Fig. 1e: Causal strength under S1 vs. structure-averaged causal strength -------------
# =========================================================================

wc_comparison <- SI_predictions %>%
  slice(1) %>%
  transmute(
    `Under S1` = S1_wc_pp,
    `Structure averaged` = SI_mean_wc
  ) %>%
  pivot_longer(
    cols = everything(),
    names_to = "estimate",
    values_to = "wc"
  ) %>%
  mutate(
    estimate = factor(
      estimate,
      levels = c(
        "Under S1",
        "Structure averaged"
      )
    )
  )


p_wc_comparison <- 
  ggplot(
    wc_comparison,
    aes(
      x = estimate,
      y = wc,
      fill = estimate,
      colour = estimate
    )
  ) +
  geom_col(
    width = 0.65,
    linewidth = 1
  ) +
  geom_text(
    aes(
      # label = sprintf("italic(w)[c] == %.2f", wc)
      label = sprintf("%.2f", wc)
    ),
    parse = TRUE,
    vjust = -0.5,
    size = 6,
    colour = "black"
  ) +
  scale_fill_manual(
    values = c(
      "Under S1" =
        scales::alpha("#0072B2", 0.6),
      "Structure averaged" =
        scales::alpha("#D55E00", 0.6)
    )
  ) +
  scale_colour_manual(
    values = c(
      "Under S1" = "#0072B2",
      "Structure averaged" = "#D55E00"
    )
  ) +
  scale_x_discrete(
    name = NULL,
    labels = c(
      "Under S1" = expression("Under "*italic(S)[1]),
      "Structure averaged" = "Structure averaged"
    )
  ) +
  scale_y_continuous(
    name = expression("Causal strength"~italic(w)[c]),
    limits = c(0, 0.45),
    # breaks = c(0, 0.25, 0.5),
    # labels = c("0", "0.25", "0.5"),
    expand = expansion(mult = c(0, 0.08))
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 18),
    # axis.title.y = element_blank(),
    axis.text.x = element_text(size = 18),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = margin(5, 12, 5, 5)
  )

p_wc_comparison

ggsave(
  filename = "figures/fig1e_wc_comparison.png",
  plot = p_wc_comparison,
  dpi = 300,
  width = 4.5,
  height = 4
)

