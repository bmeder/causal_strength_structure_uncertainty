# ============================================================
# Lu, H., Yuille, A. L., Liljeholm, M., Cheng, P. W., & Holyoak, K. J. (2008).
# Bayesian generic priors for causal learning. Psychological Review, 115, 955-984.
#
# - Figure 2, p. 959:
#     Prior distributions over w0 and w1 with sparse-and-strong priors.
#
# - Eq. 10, p. 959:
#     Generative SS prior:
#     P(w0,w1 | gen, Graph1) ∝
#       exp[-alpha*w0 - alpha*(1-w1)] +
#       exp[-alpha*(1-w0) - alpha*w1]
#
# - Eq. 11, p. 959:
#     Preventive SS prior:
#     P(w0,w1 | prev, Graph1) ∝
#       exp[-alpha*(1-w0) - alpha*(1-w1)] +
#       exp[-alpha*(1-w0) - alpha*w1]

# ============================================================

alpha <- 5
n_grid <- 201

w0 <- seq(0, 1, length.out = n_grid)
w1 <- seq(0, 1, length.out = n_grid)

prior_gen <- outer(
  w0, w1,
  FUN = function(w0, w1) {
    exp(-alpha * w0 - alpha * (1 - w1)) +
      exp(-alpha * (1 - w0) - alpha * w1)
  }
)

prior_prev <- outer(
  w0, w1,
  FUN = function(w0, w1) {
    exp(-alpha * (1 - w0) - alpha * (1 - w1)) +
      exp(-alpha * (1 - w0) - alpha * w1)
  }
)

dw <- 1 / (n_grid - 1)

prior_gen <- prior_gen / sum(prior_gen * dw * dw)
prior_prev <- prior_prev / sum(prior_prev * dw * dw)

pdf("plots/lu2008_fig2_ss_prior_grid_persp.pdf", width = 7, height = 9)

old_par <- par(no.readonly = TRUE)
par(mfrow = c(2, 1), mar = c(3, 3, 3, 1))

persp(
  x = w0,
  y = w1,
  z = prior_gen,
  theta = -35,
  phi = 25,
  expand = 0.65,
  zlim = c(0, 12),
  ticktype = "detailed",
  border = "gray35",
  shade = 0.15,
  xlab = "w0",
  ylab = "w1",
  zlab = "Density",
  main = "A. Generative SS prior"
)

persp(
  x = w0,
  y = w1,
  z = prior_prev,
  theta = -35,
  phi = 25,
  expand = 0.65,
  zlim = c(0, 12),
  ticktype = "detailed",
  border = "gray35",
  shade = 0.15,
  xlab = "w0",
  ylab = "w1",
  zlab = "Density",
  main = "B. Preventive SS prior"
)

par(old_par)
dev.off()
