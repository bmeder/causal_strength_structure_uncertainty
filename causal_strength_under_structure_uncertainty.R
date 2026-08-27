#  Housekeeping -----------------------------------------------------------

rm(list = ls())

packages <- c('tidyverse', "gridExtra")
lapply(packages, require, character.only = TRUE)

set.seed(0815)
#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

source("structure_induction_func.r") # file with functions implementing the Structure Induction Model (Meder et al., 2014; Psychological Review)
# source("diag_reasoning_func2.r") # file with functions
#source("diag_reasoning_func_linear.r") # file with functions

# Parameters
m <- 10^5

# structure priors defaults to P(S1)=P(S0)=0.5
S1_prior = 0.5

## Input: "data", a set of contingency tables (rows) with
## data[,1] = N(C=0, E=0)
## data[,2] = N(C=0, E=1)
## data[,3] = N(C=1, E=0)
## data[,4] = N(C=1, E=1)

# example data
data <- matrix(data=0, nrow=4, ncol=4)

data[1,] <- c(5, 15, 5, 15)  # P(e|c) = 0.75, P(e|¬c) = 0.75, N = 40
data[2,] <- c(15, 5, 15, 5)  # P(e|c) = 0.25, P(e|¬c) = 0.25, N = 40
data[3,] <- c(10, 10, 0, 20) # P(e|c) = 1.00, P(e|¬c) = 0.50, N = 40
data[4,] <- c(20, 0, 10, 10) # P(e|c) = 0.50, P(e|¬c) = 0.00, N = 40

data <- matrix(data=0, nrow=2, ncol=4)
data[1, ] <- c(15, 5, 5, 15)   # generative example
data[2, ] <- c(5, 15, 10, 10)    # preventive example

SI_predictions <- generate_SI_preds(data, m)


# data <- matrix(data = 0, nrow = 4, ncol = 4)
# data[1, ] <- c(5, 15, 5, 15)    # zero contingency, high effect density
# data[2, ] <- c(15, 5, 15, 5)    # zero contingency, low effect density
# data[3, ] <- c(10, 10, 0, 20)   # generative example
# data[4, ] <- c(5, 15, 15, 5)    # preventive example

# Generative S1 vs S0:

pred_gen <- generate_SI_preds(
  data = data,
  m = m,
  S1_prior = 0.5,
  model_type = "noisy_or_generative"
)

# Preventive S1 vs S0:
pred_prev <- generate_SI_preds(
  data = data,
  m = m,
  S1_prior = 0.5,
  model_type = "noisy_and_not_preventive"
)

