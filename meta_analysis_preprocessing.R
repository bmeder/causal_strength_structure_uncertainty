# Data preprocessing ------------------------------------------------------
# The starting dataset is the contingency-table compilation reported by
# Perales and Shanks (2007). It is supplemented with additional experiments
# from Buehner et al. (2003), Liljeholm and Cheng (2009), and White (2003,
# 2004). Information not consistently available in the original compilation,
# including sample size, experimental framing, cell frequencies, and mean
# causal ratings, was checked against and extracted from the primary studies.
# Variables are harmonized across datasets, preventive ratings are recoded
# onto a common signed scale, and descriptive contingency measures are added.

# Data sources ------------------------------------------------------------
# The dataset was based initially on the contingency-table compilation in
# Perales and Shanks (2007). The information was checked and supplemented
# using the original studies, and additional conditions were added from
# Buehner et al. (2003), Liljeholm and Cheng (2009), and White (2003, 2004).
#
# References:
#
# Buehner, M. J., Cheng, P. W., & Clifford, D. (2003).
#   From covariation to causation: A test of the assumption of causal power.
#   Journal of Experimental Psychology: Learning, Memory, and Cognition,
#   29, 1119–1140.
#
# Collins, D. J., & Shanks, D. R. (2006).
#   Conformity to the power PC theory of causal induction depends on type
#   of probe question. Quarterly Journal of Experimental Psychology,
#   59, 225–232.
#
# Liljeholm, M., & Cheng, P. W. (2009).
#   The influence of virtual sample size on confidence and causal-strength
#   judgments. Journal of Experimental Psychology: Learning, Memory, and
#   Cognition, 35, 157–172.
#
# Lober, K., & Shanks, D. R. (2000).
#   Is causal induction based on causal power? Critique of Cheng (1997).
#   Psychological Review, 107, 195–212.
#
# Perales, J. C., & Shanks, D. R. (2003).
#   Normative and descriptive accounts of the influence of power and
#   contingency on causal judgement. Quarterly Journal of Experimental
#   Psychology, 56A, 977–1007.
#
# Perales, J. C., & Shanks, D. R. (2004, May).
#   The cause–density effect as a tool to discriminate between causal
#   learning models. Paper presented at the Special Interest Meeting on
#   Human Contingency Learning, Lignely, Belgium.
#
# Perales, J. C., & Shanks, D. R. (2007).
#   Models of covariation-based causal judgment: A review and synthesis.
#   Psychonomic Bulletin & Review, 14, 577–596.
#
# Shanks, D. R. (2002).
#   Tests of the power PC theory of causal induction with negative
#   contingencies. Experimental Psychology, 49, 81–88.
#
# Vallée-Tourangeau, F., Murphy, R. A., Drew, S., & Baker, A. G. (1998).
#   Judging the importance of constant and variable candidate causes:
#   A test of the power PC theory. Quarterly Journal of Experimental
#   Psychology, 51A, 65–84.
#
# Wasserman, E. A., Kao, S.-F., Van Hamme, L. J., Katagiri, M., &
#   Young, M. E. (1996). Causation and association. In D. R. Shanks,
#   K. J. Holyoak, & D. L. Medin (Eds.), The psychology of learning and
#   motivation: Vol. 34. Causal learning (pp. 207–264). Academic Press.
#
# White, P. A. (2003).
#   Causal judgement as evaluation of evidence: The use of confirmatory
#   and disconfirmatory information. Quarterly Journal of Experimental
#   Psychology, 56A, 491–513.
#
# White, P. A. (2003).
#   Making causal judgments from the proportion of confirming instances:
#   The pCI rule. Journal of Experimental Psychology: Learning, Memory,
#   and Cognition, 29, 710–727.
#
# White, P. A. (2004).
#   Causal judgment from contingency information: A systematic test of
#   the pCI rule. Memory & Cognition, 32, 353–368.

# the meta-analyis by Perales and Shanks (2007; Appendix) includes key info from several studies 
# added data from multiple studies


meta_analysis_preprocessing <- function() {
  dat_Perales_Shanks_2007 <- read_csv("data/dat_Perales_Shanks_2007.csv") %>% 
    filter(!if_all(everything(), is.na)) %>%
    select(where(~ !all(is.na(.)))) %>% 
    rename(
      N_11 = a,
      N_10 = b,
      N_01 = c,
      N_00 = d
    ) %>%
    mutate(source = "Perales, J. C., & Shanks, D. R. (2007). Models of covariation-based causal judgment: A review and synthesis. Psychonomic Bulletin & Review, 14, 577-596.") %>% 
    mutate(
      reference = recode_values(
        study,
        "Perales & Shanks (2004)" ~
          "Perales, J. C., & Shanks, D. R. (2004, May). The cause–density effect as a tool to discriminate between causal learning models. Paper presented at the Special Interest Meeting on Human Contingency Learning, Lignely, Belgium.",
        "Perales & Shanks (2003)" ~
          "Perales, J. C., & Shanks, D. R. (2003). Normative and descriptive accounts of the influence of power and contingency on causal judgement. Quarterly Journal of Experimental Psychology, 56A, 977-1007.",
        "Shanks (2002)" ~
          "Shanks, D. R. (2002). Tests of the power PC theory of causal induction with negative contingencies. Experimental Psychology, 49, 81–88.",
        "Lober & Shanks (2000)" ~
          "Lober, K., & Shanks, D. R. (2000). Is causal induction based on causal power? Critique of Cheng (1997). Psychological Review, 107, 195–212.",
        "Collins & Shanks (2006)" ~
          "Collins, D. J., & Shanks, D. R. (2006). Conformity to the power PC theory of causal induction depends on type of probe question. Quarterly Journal of Experimental Psychology, 59, 225–232.",
        "Vallée-Tourangeu et al. (1998)" ~
          "Vallée-Tourangeau, F., Murphy, R. A., Drew, S., & Baker, A. G. (1998). Judging the importance of constant and variable candidate causes: A test of the power PC theory. Quarterly Journal of Experimental Psychology, 51A, 65–84.",
        "White (2003)" ~
          "White, P. A. (2003a). Causal judgement as evaluation of evidence: The use of confirmatory and disconfirmatory information. Quarterly Journal of Experimental Psychology, 56A, 491–513.",
        "Buehner, Cheng & Clifford (2003) (generative component)" ~
          "Buehner, M. J., Cheng, P. W., & Clifford, D. (2003). From covariation to causation: A test of the assumption of causal power. Journal of Experimental Psychology: Learning, Memory, & Cognition, 29, 1119–1140.",
        "Buehner, Cheng & Clifford (2003) (preventive component)" ~
          "Buehner, M. J., Cheng, P. W., & Clifford, D. (2003). From covariation to causation: A test of the assumption of causal power. Journal of Experimental Psychology: Learning, Memory, & Cognition, 29, 1119–1140.",
        "Wasserman et al. (1996)" ~
          "Wasserman, E. A., Kao, S.-F., Van-Hamme, L. J., Katagiri, M., & Young, M. E. (1996). Causation and association. In D. R. Shanks, K. J. Holyoak, & D. L. Medin (Eds.), The psychology of learning and motivation: Vol. 34. Causal learning (pp. 207–264). San Diego: Academic Press.",
        default = study
      )
    ) %>% 
    filter(!str_detect(study, "Buehner")) %>% #remove exp.1 from Buehner et al 2003 is added separately below
    filter(!str_detect(study, "White")) %>% #remove White 2003 which is added separately below
    select(
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      rating,
      reference
    ) 
  
  # Buehner, M. J., Cheng, P. W., & Clifford, D. (2003). From covariation to causation: A test of the assumption of causal power. Journal of Experimental Psychology: Learning, Memory, & Cognition, 29, 1119–1140.
  # study 1 of Buehner et al is included in the Perales and Shanks papers, Exp. 2 and 3 not
  
  dat_Buehner_Cheng_Clifford_2003 <- read_csv("data/dat_buehner_et_al_2003.csv") %>% 
    mutate(experiment = as.character(experiment)) %>% 
    # Buehner, Cheng, & Clifford, (2003) reported preventive strength on a non-negative 0–100 scale.
    # Recode preventive ratings to a common signed scale, where prevention is negative.
    mutate(
      rating = case_when(
        framing == "preventive" ~ -rating,
        TRUE ~ rating
      )) %>% 
    mutate(reference = "Buehner, M. J., Cheng, P. W., & Clifford, D. (2003). From covariation to causation: A test of the assumption of causal power. Journal of Experimental Psychology: Learning, Memory, & Cognition, 29, 1119–1140.") %>% 
    filter(experiment != 4) %>%  # np data was reported for experiment 4
    select(
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      rating,
      reference
    ) 
  
  # Liljeholm, M., & Cheng, P. W. (2009). The influence of virtual sample size on confidence and causal-strength judgments. Journal of Experimental Psychology: Learning, Memory, and Cognition, 35(1), 157-172.
  dat_Liljeholm_Cheng_2009 <- read_csv("data/dat_liljeholm_cheng_2009.csv") %>% 
    mutate(experiment = as.character(experiment)) %>% 
    # Liljeholm & Cheng (2009) reported preventive strength on a non-negative 0–100 scale.
    # Recode preventive ratings to a common signed scale, where prevention is negative.
    mutate(
      rating = case_when(
        framing == "preventive" ~ -rating,
        TRUE ~ rating
      )) %>% 
    mutate(reference = "Liljeholm, M., & Cheng, P. W. (2009). The influence of virtual sample size on confidence and causal-strength judgments. Journal of Experimental Psychology: Learning, Memory, and Cognition, 35(1), 157-172.") %>% 
    select(
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      question,
      rating,
      reference
    ) 
  
  
  # White, P. A. (2004). Causal judgment from contingency information: A systematic test of the pCI rule. Memory & Cognition, 32(3), 353-368.
  dat_White_2004 <- read_csv("data/dat_white_2004.csv") %>% 
    mutate(experiment = as.character(experiment)) %>%
    rename(
      N_11 = a,
      N_10 = b,
      N_01 = c,
      N_00 = d
    ) %>%
    select(where(~ !all(is.na(.)))) %>% 
    mutate(reference = "White, P. A. (2004). Causal judgment from contingency information: A systematic test of the pCI rule. Memory & Cognition, 32(3), 353-368.") %>% 
    select(
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      rating,
      reference
    ) 
  
  
  
  # White, P. A. (2003). Making causal judgments from the proportion of confirming instances: the pCI rule. Journal of Experimental Psychology: Learning, Memory, and Cognition, 29(4), 710-727 .
  dat_White_2003_jep <- read_csv("data/dat_white_2003_making_causal_judgments_from_proportion_confirming_instances.csv") %>% 
    mutate(study = "White (2003a)") %>% 
    mutate(experiment = as.character(experiment)) %>%
    rename(
      N_11 = a,
      N_10 = b,
      N_01 = c,
      N_00 = d
    ) %>% 
    select(where(~ !all(is.na(.)))) %>% 
    mutate(reference = "White, P. A. (2003). Making causal judgments from the proportion of confirming instances: the pCI rule. Journal of Experimental Psychology: Learning, Memory, and Cognition, 29(4), 710-727 .") %>% 
    select(
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      rating,
      reference
    ) 
  
  
  
  # White, P. A. (2003). Causal judgement as evaluation of evidence: The use of confirmatory and disconfirmatory information. The Quarterly Journal of Experimental Psychology Section A, 56(3), 491-513.
  dat_White_2003_qjep <- read_csv("data/dat_white_2003_causal_judgement_evaluation_evidence.csv") %>% 
    mutate(study = "White (2003b)") %>% 
    mutate(experiment = as.character(experiment)) %>%
    select(where(~ !all(is.na(.)))) %>% 
    rename(
      N_11 = a,
      N_10 = b,
      N_01 = c,
      N_00 = d
    ) %>% 
    mutate(reference = "White, P. A. (2003). Causal judgement as evaluation of evidence: The use of confirmatory and disconfirmatory information. The Quarterly Journal of Experimental Psychology Section A, 56(3), 491-513.") %>% 
    select(
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      rating,
      reference
    ) 
  
  
  
  # merge data
  df_data <- 
  bind_rows(
    dat_Perales_Shanks_2007,
    dat_Buehner_Cheng_Clifford_2003,
    dat_Liljeholm_Cheng_2009,
    dat_White_2003_qjep,
    dat_White_2003_jep,
    dat_White_2004
  ) %>%
    mutate(id = factor(row_number()),
           # maximum likelihood estimates of conditional and unconditional probabilities
           pC = (N_11 + N_10) / (N_11 + N_10 + N_01 + N_00), # base rate of cause C, P(c)
           pEC  = N_11 / (N_11 + N_10), # likelihood of the effect if the cause is present
           pEnoC  = N_01 / (N_01 + N_00)) %>% # likelihood of the effect if the cause is absent 
    mutate(# maximum likelihood estimate of contingency delta P
      deltaP = pEC - pEnoC) %>% 
    mutate(# maximum likelihood estimate of causal power (with different calculation depending on whether the contingency is positive or negative)
      power = ifelse(
        is.nan(deltaP), NaN,
        ifelse(
          deltaP == 0, 0,
          ifelse(
            deltaP > 0,
            deltaP / (1 - pEnoC),
            -deltaP / pEnoC
          )
        )
      ),
    ) %>% 
    mutate(sample_direction = case_when(
      deltaP > 0 ~ "generative",
      deltaP == 0 ~ "zero",
      deltaP < 0 ~ "preventive" )) %>% 
    select(
      id, 
      study,
      experiment,
      n,
      framing,
      N_11,
      N_10,
      N_01,
      N_00,
      pEC,
      pEnoC,
      deltaP,
      power,
      sample_direction,
      rating,
      reference
    ) 
  
  write_csv(df_data, "data/studies_model_comparison.csv")
  
  return(df_data)
}