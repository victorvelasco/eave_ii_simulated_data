# ==============================================================================
# Treatment model
# ==============================================================================
# Fit a binomial logistic regression to estimate the probability of receiving
# AZ vaccine (vs remaining unvaccinated) as a function of age group and sex.
# Predictors: SexMale and AgeGroup* columns from demographics_coded.
tmt_model <- glm(cbind(AZ, Unvaccinated) ~ . - PB - total,
                 family = "binomial",
                 data   = demographics_coded)

if (!tmt_model$converged) {
  warning("derive_parameters.R: treatment model glm did not converge.")
}

# gammas is a list with one element per time point (K+1 elements for K visits).
# Each element is the coefficient vector from the treatment model at that visit.
# For K=0 (single visit) this is a length-1 list.
gammas0 <- coef(tmt_model)
gammas  <- list(gammas0)

# ==============================================================================
# Outcome model parameters
# ==============================================================================
# betas is a named list, one entry per outcome.
# Each entry is a list with one element per time point (length K+1).
# Each element is a numeric vector: c(intercept, log_OR_treatment).
#   intercept   = qlogis(baseline risk) — set to a 0.5% baseline event rate.
#   log_OR      = log odds ratio for AZ vaccination vs unvaccinated.
#
# Source: Simpson et al
betas <- list(
  Thrombocytopenia = list(c(qlogis(0.005), log(1.42))),
  ITP              = list(c(qlogis(0.005), log(5.77))),
  Venous           = list(c(qlogis(0.005), log(1.03))),
  Arterial         = list(c(qlogis(0.005), log(1.22))),
  Hemorrhagic      = list(c(qlogis(0.005), log(1.48)))
)

# ==============================================================================
# Archive: alternative specifications (not currently used)
# ==============================================================================
# --- Multi-visit treatment models (time-window approach, K > 0) ---
# Fit separate models for each 7-day post-vaccination window.
# Kept for reference; superseded by the single-visit model above.
#
# visits <- c("6d", "13d", "20d", "27d")
# tmt_models <- list()
# treatments_coded[is.na(treatments_coded)] <- 1
#
# tmt_models[[1]] <- glm(cbind(`AZ_0-6`,
#                              `AZ_7-13`+`AZ_14-20`+`AZ_21-27`+`AZ_27+`+Unvaccinated)
#                        ~ . - ..., family = "binomial", data = treatments_coded_subset)
# tmt_models[[2]] <- glm(cbind(`AZ_7-13`,
#                              `AZ_14-20`+`AZ_21-27`+`AZ_27+`+Unvaccinated)
#                        ~ . - ..., family = "binomial", data = treatments_coded_subset)
# tmt_models[[3]] <- glm(cbind(`AZ_14-20`, `AZ_21-27`+`AZ_27+`+Unvaccinated)
#                        ~ . - ..., family = "binomial", data = treatments_coded_subset)
# tmt_models[[4]] <- glm(cbind(`AZ_21-27`, `AZ_27+`+Unvaccinated)
#                        ~ . - ..., family = "binomial", data = treatments_coded_subset)
# gammas <- lapply(tmt_models, coef)

# --- Alternative beta values (higher effect sizes) ---
# betas <- list(
#   c(qlogis(0.005), log(3.43)),
#   c(qlogis(0.005), log(4.60)),
#   c(qlogis(0.005), log(7.81)),
#   c(qlogis(0.005), log(14.07))
# )
