library(tidyverse)
source("script/aux_functions.R")

demographics <- read.csv("data/vaccination_by_sex_and_age.csv", stringsAsFactors = TRUE) |>
  mutate(total = AZ + PB + Unvaccinated)

demographics_coded <- cbind(
  encode_binary(demographics$sex, name = "Sex"),
  encode_ordinal(demographics$age_group, name = "AgeGroup"),
  AZ = demographics$AZ, PB = demographics$PB,
  Unvaccinated = demographics$Unvaccinated, total = demographics$total
) |> as.data.frame()

K <- 0
thetas    <- c(2, rep(1, nlevels(demographics$age_group) - 1))
a_k_values <- c("0", "1")
source("script/derive_parameters.R")   # defines betas, gammas

outcome <- "Venous"; l <- 1; rho <- rep(-0.1, K + 1)
nsim    <- 1000

set.seed(42)
source("script/estimate_cdf.R")
rlist <- lapply(seq_len(nsim), function(isim) {
  source("script/simulate_data.R", local = TRUE)
  source("script/results.R",       local = TRUE)
  as.data.frame(cbind(outcome = outcome, rho = rho[1], isim = isim, param.estimates))
})
new_rows <- do.call(rbind, rlist)
new_rows <- cbind(parameter = rep(c("a_0", "a_1"), nrow(new_rows) / 2), new_rows)

results_old <- readRDS("results.RDS")
results     <- rbind(results_old, new_rows)
for (j in 3:9) results[, j] <- as.numeric(results[, j])
saveRDS(results, "results.RDS")
message("Done. results.RDS now has ", nrow(results), " rows.")
