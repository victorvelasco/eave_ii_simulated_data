library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

setwd("/Users/vvelasco/bts3/eave_ii_simulated_data/")

source("script/aux_functions.R")
# Load demographics table
# demographics <- read.csv("data/vaccination_by_sex_and_age.csv", stringsAsFactors = TRUE)
demographics <- read.csv("data/vaccination_by_sex_and_age.csv", stringsAsFactors = TRUE)
treatments <- read.csv("data/events_by_sex_age_type_with_unvaxxed.csv", stringsAsFactors = TRUE)


### demographics <- demographics |>
###   filter(endpoint_group == "any_haem") |> select(-endpoint_group)

treatments <- treatments |>
  pivot_wider(names_from = vacc_status, values_from = event_count) |>
  relocate(
    endpoint_group, age_group, sex,
    `AZ_0-6`, `AZ_7-13`, `AZ_14-20`, `AZ_21-27`, `AZ_27+`,
    `PB_0-6`, `PB_7-13`, `PB_14-20`, `PB_21-27`, `PB_27+`,
    Unvaccinated
  )
treatments_coded <- cbind(
  endpoint_group = treatments$endpoint_group,
  encode_binary(treatments$sex, name = "Sex"),
  encode_ordinal(treatments$age_group, name = "AgeGroup"),
  as.data.frame(treatments)[, -c(1:3)]
)

treatments[is.na(treatments)] <- 1
demographics <- demographics |>
  mutate(total = AZ + PB + Unvaccinated)

demographics_coded <- cbind(
  encode_binary(demographics$sex, name = "Sex"),
  encode_ordinal(demographics$age_group, name = "AgeGroup"),
  AZ = demographics$AZ,
  PB = demographics$PB,
  Unvaccinated = demographics$Unvaccinated,
  total = demographics$total
)
demographics_coded <- as.data.frame(demographics_coded)

K <- 0 # Number of visits minus one

# Parameter determining the strength of confounding (five scenarios)
rhos <- c(-0.1, -0.3, -0.5, -0.7, -0.9)

# Parameters for causal quantity of interest
# Parameters for allocation of treatment to individuals
source("script/derive_parameters.R")

# Parameters for confounding mechanism
thetas <- c(
  1,
  rep(0.1, nlevels(demographics$age_group)-1)
)

# Possible values of treatment allocation vector
##### a_k_values <- c(
#####    "0", "1",
#####    "00", "01", "10", "11",
#####    "000", "001", "010", "011", "100", "101", "110", "111" ,
#####    "0000", "0001", "0010", "0011", "0100", "0101", "0110", "0111" ,
#####    "1000", "1001", "1010", "1011", "1100", "1101", "1110", "1111"
##### )
##### a_k_values <- c(
#####    "0", "1",
#####    "00", "01", "10",
#####    "000", "001", "010", "100",
#####    "0000", "0001", "0010", "0100",
#####    "1000"
##### )

a_k_values <- c("0", "1")
nsim <- 100

set.seed(1)
outcomes <- names(betas)
models <- data.frame(outcome = rep(outcomes, each = length(rhos)), rho = rep(1:length(rhos), length(outcomes)))

Sys.time()
results_list <- mclapply(
  1:nrow(models), 
  function(r) {
    outcome <- models[r, "outcome"]
    l <- models[r, "rho"]
    rho <- rep(rhos[l], K+1)
    rlist <- list()
    for (isim in 1:nsim) {
      cat("isim = ", isim, "\n")
      source("script/estimate_cdf.R")
      source("script/simulate_data.R")
      source("script/results.R")
      rlist[[isim]] <- cbind(outcome = outcome, rho = rhos[l], isim = isim, param.estimates)
      rlist[[isim]] <- as.data.frame(rlist[[isim]])
    }
    
    return(do.call(rbind, rlist))
  },
  mc.cores = 12
)
results <- do.call(rbind, results_list)
results <- cbind(parameter = rep(c("a_0", "a_1"), nrow(results)/2), results)
rownames(results)<-NULL
results[, 3] <- as.numeric(results[, 3])
results[, 4] <- as.numeric(results[, 4])
results[, 5] <- as.numeric(results[, 5])
results[, 6] <- as.numeric(results[, 6])
results[, 7] <- as.numeric(results[, 7])
results[, 8] <- as.numeric(results[, 8])
results[, 9] <- as.numeric(results[, 9])
Sys.time()
saveRDS(results, "results.RDS")
source("script/create_tables.R")

outcomes <- names(betas)
for (outcome in outcomes) {
  cat("outcome = ", outcome, "\n")
  for (l in 1:length(rhos)) {
    cat("rho = ", rhos[l], "\n")
    rho <- rep(rhos[l], K+1)
    for (isim in 1:nsim) {
      cat("isim = ", isim, "\n")
      source("script/estimate_cdf.R")
      source("script/simulate_data.R")
      source("script/results.R")
      results_list[[i]] <- cbind(outcome = outcome, rho = rhos[l], isim = isim, param.estimates)
      results_list[[i]] <- as.data.frame(results_list[[i]])
      i <- i + 1
    }

    write_csv(as.data.frame(b), paste0("simulations/", outcome, "/", l, "/confounders.csv"))
    write_csv(as.data.frame(a), paste0("simulations/", outcome, "/", l, "/treatment.csv"))
    write_csv(as.data.frame(y), paste0("simulations/", outcome, "/", l, "/outcome.csv"))
  }
}
results <- do.call(rbind, results_list)
results <- cbind(parameter = rep(c("a_0", "a_1"), nrow(results)/2), results)
rownames(results)<-NULL
results[, 3] <- as.numeric(results[, 3])
results[, 4] <- as.numeric(results[, 4])
results[, 5] <- as.numeric(results[, 5])
results[, 6] <- as.numeric(results[, 6])
results[, 7] <- as.numeric(results[, 7])
results[, 8] <- as.numeric(results[, 8])
results[, 9] <- as.numeric(results[, 9])
saveRDS(results, "results.RDS")
source("script/create_tables.R")

