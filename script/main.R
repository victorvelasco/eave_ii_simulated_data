library(tidyverse)
library(parallel)
library(foreach)
library(iterators)
library(doParallel)

here::i_am("script/main.R")

source("script/aux_functions.R")
# Load demographics table. Sent on 8 August 2025
demographics <- read.csv("data/vaccination_by_sex_and_age.csv", stringsAsFactors = TRUE)
# Load treatments table. Extracted on 8 November 2025
treatments <- read.csv("data/events_by_sex_age_type_with_unvaxxed.csv", stringsAsFactors = TRUE)

# Treatments to wide format
treatments <- treatments |>
  pivot_wider(names_from = vacc_status, values_from = event_count) |>
  relocate(
    endpoint_group, age_group, sex,
    `AZ_0-6`, `AZ_7-13`, `AZ_14-20`, `AZ_21-27`, `AZ_27+`,
    `PB_0-6`, `PB_7-13`, `PB_14-20`, `PB_21-27`, `PB_27+`,
    Unvaccinated
  )

# Code covariates with binary variables
treatments_coded <- cbind(
  endpoint_group = treatments$endpoint_group,
  encode_binary(treatments$sex, name = "Sex"),
  encode_ordinal(treatments$age_group, name = "AgeGroup"),
  as.data.frame(treatments)[, -c(1:3)]
)

# Some cells are set to NA because the observed number is between 1 and 5. Replace NAs with 1s
treatments[is.na(treatments)] <- 1

# In demographics table, work out total population per strata
demographics <- demographics |>
  mutate(total = AZ + PB + Unvaccinated)

source("script/print_demographics.R")

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
  2,
  rep(1, nlevels(demographics$age_group)-1)
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
nsim <- 1000

outcomes <- names(betas)
models <- data.frame(
  outcome = rep(outcomes, each = length(rhos)),
  rho     = rep(seq_along(rhos), length(outcomes))
)

##### ---- Sequential for loop (kept for reference / debugging) ----
##### i <- 1
##### results_list <- list()
##### Sys.time()
##### for (outcome in outcomes) {
#####   cat("outcome = ", outcome, "\n")
#####   for (l in 1:length(rhos)) {
#####     cat("rho = ", rhos[l], "\n")
#####     rho <- rep(rhos[l], K+1)
#####     for (isim in 1:nsim) {
#####       cat("isim = ", isim, "\n")
#####       source("script/estimate_cdf.R")
#####       source("script/simulate_data.R")
#####       source("script/results.R")
#####       results_list[[i]] <- cbind(outcome = outcome, rho = rhos[l], isim = isim, param.estimates)
#####       results_list[[i]] <- as.data.frame(results_list[[i]])
#####       i <- i + 1
#####     }
#####     dir.create(paste0("simulations/", outcome, "/", l), recursive = TRUE, showWarnings = FALSE)
#####     write_csv(as.data.frame(b), paste0("simulations/", outcome, "/", l, "/confounders.csv"))
#####     write_csv(as.data.frame(a), paste0("simulations/", outcome, "/", l, "/treatment.csv"))
#####     write_csv(as.data.frame(y), paste0("simulations/", outcome, "/", l, "/outcome.csv"))
#####   }
##### }
##### Sys.time()
##### results <- do.call(rbind, results_list)
##### results <- cbind(parameter = rep(c("a_0", "a_1"), nrow(results)/2), results)
##### rownames(results) <- NULL
##### for (j in 3:9) results[, j] <- as.numeric(results[, j])
##### saveRDS(results, "results.RDS")
##### source("script/create_tables.R")

##### ---- Parallel version ----
library(NoSleepR)

# Use L'Ecuyer-CMRG so each forked worker gets an independent RNG stream
RNGkind("L'Ecuyer-CMRG")
set.seed(1)

Sys.time()
nosleep_on()
results_list <- mclapply(
  seq_len(nrow(models)),
  function(r) {
    tryCatch({
      outcome <- models[r, "outcome"]
      l       <- models[r, "rho"]
      rho     <- rep(rhos[l], K + 1)

      # estimate_cdf depends only on (outcome, rho): compute once per model,
      # not once per simulation replicate
      source("script/estimate_cdf.R", local = TRUE)

      rlist <- list()
      for (isim in seq_len(nsim)) {
        source("script/simulate_data.R", local = TRUE)
        source("script/results.R",       local = TRUE)
        rlist[[isim]] <- as.data.frame(
          cbind(outcome = outcome, rho = rhos[l], isim = isim, param.estimates)
        )
      }
      do.call(rbind, rlist)
    }, error = function(e) {
      message("Error in r=", r,
              " (outcome=", models[r, "outcome"],
              ", rho=", rhos[models[r, "rho"]], "): ", e$message)
      NULL
    })
  },
  mc.cores = 14
)
Sys.time()
nosleep_off()

failed <- vapply(results_list, is.null, logical(1))
if (any(failed)) warning("Failed models: ", paste(which(failed), collapse = ", "))

results <- do.call(rbind, results_list[!failed])
results <- cbind(parameter = rep(c("a_0", "a_1"), nrow(results) / 2), results)
rownames(results) <- NULL
for (j in 3:9) results[, j] <- as.numeric(results[, j])

saveRDS(results, "results.RDS")
source("script/create_tables.R")


