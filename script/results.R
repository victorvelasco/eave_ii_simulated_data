### library(tidyverse)
### library(parallel)
### library(foreach)
### library(iterators)
### library(doParallel)
### 
### setwd("/Users/vvelasco/bts3/eave_ii_simulated_data/")
### K <- 0 # Number of visits minus one
### 
### # Parameter determining the strength of confounding (five scenarios)
### rhos <- c(-0.1, -0.3, -0.5, -0.7, -0.9)
### 
### source("script/aux_functions.R")
### # Parameters for causal quantity of interest
### source("script/derive_parameters.R")

### for (outcome in names(betas)) {
###   for (l in 1:length(rhos)) {
model.outcome.unweighted <- list()

### b <- read_csv(paste0("simulations/", outcome, "/", l, "/confounders.csv")) |>
###   mutate_all(factor) |>
###   as.data.frame()
### b_coded <- cbind(
###   encode_binary(b$sex, name = "Sex"),
###   encode_ordinal(b$age_group, name = "AgeGroup")
### ) 
### a <- read_csv(paste0("simulations/", outcome, "/", l, "/treatment.csv")) |> as.data.frame() |> as.matrix()
### y <- read_csv(paste0("simulations/", outcome, "/", l, "/outcome.csv")) |> as.data.frame() |> as.matrix()

failure_times <- rowSums(y)
for (k in 1:(K+1)) {
  if (k == 1) {
    aSum <- matrix(a[, 1:k] , ncol = 1)
  } else {
    aSum <- rowSums(a[, 1:k])
  }
  model.outcome.unweighted[[k]] <- glm(
    as.factor(y[,k+1]==0) ~ aSum, 
    family = binomial, 
    subset = which(failure_times>=k)
  ) 
}

N <- nrow(b)
w <- matrix(1, N, K+1)
model.outcome.weighted <- list()
for (k in 1:(K+1)) {
  if (k == 1) {
    fit.numer <- glm(a[, 1] ~ 1, family = binomial)
    fit.denom <- glm(a[, 1] ~ b_coded, family = binomial)
  } else {
    if (k == 2) {
      aSum <- matrix(a[, 1:(k-1)] , ncol = 1)
    } else {
      aSum <- rowSums(a[, 1:(k-1)])
    }
    fit.numer <- glm(a[, k] ~ aSum, family = binomial, subset = which(failure_times>=k))
    fit.denom <- glm(a[, k] ~ aSum + b_coded, family = binomial, subset = which(failure_times>=k))
  }
  
  phat.denom <- a[failure_times>=k, k] * fitted(fit.denom) + (1-a[failure_times>=k, k]) * (1- fitted(fit.denom))
  phat.numer <- a[failure_times>=k, k] * fitted(fit.numer) + (1-a[failure_times>=k, k]) * (1- fitted(fit.numer))
  if (k == 1) {
    w[, k] <- phat.numer / phat.denom
  } else {
    w[failure_times>=k, k] <- w[failure_times>=k,k-1] * phat.numer / phat.denom
  }
  if (k == 1) {
    aSum <- matrix(a[, 1], ncol = 1)
  } else {
    aSum <- rowSums(a[, 1:k])
  }
  model.outcome.weighted[[k]] <- glm(
    as.factor(y[,k+1]==0) ~ aSum, 
    family = binomial, 
    weights = w[, k], 
    subset = which(failure_times>=k)
  )
}

param.estimates <- do.call(rbind, lapply(1:(K+1), 
  function(i) cbind(
    betas[[outcome]][[i]],
    #### betas[[i]],
    coefficients(summary(model.outcome.unweighted[[i]]))[,1:2],
    coefficients(summary(model.outcome.weighted[[i]]))[,1:2]
  ) |> round(4)
))

colnames(param.estimates)[1] <- "True value"
colnames(param.estimates)[2] <- "Estimate (unweighted)"
colnames(param.estimates)[3] <- "SD (unweighted)"
colnames(param.estimates)[4] <- "Estimate (weighted)"
colnames(param.estimates)[5] <- "SD (weighted)"
write_csv(as.data.frame(param.estimates), paste0("output/", outcome, "/", l, "/results.csv"))
###   }
### }
