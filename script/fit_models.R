failure_times <- rowSums(y)
# First compute weights
w <- matrix(1, N, K+1)
N <- nrow(y)
K <- ncol(a)
i <- 1

for (k in 1:(K+1)) {
  if (k == 1) {
    fit.numer <- glm(a[, 1] ~ 1, family = binomial)
    fit.denom <- glm(a[, 1] ~ b_coded, family = binomial)
  } else {
    fit.numer <- glm(a[, k] ~ a[,(1:(k-1))], family = binomial, subset = which(failure_times>=k))
    fit.denom <- glm(a[, k] ~ a[,(1:(k-1))] + b_coded, family = binomial, subset = which(failure_times>=k))
  }
  
  phat.denom <- a[failure_times>=k, k] * fitted(fit.denom) + (1-a[failure_times>=k, k]) * (1- fitted(fit.denom))
  phat.numer <- a[failure_times>=k, k] * fitted(fit.numer) + (1-a[failure_times>=k, k]) * (1- fitted(fit.numer))
  if (k == 1) {
    w[, k] <- phat.numer / phat.denom
  } else {
    w[failure_times>=k, k] <- w[failure_times>=k,k-1] * phat.numer / phat.denom
  }
}
