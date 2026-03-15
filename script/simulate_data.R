# Initialise matrices to store variables
N <- 1e6
a <- matrix(NA, nrow = N, ncol = K+1); colnames(a) <- paste0("A", 0:K)
y <- cbind(
  matrix(1,  nrow = N, ncol = 1),
  matrix(NA, nrow = N, ncol = K+1)
); colnames(y) <- paste0("Y", 0:(K+1))

# Step 1. Sample B ~ p(B). Set k = 0
b <- demographics[sample(1:nrow(demographics), size = N, replace = TRUE, prob = demographics$n), -ncol(demographics)]
b <- b[, c("sex", "age_group")]
# Variables age, SCSIMD5 and NumComorbidities are coded as ordinal predictors
b_coded <- cbind(
  encode_binary(b$sex, name = "Sex"),
  encode_ordinal(b$age_group, name = "AgeGroup")
)

# In the first iteration, all subjects are alive (Y0 = 1 for every subject)
failed <- rep(FALSE, N)
survived <- rep(TRUE, N)
failures <- rep(K+1, N)
Nsurvived <- N
# For every visit k=0,...,K, we iterate Steps 2-8 of the algorithm
for (k in 0:K) {
  # (Skip) Step 2. Sample time-varying covariates
  # Step 3. Sample A_k ~ p(A_k | B, \bar{A}_{k-1}, Y_k = 1)
  if (k == 0) {
    #### eligible <- survived
    #### Neligible <- sum(eligible)
    linpred.a <- as.vector(cbind(rep(1, Nsurvived), b_coded[which(survived), ]) %*% gammas[[k+1]])
  } else {
    #### eligible <- survived & rowSums(a[, 1:(k+1)], na.rm = TRUE) == 0
    #### Neligible <- sum(eligible)
    #### linpred.a <- as.vector(cbind(rep(1, Neligible), b_coded[which(eligible), ]) %*% gammas[[k+1]])
    linpred.a <- as.vector(cbind(rep(1, Neligible), b_coded[which(survived), ]) %*% gammas[[k+1]])
  }
  #### a[eligible, k+1] <- as.integer(runif(Neligible) < expit(linpred.a))
  a[survived, k+1] <- as.integer(runif(Nsurvived) < expit(linpred.a))
  #### if (any(survived & !eligible)) {
  ####   a[survived & !eligible, k+1] <- 0
  #### }
  a[!survived, k+1] <- NA

  # Step 4. Calculate H_k
  h <- as.vector(b_coded[survived, ] %*% thetas)

  # Step 5. Calculate U_H, Z_H. Use risk score as in Appendix H
  if (k == 0) {
    a_bar <- as.character(a[,1]) # This is the history of allocated treatments, to date
  } else {
    a_bar <- apply(a[survived,1:(k+1)], 1, function(x) paste(x, collapse = ""))
  }
  u_h <- rep(0, Nsurvived)
  u_h_min <- rep(0, Nsurvived)
  u_h_max <- rep(0, Nsurvived)
  for (aa in unique(a_bar)) {
    u_h_min[a_bar == aa] <- fhat[[aa]](h[a_bar == aa] - 1)
    u_h_max[a_bar == aa] <- fhat[[aa]](h[a_bar == aa])
  }
  u_h <- runif(Nsurvived, u_h_min, u_h_max)
  z_h <- qnorm(u_h)

  # Step 6. Sample Z_Y  and calculate U_Y
  z_y <- rnorm(Nsurvived, rho[k+1]*z_h, sqrt(1-rho[k+1]^2))
  u_y <- pnorm(z_y)

  # Step 7. Compute Y_{k+1}
  if (k == 0) {
    linpred.y <- as.vector(cbind(matrix(1, nrow = Nsurvived, ncol = 1), a[, 1]) %*% betas[[outcome]][[1]])
    #### linpred.y <- as.vector(cbind(matrix(1, nrow = Nsurvived, ncol = 1), a[, 1]) %*% betas[[1]])
  } else {
    linpred.y <- as.vector(cbind(matrix(1, nrow = Nsurvived, ncol = 1), a[survived, 1:(k+1)]) %*% betas[[outcome]][[k+1]])
    #### linpred.y <- as.vector(cbind(matrix(1, nrow = Nsurvived, ncol = 1), rowSums(a[survived, 1:(k+1)])) %*% betas[[k+1]])
  }
  failed[failed == FALSE] <- u_y < expit(linpred.y)
  y[, k+2] <- ifelse(failed, 0, 1)
  survived <- !failed
  Nsurvived <- sum(survived)

  # Step 8. Set k = k + 1 and return to Step 2
}
