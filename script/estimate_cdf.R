
m <- 1e6 # Sample size for estimating the CDF

fhat <- list() # List to store empirical CDFs

for (tmt_allocation in a_k_values) {
  # Step 1. For j=1,...,m, simulate B_j ~ p(B_j). Set k = 0.
  b <- demographics[sample(1:nrow(demographics), size = m, replace = TRUE, prob = demographics$n), -ncol(demographics)]
  
  b <- b[, c("sex", "age_group")]
  rownames(b) <- NULL
  # Categorical and ordinal encoding of demographic covariates, B
  b_coded <- cbind(
    encode_binary(b$sex, name = "Sex"),
    encode_ordinal(b$age_group, name = "AgeGroup")
  )
  for (k in 0:nchar(tmt_allocation)) {
    cat("Tmt allocation = ", tmt_allocation, " k = ", k, "\n")
    # (Skip) Step 2. Sample time-varying covariates
    # Step 3. Compute H_k, \hat{F}_H.
    h <- as.vector(b_coded %*% thetas)
    if (k == nchar(tmt_allocation)) {
      fhat[[tmt_allocation]] <- ecdf(h)
    }
    
    # Step 4. If k = K: Stop
    if (k == nchar(tmt_allocation)) break
    
    # Step 5. Compute R_k, ranking of H_k
    r <- rank(h)
    
    # Step 6. Computer U_H
    u_h <- (r - runif(m))/m
    
    # Step 7. Compute Z_H
    z_h <- qnorm(u_h)
    
    # Step 8. Sample Z_Y, calculate U_Y
    z_y <- rnorm(m, rho[k+1]*z_h, sqrt(1-rho[k+1]^2))
    u_y <- pnorm(z_y)
    
    # Step 9. Compute Y_{k+1}
    #### failed <- u_y < expit(sum(as.vector(c(1, as.numeric(strsplit(substr(tmt_allocation, 1, k+1), "")[[1]]))) * betas[[outcome]][[k+1]]))
    failed <- u_y < expit(
      sum(as.vector(c(1, sum(as.numeric(strsplit(substr(tmt_allocation, 1, k+1), "")[[1]])))) * betas[[outcome]][[k+1]])
    )
    
    # Step 10. Replace failed individuals with individuals that have survived
    random_replacements <- sample(which(!failed), size = sum(failed), replace = TRUE)
    b[which(failed), ] <- b[random_replacements, ]
    b_coded[which(failed), ] <- b_coded[random_replacements, ]
    
    # Step 11. Set k = k + 1 and return to step 2
  }
}
