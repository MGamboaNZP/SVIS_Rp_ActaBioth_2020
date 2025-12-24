############################################################
# Probability mass function of Rp in a stochastic SIV model
#
# Based on:
# Gamboa, M., & Lopez-Herrero, M. J. (2020).
# Measuring Infection Transmission in a Stochastic SIV Model
# with Infection Reintroduction and Imperfect Vaccine.
# Acta Biotheoretica, 68(4), 395–420.
# DOI: 10.1007/s10441-019-09373-9
#
# This function computes the probability mass function
# of the number of secondary infections Rp, conditional
# on the number of vaccinated individuals.
############################################################

compute_Rp_distribution <- function(
  N,
  v0,
  beta,
  gamma,
  q,
  xi
) {
  
  ##########################################################
  # Auxiliary rate function
  ##########################################################
  
  fq_vi <- function(N, i, v, beta, xi, gamma, q) {
    (N - i - v) * (beta * i / N + xi) +
      (beta * q * i * v / N) +
      (q * v * xi) +
      (gamma * i)
  }
  
  ##########################################################
  # Storage for results
  ##########################################################
  
  Rp_distributions <- vector("list", v0 + 1)
  
  ##########################################################
  # Loop over number of vaccinated individuals
  ##########################################################
  
  for (v in 0:v0) {
    
    ########################################################
    # Build diagonal matrices
    ########################################################
    
    DQ <- diag(sapply(1:(N - v), function(i)
      fq_vi(N, i, v, beta, xi, gamma, q)))
    
    D_tilde <- diag(sapply(1:(N - v), function(i)
      (beta * i / N + xi) * (N - v - i)))
    
    D_v <- diag(sapply(1:(N - v), function(i)
      (beta * i / N + xi) * q * v))
    
    ########################################################
    # Initialization: k = 0
    ########################################################
    
    z_vk <- matrix(0, nrow = N - v, ncol = N)
    
    for (i in 1:(N - v)) {
      z_vk[i, 1] <- (i * gamma) /
        fq_vi(N, i, v, beta, xi, gamma, q)
    }
    
    ########################################################
    # Recursion over k
    ########################################################
    
    for (k in 1:(N - 1)) {
      
      # Shifted vector (tilde)
      z_tilde_prev <- c(z_vk[2:(N - v), k], 0)
      
      # Shifted vector for v - 1 (hat)
      if (v > 0) {
        z_hat_prev <- zv_prev[2:(N - v + 1), k]
      } else {
        z_hat_prev <- rep(0, N - v)
      }
      
      rhs <- D_tilde %*% z_tilde_prev +
        D_v %*% z_hat_prev
      
      z_vk[, k + 1] <- solve(DQ, rhs)
    }
    
    ########################################################
    # Store and update
    ########################################################
    
    Rp_distributions[[v + 1]] <- z_vk
    zv_prev <- z_vk
  }
  
  ##########################################################
  # Output
  ##########################################################
  
  return(Rp_distributions)
}


