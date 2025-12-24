# Rp Distribution in a Stochastic SIV Model

This repository provides the R implementation used to compute the probability
mass function of the number of secondary infections (Rp) in a stochastic
Susceptible–Infected–Vaccinated (SIV) model with infection reintroduction and
imperfect vaccination.

The code implements the recursive numerical algorithm described in the
following peer-reviewed publication:

> **Gamboa, M., & Lopez-Herrero, M. J. (2020).**  
> Measuring Infection Transmission in a Stochastic SIV Model with Infection
> Reintroduction and Imperfect Vaccine.  
> *Acta Biotheoretica*, 68(4), 395–420.  
> https://doi.org/10.1007/s10441-019-09373-9

---

## Repository structure

- `compute_Rp_distribution.R`: R script implementing the recursive algorithm
  to compute the probability mass function of the number of secondary infections (Rp).
## Usage

```r
source("R/compute_Rp_distribution.R")

Rp <- compute_Rp_distribution(
  N = 100,
  v0 = 50,
  beta = 10,
  gamma = 1,
  q = 0.03,
  xi = 0.01
)


Rp_v0_i0 <- Rp[[v0+1]][i0, ]
k_values <- 0:(length(Rp_v0_i0) - 1)
Expected_value<-sum(k_values * Rp_v0_i0)

# mass probability function: Rp_v0_i0
# support: k_values
# expected value: Expected_value