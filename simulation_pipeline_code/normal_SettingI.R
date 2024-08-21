#####################################################
## Load packages, set seed
#####################################################

library(fusQIF)
seed <- 1

#####################################################
## Set up simulation parameters
#####################################################

N <- 10000
M <- 1000
family <- "gaussian"
corstr <- "AR1"
subject_indicator <- c(rep(1,N/2),rep(2,N/2))
set.seed(500)
response_indicator <- sort(sample(1:20, 1000, replace=TRUE))
lambda_seq <- seq(0,2.5,0.05)

comb_scheme <- c(1,1,2,3,3,4,4,4,5,5,1,1,2,3,3,4,4,4,5,5,1,2,2,3,3,4,4,4,5,5,1,2,2,3,3,4,4,4,5,5)
theta_1 <- c(-4, 1, -2)
theta_2 <- c(4, -1, 2)
theta_3 <- c(0.8, 0.2, 0.6)
theta_4 <- c(1,-2,3)
theta_5 <- c(-1,2,-3)
intercepts <- list()
betas <- list()
J <- length(unique(response_indicator))
K <- length(unique(subject_indicator))
for(k in 1:K){
  intercepts[[k]] <- vector("numeric",M)
  betas[[k]] <- matrix(0, 2, M)
  intercepts[[k]][which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))] <- theta_1[1]
  intercepts[[k]][which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))] <- theta_2[1]
  intercepts[[k]][which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))] <- theta_3[1]
  intercepts[[k]][which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==4))] <- theta_4[1]
  intercepts[[k]][which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==5))] <- theta_5[1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))] <- theta_1[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))] <- theta_2[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))] <- theta_3[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==4))] <- theta_4[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==5))] <- theta_5[-1]
  betas[[k]] <- t(betas[[k]])
}

#####################################################
## Simulate data
#####################################################

data <- dataset.normal.X2.selective(s=seed, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas)

#####################################################
## Run fusion
#####################################################

time <- system.time(results <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                                         family, corstr, lambda_seq))




