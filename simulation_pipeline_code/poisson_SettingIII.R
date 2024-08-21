#####################################################
## Load packages, set seed
#####################################################

library(fusQIF)
seed <- 1

#####################################################
## Set up simulation parameters
#####################################################

N <- 700
M <- 15103
family <- "poisson"
corstr <- "AR1"
subject_indicator <- c(rep(1,550),rep(2,150))
set.seed(500)
response_indicator <- c(rep(1,149), rep(2,3035), rep(3,887), rep(4,2105), rep(5,372), 
                        rep(6,1130), rep(7,95), rep(8,388), rep(9,263), rep(10,230), 
                        rep(11,3605), rep(12,1187), rep(13,1149), rep(14,115), rep(15,393))
lambda_seq <- seq(0,0.9,0.06)

comb_scheme <- c(1,1,2,2,2,2,1,1,3,1,1,1,2,1,1,1,1,1,1,2,4,2,1,1,1,1,1,2,1,1)
theta_1 <- c(-4, 1, -2)/5
theta_2 <- c(4, -1, 2)/5
theta_3 <- c(0.8, 0.2, 0.6)
theta_4 <- c(1,-2,3)/5
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
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))] <- theta_1[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))] <- theta_2[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))] <- theta_3[-1]
  betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==4))] <- theta_4[-1]
  betas[[k]] <- t(betas[[k]])
}

#####################################################
## Simulate data
#####################################################

data <- dataset.poisson.X2.selective.large(s=seed, N, M, M_2=2^14, response_indicator, subject_indicator, comb_scheme, intercepts, betas)

#####################################################
## Run fusion
#####################################################

time <- system.time(results <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                                         family, corstr, lambda_seq))


