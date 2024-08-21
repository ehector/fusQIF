#####################################################
## Load packages, set seed
#####################################################

library(fusQIF)
seed <- 1

#####################################################
## Set up simulation parameters
#####################################################

N <- 5000
M <- 500
family <- "binomial"
corstr <- "AR1"
subject_indicator <- c(rep(1,N/2),rep(2,N/2))
set.seed(500)
response_indicator <- sort(sample(1:10, 500, replace=TRUE))
lambda_seq <- seq(0,2.5,0.05)

comb_scheme <- c(c(1,1,2,3,3,4,4,4,5,5), c(1,2,2,3,3,4,4,4,5,5))
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
  intercepts[[k]] <- c(
    rep(theta_1[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))), 
    rep(theta_2[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))), 
    rep(theta_3[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))),
    rep(theta_4[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==4))),
    rep(theta_5[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==5))))
  betas[[k]] <- rbind(
    t(matrix(rep(theta_1[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))), nrow=2)),
    t(matrix(rep(theta_2[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))), nrow=2)),
    t(matrix(rep(theta_3[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))), nrow=2)),
    t(matrix(rep(theta_4[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==4))), nrow=2)),
    t(matrix(rep(theta_5[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==5))), nrow=2)) ) 
}

#####################################################
## Simulate data
#####################################################

data <- dataset.binomial.X2.selective(s=seed, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas)

#####################################################
## Run fusion
#####################################################

time <- system.time(results <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                                         family, corstr, lambda_seq))

