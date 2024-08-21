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
family <- "poisson"
corstr <- "AR1"
subject_indicator <- rep(1,N)
set.seed(500)
response_indicator <- sort(sample(1:25, 1000, replace=TRUE))
lambda_seq <- c(seq(0,0.5,0.025), seq(0.6,2,0.05))

comb_scheme <- rep(1,25)
theta_1 <- c(0.1, -0.3, -0.6)
intercepts <- list()
betas <- list()
J <- length(unique(response_indicator))
K <- length(unique(subject_indicator))
for(k in 1:K){
  intercepts[[k]] <- c(
    rep(theta_1[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))))
  betas[[k]] <- rbind(
    t(matrix(rep(theta_1[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))), nrow=2)) ) 
}

#####################################################
## Simulate data
#####################################################

data <- dataset.poisson.X2.selective(s=seed, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas)

#####################################################
## Run fusion
#####################################################

time <- system.time(results <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                                         family, corstr, lambda_seq))


