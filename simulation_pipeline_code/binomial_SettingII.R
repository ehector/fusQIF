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
set.seed(123)
response_indicator <- sort(sample(1:5, 500, replace=TRUE))
lambda_seq <- seq(0,1,0.05)

comb_scheme <- 1:10
ints <- c(2,3.5,-2.5,-3.25,-1.75,-1,-0.25,0.5,1.25,-4)
beta1 <- c(1.25,-4,3.5,-3.25,-0.25,2,2.75,0.5,-1,-2.5)
beta2 <- c(-1,-3.25,0.5,2,-4,1.25,3.5,-0.25,-1.75,2.75)
theta <- list()
theta[[1]] <- c(ints[1], beta1[1], beta2[1])
theta[[2]] <- c(ints[2], beta1[2], beta2[2])
theta[[3]] <- c(ints[3], beta1[3], beta2[3])
theta[[4]] <- c(ints[4], beta1[4], beta2[4])
theta[[5]] <- c(ints[5], beta1[5], beta2[5])
theta[[6]] <- c(ints[6], beta1[6], beta2[6])
theta[[7]] <- c(ints[7], beta1[7], beta2[7])
theta[[8]] <- c(ints[8], beta1[8], beta2[8])
theta[[9]] <- c(ints[9], beta1[9], beta2[9])
theta[[10]] <- c(ints[10], beta1[10], beta2[10])
intercepts <- list()
betas <- list()
J <- length(unique(response_indicator))
K <- length(unique(subject_indicator))
c <- 1
for(k in 1:K){
  intercepts[[k]] <- vector("numeric",M)
  betas[[k]] <- matrix(0, 2, M)
  for(j in 1:J){
    intercepts[[k]][which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==c))] <- theta[[c]][1]
    betas[[k]][,which(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==c))] <- theta[[c]][-1]
    c <- c+1
  }
  betas[[k]] <- t(betas[[k]]) 
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

save.image("binomial_II_result.RData")
