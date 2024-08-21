\name{fusQIF-package}
\alias{fusQIF-package}
\alias{fusQIF}
\docType{package}
\title{\packageTitle{fusQIF}}
\description{\packageDescription{fusQIF}}
\details{

  THIS SOURCE CODE IS SUPPLIED "AS IS"" WITHOUT WARRANTY OF ANY KIND, AND ITS AUTHOR DISCLAIMs ANY AND ALL WARRANTIES, INCLUDING BUT NOT LIMITED TO ANY IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, AND ANY WARRANTIES OR NON INFRINGEMENT. THE USER ASSUMES ALL LIABILITY AND RESPONSIBILITY FOR USE OF THIS SOURCE CODE, AND THE AUTHOR WILL NOT BE LIABLE FOR DAMAGES OF ANY KIND RESULTING FROM ITS USE. Without limiting the generality of the foregoing, the author does not warrant that the Source Code will be error-free, will operate without interruption, or will meet the needs of the user.
  
  The DESCRIPTION file: \packageDESCRIPTION{fusQIF}
  \packageIndices{fusQIF}

}
\author{
  \packageAuthor{fusQIF}
  
  Maintainer: \packageMaintainer{fusQIF}
}
\references{

  QIF sub-routine implementation is courtesy of Dr. Lan Luo.
  
  The posdef.matrix function was written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708.
  
  Emily C. Hector and Peter X.-K. Song (2020). Doubly distributed supervised learning and inference with high-dimensional correlated outcomes. Journal of Machine Learning Research, 21:1–35.
  
  Emily C. Hector and Peter X.-K. Song (2020). A distributed and integrated method of moments for high-dimensional correlated data analysis. Journal of the American Statistical Association, 116:805-818.
  
}

\keyword{package}
\seealso{
  
}
\examples{

N <- 10000
M <- 500
subject_indicator <- c(rep(1,N/2),rep(2,N/2))
set.seed(500)
response_indicator <- c(rep(1,45),rep(2,56),rep(3,46),rep(4,42),rep(5,55),rep(6,50),rep(7,47),rep(8,50),rep(9,59),rep(10,50))

comb_scheme <- c(c(rep(1,4),rep(2,3),rep(3,3)), c(rep(1,4),rep(2,3),rep(3,3)))
theta_1 <- c(-0.4, 0.1, -0.2)
theta_2 <- c(0.1, -0.3, -0.6)
theta_3 <- c(-0.8, 0.2, 0.4)
intercepts <- list()
betas <- list()
J <- length(unique(response_indicator))
K <- length(unique(subject_indicator))
for(k in 1:K){
  intercepts[[k]] <- c(
    rep(theta_1[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))), 
    rep(theta_2[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))), 
    rep(theta_3[1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))))
  betas[[k]] <- rbind(
    t(matrix(rep(theta_1[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==1))), nrow=2)),
    t(matrix(rep(theta_2[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==2))), nrow=2)),
    t(matrix(rep(theta_3[-1], sum(response_indicator %in% which(comb_scheme[((k-1)*J+1):(k*J)]==3))), nrow=2)) ) 
}

#####################################
## Gaussian Example
#####################################

# Generate a dataset with Gaussian outcome
data <- dataset.normal.X2.selective(s=1, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas)
# Fusion analysis
results.gaussian <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                              family="gaussian", corstr="AR1", lambda_seq=seq(0,2.5,0.05), rho=1)
                         
#####################################
## Binomial Example
#####################################

# Generate a dataset with Bernoulli outcome
data <- dataset.binomial.X2.selective(s=1, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas)
# Fusion analysis
results.binomial <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                              family="binomial", corstr="AR1", lambda_seq=seq(0,2.5,0.05), rho=1)
                    
#####################################
## Poisson Example
#####################################

# Generate a dataset with Poisson outcome
data <- dataset.poisson.X2.selective(s=1, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas)
# Fusion analysis
results.poisson <- qiffusion(response~X1+X2, data, id="id", response_indicator, subject_indicator,
                              family="poisson", corstr="AR1", lambda_seq=c(seq(0,0.5,0.025), seq(0.6,2,0.05)), rho=1)

}
