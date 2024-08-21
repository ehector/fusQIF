posdef.matrix <- function(n, ev = runif(n, 0, 10), seed=500){
  ## Function written by Ravi Varadhan: https://stat.ethz.ch/pipermail/r-help/2008-February/153708
  ## Generating a random positive-definite matrix with user-specified positive eigenvalues
  ## If eigenvalues are not specified, they are generated from a uniform distribution
  set.seed(seed)
  Z <- matrix(ncol=n, rnorm(n^2))
  decomp <- qr(Z)
  Q <- qr.Q(decomp) 
  R <- qr.R(decomp)
  d <- diag(R)
  ph <- d / abs(d)
  O <- Q %*% diag(ph)
  Z <- t(O) %*% diag(ev) %*% O
  return(Z)
}

dataset.normal.X2.selective <- function(s, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas){
  family <- "gaussian"
  set.seed(s)
  sd <- 4.0
  r <- 0.8
  A <- sd^2 * r^abs(outer(1:max(table(response_indicator)), 1:max(table(response_indicator)) , "-"))
  
  S <- posdef.matrix(length(unique(response_indicator)), seed=s)
  
  if (length(unique(table(response_indicator)))==1){
    Sigma <- kronecker(S,A)
  } else{
    Sigma <- kronecker(S,A)
    for(i in 1:length(unique(response_indicator))){
      if(table(response_indicator)[i] < max(table(response_indicator))){
        start <- sum(table(response_indicator)[1 : i])+1
        end <- start + max(table(response_indicator))-table(response_indicator)[i]-1
        Sigma <- Sigma[-c(start:end), -c(start:end)]
      }
    }
  }
  
  ## Generate N replicates of an M-dimensional Normal vector with mean beta_0+beta_1*x
  epsilon <- do.call(cbind, lapply(unique(subject_indicator), 
                                   function(x) MASS::mvrnorm(N/length(unique(subject_indicator)), rep(0,M), Sigma*runif(1,0,1), tol=1e-6) ))
  X_1 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s))
  X_2 <- MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s*2000))
  new_eps <- do.call(rbind, lapply(unique(subject_indicator), function(x){ epsilon[,((x-1)*M+1):(((x-1)*M+M))] }))
  Y <- matrix(0, M, N)
  for(k in 1:length(unique(subject_indicator))){
    Y[,min(which(subject_indicator==k)):max(which(subject_indicator==k))] <- 
      t(new_eps[min(which(subject_indicator==k)):max(which(subject_indicator==k)),]) + 
      (matrix(rep(intercepts[[k]],length(min(which(subject_indicator==k)):max(which(subject_indicator==k)))), M)) +
      apply(X_1[min(which(subject_indicator==k)):max(which(subject_indicator==k)),], 1, function(x) betas[[k]][,1]*x) + 
      apply(X_2[min(which(subject_indicator==k)):max(which(subject_indicator==k)),], 1, function(x) betas[[k]][,2]*x)
  }
  
  data <- data.frame(id=sort(rep(1:N,M)), c(t(X_1)), c(t(X_2)), c(Y))
  colnames(data) <- c("id","X1","X2","response")
  return(data)
}

dataset.binomial.X2.selective <- function(s, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas){
  family <- "binomial"
  set.seed(s)
  X_1 <- t(MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s)))
  X_2 <- t(MASS::mvrnorm(N, rep(0,M), posdef.matrix(M, seed=s*2000)))
  if(length(unique(subject_indicator))==1){
    r <- 0.8
    A <- r^abs(outer(1:max(table(response_indicator)), 1:max(table(response_indicator)) , "-"))
    S <- runif(1,0,1)^abs(outer(1:length(unique(response_indicator)), 1:length(unique(response_indicator)) , "-"))
    
    if (length(unique(table(response_indicator)))==1){
      Sigma <- kronecker(S,A)
    } else{
      Sigma <- kronecker(S,A)
      for(i in 1:length(unique(response_indicator))){
        if(table(response_indicator)[i] < max(table(response_indicator))){
          start <- sum(table(response_indicator)[1 : i])+1
          end <- start + max(table(response_indicator))-table(response_indicator)[i]-1
          Sigma <- Sigma[-c(start:end), -c(start:end)]
        }
      }
    }
    Y <- matrix(0, N, M)
    for(n in 1:N){
      X_1n <- X_1[,n]
      X_2n <- X_2[,n]
      Y[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=intercepts, betas=betas, xformula= ~X_1n + X_2n, link="logit", cor.matrix=Sigma)$simdata$y
    }
    Y[((k-1)*M*n_k+1):(k*M*n_k)] <- Yk
  } else {
    Y <- vector("numeric",N*M)
    for(k in 1:(length(unique(subject_indicator)))){
      r <- runif(1,-1,1)
      A <- r^abs(outer(1:max(table(response_indicator)), 1:max(table(response_indicator)) , "-"))
      S <- runif(1,0,1)^abs(outer(1:length(unique(response_indicator)), 1:length(unique(response_indicator)) , "-"))
      
      if (length(unique(table(response_indicator)))==1){
        Sigma <- kronecker(S,A)
      } else{
        Sigma <- kronecker(S,A)
        for(i in 1:length(unique(response_indicator))){
          if(table(response_indicator)[i] < max(table(response_indicator))){
            start <- sum(table(response_indicator)[1 : i])+1
            end <- start + max(table(response_indicator))-table(response_indicator)[i]-1
            Sigma <- Sigma[-c(start:end), -c(start:end)]
          }
        }
      }
      n_k <- sum(subject_indicator==k)
      
      Yk <- matrix(0, n_k, M)
      for(n in 1:n_k){
        X_1n <- X_1[,which(subject_indicator==k)[n]]
        X_2n <- X_2[,which(subject_indicator==k)[n]]
        Yk[n,] <- SimCorMultRes::rbin(clsize=M, intercepts=intercepts[[k]], betas=betas[[k]], xformula= ~X_1n + X_2n, link="logit", cor.matrix=Sigma)$simdata$y
      }
      Y[((k-1)*M*n_k+1):(k*M*n_k)] <- c(t(Yk))
    } 
  }
  data <- data.frame(id=sort(rep(1:N,M)), c(X_1), c(X_2), Y)
  colnames(data) <- c("id","X1","X2","response")
  return(data)
}

dataset.poisson.X2.selective <- function(s, N, M, response_indicator, subject_indicator, comb_scheme, intercepts, betas){
  family <- "poisson"
  set.seed(s)
  sd <- 0.4
  r <- 0.8
  A <- sd^2 * r^abs(outer(1:max(table(response_indicator)), 1:max(table(response_indicator)) , "-"))
  
  S <- posdef.matrix(length(unique(response_indicator)), seed=s)
  
  if (length(unique(table(response_indicator)))==1){
    Sigma <- kronecker(S,A)
  } else{
    Sigma <- kronecker(S,A)
    for(i in 1:length(unique(response_indicator))){
      if(table(response_indicator)[i] < max(table(response_indicator))){
        start <- sum(table(response_indicator)[1 : i])+1
        end <- start + max(table(response_indicator))-table(response_indicator)[i]-1
        Sigma <- Sigma[-c(start:end), -c(start:end)]
      }
    }
  }
  
  Sigma_list <- as.matrix(Matrix::bdiag(lapply(unique(subject_indicator), function(x) { 
    Sigma
  })))
  epsilon <- MASS::mvrnorm(N/length(unique(subject_indicator)), rep(0,M*length(unique(subject_indicator))), Sigma_list, tol = 1e-6)
  covariates <- MASS::mvrnorm(N, c(rep(0,M), rep(1,M)), as.matrix(Matrix::bdiag(posdef.matrix(M, seed=s), posdef.matrix(M, seed=s*2000))))
  X_1 <- covariates[,1:M]
  X_2 <- covariates[,(M+1):(2*M)]
  new_eps <- do.call(rbind, lapply(unique(subject_indicator), function(x){ epsilon[,((x-1)*M+1):(((x-1)*M+M))] }))
  error <- t(epsilon)
  Xbeta <- matrix(0, M, N)
  for(k in 1:length(unique(subject_indicator))){
    Xbeta[,min(which(subject_indicator==k)):max(which(subject_indicator==k))] <-
      (matrix(rep(intercepts[[k]],length(min(which(subject_indicator==k)):max(which(subject_indicator==k)))), M)) +
      apply(X_1[min(which(subject_indicator==k)):max(which(subject_indicator==k)),], 1, function(x) betas[[k]][,1]*x) + 
      apply(X_2[min(which(subject_indicator==k)):max(which(subject_indicator==k)),], 1, function(x) betas[[k]][,2]*x)
  }
  sample <- list()
  sample[[1]] <- qpois(pnorm(error, sd = sqrt(diag(Sigma))), exp(Xbeta))
  sample[[2]] <- X_1
  sample[[3]] <- X_2
  
  data_short <- cbind(id=seq(1,N), X1=sample[[2]], X2=sample[[3]])
  colnames(data_short)[2:(M+1)] <- paste("X1", seq(1, M, 1), sep="_")
  colnames(data_short)[(M+2):(2*M+1)] <- paste("X2", seq(1, M, 1), sep="_")
  data <- reshape(as.data.frame(data_short), direction="long", varying=c(2:(2*M+1)), sep="_")
  data <- data[order(data$id),-which(names(data)=="time")]
  data <- cbind(data, response=c(sample[[1]]))
  return(data)
}
