qiffusion <- function(x, ...) UseMethod("qiffus")

print.qiffus <- function(x, ...){
  cat("Call:\n")
  print(x$call)
  cat("\nCoefficients:\n")
  print(x$coefficients)
  cat("\nVariance:\n")
  print(x$vcov)
}

summary.qiffus <- function(object, ...){
  se <- sqrt(diag(object$vcov))
  zval <- coef(object) / se
  TAB <- cbind(Estimate = coef(object),
               StdErr = se,
               z.value = zval,
               p.value = 2*pnorm(-abs(zval)))
  res <- list(call=object$call,
              coefficients=TAB)
  class(res) <- "summary.qiffus"
  return(res) 
}

print.summary.qiffus <- function(x, ...)
{
  cat("Call:\n")
  print(x$call)
  cat("\n")
  printCoefmat(x$coefficients, P.values=TRUE, has.Pvalue=TRUE)
}

qiffusion <- function(formula, data, id=id, response_indicator=NULL, subject_indicator=NULL, family, corstr, lambda_seq=NULL, rho=1, ...){
  cl <- match.call()
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mf[,"(id)"] <- data[,id]
  mt <- attr(mf, "terms")
  response <- data.frame(id=mf[,"(id)"], model.response(mf, "numeric"))
  covariates <- data.frame(id=mf[,"(id)"], model.matrix(mt, mf))
  colnames(covariates)[match("X.Intercept.",colnames(covariates))] <- "(Intercept)"
  
  if (is.null(response_indicator)){
    warning(gettextf("Response indicator is null. Using one response block"), domain = NA)
    response_indicator <- as.vector(rep(1, table(response[,"id"])[1]))
  }
  if (is.null(subject_indicator)){
    warning(gettextf("Subject indicator is null. Using one subject block"), domain=NA)
    subject_indicator <- as.vector(rep(1, length(unique(response[,"id"]))))
  }
  if (is.null(lambda_seq)) lambda_seq <- seq(0, 10, 0.2)
  if (is.empty.model(mt)) stop("Model is empty")
  if (family != "gaussian" & family!="binomial" & family!="poisson") {
    warning(gettextf("family = '%s' is not supported. Using gaussian", family), domain = NA)
    family <- "gaussian"
  }
  if (corstr != "AR1" & corstr != "CS" & corstr != "independence" & corstr != "CS+AR1") stop("corstr must be supported type")
  
  response <- cbind(matrix(response[,-which(colnames(response) == "id")], length(response_indicator), length(subject_indicator)), response_indicator)
  
  output <- qiffusion.compute.mean(response, covariates, response_indicator, subject_indicator, family, corstr, lambda_seq, rho)
  
  names(output$coefficients) <- c(sapply(unique(output$comb_scheme), function(x) paste(colnames(covariates)[-1], x,sep="-")))
  colnames(output$vcov) <- names(output$coefficients)
  rownames(output$vcov) <-names(output$coefficients)
  output$response <- data.frame(id=mf[,"(id)"], model.response(mf, "numeric"))
  output$covariates <- data.frame(id=mf[,"(id)"], model.matrix(mt, mf))
  colnames(output$response)[1] <- id
  output$id <- id
  
  output <- c(output, list(call=cl, formula=formula))
  class(output) <- "qiffus"
  return(output)
}

qiffusion.compute.mean <- function(response, covariates, response_indicator, subject_indicator, family, corstr, lambda_seq, rho, ...){
  time1 <- proc.time()
  output <- list()
  
  colnames(response) <- c(subject_indicator,"response_indicator")
  
  J <- length(unique(response_indicator))
  K <- length(unique(subject_indicator))
  N <- length(subject_indicator)
  M <- length(response_indicator)
  p <- dim(covariates)[2]-1
  
  if (corstr=="CS" | corstr=="AR1") d <- 2
  if (corstr=="independence") d <- 1
  if (corstr=="CS+AR1") d <- 3
  
  psi_g <- matrix(0, N, J*K*(p*d))
  p_1 <- p*d
  MCLE_mean <- list()
  print("Computing block coefficients.", quote=FALSE)
  time1 <- proc.time()-time1
  time_jk <- list()
  for (k in 1:K) {
    MCLE_mean[[k]] <- list()
    time_jk[[k]] <- list()
    for (j in 1:J) {
      time_j <- proc.time()
      block_y <- matrix(response[response[,(colnames(response)=="response_indicator")]==j, which(colnames(response)==k)], 
                        nrow=sum(response[,(colnames(response)=="response_indicator")]==j))
      ids_block_x <- covariates[covariates[,"id"] %in% unique(covariates[,"id"])[subject_indicator==k],]
      block_x <- as.data.frame(ids_block_x[rep(response_indicator, length(unique(ids_block_x[,"id"])))==j,])
      id <- block_x[,which(colnames(block_x)=="id")]
      block_x <- as.matrix(block_x[, -which(colnames(block_x)=="id")])
      m <- dim(block_y)[1]
      n <- dim(block_y)[2]
      init_betas <- coef(glm(c(block_y) ~ 0 + block_x, family=family))
      qif_block_fit <- increQIF_sub(block_x, c(block_y), nobs=rep(m,n), family, corstr, init_betas, maxit=10000, tol=1e-6)
      MCLE_jk <- as.vector(qif_block_fit$beta_sub)
      MCLE_mean[[k]][[j]] <- MCLE_jk
      time_jk[[k]][[j]] <- proc.time()-time_j
    }
    time_jk[[k]] <- do.call(rbind, time_jk[[k]])
  }
  time_list <- do.call(rbind, time_jk)
  time_max <- time_list[which.max(time_list[,1]),]
  time2 <- proc.time()
  
  print("Fusing coefficients.", quote=FALSE)
  Delta <- delta_constructor(p,J,K)
  init_betas <- unlist(MCLE_mean)
  init_gamma <- Delta %*% init_betas
  init_t <- vector("numeric", (K*(K-1)/2*J*J+K*J*(J-1)/2)*p)
  
  nobs <- lapply(1:K, function(k) {
    lapply(1:J, function(j) rep(table(response_indicator)[j], table(subject_indicator)[k])) })
  
  tol_1 <- 1e-6
  tol_2 <- 1e-6
  maxit_1 <- 10000
  maxit_2 <- 10000
  X <- X_constructor(as.matrix(covariates[,-match("id",colnames(covariates))]), M, N, p)
  n_k_list <- as.vector(table(subject_indicator))
  m_j_list <- as.vector(table(response_indicator))
  M <- dim(X)[1]
  N <- dim(X)[2]
  p <- dim(X)[3]
  
  opt_list <- list()
  stop <- FALSE
  l <- 1
  while(l<=length(lambda_seq) & !stop){
    opt_list[[l]] <- tryCatch(
      try_doADMM_MCP(X, y=response[,-dim(response)[2]], nobs, family, corstr, init_betas, init_gamma, init_t, J, K, M, N, p,
                                      n_k=n_k_list, m_j=m_j_list, lambda=lambda_seq[l], rho=rho, delta=3, 
                                      tol_1=tol_1, tol_2, maxit_1=maxit_1, maxit_2=maxit_2),
      error = function(c) list(message="error")
    )
    if(!opt_list[[l]]$convergence) {
      opt_list[[l]]$message <- "failed to converge"
      warning(paste0("Iteration ", l, " failed to converge. Increase maximum iterations."), domain = NA)
    }
    if(opt_list[[l]]$message==""){
      diff_pattern <- matrix(0, length(opt_list[[l]]$gamma),5)
      it <- 1
      for (k in 1:K){
        for(j in 1:(J-1)){
          for(j_p in (j+1):J){
            for(q in 1:p){
              diff_pattern[it,] <- c(j,k,j_p,k,q)
              it <- it + 1 
            }
          }
        }
        if(k!=K){
          for(k_p in (k+1):K){
            for(j in 1:J){
              for(j_p in 1:J){
                for(q in 1:p){
                  diff_pattern[it,] <- c(j,k,j_p,k_p,q)
                  it <- it + 1  
                }
              }
            }
          } 
        }
      }
      comb_pattern <- unique(diff_pattern[opt_list[[l]]$gamma==0,1:4,drop=FALSE])[apply(unique(diff_pattern[opt_list[[l]]$gamma==0,1:4,drop=FALSE]), 1, function(x) {
        sum(apply(diff_pattern[opt_list[[l]]$gamma==0, 1:4,drop=FALSE], 1, function(y) sum(y==x))==4)==p
      }),1:4,drop=FALSE]
      comb_scheme <- c(1:(J*K))
      if(dim(comb_pattern)[1]!=0){
        for(it in 1:(dim(comb_pattern)[1])){
          j <- comb_pattern[it,3]
          k <- comb_pattern[it,4]
          comb_scheme[J*(k-1)+j] <- comb_scheme[J*(comb_pattern[it,2]-1) + comb_pattern[it,1]]
        }
      }
      opt_list[[l]]$comb_scheme <- as.numeric(factor(comb_scheme, levels=unique(comb_scheme), ordered=TRUE))
      
      opt_list[[l]]$BIC <- N*diqif_eval(X, y=response[,-dim(response)[2]], nobs, family, corstr, opt_list[[l]]$beta, J, K, M, N, p, 
                                        n_k=n_k_list, m_j=m_j_list, rho=1, gamma=opt_list[[l]]$gamma, t=opt_list[[l]]$t)$Obj -
        (J*K*d-length(unique(opt_list[[l]]$comb_scheme)))*p*log(N)
      
      init_betas <- opt_list[[l]]$beta
      init_gamma <- Delta %*% init_betas 
      # if(family=="binomial" & l!=1) if(opt_list[[l]]$BIC>opt_list[[l-1]]$BIC) stop <- TRUE 
    }
    if(length(unique(opt_list[[l]]$comb_scheme))==1) stop <- TRUE
    if(opt_list[[l]]$message!="") stop <- TRUE
    l <- l+1
  }
  best_l <- which.min(lapply(opt_list[unlist(lapply(opt_list, function(x) x$message==""))], function(x) x$BIC))
  comb_scheme <- opt_list[[best_l]]$comb_scheme
  final_opt <- diqif_eval(X, y=response[,-dim(response)[2]], nobs, family, corstr, beta=opt_list[[best_l]]$beta, J, K, M, N, p, n_k=n_k_list, m_j=m_j_list, rho=1, 
                                   gamma=opt_list[[best_l]]$gamma, t=opt_list[[best_l]]$t)
  psi_g <- final_opt$gi_list
  S <- t(do.call(rbind, lapply(1:(K*J), function(x) final_opt$G_sub[((x-1)*p_1+1):((x-1)*p_1+p_1),((x-1)*p+1):((x-1)*p+p)])))
  
  unordered_psi_g <- psi_g
  unordered_S <- S
  
  psi_g <- unordered_psi_g[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
  S_all_reordered <- S[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
  
  print("Computing combined estimate.", quote=FALSE)
  rank <- qr(psi_g)$rank
  if(dim(psi_g)[2] != rank){
    psi_g_pca <- prcomp(psi_g, center = FALSE,scale. = FALSE)
    old_psi_g <- psi_g
    old_J <- J
    old_K <- K
    old_S <- S
    ta <- cumsum(table(comb_scheme))
    psi_g <- psi_g_pca$x[,1:rank]
    S_matrix <- (as.matrix(Matrix::bdiag(lapply(sort(unique(comb_scheme)), function(c) {
      S_all_reordered[,sapply(comb_scheme, function(x) rep(x,p_1))[order(sapply(comb_scheme, function(x) rep(x,p_1)))]==c]
    }))) %*% psi_g_pca$rotation )[,1:rank]
    S_MCLE_mean_matrix <- (t(final_opt$G_sub%*%opt_list[[best_l]]$beta)[,order(sapply(comb_scheme, function(x) rep(x,p_1)))] %*% 
        psi_g_pca$rotation)[,1:rank]
  } else {
    S_matrix <- as.matrix(Matrix::bdiag(lapply(sort(unique(comb_scheme)), function(c) {
      S_all_reordered[,sapply(comb_scheme, function(x) rep(x,p_1))[order(sapply(comb_scheme, function(x) rep(x,p_1)))]==c]
    })))
    S_MCLE_mean_matrix <- t(final_opt$G_sub%*%opt_list[[best_l]]$beta)[,order(sapply(comb_scheme, function(x) rep(x,p_1)))]
  }
  
  V_psi <- t(psi_g)%*%psi_g/N
  output$W <- solve(V_psi)
  
  scale <- S_matrix %*% output$W %*% t(S_matrix)
  main <- S_matrix %*% output$W %*% matrix(S_MCLE_mean_matrix, ncol=1)
  output$coefficients <- as.vector(solve(scale)%*%main)
  output$vcov <- solve(scale)*N
  
  time2 <- proc.time()-time2
  output$comb_scheme <- comb_scheme
  output$family <- family
  output$corstr <- corstr
  output$lambda_seq <- lambda_seq
  output$response_indicator <- response_indicator
  output$subject_indicator <- subject_indicator
  output$J <- J
  output$K <- K
  output$MCLE <- list()
  output$MCLE$beta <- MCLE_mean
  output$opt_list <- opt_list
  output$time <- time1+time2+time_max
  
  return(output)
}
