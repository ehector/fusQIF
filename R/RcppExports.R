# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

delta_constructor <- function(p, J, K) {
    .Call(`_fusQIF_delta_constructor`, p, J, K)
}

X_constructor <- function(covariates, M, N, p) {
    .Call(`_fusQIF_X_constructor`, covariates, M, N, p)
}

increQIF_sub <- function(X, y, nobs, family, corstr, beta_old, maxit, tol) {
    .Call(`_fusQIF_increQIF_sub`, X, y, nobs, family, corstr, beta_old, maxit, tol)
}

diqif_min <- function(X, y, nobs, family, corstr, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit, tol) {
    .Call(`_fusQIF_diqif_min`, X, y, nobs, family, corstr, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit, tol)
}

diqif_eval <- function(X, y, nobs, family, corstr, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t) {
    .Call(`_fusQIF_diqif_eval`, X, y, nobs, family, corstr, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t)
}

doADMM_MCP <- function(X, y, nobs, family, corstr, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2) {
    .Call(`_fusQIF_doADMM_MCP`, X, y, nobs, family, corstr, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2)
}

try_doADMM_MCP <- function(X, y, nobs, family, corstr, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2) {
    .Call(`_fusQIF_try_doADMM_MCP`, X, y, nobs, family, corstr, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2)
}

linear_min <- function(X, y, nobs, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t) {
    .Call(`_fusQIF_linear_min`, X, y, nobs, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t)
}

binomial_min <- function(X, y, nobs, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit, tol) {
    .Call(`_fusQIF_binomial_min`, X, y, nobs, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit, tol)
}

poisson_min <- function(X, y, nobs, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit, tol) {
    .Call(`_fusQIF_poisson_min`, X, y, nobs, beta_old, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit, tol)
}

linear_eval <- function(X, y, nobs, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t) {
    .Call(`_fusQIF_linear_eval`, X, y, nobs, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t)
}

binomial_eval <- function(X, y, nobs, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t) {
    .Call(`_fusQIF_binomial_eval`, X, y, nobs, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t)
}

poisson_eval <- function(X, y, nobs, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t) {
    .Call(`_fusQIF_poisson_eval`, X, y, nobs, beta, J, K, M, N, p, n_k, m_j, rho, gamma, t)
}

indep_doADMM_MCP <- function(X, y, nobs, family, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2) {
    .Call(`_fusQIF_indep_doADMM_MCP`, X, y, nobs, family, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2)
}

try_indep_doADMM_MCP <- function(X, y, nobs, family, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2) {
    .Call(`_fusQIF_try_indep_doADMM_MCP`, X, y, nobs, family, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, tol_1, tol_2, maxit_1, maxit_2)
}

