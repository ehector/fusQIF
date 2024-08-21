// [[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include <math.h>
#include <RcppArmadillo.h>
#include <sstream>
using namespace Rcpp; 
using namespace arma;

// [[Rcpp::export]]
arma::mat delta_constructor(const int& p, const int& J, const int& K)
{
  arma::mat Delta=zeros<mat>((K*(K-1)/2*J*J+K*J*(J-1)/2)*p, J*K*p);
  int it = 0;
  
  for (int k=0; k<K; k++){
    for(int j=0; j<J-1; j++){
      for(int j_p=j+1; j_p<J; j_p++){
        for(int q=0; q<p; q++){
          Delta(it, (k*J+j)*p+q) += 1;
          Delta(it, (k*J+j_p)*p+q) += -1;
          it += 1; 
        }
      }
    }
    if(k!=(K-1)){
      for(int k_p=k+1; k_p<K; k_p++){
        for(int j=0; j<J; j++){
          for(int j_p=0; j_p<J; j_p++){
            for(int q=0; q<p; q++){
              Delta(it, (k*J+j)*p+q) += 1;
              Delta(it, (k_p*J+j_p)*p+q) += -1;
              it += 1;  
            }
          }
        }
      } 
    }
  }
  return Delta;
}

// [[Rcpp::export]]
arma::cube X_constructor(const arma::mat& covariates, const int& M, const int& N, const int& p){
  arma::cube X(M, N, p);
  for(int q=0; q<p; q++){
    X.slice(q) = reshape(covariates.col(q), M, N);
  }
  return X;
}

//[[Rcpp::export]]
List increQIF_sub(arma::mat X, arma::vec y, arma::vec nobs, String family, String corstr,
                  arma::vec beta_old, int maxit, double tol){
  
  int niter2= 0;
  bool stop_flag2= FALSE;
  bool converged2= FALSE;
  
  int N = X.n_rows;
  int p = X.n_cols;
  int n = nobs.n_elem;
  
  arma::vec gb_sub;
  arma::mat Gb_sub;
  arma::mat Cb_sub;
  
  arma::vec beta_sub = beta_old;
  
  arma::vec mu;
  arma::vec vu;
  arma::mat M0;
  arma::mat M1;
  arma::mat M2;
  
  arma::mat qif1dev_sub;
  arma::mat qif2dev_sub;
  
  
  while (!stop_flag2){
    niter2 += 1;
    if(corstr=="independence"){
      gb_sub = zeros<vec>(p );
      Gb_sub = zeros<mat>(p , p);
      Cb_sub = zeros<mat>(p , p);
    }else if(corstr=="CS+AR1"){
      gb_sub = zeros<vec>(p * 3);
      Gb_sub = zeros<mat>(p * 3, p);
      Cb_sub = zeros<mat>(p * 3, p * 3);
    }else{
      gb_sub = zeros<vec>(p * 2);
      Gb_sub = zeros<mat>(p * 2, p);
      Cb_sub = zeros<mat>(p * 2, p * 2);
    }
    
    arma::vec eta2 = X * beta_sub;
    
    if(family == "gaussian") {
      mu = eta2 ;
      vu = ones<vec>(N) ;
    } else if(family == "binomial") {
      mu = exp(eta2)/( 1 + exp(eta2) ) ;
      vu = exp(eta2)/(pow( (1 + exp(eta2)), 2 )) ;
    } else if(family == "poisson") {
      mu = exp(eta2) ;
      vu = exp(eta2) ;
    } else{Rcpp::stop("Unknown distribution family\n");}
    
    int loc1 = 0;
    int loc2 = -1;
    for (int i = 0; i < n; i++){
      loc1 = loc2 + 1 ;
      loc2 = loc1 + nobs[i] -1;
      
      arma::mat Xi = X.rows(loc1, loc2) ;
      arma::vec yi = y.subvec(loc1, loc2);
      arma::vec mui = mu.subvec(loc1, loc2);
      arma::vec vui = vu.subvec(loc1, loc2) ;
      
      int mi = nobs[i] ;
      
      arma::mat Ai_half = diagmat(sqrt(vui)) ;
      arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
      
      if (corstr == "independence") {
        M0 = eye<mat>(mi, mi);
        
      } else if (corstr == "CS") {
        M0 = eye<mat>(mi, mi);
        M1 = ones(mi, mi) - M0;  
      } else if (corstr == "AR1") {
        M0 = eye<mat>(mi, mi);
        M1 = zeros<mat>(mi, mi);
        for (int j = 0; j < mi; j++ ){
          for (int k = 0; k < mi; k++){
            if(abs(j-k)==1){ M1(j, k) = 1; }
          }
        }           
      } else if (corstr == "unstructured"){
        int m = N/n;
        M0 = eye<mat>(m, m);
        M1 = zeros<mat>(m, m);
        
        int loc3 = 0;
        int loc4 = -1;
        for (int i = 0; i < n; i++){
          loc3 = loc4 +1;
          loc4 = loc3 + nobs[i] -1;
          arma::mat Xi = X.rows(loc3, loc4);
          arma::vec yi = y.subvec(loc3, loc4);
          arma::vec mui = mu.subvec(loc3, loc4);
          M1 += (yi - mui) * (yi - mui).t(); 
        }
        M1 = M1 / n;
      }else if(corstr=="CS+AR1"){
        M0 = eye<mat>(mi, mi);
        M1 = ones(mi, mi) - M0;  
        M2 = zeros<mat>(mi, mi);
        for (int j = 0; j < mi; j++ ){
          for (int k = 0; k < mi; k++){
            if(abs(j-k)==1){ M2(j, k) = 1; }
          }
        }     
      }else {Rcpp::stop("Unknown correlation structure\n");}   
      
      if(corstr=="independence"){
        vec gi_sub = Xi.t() * (yi - mui);
        gb_sub += gi_sub;
        Gb_sub += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
        Cb_sub += gi_sub * gi_sub.t();
        
      }else if(corstr == "CS+AR1"){
        vec gi_sub = join_cols( join_cols( Xi.t() * (yi - mui), 
                                           Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
                                           Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
        gb_sub += gi_sub;
        Gb_sub += join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                       Xi.t() * Ai_half * M1 * Ai_half * Xi),
                                       Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
        Cb_sub += gi_sub * gi_sub.t();
        
      }else{
        arma::vec gi_sub = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
        gb_sub += gi_sub;
        Gb_sub += join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                            Xi.t() * Ai_half * M1 * Ai_half * Xi) ;
        Cb_sub += gi_sub * gi_sub.t();
        
      }
    }       
    
    qif1dev_sub = Gb_sub.t() * pinv(Cb_sub) * gb_sub ;
    
    qif2dev_sub = Gb_sub.t() * pinv(Cb_sub) * Gb_sub ;
    
    vec d_beta_sub = solve(qif2dev_sub, qif1dev_sub);
    
    double df_beta_sub = as_scalar(qif1dev_sub.t() * d_beta_sub);
    
    beta_sub += d_beta_sub;
    
    if(fabs(df_beta_sub) < tol) {converged2 = TRUE; stop_flag2 = TRUE;}
    if (niter2 > maxit) {stop_flag2 = TRUE;}
  }
  
  return List::create(Named("beta_sub") = beta_sub);
}

//[[Rcpp::export]]
List diqif_min(const arma::cube& X, const arma::mat& y, const List& nobs, const String& family, const String& corstr,
               const arma::vec& beta_old, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
               const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t, const int& maxit, const double& tol){
  
  int niter2= 0;
  bool stop_flag2= FALSE;
  bool converged2= FALSE;
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::vec gb_sub;
  arma::mat Gb_sub;
  arma::mat Cb_sub;
  arma::mat gi_list;
  
  //initialization by the beta estimated from previous data
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_sub = beta_old;
  arma::vec beta_jk;
  
  arma::vec mu;
  arma::vec vu;
  arma::mat M0;
  arma::mat M1;
  arma::mat M2;
  
  arma::mat Delta = delta_constructor(p, J, K);
  arma::mat qif1dev_sub;
  arma::mat qif2dev_sub;
  
  while (!stop_flag2){
    niter2 += 1;
    //reset gb to 0 after each iteration
    if(corstr=="independence"){
      gb_sub = zeros<vec>(p*J*K);
      Gb_sub = zeros<mat>(p*J*K, p*J*K);
      Cb_sub = zeros<mat>(p*J*K, p*J*K);
      gi_list = zeros<mat>(N, J*K*p);
    }else if(corstr=="CS+AR1"){
      gb_sub = zeros<vec>(p*3*J*K);
      Gb_sub = zeros<mat>(p*3*J*K, p*J*K);
      Cb_sub = zeros<mat>(p*3*J*K, p*3*J*K);
      gi_list = zeros<mat>(N, p*3*J*K);
    }else{
      gb_sub = zeros<vec>(p*2*J*K);
      Gb_sub = zeros<mat>(p*2*J*K, p*J*K);
      Cb_sub = zeros<mat>(p*2*J*K, p*2*J*K);
      gi_list = zeros<mat>(N, p*2*J*K);
    }
    
    n_k_cumsum = 0;
    
    for (int k=0; k<K; k++){
      m_j_cumsum = 0;
      for(int j=0; j<J; j++){
        //rows and columns of Y and X have to be sorted
        
        block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
        x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
        block_x = x_p.slice(0).as_col();
        for(int q=1; q<p; q++){
          block_x = join_rows(block_x, x_p.slice(q).as_col());
        }
        
        beta_jk = beta_sub.subvec((k*J+j)*p, (k*J+j)*p+p-1);
        nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
        int n = nobs_jk.n_elem;
        
        // update gb_new with beta_new over iterations
        arma::vec eta2 = block_x * beta_jk;
        
        if(family == "gaussian") {
          mu = eta2 ;
          vu = ones<vec>(sum(nobs_jk)) ;
        } else if(family == "binomial") {
          mu = exp(eta2)/( 1 + exp(eta2) ) ;
          vu = exp(eta2)/(pow( (1 + exp(eta2)), 2 )) ;
        } else if(family == "poisson") {
          mu = exp(eta2) ;
          vu = exp(eta2) ;
        } else{Rcpp::stop("Unknown distribution family\n");}

        int loc1 = 0;
        int loc2 = -1;
        for (int i = 0; i < n; i++){
          loc1 = loc2 + 1 ;
          loc2 = loc1 + nobs_jk[i] -1;
          arma::mat Xi = block_x.rows(loc1, loc2) ;
          arma::vec yi = block_y.subvec(loc1, loc2);
          arma::vec mui = mu.subvec(loc1, loc2);
          arma::vec vui = vu.subvec(loc1, loc2) ;
          int mi = nobs_jk[i] ;
          
          arma::mat Ai_half = diagmat(sqrt(vui)) ;
          arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
          
          if (corstr == "independence") {
            M0 = eye<mat>(mi, mi);
          } else if (corstr == "CS") {
            M0 = eye<mat>(mi, mi);
            // M1 is a matrix with 1 on off-diagonal elements
            M1 = ones(mi, mi) - M0;  
          } else if (corstr == "AR1") {
            M0 = eye<mat>(mi, mi);
            M1 = zeros<mat>(mi, mi);
            for (int j_p = 0; j_p < mi; j_p++ ){
              for (int k_p = 0; k_p < mi; k_p++){
                if(std::abs(j_p-k_p)==1){ M1(j_p, k_p) = 1; }
              }
            }           
          } else if (corstr == "unstructured"){
            int m = N/n;
            M0 = eye<mat>(m, m);
            M1 = zeros<mat>(m, m);
            int loc3 = 0;
            int loc4 = -1;
            for (int i_p = 0; i_p < n; i_p++){
              loc3 = loc4 +1;
              loc4 = loc3 + nobs_jk[i] -1;
              arma::mat Xi = block_x.rows(loc3, loc4);
              arma::vec yi = block_y.subvec(loc3, loc4);
              arma::vec mui = mu.subvec(loc3, loc4);
              M1 += (yi - mui) * (yi - mui).t(); 
            }
            M1 = M1 / n;
          }else if(corstr=="CS+AR1"){
            M0 = eye<mat>(mi, mi);
            M1 = ones(mi, mi) - M0;  
            M2 = zeros<mat>(mi, mi);
            for (int j_p = 0; j_p < mi; j_p++ ){
              for (int k_p = 0; k_p < mi; k_p++){
                if(std::abs(j_p-k_p)==1){ M2(j_p, k_p) = 1; }
              }
            }     
          }else {Rcpp::stop("Unknown correlation structure\n");}   
          
          if(corstr=="independence"){
            arma::vec gi_sub = Xi.t() * (yi - mui);
            gi_list(n_k_cumsum+ i, span((k*J+j)*p, (k*J+j)*p+p-1)) = gi_sub.t();
            gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1)) += gi_sub;
            Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
            //Cb_sub += gi_sub * gi_sub.t();
            
          }else if(corstr == "CS+AR1"){
            arma::vec gi_sub = join_cols( join_cols( Xi.t() * (yi - mui), 
                                                     Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
                                                     Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
            gi_list(n_k_cumsum+ i, span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1)) = gi_sub.t();
            gb_sub(span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1)) += gi_sub;
            Gb_sub(span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += 
              join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                   Xi.t() * Ai_half * M1 * Ai_half * Xi),
                                   Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
            //Cb_sub += gi_sub * gi_sub.t();
            
          }else{
            arma::vec gi_sub = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
            gi_list(n_k_cumsum+i, span((k*J+j)*p*2, (k*J+j)*p*2+p*2-1)) = gi_sub.t();
            gb_sub(span((k*J+j)*p*2, (k*J+j)*p*2+p*2-1)) += gi_sub;
            Gb_sub(span((k*J+j)*p*2, (k*J+j)*p*2+p*2-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += 
              join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                        Xi.t() * Ai_half * M1 * Ai_half * Xi);
            //Cb_sub += gi_sub * gi_sub.t();
          }
        } 
        m_j_cumsum += m_j[j];
      }
      n_k_cumsum += n_k[k];
    }
    
    Cb_sub = gi_list.t() * gi_list;
    //qif1dev_sub = 2*Gb_sub.t() * pinv(Cb_sub) * gb_sub ;
    //qif2dev_sub = 2*Gb_sub.t() * pinv(Cb_sub) * Gb_sub ;
    
    //Now put in the penalties
    //Rcout << "pinv() for Cb_sub" << std::endl;
    arma::mat obj1_dev = Gb_sub.t() * pinv(Cb_sub) * gb_sub/N - rho * Delta.t()*(Delta*beta_sub - gamma + t/rho);
    arma::mat obj2_dev = Gb_sub.t() * pinv(Cb_sub) * Gb_sub/N + rho * Delta.t()*Delta;
    
    vec d_beta_sub = solve(obj2_dev, obj1_dev);
    //Rcout << "pinv() for iteration " << niter2 << std::endl;
    //vec d_beta_sub = pinv(obj2_dev) * obj1_dev;
    double df_beta_sub = as_scalar(obj1_dev.t() * d_beta_sub);
    
    beta_sub += d_beta_sub;
    if(std::abs(df_beta_sub) < tol) {converged2 = TRUE; stop_flag2 = TRUE;}
    if (niter2 > maxit) {stop_flag2 = TRUE;}
  }
  
  //if (converged2==FALSE) {Rcpp::stop("algorithm for single data batch reached 'maxit' but did not converge\n"); }
  
  mat beta_cov = pinv(Gb_sub.t() * pinv(Cb_sub) * Gb_sub);
  //vec beta_sd = sqrt(beta_cov.diag());
  //double obj = as_scalar(gb_sub.t() * pinv(Cb_sub) * gb_sub) + penalty(rho, beta_sub, gamma, t);
  
  return List::create(Named("beta_sub") = beta_sub,
                      Named("beta_cov") = beta_cov,
                      //Named("Obj") = obj,
                      Named("g_sub") = gb_sub,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub,
                      Named("C_sub") = Cb_sub,
                      Named("obj1_dev_1") = Gb_sub.t() * pinv(Cb_sub) * gb_sub,
                      Named("obj1_dev_2") = rho * Delta.t()*(Delta*beta_sub - gamma -t/rho),
                      Named("obj2_dev_1") = Gb_sub.t() * pinv(Cb_sub) * Gb_sub,
                      Named("obj2_dev_2") = rho * Delta.t()*Delta
  );
}

//[[Rcpp::export]]
List diqif_eval(const arma::cube& X, const arma::mat& y, const List& nobs, const String& family, const String& corstr,
                const arma::vec& beta, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t){
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::vec gb_sub;
  arma::mat Gb_sub;
  arma::mat Cb_sub;
  arma::mat gi_list;
  if(corstr=="independence"){
    arma::vec gb_sub = zeros<vec>(p*J*K);
    Gb_sub = zeros<mat>(p*J*K, p*J*K);
    Cb_sub = zeros<mat>(p*J*K, p*J*K);
    gi_list = zeros<mat>(N, J*K*p);
  }else if(corstr=="CS+AR1"){
    gb_sub = zeros<vec>(p*3*J*K);
    Gb_sub = zeros<mat>(p*3*J*K, p*J*K);
    Cb_sub = zeros<mat>(p*3*J*K, p*3*J*K);
    gi_list = zeros<mat>(N, p*3*J*K);
  }else{
    gb_sub = zeros<vec>(p*2*J*K);
    Gb_sub = zeros<mat>(p*2*J*K, p*J*K);
    Cb_sub = zeros<mat>(p*2*J*K, p*2*J*K);
    gi_list = zeros<mat>(N, p*2*J*K);
  }
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_jk;
  
  arma::vec mu;
  arma::vec vu;
  arma::mat M0;
  arma::mat M1;
  arma::mat M2;
  
  n_k_cumsum = 0;
  
  for (int k=0; k<K; k++){
    m_j_cumsum = 0;
    for(int j=0; j<J; j++){
      //rows and columns of Y and X have to be sorted
      
      block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
      x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
      block_x = x_p.slice(0).as_col();
      for(int q=1; q<p; q++){
        block_x = join_rows(block_x, x_p.slice(q).as_col());
      }
      
      beta_jk = beta.subvec((k*J+j)*p, (k*J+j)*p+p-1);
      nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
      int n = nobs_jk.n_elem;
      
      // update gb_new with beta_new over iterations
      arma::vec eta2 = block_x * beta_jk;
      
      if(family == "gaussian") {
        mu = eta2 ;
        vu = ones<vec>(sum(nobs_jk)) ;
      } else if(family == "binomial") {
        mu = exp(eta2)/( 1 + exp(eta2) ) ;
        vu = exp(eta2)/(pow( (1 + exp(eta2)), 2 )) ;
      } else if(family == "poisson") {
        mu = exp(eta2) ;
        vu = exp(eta2) ;
      } else{Rcpp::stop("Unknown distribution family\n");}
      
      int loc1 = 0;
      int loc2 = -1;
      for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs_jk[i] -1;
        arma::mat Xi = block_x.rows(loc1, loc2) ;
        arma::vec yi = block_y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        arma::vec vui = vu.subvec(loc1, loc2) ;
        int mi = nobs_jk[i] ;
        
        arma::mat Ai_half = diagmat(sqrt(vui)) ;
        arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
        
        if (corstr == "independence") {
          M0 = eye<mat>(mi, mi);
        } else if (corstr == "CS") {
          M0 = eye<mat>(mi, mi);
          // M1 is a matrix with 1 on off-diagonal elements
          M1 = ones(mi, mi) - M0;  
        } else if (corstr == "AR1") {
          M0 = eye<mat>(mi, mi);
          M1 = zeros<mat>(mi, mi);
          for (int j_p = 0; j_p < mi; j_p++ ){
            for (int k_p = 0; k_p < mi; k_p++){
              if(std::abs(j_p-k_p)==1){ M1(j_p, k_p) = 1; }
            }
          }           
        } else if (corstr == "unstructured"){
          int m = N/n;
          M0 = eye<mat>(m, m);
          M1 = zeros<mat>(m, m);
          int loc3 = 0;
          int loc4 = -1;
          for (int i_p = 0; i_p < n; i_p++){
            loc3 = loc4 +1;
            loc4 = loc3 + nobs_jk[i] -1;
            arma::mat Xi = block_x.rows(loc3, loc4);
            arma::vec yi = block_y.subvec(loc3, loc4);
            arma::vec mui = mu.subvec(loc3, loc4);
            M1 += (yi - mui) * (yi - mui).t(); 
          }
          M1 = M1 / n;
        }else if(corstr=="CS+AR1"){
          M0 = eye<mat>(mi, mi);
          M1 = ones(mi, mi) - M0;  
          M2 = zeros<mat>(mi, mi);
          for (int j_p = 0; j_p < mi; j_p++ ){
            for (int k_p = 0; k_p < mi; k_p++){
              if(std::abs(j_p-k_p)==1){ M2(j_p, k_p) = 1; }
            }
          }     
        }else {Rcpp::stop("Unknown correlation structure\n");}   
        
        if(corstr=="independence"){
          arma::vec gi_sub = Xi.t() * (yi - mui);
          gi_list(n_k_cumsum+ i, span((k*J+j)*p, (k*J+j)*p+p-1)) = gi_sub.t();
          gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1)) += gi_sub;
          Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Ai_half * M0 * Ai_half * Xi ;
        }else if(corstr == "CS+AR1"){
          arma::vec gi_sub = join_cols( join_cols( Xi.t() * (yi - mui), 
                                                   Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui)),
                                                   Xi.t() * Ai_half * M2 * Ai_inv_half * (yi-mui) );
          gi_list(n_k_cumsum+ i, span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1)) = gi_sub.t();
          gb_sub(span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1)) += gi_sub;
          Gb_sub(span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1), span((k*J+j)*p*3, (k*J+j)*p*3+p*3-1)) += 
            join_cols( join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                                 Xi.t() * Ai_half * M1 * Ai_half * Xi),
                                 Xi.t() * Ai_half * M2 * Ai_half * Xi) ;
          
        }else{
          arma::vec gi_sub = join_cols(Xi.t() * (yi - mui), Xi.t() * Ai_half * M1 * Ai_inv_half * (yi-mui) );
          gi_list(n_k_cumsum+i, span((k*J+j)*p*2, (k*J+j)*p*2+p*2-1)) = gi_sub.t();
          gb_sub(span((k*J+j)*p*2, (k*J+j)*p*2+p*2-1)) += gi_sub;
          Gb_sub(span((k*J+j)*p*2, (k*J+j)*p*2+p*2-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += 
            join_cols(Xi.t() * Ai_half * M0 * Ai_half * Xi, 
                      Xi.t() * Ai_half * M1 * Ai_half * Xi);
        }
      } 
      m_j_cumsum += m_j[j];
    }
    n_k_cumsum += n_k[k];
  }
  
  Cb_sub = gi_list.t() * gi_list;
  
  double obj = as_scalar(gb_sub.t() * pinv(Cb_sub) * gb_sub) /N;
  
  return List::create(Named("Obj") = obj,
                      Named("g_sub") = gb_sub,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub,
                      Named("C_sub") = Cb_sub
  );
}

// [[Rcpp::export]]
List doADMM_MCP(const arma::cube& X, const arma::mat& y, const List& nobs, const String& family, const String& corstr,
                const arma::vec& init_betas, const arma::vec& init_gamma, const arma::vec& init_t, const int& J, const int& K, const int& M, 
                const int& N, const int& p, const arma::vec& n_k, const arma::vec& m_j, const double& lambda, const double& rho, 
                const double& delta, const double& tol_1, const double& tol_2, const int& maxit_1, const int& maxit_2){
  std::ostringstream __my_cerr;
  arma::set_cerr_stream(__my_cerr);
  bool convergence = FALSE;
  int r = 0;
  arma::vec old_beta;
  arma::vec beta = init_betas;
  arma::vec gamma = init_gamma;
  arma::vec gamma_unpen;
  arma::vec t = init_t;
  arma::mat Delta = delta_constructor(p, J, K);
  arma::mat Delta_beta;
  double gamma_diff;
  double first_diff;
  double residual;
  
  while((!convergence) & (r <= maxit_1)){
    old_beta = beta;
    
    // First update beta
    List minimized = diqif_min(X, y, nobs, family, corstr, old_beta, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit_2, tol_2);
    arma::vec minimized_beta = minimized[0];
    beta = minimized_beta;
    Delta_beta = Delta * beta;
    
    // Second update gamma
    if(r==0) { gamma_diff = max(abs(Delta_beta-gamma)); }
    gamma_unpen = Delta_beta + t/rho;
    gamma = zeros<vec>((K*(K-1)/2*J*J+K*J*(J-1)/2)*p);
    for(int t=0; t<((K*(K-1)/2*J*J+K*J*(J-1)/2)); t++){
      double sum_abs = sum(abs( gamma_unpen( span(t*p, t*p+p-1) ) ));
      if(sum_abs <= (delta*lambda)){
        for(int q=0; q<p; q++){
          double shrink = std::abs(gamma_unpen[t*p+q])-lambda/rho;
          gamma[t*p+q] = ( sign(gamma_unpen[t*p+q]) * std::max(shrink, 0.0) ) /  ( 1 - 1/ ( delta*rho ) );
        }
      } else {
        gamma( span(t*p, t*p+p-1) ) = gamma_unpen( span(t*p, t*p+p-1) );
      }
    }
    
    // Third update t
    t = t + rho*(Delta_beta - gamma);
    
    residual = sum(pow(Delta_beta - gamma, 2.0));
    if(residual < tol_1){
      convergence = TRUE;
      } else {
        r = r+1;
        }
    }
  
  std::string text(__my_cerr.str());
  
  return List::create(Named("convergence") = convergence,
                      Named("beta") = beta,
                      Named("gamma") = gamma,
                      Named("t") = t,
                      Named("lambda") = lambda,
                      Named("rho") = rho,
                      Named("r") = r,
                      Named("first_diff") = first_diff,
                      Named("residual") = residual,
                      Named("message") = text
  );
}

// [[Rcpp::export]]
List try_doADMM_MCP(const arma::cube& X, const arma::mat& y, const List& nobs, const String& family, const String& corstr,
                    const arma::vec& init_betas, const arma::vec& init_gamma, const arma::vec& init_t, const int& J, const int& K, 
                    const int& M, const int& N, const int& p, const arma::vec& n_k, const arma::vec& m_j, const double& lambda, 
                    const double& rho, const double& delta, const double& tol_1, const double& tol_2, const int& maxit_1, const int& maxit_2){
  List ADMM;
  try{
    ADMM = doADMM_MCP(X, y, nobs, family, corstr, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, 
                      tol_1, tol_2, maxit_1, maxit_2);
  } 
  catch (const std::exception &__ex) {
    return(List::create(Named("message") = std::string("error")));
  }
  catch (...) {
    return(List::create(Named("message") = std::string("error")));
  }
  return(ADMM);
}

//[[Rcpp::export]]
List linear_min(const arma::cube& X, const arma::mat& y, const List& nobs, 
                const arma::vec& beta_old, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t){
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::vec gb_sub;
  arma::mat Gb_sub;
  arma::mat Cb_sub;
  arma::mat gi_list;
  
  //initialization by the beta estimated from previous data
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_sub = beta_old;
  arma::vec beta_jk;
  
  arma::vec mu;
  arma::vec vu;
  
  arma::mat Delta = delta_constructor(p, J, K);
  
  //reset gb to 0 after each iteration
  gb_sub = zeros<vec>(p*J*K);
  Gb_sub = zeros<mat>(p*J*K, p*J*K);
  gi_list = zeros<mat>(N, J*K*p);
  
  n_k_cumsum = 0;
  
  for (int k=0; k<K; k++){
    m_j_cumsum = 0;
    for(int j=0; j<J; j++){
      //rows and columns of Y and X have to be sorted
      
      block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
      x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
      block_x = x_p.slice(0).as_col();
      for(int q=1; q<p; q++){
        block_x = join_rows(block_x, x_p.slice(q).as_col());
      }
      nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
      int n = nobs_jk.n_elem;
      
      int loc1 = 0;
      int loc2 = -1;
      for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs_jk[i] -1;
        arma::mat Xi = block_x.rows(loc1, loc2) ;
        arma::vec yi = block_y.subvec(loc1, loc2);

        arma::vec gi_sub = Xi.t() * yi;
        gi_list(n_k_cumsum+ i, span((k*J+j)*p, (k*J+j)*p+p-1)) = gi_sub.t();
        gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1)) += gi_sub;
        Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Xi;
      } 
      m_j_cumsum += m_j[j];
    }
    n_k_cumsum += n_k[k];
  }
  
  //Now put in the penalties
  arma::vec obj1_dev = gb_sub/N + rho*Delta.t()*gamma - Delta.t()*t;
  arma::mat obj2_dev = rho*Delta.t()*Delta + Gb_sub/N;
  
  beta_sub = inv(obj2_dev)*obj1_dev;
  
  return List::create(Named("beta_sub") = beta_sub,
                      Named("g_sub") = gb_sub,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub
  );
}

//[[Rcpp::export]]
List binomial_min(const arma::cube& X, const arma::mat& y, const List& nobs, 
                const arma::vec& beta_old, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t, const int& maxit, const double& tol){
  int niter2= 0;
  bool stop_flag2= FALSE;
  bool converged2= FALSE;
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::vec gb_sub;
  arma::mat Gb_sub;
  arma::mat Cb_sub;
  arma::mat gi_list;
  
  //initialization by the beta estimated from previous data
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_sub = beta_old;
  arma::vec beta_jk;
  
  arma::vec mu;
  arma::vec vu;
  arma::mat M0;
  arma::mat M1;
  arma::mat M2;
  
  arma::mat Delta = delta_constructor(p, J, K);
  arma::mat qif1dev_sub;
  arma::mat qif2dev_sub;
  
  while (!stop_flag2){
    niter2 += 1;
    //reset gb to 0 after each iteration
    gb_sub = zeros<vec>(p*J*K);
    Gb_sub = zeros<mat>(p*J*K, p*J*K);
    Cb_sub = zeros<mat>(p*J*K, p*J*K);
    gi_list = zeros<mat>(N, J*K*p);
    
    n_k_cumsum = 0;
    
    for (int k=0; k<K; k++){
      m_j_cumsum = 0;
      for(int j=0; j<J; j++){
        //rows and columns of Y and X have to be sorted
        
        block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
        x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
        block_x = x_p.slice(0).as_col();
        for(int q=1; q<p; q++){
          block_x = join_rows(block_x, x_p.slice(q).as_col());
        }
        
        beta_jk = beta_sub.subvec((k*J+j)*p, (k*J+j)*p+p-1);
        nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
        int n = nobs_jk.n_elem;
        
        // update gb_new with beta_new over iterations
        arma::vec eta2 = block_x * beta_jk;
        
        mu = exp(eta2)/( 1 + exp(eta2) ) ;
        vu = exp(eta2)/(pow( (1 + exp(eta2)), 2 )) ;
        
        int loc1 = 0;
        int loc2 = -1;
        for (int i = 0; i < n; i++){
          loc1 = loc2 + 1 ;
          loc2 = loc1 + nobs_jk[i] -1;
          arma::mat Xi = block_x.rows(loc1, loc2) ;
          arma::vec yi = block_y.subvec(loc1, loc2);
          arma::vec mui = mu.subvec(loc1, loc2);
          arma::vec vui = vu.subvec(loc1, loc2) ;
          
          arma::mat Ai_half = diagmat(sqrt(vui)) ;
          arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
          
          arma::vec gi_sub = Xi.t() * (yi - mui);
          gi_list(n_k_cumsum+ i, span((k*J+j)*p, (k*J+j)*p+p-1)) = gi_sub.t();
          gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1)) += gi_sub;
          Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Ai_half * Ai_half * Xi ;
        } 
        m_j_cumsum += m_j[j];
      }
      n_k_cumsum += n_k[k];
    }
    
    //Now put in the penalties
    arma::vec obj1_dev = gb_sub/N - rho * Delta.t()*(Delta*beta_sub - gamma + t/rho);
    arma::mat obj2_dev = Gb_sub/N + rho * Delta.t()*Delta;
    
    vec d_beta_sub = solve(obj2_dev, obj1_dev);
    double df_beta_sub = as_scalar(obj1_dev.t() * d_beta_sub);
    
    beta_sub += d_beta_sub;
    if(std::abs(df_beta_sub) < tol) {converged2 = TRUE; stop_flag2 = TRUE;}
    if (niter2 > maxit) {stop_flag2 = TRUE;}
  }

  return List::create(Named("beta_sub") = beta_sub,
                      Named("g_sub") = gb_sub,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub
  );
}

//[[Rcpp::export]]
List poisson_min(const arma::cube& X, const arma::mat& y, const List& nobs, 
                  const arma::vec& beta_old, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                  const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t, const int& maxit, const double& tol){
  int niter2= 0;
  bool stop_flag2= FALSE;
  bool converged2= FALSE;
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::vec gb_sub;
  arma::mat Gb_sub;
  arma::mat Cb_sub;
  arma::mat gi_list;
  
  //initialization by the beta estimated from previous data
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_sub = beta_old;
  arma::vec beta_jk;
  
  arma::vec mu;
  arma::vec vu;
  arma::mat M0;
  arma::mat M1;
  arma::mat M2;
  
  arma::mat Delta = delta_constructor(p, J, K);
  arma::mat qif1dev_sub;
  arma::mat qif2dev_sub;
  
  while (!stop_flag2){
    niter2 += 1;
    //reset gb to 0 after each iteration
    gb_sub = zeros<vec>(p*J*K);
    Gb_sub = zeros<mat>(p*J*K, p*J*K);
    Cb_sub = zeros<mat>(p*J*K, p*J*K);
    gi_list = zeros<mat>(N, J*K*p);
    
    n_k_cumsum = 0;
    
    for (int k=0; k<K; k++){
      m_j_cumsum = 0;
      for(int j=0; j<J; j++){
        //rows and columns of Y and X have to be sorted
        
        block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
        x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
        block_x = x_p.slice(0).as_col();
        for(int q=1; q<p; q++){
          block_x = join_rows(block_x, x_p.slice(q).as_col());
        }
        
        beta_jk = beta_sub.subvec((k*J+j)*p, (k*J+j)*p+p-1);
        nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
        int n = nobs_jk.n_elem;
        
        // update gb_new with beta_new over iterations
        arma::vec eta2 = block_x * beta_jk;
        
        mu = exp(eta2) ;
        vu = exp(eta2) ;
        
        int loc1 = 0;
        int loc2 = -1;
        for (int i = 0; i < n; i++){
          loc1 = loc2 + 1 ;
          loc2 = loc1 + nobs_jk[i] -1;
          arma::mat Xi = block_x.rows(loc1, loc2) ;
          arma::vec yi = block_y.subvec(loc1, loc2);
          arma::vec mui = mu.subvec(loc1, loc2);
          arma::vec vui = vu.subvec(loc1, loc2) ;
          
          arma::mat Ai_half = diagmat(sqrt(vui)) ;
          arma::mat Ai_inv_half = diagmat(1/sqrt(vui)) ;
          
          arma::vec gi_sub = Xi.t() * (yi - mui);
          gi_list(n_k_cumsum+ i, span((k*J+j)*p, (k*J+j)*p+p-1)) = gi_sub.t();
          gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1)) += gi_sub;
          Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Ai_half * Ai_half * Xi ;
        } 
        m_j_cumsum += m_j[j];
      }
      n_k_cumsum += n_k[k];
    }
    
    //Now put in the penalties
    arma::vec obj1_dev = gb_sub/N - rho * Delta.t()*(Delta*beta_sub - gamma + t/rho);
    arma::mat obj2_dev = Gb_sub/N + rho * Delta.t()*Delta;
    
    vec d_beta_sub = solve(obj2_dev, obj1_dev);
    double df_beta_sub = as_scalar(obj1_dev.t() * d_beta_sub);
    
    beta_sub += d_beta_sub;
    if(std::abs(df_beta_sub) < tol) {converged2 = TRUE; stop_flag2 = TRUE;}
    if (niter2 > maxit) {stop_flag2 = TRUE;}
  }
  
  return List::create(Named("beta_sub") = beta_sub,
                      Named("g_sub") = gb_sub,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub
  );
}

//[[Rcpp::export]]
List linear_eval(const arma::cube& X, const arma::mat& y, const List& nobs, 
                const arma::vec& beta, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t){
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_jk;
  arma::vec gi_list = zeros<vec>(N*M);
  arma::mat Gb_sub = zeros<mat>(p*J*K, p*J*K);
  
  arma::vec mu;
  arma::vec vu;
  
  n_k_cumsum = 0;
  double ll =0;
  
  for (int k=0; k<K; k++){
    m_j_cumsum = 0;
    for(int j=0; j<J; j++){
      //rows and columns of Y and X have to be sorted
      
      block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
      x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
      block_x = x_p.slice(0).as_col();
      for(int q=1; q<p; q++){
        block_x = join_rows(block_x, x_p.slice(q).as_col());
      }
      
      beta_jk = beta.subvec((k*J+j)*p, (k*J+j)*p+p-1);
      nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
      int n = nobs_jk.n_elem;
      
      // update gb_new with beta_new over iterations
      arma::vec eta2 = block_x * beta_jk;
      
      mu = eta2 ;
      
      int loc1 = 0;
      int loc2 = -1;
      for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs_jk[i] -1;
        arma::mat Xi = block_x.rows(loc1, loc2) ;
        arma::vec yi = block_y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        gi_list(span(loc1, loc2)) = yi - mui;
        Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Xi;
        
        ll += as_scalar((yi - mui).t()*(yi - mui));
      } 
      m_j_cumsum += m_j[j];
    }
    n_k_cumsum += n_k[k];
  }
  
  return List::create(Named("ll") = log(ll/N),
               Named("gi_list") = gi_list,
               Named("G_sub") = Gb_sub
  );
}

//[[Rcpp::export]]
List binomial_eval(const arma::cube& X, const arma::mat& y, const List& nobs, 
                 const arma::vec& beta, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                 const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t){
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_jk;
  arma::vec gi_list = zeros<vec>(N*M);
  arma::mat Gb_sub = zeros<mat>(p*J*K, p*J*K);
  
  arma::vec mu;
  arma::vec vu;
  
  n_k_cumsum = 0;
  double ll =0;
  
  for (int k=0; k<K; k++){
    m_j_cumsum = 0;
    for(int j=0; j<J; j++){
      //rows and columns of Y and X have to be sorted
      
      block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
      x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
      block_x = x_p.slice(0).as_col();
      for(int q=1; q<p; q++){
        block_x = join_rows(block_x, x_p.slice(q).as_col());
      }
      
      beta_jk = beta.subvec((k*J+j)*p, (k*J+j)*p+p-1);
      nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
      int n = nobs_jk.n_elem;
      
      // update gb_new with beta_new over iterations
      arma::vec eta2 = block_x * beta_jk;
      
      mu = exp(eta2)/( 1 + exp(eta2) ) ;
      
      int loc1 = 0;
      int loc2 = -1;
      for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs_jk[i] -1;
        arma::mat Xi = block_x.rows(loc1, loc2) ;
        arma::vec yi = block_y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        gi_list(span(loc1, loc2)) = yi - mui;
        Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Xi;
        
        for(int r=0; r<yi.size(); r++){
          if(yi(r)==1) {ll += yi(r)*log(mui(r));}
          if(yi(r)==0) {ll += (1-yi(r))*log(1-mui(r));}
        }
      } 
      m_j_cumsum += m_j[j];
    }
    n_k_cumsum += n_k[k];
  }
  
  return List::create(Named("ll") = ll/N,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub
  );
}

//[[Rcpp::export]]
List poisson_eval(const arma::cube& X, const arma::mat& y, const List& nobs, 
                   const arma::vec& beta, const int& J, const int& K, const int& M, const int& N, const int& p, const arma::vec& n_k, 
                   const arma::vec& m_j, const double& rho, const arma::vec& gamma, const arma::vec& t){
  
  int n_k_cumsum;
  int m_j_cumsum;
  
  arma::cube x_p;
  arma::vec block_y;
  arma::mat block_x;
  arma::vec nobs_jk;
  arma::vec beta_jk;
  arma::vec gi_list = zeros<vec>(N*M);
  arma::mat Gb_sub = zeros<mat>(p*J*K, p*J*K);
  
  arma::vec mu;
  arma::vec vu;
  
  n_k_cumsum = 0;
  double ll =0;
  
  for (int k=0; k<K; k++){
    m_j_cumsum = 0;
    for(int j=0; j<J; j++){
      //rows and columns of Y and X have to be sorted
      
      block_y = (y.submat(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1)).as_col();
      x_p = (X.tube(m_j_cumsum , n_k_cumsum, (m_j_cumsum+m_j[j])-1, (n_k_cumsum+n_k[k])-1));
      block_x = x_p.slice(0).as_col();
      for(int q=1; q<p; q++){
        block_x = join_rows(block_x, x_p.slice(q).as_col());
      }
      
      beta_jk = beta.subvec((k*J+j)*p, (k*J+j)*p+p-1);
      nobs_jk = Rcpp::as<arma::vec>( Rcpp::as<Rcpp::List>(nobs[k])[j] );
      int n = nobs_jk.n_elem;
      
      // update gb_new with beta_new over iterations
      arma::vec eta2 = block_x * beta_jk;
      
      mu = exp(eta2);
      
      int loc1 = 0;
      int loc2 = -1;
      for (int i = 0; i < n; i++){
        loc1 = loc2 + 1 ;
        loc2 = loc1 + nobs_jk[i] -1;
        arma::mat Xi = block_x.rows(loc1, loc2) ;
        arma::vec yi = block_y.subvec(loc1, loc2);
        arma::vec mui = mu.subvec(loc1, loc2);
        gi_list(span(loc1, loc2)) = yi - mui;
        Gb_sub(span((k*J+j)*p, (k*J+j)*p+p-1), span((k*J+j)*p, (k*J+j)*p+p-1)) += Xi.t() * Xi;
        
        for(int r=0; r<yi.size(); r++){
          ll += -mui(r) + log(mui(r))*yi(r);
        }
      } 
      m_j_cumsum += m_j[j];
    }
    n_k_cumsum += n_k[k];
  }
  
  return List::create(Named("ll") = ll/N,
                      Named("gi_list") = gi_list,
                      Named("G_sub") = Gb_sub
  );
}

// [[Rcpp::export]]
List indep_doADMM_MCP(const arma::cube& X, const arma::mat& y, const List& nobs, const String& family,
                       const arma::vec& init_betas, const arma::vec& init_gamma, const arma::vec& init_t, const int& J, const int& K, const int& M, 
                       const int& N, const int& p, const arma::vec& n_k, const arma::vec& m_j, const double& lambda, const double& rho, 
                       const double& delta, const double& tol_1, const double& tol_2, const int& maxit_1, const int& maxit_2){
  std::ostringstream __my_cerr;
  arma::set_cerr_stream(__my_cerr);
  bool convergence = FALSE;
  int r = 0;
  arma::vec old_beta;
  arma::vec beta = init_betas;
  arma::vec gamma = init_gamma;
  arma::vec gamma_unpen;
  arma::vec t = init_t;
  arma::mat Delta = delta_constructor(p, J, K);
  arma::mat Delta_beta;
  double gamma_diff;
  double first_diff;
  double residual;
  
  while((!convergence) & (r <= maxit_1)){
    old_beta = beta;
    List minimized;
    // First update beta
    if(family=="gaussian") minimized = linear_min(X, y, nobs, old_beta, J, K, M, N, p, n_k, m_j, rho, gamma, t);
    if(family=="binomial") minimized = binomial_min(X, y, nobs, old_beta, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit_2, tol_2);
    if(family=="poisson") minimized = poisson_min(X, y, nobs, old_beta, J, K, M, N, p, n_k, m_j, rho, gamma, t, maxit_2, tol_2);
    arma::vec minimized_beta = minimized[0];
    beta = minimized_beta;
    Delta_beta = Delta * beta;
    
    // Second update gamma
    if(r==0) { gamma_diff = max(abs(Delta_beta-gamma)); }
    gamma_unpen = Delta_beta + t/rho;
    gamma = zeros<vec>((K*(K-1)/2*J*J+K*J*(J-1)/2)*p);
    for(int t=0; t<((K*(K-1)/2*J*J+K*J*(J-1)/2)); t++){
      double sum_abs = sum(abs( gamma_unpen( span(t*p, t*p+p-1) ) ));
      if(sum_abs <= (delta*lambda)){
        for(int q=0; q<p; q++){
          double shrink = std::abs(gamma_unpen[t*p+q])-lambda/rho;
          gamma[t*p+q] = ( sign(gamma_unpen[t*p+q]) * std::max(shrink, 0.0) ) /  ( 1 - 1/ ( delta*rho ) );
        }
      } else {
        gamma( span(t*p, t*p+p-1) ) = gamma_unpen( span(t*p, t*p+p-1) );
      }
    }
    
    // Third update t
    t = t + rho*(Delta_beta - gamma);
    
    residual = sum(pow(Delta_beta - gamma, 2.0));
    if(residual < tol_1){
      convergence = TRUE;
    } else {
      r = r+1;
    }
  }
  
  std::string text(__my_cerr.str());
  
  return List::create(Named("convergence") = convergence,
                      Named("beta") = beta,
                      Named("gamma") = gamma,
                      Named("t") = t,
                      Named("lambda") = lambda,
                      Named("rho") = rho,
                      Named("r") = r,
                      Named("first_diff") = first_diff,
                      Named("residual") = residual,
                      Named("message") = text
  );
}

// [[Rcpp::export]]
List try_indep_doADMM_MCP(const arma::cube& X, const arma::mat& y, const List& nobs, const String& family,
                           const arma::vec& init_betas, const arma::vec& init_gamma, const arma::vec& init_t, const int& J, const int& K, 
                           const int& M, const int& N, const int& p, const arma::vec& n_k, const arma::vec& m_j, const double& lambda, 
                           const double& rho, const double& delta, const double& tol_1, const double& tol_2, const int& maxit_1, const int& maxit_2){
  List ADMM;
  try{
    ADMM = indep_doADMM_MCP(X, y, nobs, family, init_betas, init_gamma, init_t, J, K, M, N, p, n_k, m_j, lambda, rho, delta, 
                             tol_1, tol_2, maxit_1, maxit_2);
  } 
  catch (const std::exception &__ex) {
    return(List::create(Named("message") = std::string("error")));
  }
  catch (...) {
    return(List::create(Named("message") = std::string("error")));
  }
  return(ADMM);
}
