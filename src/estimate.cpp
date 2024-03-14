#include <RcppArmadillo.h>
#include "auxiliary.h"

// [[Rcpp::depends(RcppArmadillo)]]

arma::vec lt_to_vec(const arma::mat& L){
  int r = L.n_rows;
  int k = 0;
  arma::vec v = arma::zeros<arma::vec>(r*(r+1)/2);
  for(int i = 1; i <= r; i++) {
    for(int j = 0; j < i; j++) {
      v(k) = L(i-1, j);
      k += 1;
    }
  }
  return v;
}



arma::mat vec_to_lt(const arma::vec& v, const int& r){
  int k = 0;
  arma::mat L = arma::zeros<arma::mat>(r, r);
  for(int i = 1; i <= r; i++) {
    for(int j = 0; j < i; j++) {
      L(i-1, j) = v(k);
      k += 1;
    }
  }
  return L;
}


// [[Rcpp::export]]
arma::mat minimizeQ(const arma::mat& Lt, const arma::mat& Lr, const Rcpp::List& auxmat, const double lambda,
                     const arma::vec& weight, const arma::mat& L0, const int& maxIt) {
  //Rcpp::List outList;
  Rcpp::List Hessian;
  arma::mat B = Rcpp::as<arma::mat>(auxmat["B"]);
  arma::mat U = Rcpp::as<arma::mat>(auxmat["U"]);
  arma::mat V = Rcpp::as<arma::mat>(auxmat["V"]);
  arma::mat W = Rcpp::as<arma::mat>(auxmat["W"]);
  int r = U.n_rows;
  arma::mat H = arma::zeros<arma::mat>(r*(r+1)/2, r*(r+1)/2);
  arma::mat G = arma::zeros<arma::mat>(r, r);
  arma::vec g = arma::zeros<arma::vec>(r*(r+1)/2);
  arma::vec eta = arma::zeros<arma::vec>(r*(r+1)/2);
  arma::mat Eta = arma::zeros<arma::mat>(r, r);
  
  arma::mat S = arma::zeros<arma::mat>(r, r);
  arma::mat L = L0;
  for(int k = 0; k < maxIt; k++) {
    Hessian = RhessianQ(Lt, Lr, S, B, U, V, W, lambda, weight, L);
    H = Rcpp::as<arma::mat>(Hessian["H"]);
    G = Rcpp::as<arma::mat>(Hessian["G"]);
    g = lt_to_vec(G);
    eta = solve(H, -g);
    Eta = vec_to_lt(eta, r);
    L = expM(L, Eta);
  }
  //outList["L"] = L;
  //outList["Eta"] = Eta;
  //outList["H"] = H;
  return L;
}


// [[Rcpp::export]]
arma::mat minimizeL(const arma::mat& Lt, const arma::mat& Lr, const Rcpp::List& auxmat, const double& lambda,
                     const arma::vec& weight) {
  //Rcpp::List outList;
  int n = Lt.n_rows;
  int m = Lt.n_cols;
  arma::mat B = Rcpp::as<arma::mat>(auxmat["B"]);
  arma::mat U = Rcpp::as<arma::mat>(auxmat["U"]);
  arma::mat V = Rcpp::as<arma::mat>(auxmat["V"]);
  arma::mat W = Rcpp::as<arma::mat>(auxmat["W"]);
  int r = U.n_rows;
  arma::mat E = arma::zeros<arma::vec>(r*r);
  arma::mat E1 = arma::zeros<arma::vec>(r*(r+1)/2);
  arma::mat D = arma::zeros<arma::mat>(r*r, r*r);
  arma::mat D0 = arma::zeros<arma::mat>(r*(r+1)/2, r*r);
  arma::mat D1 = arma::zeros<arma::mat>(r*(r+1)/2, r*(r+1)/2);
  for(int i = 0; i < n; i++) {
    const arma::mat Ri = Lr.cols(i*m,(i+1)*m-1);
    const arma::mat Bi = B.cols(i*r,(i+1)*r-1);
    // Term I
    E += - weight(i) * arma::kron(Bi.t(), Bi.t()) * Ri.as_col();
    // Term II
    D += weight(i) * arma::kron(Bi.t() * Bi, Bi.t() * Bi);
    // Term III and Term IV
    for(int j = 0; j < m; j++) {
      arma::rowvec bij = Bi.row(j);
      D += - weight(i) * arma::kron(bij.t() * bij, bij.t() * bij);
      E += weight(i) * Ri(j,j) * arma::kron(bij.t(), bij.t());
    }
  }
  // penalty term
  D += lambda * (arma::kron(W, U) + arma::kron(V, V));
  D = 2*D;
  E = 2*E;
  // end of compution of DE
  int k = 0;
  for(int ii = 1; ii <= r; ii++) {
    for(int jj = 0; jj < ii; jj++) {
      int k1 = r * jj + (ii-1);
      int k2 = r * (ii-1) + jj;
      if(ii == (jj+1)) {
        D0.row(k) = D.row(k1);
      }else{
        D0.row(k) = D.row(k1) + D.row(k2);
      }
      k += 1;
    } 
  }
  k = 0;
  for(int ii = 1; ii <= r; ii++) {
    for(int jj = 0; jj < ii; jj++) {
      int k1 = r * jj + (ii-1);
      int k2 = r * (ii-1) + jj;
      if(ii == (jj+1)){
        D1.col(k) = D0.col(k1);
        E1(k) = E(k1);
      }else{
        D1.col(k) = D0.col(k1) + D0.col(k2);
        E1(k) = E(k1) + E(k2);
      }
      k += 1;
    }
  }
  arma::vec c0 = solve(D1, -E1);
  arma::mat C0 = vec_to_lt(c0, r);
  C0 = C0 + C0.t();
  C0 = C0 - diagmat(C0)/2;
  //outList["C0"] = C0;
  //outList["D"] = D1;
  //outList["E"] = E1;
  return C0;
}


// [[Rcpp::export]]
arma::mat C_hat(const arma::mat& Lt, const arma::mat& Lr, const Rcpp::List& auxmat, const double lambda,
                const arma::vec& weight, const int& maxIt) {
  arma::mat U = Rcpp::as<arma::mat>(auxmat["U"]);
  int r = U.n_rows;
  arma::mat Chat = minimizeL(Lt, Lr, auxmat, lambda, weight);
  arma::vec eigval = arma::eig_sym(Chat);
  double rho = min(eigval);
  if(rho < 0) {
    arma::mat tmp = Chat - (rho * 1.1) * diagmat(arma::ones<arma::vec>(r));
    tmp = chol(tmp).t();
    arma::mat Lmin = minimizeQ(Lt, Lr, auxmat, lambda, weight, tmp, maxIt);
    Chat = Lmin * Lmin.t();
  }
  return Chat;
}





// [[Rcpp::export]]
Rcpp::List cov_basis(const arma::mat& Lt, const arma::mat& Ly, const arma::mat& Lr, const int& r, const double& lambda,
                     const double& ext, const int& maxIt){
  Rcpp::List outList;
  int n = Lt.n_rows;
  int m = Lt.n_cols;
  arma::vec weight = 1.0/(n*(m*(m-1))) * arma::ones<arma::vec>(n);
  arma::vec delta_vec = arma::zeros<arma::vec>(n);
  arma::vec domain = arma::zeros<arma::vec>(2);
  arma::mat Lp = arma::zeros<arma::mat>(m, m*n);
  arma::mat Ls = arma::zeros<arma::mat>(n, m);
  double min_domain = Lt.min();
  double max_domain = Lt.max();
  domain(0) = min_domain;
  domain(1) = max_domain;
  if (ext > 0) {
    Ls = ext + (1-ext*2) * (Lt - min_domain) / (max_domain - min_domain);
  }else {
    Ls = (Lt - min_domain) / (max_domain - min_domain);
  }
  Rcpp::List auxmat = auxiliary_mat(r, Ls);
  arma::mat C = C_hat(Ls, Lr, auxmat, lambda, weight, maxIt);
  for(int i = 0; i < n; ++i) {
    arma::vec ls = Ls.row(i).t();
    arma::mat basis_mat = evaluate_basis(r, domain, ls);
    Lp.cols(i*m, (i+1)*m-1) = basis_mat * C * basis_mat.t();
  }
  double error = arma::norm(Lr - Lp, "fro");
  outList["C"] = C;
  outList["error"] = error;
  //outList["Lr"] = Lr;
  //outList["Lp"] = Lp;
  //outList["Ls"] = Ls;
  //outList["domain"] = domain;
  return outList;
}


// // [[Rcpp::export]]
// Rcpp::List error_gof(const int& s, const int& e, const Rcpp::List& Lt, const Rcpp::List& Ly, const int& r,
//                  const double& lambda, const double& ext, const int& maxIt){
//   Rcpp::List outList;
//   Rcpp::List Lt_new(e-s+1);
//   Rcpp::List Ly_new(e-s+1);
//   for(int i = 0; i < e-s+1; ++i) {
//     Lt_new[i] = Lt[i+s-1];
//     Ly_new[i] = Ly[i+s-1];
//   }
//   Rcpp::List cov_obj_new = cov_basis(Lt_new, Ly_new, r, lambda, ext, maxIt);
//   Rcpp::List Lr_new = cov_obj_new["Lr"];
//   Rcpp::List Lp_new = cov_obj_new["Lp"];
//   arma::mat C_new = cov_obj_new["C"];
//   double error = 0.0;
//   int sub_n = Lr_new.size();
//   for(int i = 0; i < sub_n; ++i) {
//     error += arma::norm(Rcpp::as<arma::mat>(Lr_new[i]) - Rcpp::as<arma::mat>(Lp_new[i]), "fro");
//   }
//   outList["C"] = C_new;
//   outList["error"] = error;
//   return outList;
// }

