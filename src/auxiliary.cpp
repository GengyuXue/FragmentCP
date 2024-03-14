#include <RcppArmadillo.h>

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat expM(const arma::mat& X, const arma::mat& S) {
  arma::mat R = X;
  R += S;
  R.diag().zeros();
  arma::vec D = X.diag();
  arma::vec Z = S.diag() / D;
  Z = exp(Z);
  D %= Z;
  R.diag() = D;
  return R;
}


// [[Rcpp::export]]
arma::mat RgradQ(const arma::mat& Lt,
                 const arma::mat& Lr,
                 const arma::mat& L,
                 const arma::mat& B,
                 const arma::mat& U,
                 const arma::mat& V,
                 const arma::mat& W,
                 const double lam,
                 const arma::vec& weight) {
  int n = Lt.n_rows;
  int p = L.n_rows;
  int m = Lt.n_cols;
  const arma::mat C = L * L.t();
  arma::mat H = arma::zeros<arma::mat>(p, p);
  for(int i = 0; i < n; i++) {
    const arma::mat Ri = Lr.cols(i*m,(i+1)*m-1);
    const arma::mat Bi = B.cols(i*p,(i+1)*p-1);
    arma::mat tmp = -4 * (Bi.t() * Ri * Bi * L);
    tmp += 4 * (Bi.t() * Bi * C * Bi.t() * Bi * L);
    for(int j = 0; j < m; j++) {
      const arma::rowvec bij = Bi.row(j); // Row vector
      tmp -= 4 * (bij.t() * bij * C * bij.t() * bij * L);
      tmp += 4 * (Ri(j,j) * bij.t() * bij * L);
    }
    tmp *= weight(i);
    H += tmp;
  }
  H += 2 * lam * (U * C * W * L + W * C * U * L) + 4 * lam * V * C * V * L;
  return H;
}



// [[Rcpp::export]]
Rcpp::List RhessianQ(const arma::mat& Lt, const arma::mat& Lr, const arma::mat& S, const arma::mat& B,
                     const arma::mat& U, const arma::mat& V, const arma::mat& W, const double lam, const arma::vec& weight, const arma::mat& X) {
  int n = Lt.n_rows;
  int m = Lt.n_cols;
  arma::mat L = expM(X, S);
  int p = L.n_rows;
  arma::mat C = L * L.t();
  int d = p * (p + 1) / 2;
  arma::mat H = arma::zeros<arma::mat>(d, d);
  arma::mat G = RgradQ(Lt, Lr, L, B, U, V, W, lam, weight);
  
  double tmp1, tmp2, tmp3;
  arma::mat Eim;
  int j, k;
  double coef;
  
  Rcpp::List outList;
  
  for(int i = 0; i < n; i++) {
    const arma::mat Ri = Lr.cols(i*m,(i+1)*m-1);
    const arma::mat Bi = B.cols(i*p,(i+1)*p-1);
    arma::mat E1 = Bi.t() * Ri * Bi;
    arma::mat E2 = Bi.t() * Bi;
    for(int j1 = 0; j1 < p; j1++) {
      for(int j2 = 0; j2 <= j1; j2++) {
        j = j1 * (j1 + 1) / 2 + j2;
        for(int k1 = 0; k1 < p; k1++) {
          for(int k2 = 0; k2 <= k1; k2++) {
            k = k1 * (k1 + 1) / 2 + k2;
            coef = 1;
            if(j1 == j2 && k1 == k2) coef = exp(S(j1, j2) / X(j1, j2)) * exp(S(k1, k2) / X(k1, k2));
            else if(j1 == j2) coef = exp(S(j1, j2) / X(j1, j2));
            else if(k1 == k2) coef = exp(S(k1, k2) / X(k1, k2));
            coef *= weight(i);
            
            tmp1 = arma::as_scalar(L.col(k2).t() * E2 * L.col(j2));
            tmp2 = arma::as_scalar(E2.row(j1) * L.col(k2)) * arma::as_scalar(E2.row(k1) * L.col(j2));
            
            if(j2 == k2) {
              // Term I
              H(j, k) -= coef * E1(j1, k1);
              
              // Term II
              tmp3 = arma::as_scalar(E2.row(j1) * C * E2.col(k1));
              H(j, k) += coef * (E2(j1, k1) * tmp1 + tmp2 + tmp3);
            } else {
              // Term I is zero
              // Term II
              H(j, k) += coef * (E2(j1, k1) * tmp1 + tmp2);
            }
            
            for(int mm = 0; mm < m; mm++) {
              arma::mat Eim = Bi.row(mm).t() * Bi.row(mm);
              double tmp1 = arma::as_scalar(L.col(k2).t() * Eim * L.col(j2)); 
              double tmp2 = arma::as_scalar(Eim.row(j1) * L.col(k2)) * arma::as_scalar(Eim.row(k1) * L.col(j2)); 
              
              if(j2 == k2) {
                // Term III
                double tmp3 = arma::as_scalar(Eim.row(j1) * C * Eim.col(k1)); // Dot product and multiplication
                H(j,k) -= coef * (Eim(j1,k1) * tmp1 + tmp2 + tmp3);
                
                // Term IV
                H(j,k) += coef * Ri(mm,mm) * Eim(j1,k1);
              } else {
                // Term III
                H(j,k) -= coef * (Eim(j1,k1) * tmp1 + tmp2);
                // Term IV is zero
              }
            }
          }
        }
      }
    }
  }
  H = 4 * H;
  for(int j = 0; j < p; j++) {
    int jj = (j+2)*(j+1)/2 - 1;
    H(jj, jj) += G(j, j) * std::exp(S(j, j) / X(j, j)) / X(j, j);
  }
  for(int j1 = 0; j1 < p; j1++) {
    for(int j2 = 0; j2 <= j1; j2++) {
      int j = j1*(j1+1)/2 + j2;
      for(int k1 = 0; k1 < p; k1++) {
        for(int k2 = 0; k2 <= k1; k2++) {
          int k = k1*(k1+1)/2 + k2;
          if(j1==j2) coef = std::exp(S(j1,j2)/X(j1,j2));
          else coef = 1;
          if(j2 == k2) {
            //Term I
            tmp1 = arma::as_scalar(L.col(k2).t() * W * L.col(j2));
            tmp1 = U(j1, k1) * tmp1;
            tmp2 = arma::as_scalar(U.row(j1) * L.col(k2)) * arma::as_scalar(W.row(k1) * L.col(j2));
            tmp3 = arma::as_scalar(U.row(j1) * C * W.col(k1));
            H(j, k) += 2 * lam * coef * (tmp1 + tmp2 + tmp3);
            //Term II
            tmp1 = arma::as_scalar(L.col(k2).t() * U * L.col(j2));
            tmp1 = W(j1, k1) * tmp1;
            tmp2 = arma::as_scalar(W.row(j1) * L.col(k2)) * arma::as_scalar(U.row(k1) * L.col(j2));
            tmp3 = arma::as_scalar(W.row(j1) * C * U.col(k1));
            H(j, k) += 2 * lam * coef * (tmp1 + tmp2 + tmp3);
            //Term III
            tmp1 = arma::as_scalar(L.col(k2).t() * V * L.col(j2));
            tmp1 = V(j1, k1) * tmp1;
            tmp2 = arma::as_scalar(V.row(j1) * L.col(k2)) * arma::as_scalar(V.row(k1) * L.col(j2));
            tmp3 = arma::as_scalar(V.row(j1) * C * V.col(k1));
            H(j, k) += 4 * lam * coef * (tmp1 + tmp2 + tmp3);
          } else {
            //Term I
            tmp1 = arma::as_scalar(L.col(k2).t() * W * L.col(j2));
            tmp1 = U(j1, k1) * tmp1;
            tmp2 = arma::as_scalar(U.row(j1) * L.col(k2)) * arma::as_scalar(W.row(k1) * L.col(j2));
            H(j, k) += 2 * lam * coef * (tmp1 + tmp2);
            //Term II
            tmp1 = arma::as_scalar(L.col(k2).t() * U * L.col(j2));
            tmp1 = W(j1,k1) * tmp1;
            tmp2 = arma::as_scalar(W.row(j1) * L.col(k2)) * arma::as_scalar(U.row(k1) * L.col(j2));
            H(j, k) += 2 * lam * coef * (tmp1 + tmp2);
            //Term III
            tmp1 = arma::as_scalar(L.col(k2).t() * V * L.col(j2));
            tmp1 = V(j1, k1) * tmp1;
            tmp2 = arma::as_scalar(V.row(j1) * L.col(k2)) * arma::as_scalar(V.row(k1) * L.col(j2));
            H(j, k) += 4 * lam * coef * (tmp1 + tmp2);
          }
        }
      }
    }
  }
  outList["G"] = G;
  outList["H"] = H;
  return outList;
}


// [[Rcpp::export]]
arma::mat deriv_fourier(const int& K, const arma::vec& domain, const arma::vec& grid, const int& order) {
  int m = grid.size();
  arma::mat V = arma::zeros<arma::mat>(m, K);
  double s;
  double L = domain(1) - domain(0);
  double a = domain(0);
  for(int k = 0; k < K; k++) {
    if(k == 0){
      if(order == 0) {
        V.col(k) = arma::ones<arma::vec>(m);
      }
    }else if((k+1) % 2 == 0) {
      s = 1;
      if((order % 4 == 1) || (order % 4 == 2)) {
        s = -1;
      } 
      if(order % 2 == 1) {
        V.col(k) = s * sqrt(2/L) * pow(((k+1)*arma::datum::pi/L), order) * sin((k+1)*arma::datum::pi*(grid-a)/L);
      }else {
        V.col(k) = s * sqrt(2/L) * pow(((k+1)*arma::datum::pi/L), order) * cos((k+1)*arma::datum::pi*(grid-a)/L);
      }
    }else{
      s = 1;
      if(order % 4 == 2 || order % 4 == 3) {
        s = -1;
      }
      if(order % 2 == 1) {
        V.col(k) = s * sqrt(2/L) * pow((k*arma::datum::pi/L), order) * cos(k*arma::datum::pi*(grid-a)/L);
      }else {
        V.col(k) = s * sqrt(2/L) * pow((k*arma::datum::pi/L), order) * sin(k*arma::datum::pi*(grid-a)/L);
      }
    }
  }
  return V;
}


// [[Rcpp::export]]
arma::mat evaluate_basis(const int& K, const arma::vec& domain, const arma::vec& grid) {
  int m = grid.size();
  double a = domain(0);
  double b = domain(1);
  double x = 0;
  double y = 1;
  double alpha = (b*x-a*y)/(x-y);
  double beta = (a-b)/(x-y);
  arma::vec pts = (grid - alpha) / beta;
  arma::mat res = arma::ones<arma::mat>(m, K);
  for(int k = 1; k < K; k++) {
    if ((k+1) % 2 == 0) {
      res.col(k) = sqrt(2) * cos((k+1) * arma::datum::pi * pts);
    }else {
      res.col(k) = sqrt(2) * sin(k * arma::datum::pi * pts);
    }
  }
  res = res/sqrt(beta);
  return res;
}


// [[Rcpp::export]]
Rcpp::List auxiliary_mat(const int& r, const arma::mat& Lt){
  Rcpp::List outList;
  int n = Lt.n_rows;
  int m = Lt.n_cols;
  arma::vec domain = {0,1};
  arma::vec lt;
  arma::mat B = arma::zeros<arma::mat>(m, r*n);
  for(int i = 0; i < n; i++) {
    lt = Lt.row(i).t();
    B.cols(i*r, (i+1)*r-1) = evaluate_basis(r, domain, lt);
  }
  
  int length = 1000;
  arma::vec pts = arma::linspace(0.5/length, 1-0.5/length, length);
  arma::mat W = deriv_fourier(r, domain, pts, 2);
  W = W.t() * W / length;
  arma::mat V = deriv_fourier(r, domain, pts, 1);
  V = V.t() * V / length;
  arma::mat U = deriv_fourier(r, domain, pts, 0);
  U = U.t() * U / length;
  outList["B"] = B;
  outList["W"] = W;
  outList["V"] = V;
  outList["U"] = U;
  return outList;
}
