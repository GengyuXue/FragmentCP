#ifndef AUXILIARY_H
#define AUXILIARY_H

arma::mat expM(const arma::mat& X, const arma::mat& S);
arma::mat RgradQ(const arma::mat& Lt,
                 const arma::mat& Lr,
                 const arma::mat& L,
                 const arma::mat& B,
                 const arma::mat& U,
                 const arma::mat& V,
                 const arma::mat& W,
                 const double lam,
                 const arma::vec& weight);
Rcpp::List RhessianQ(const arma::mat& Lt, const arma::mat& Lr, const arma::mat& S, const arma::mat& B,
                     const arma::mat& U, const arma::mat& V, const arma::mat& W, const double lam, const arma::vec& weight, const arma::mat& X);
arma::mat deriv_fourier(const int& K, const arma::vec& domain, const arma::vec& grid, const int& order);
arma::mat evaluate_basis(const int& K, const arma::vec& domain, const arma::vec& grid);
Rcpp::List auxiliary_mat(const int& r, const arma::mat& Lt);

#endif
