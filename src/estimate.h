#ifndef ESTIMATE_H
#define ESTIMATE_H

Rcpp::List minimizeQ(const arma::mat& Lt, const arma::mat& Lr, const Rcpp::List& auxmat, const double lambda,
                     const arma::vec& weight, const arma::mat& L0, const int& maxIt);
Rcpp::List minimizeL(const arma::mat& Lt, const arma::mat& Lr, const Rcpp::List& auxmat, const double& lambda,
                     const arma::vec& weight);
arma::mat C_hat(const arma::mat& Lt, const arma::mat& Lr, const Rcpp::List& auxmat, const double lambda,
                const arma::vec& weight, const int& maxIt);
Rcpp::List cov_basis(const arma::mat& Lt, const arma::mat& Ly, const arma::mat& Lr, const int& r, const double& lambda,
                     const double& ext, const int& maxIt);

#endif
