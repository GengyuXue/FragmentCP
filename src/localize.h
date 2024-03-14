#ifndef LOCALIZE_H
#define LOCALIZE_H

Rcpp::List DP_fragment(const arma::mat& Lt, const arma::mat& Ly, const arma::mat& Lr, const int& r, const double& lambda, const double& xi,
                       const double& ext, const int& maxIt, const int& Delta);
                       
#endif
