#include <RcppArmadillo.h>
#include "estimate.h"

// [[Rcpp::export]]
Rcpp::List DP_fragment(const arma::mat& Lt, const arma::mat& Ly, const arma::mat& Lr, const int& r, const double& lambda, const double& xi,
                       const double& ext, const int& maxIt, const int& Delta){
  Rcpp::List outList;
  int n = Lt.n_rows;
  int m = Lt.n_cols;
  //Rcpp::List C_temp1(n);
  //Rcpp::List C_temp2(n);
  //Rcpp::List C_final(n);
  double b = 0.0;
  arma::vec bestvalue = arma::zeros<arma::vec>(n+1);
  arma::vec partition = arma::zeros<arma::vec>(n+1);
  bestvalue(0) = -xi;
  for(int i = 1; i < n+1; ++i){
    Rcpp::Rcout << "i is" << std::endl << i << std::endl;
    bestvalue(i) = R_PosInf;
    for(int l = 1; l < i+1; ++l){
      //Rcpp::Rcout << "l is" << std::endl << l << std::endl;
      if ((i-l+1) > Delta){
        arma::mat Lt_new = Lt.rows((l-1), (i-1));
        arma::mat Ly_new = Ly.rows((l-1), (i-1));
        arma::mat Lr_new = Lr.cols((l-1)*m, i*m-1);
        Rcpp::List fit = cov_basis(Lt_new, Ly_new, Lr_new, r, lambda, ext, maxIt);
        //Rcpp::List fit = error_gof(l, i, Lt, Ly, r, lambda, ext, maxIt);
        //C_temp1[l-1] = Rcpp::as<arma::mat>(fit["C"]);
        //b = bestvalue(l-1) + xi + Rcpp::as<double>(fit["error"]);
        b = bestvalue(l-1) + xi + Rcpp::as<double>(fit["error"]);
      }else{
        b = bestvalue(l-1) + xi; 
      }
      if (b < bestvalue(i)){
        bestvalue(i) = b;
        partition(i) = l-1;
        //C_temp2[i-1] = C_temp1[l-1];
        //Rcpp::Rcout << partition(i) << std::endl;
      }
    }
  }
  int R = n;
  int L = partition(R);
  while(R > 0){
    //for(int t = L; t < R; ++t){
    //  C_final[t] = C_temp2[R-1];
    //}
    R = L;
    L = partition(R);
  }
  //outList["C_list"] = C_final;
  outList["partition"] = partition.subvec(1,n);
  outList["bestvalue"] = bestvalue.subvec(1,n);
  return outList;
}