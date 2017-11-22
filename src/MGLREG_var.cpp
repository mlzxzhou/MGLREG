# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
// [[Rcpp::export]]
double MGLREG_var(arma::mat cov, arma::mat prob_real,arma::mat prob_sam, arma::rowvec gamma1, arma::rowvec gamma2,
                 arma::mat trait, arma::vec err, int dim_cov, int num_sub) {
  arma::mat gra = arma::zeros(5*dim_cov,5*dim_cov);
  arma::mat var_gee = arma::zeros(5*dim_cov,5*dim_cov);

  for (int i = 0; i < num_sub; i++) {
    arma::mat Xi = cov.row(i);

    double p1_i = prob_real(i,1);
    double p2_i = prob_real(i,2);

    double pi1 = prob_sam(i,1);
    double pi2 = prob_sam(i,2);

    arma::mat beta_dev = Xi.t();
    arma::mat gamma1_dev = (trait(i,1)-p1_i)*beta_dev;
    arma::mat gamma2_dev = (trait(i,2)-p2_i)*beta_dev;
    arma::mat Ai = arma::join_cols(arma::join_cols(beta_dev,gamma1_dev),gamma2_dev);

    arma::mat psi1_dev = beta_dev * (gamma1 * beta_dev * p1_i*(1-p1_i)-gamma2 * beta_dev*(p1_i*p2_i));
    arma::mat psi2_dev = beta_dev * (-gamma1 * beta_dev * p1_i*p2_i+gamma2 * beta_dev*(p2_i*(1-p2_i)));
    arma::mat Bi = arma::join_cols(psi1_dev,psi2_dev);

    arma::mat Di_11(dim_cov,dim_cov);
    Di_11.fill(pi1*(1-pi1));
    arma::mat Di_12(dim_cov,dim_cov);
    Di_12.fill(-pi1*pi2);
    arma::mat Di_21(dim_cov,dim_cov);
    Di_21.fill(-pi1*pi2);
    arma::mat Di_22(dim_cov,dim_cov);
    Di_22.fill(pi2*(1-pi2));
    arma::mat Di = arma::join_cols(arma::join_rows(Di_11 % (beta_dev * beta_dev.t()),Di_12 % (beta_dev * beta_dev.t())),
                                   arma::join_rows(Di_21 % (beta_dev * beta_dev.t()),Di_22 % (beta_dev * beta_dev.t())));

    arma::mat Ci1_11(dim_cov,dim_cov);
    Ci1_11.fill(-p1_i*(1-p1_i));
    arma::mat Ci1_12(dim_cov,dim_cov);
    Ci1_12.fill(p1_i*p2_i);
    arma::mat Ci1_21(dim_cov,dim_cov);
    Ci1_21.fill(p1_i*p2_i);
    arma::mat Ci1_22(dim_cov,dim_cov);
    Ci1_22.fill(-p2_i*(1-p2_i));
    arma::mat Ci_1 = arma::join_cols(arma::join_rows(Ci1_11 % (beta_dev * beta_dev.t()),Ci1_12 % (beta_dev * beta_dev.t())),
                                     arma::join_rows(Ci1_21 % (beta_dev * beta_dev.t()),Ci1_22 % (beta_dev * beta_dev.t())));
    arma::mat Ci = arma::join_cols(arma::zeros(beta_dev.n_rows,Di.n_cols),Ci_1);
    double err_1 = err(i);

    arma::mat dev_11 = -Ai * Ai.t();
    arma::mat temp(Ci.n_rows,Ci.n_cols);
    temp.fill(err_1);
    arma::mat dev_12 = -Ai * Bi.t() + Ci % temp;
    arma::mat dev_21 = arma::zeros(dev_12.n_cols,dev_12.n_rows);
    arma::mat dev_22 = Di;

    gra = gra + arma::join_cols(arma::join_rows(dev_11,dev_12),arma::join_rows(dev_21,dev_22));

    arma::mat U_i = Ai * err_1;
    arma::mat psi1_dev_1 = beta_dev * (trait(i,1)-pi1);
    arma::mat psi2_dev_1 = beta_dev * (trait(i,2)-pi2);
    arma::mat L_i = arma::join_cols(psi1_dev_1,psi2_dev_1);
    arma::mat R_i = arma::join_cols(U_i,L_i);
    var_gee = var_gee + R_i * R_i.t();
  }
  arma::mat var_all = inv(gra) * var_gee * (inv(gra).t());
  double var_beta = var_all(1,1);
  return(var_beta);
}

// [[Rcpp::export]]
double MGLREG_var2(arma::mat X, arma::mat X_1, arma::mat prob_real,arma::mat prob_sam, arma::rowvec gamma1, arma::rowvec gamma2, arma::mat trait, arma::vec err, int dim_cov, int dim_cov_1, int num_sub) {
  arma::mat gra = arma::zeros(dim_cov*3+dim_cov_1*2,dim_cov*3+dim_cov_1*2);
  arma::mat var_gee = arma::zeros(dim_cov*3+dim_cov_1*2,dim_cov*3+dim_cov_1*2);

  for (int i = 0; i < num_sub; i++) {
    arma::mat Xi = X.row(i);
    arma::mat Xi_1 = X_1.row(i);

    double p1_i = prob_real(i,1);
    double p2_i = prob_real(i,2);

    double pi1 = prob_sam(i,1);
    double pi2 = prob_sam(i,2);

    arma::mat beta_dev = Xi.t();
    arma::mat gamma1_dev = (trait(i,1)-p1_i)*beta_dev;
    arma::mat gamma2_dev = (trait(i,2)-p2_i)*beta_dev;
    arma::mat Ai = arma::join_cols(arma::join_cols(beta_dev,gamma1_dev),gamma2_dev);

    arma::mat beta_dev_1 = Xi_1.t();
    arma::mat psi1_dev = beta_dev_1 * (gamma1 * beta_dev * p1_i*(1-p1_i)-gamma2 * beta_dev*(p1_i*p2_i));
    arma::mat psi2_dev = beta_dev_1 * (-gamma1 * beta_dev * p1_i*p2_i+gamma2 * beta_dev*(p2_i*(1-p2_i)));
    arma::mat Bi = arma::join_cols(psi1_dev,psi2_dev);

    arma::mat Di_11(dim_cov_1,dim_cov_1);
    Di_11.fill(pi1*(1-pi1));
    arma::mat Di_12(dim_cov_1,dim_cov_1);
    Di_12.fill(-pi1*pi2);
    arma::mat Di_21(dim_cov_1,dim_cov_1);
    Di_21.fill(-pi1*pi2);
    arma::mat Di_22(dim_cov_1,dim_cov_1);
    Di_22.fill(pi2*(1-pi2));
    arma::mat Di = arma::join_cols(arma::join_rows(Di_11 % (beta_dev_1 * beta_dev_1.t()),Di_12 % (beta_dev_1 * beta_dev_1.t())),arma::join_rows(Di_21 % (beta_dev_1 * beta_dev_1.t()),Di_22 % (beta_dev_1 * beta_dev_1.t())));

    arma::mat Ci1_11(dim_cov,dim_cov_1);
    Ci1_11.fill(-p1_i*(1-p1_i));
    arma::mat Ci1_12(dim_cov,dim_cov_1);
    Ci1_12.fill(p1_i*p2_i);
    arma::mat Ci1_21(dim_cov,dim_cov_1);
    Ci1_21.fill(p1_i*p2_i);
    arma::mat Ci1_22(dim_cov,dim_cov_1);
    Ci1_22.fill(-p2_i*(1-p2_i));
    arma::mat Ci_1 = arma::join_cols(arma::join_rows(Ci1_11 % (beta_dev * beta_dev_1.t()),Ci1_12 % (beta_dev * beta_dev_1.t())),arma::join_rows(Ci1_21 % (beta_dev * beta_dev_1.t()),Ci1_22 % (beta_dev * beta_dev_1.t())));
    arma::mat Ci = arma::join_cols(arma::zeros(beta_dev.n_rows,Di.n_cols),Ci_1);
    double err_1 = err(i);

    arma::mat dev_11 = -Ai * Ai.t();
    arma::mat temp(Ci.n_rows,Ci.n_cols);
    temp.fill(err_1);
    arma::mat dev_12 = -Ai * Bi.t() + Ci % temp;
    arma::mat dev_21 = arma::zeros(dev_12.n_cols,dev_12.n_rows);
    arma::mat dev_22 = Di;
      
    gra = gra + arma::join_cols(arma::join_rows(dev_11,dev_12),arma::join_rows(dev_21,dev_22));
      
    arma::mat U_i = Ai * err_1;
    arma::mat psi1_dev_1 = beta_dev_1 * (trait(i,1)-pi1);
    arma::mat psi2_dev_1 = beta_dev_1 * (trait(i,2)-pi2);
    arma::mat L_i = arma::join_cols(psi1_dev_1,psi2_dev_1);
    arma::mat R_i = arma::join_cols(U_i,L_i);
    var_gee = var_gee + R_i * R_i.t();
  }
  arma::mat var_all = inv(gra) * var_gee * (inv(gra).t());
  double sd_beta = sqrt(var_all(1,1));
  return(sd_beta);
}