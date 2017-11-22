#' MGLReg on multiple-group sample (two-phase).
#'
#' \code{STACPreg3} This function is to deal with the two-phase case discussed in the paper by running the unbiased linear regression model MGLReg. The output includes the estimation, asymptotical variance and the corresponding p-value of the genetic effect.
#'
#' @param y Response. Currently, we only consider continuous variables. A vector with length n ( n is number of observations).
#'
#' @param G Genetic covariate, valued 0,1,2 according to the allele numbers for a certain SNP. A vector with length n (n is number of observations).
#'
#' @param C Environmental Covariates such as age, sex and first few PCs to adjust for strafication. An n by q matrix (n: number of observations, q: number of covariates).
#'
#' @param D Primary status. Valued 0,1,2 in the three-subgroup case, representing the belongs of obervation to a certain subgroup.
#'
#' @param p Phase status, value 0 for the first phase, 1 for the second phase. We assume that all the subjects are sorted by phase number.
#'
#' @param p_til Subgroup prevalence in the whole population. A vector with length 3 in the three-subgroup case.
#'
#' @return A list object.
#'
#' @return coef: Estimation of coefficient of genetic covariate G in the regression model given by our method.
#'
#'         pvalue : p-value of the corresponding genetic effect
#'
#' @author Fan Zhou
#'
#' @references
#'
#' Fan Zhou, Haibo Zhou, Tengfei Li, Hongtu Zhu. "Analysis of Secondary Phenotypes in Multiple-group Association Studies"
#'
#' @export

MGLReg2 <- function(y, G, C, D, p, p_til){
  num_sub <- length(y)
  X <- cbind(rep(1,num_sub),G,C,p)
  X=as.matrix(X)
  dim_cov <- dim(X)[2]
  X_1 <- X
  dim_cov_1 <- dim(X_1)[2]
  d_0 <- which(D==0)
  d_1 <- which(D==1)
  d_2 <- which(D==2)
  trait <- matrix(0,num_sub,3)
  trait[d_0,1] <- 1
  trait[d_1,2] <- 1
  trait[d_2,3] <- 1
  phase_1=which(p==0)
  phase_2=which(p==1)
  num_p1=length(phase_1)
  num_p2=length(phase_2)
  pi_til <- as.numeric(table(D)/length(D))
  pi_til_1=as.numeric(table(D[phase_1])/length(D[phase_1]))
  pi_til_2=as.numeric(table(D[phase_2])/length(D[phase_2]))
  prob_adjust_1 <- rep(1,3)  #The weights to adjust the prob into the prob in whole population
  prob_adjust_1[2] <- (p_til[1]/p_til[2])*(pi_til_1[2]/pi_til_1[1])
  prob_adjust_1[3] <- (p_til[1]/p_til[3])*(pi_til_1[3]/pi_til_1[1])
  prob_adjust_2 <- rep(1,3)  #The weights to adjust the prob into the prob in whole population
  prob_adjust_2[2] <- (p_til[1]/p_til[2])*(pi_til_2[2]/pi_til_2[1])
  prob_adjust_2[3] <- (p_til[1]/p_til[3])*(pi_til_2[3]/pi_til_2[1])
  eta_1=eta_2=rep(0,2)
  eta_1[1]=log(prob_adjust_1[2])
  eta_1[2]=log(prob_adjust_1[3])
  eta_2[1]=log(prob_adjust_2[2])
  eta_2[2]=log(prob_adjust_2[3])

  ##First step: calculate P(D|X)##
  #construct the vector of response D#
  trait_1 <- factor(D,levels=c("0","1","2"))
  prob_fit_DXG_1 <- multinom(trait_1~X_1+0,trace=FALSE,Hess = TRUE)
  out_1 <- summary(prob_fit_DXG_1)
  psi1 <- out_1$coefficients[1,]
  psi2 <- out_1$coefficients[2,]
  psi=rbind(t(t(psi1)),t(t(psi2)))
  psi[1]=psi[1]-eta_1[1]
  psi[dim_cov_1+1]=psi[dim_cov_1+1]-eta_1[2]
  psi[dim_cov_1]=psi[dim_cov_1]+eta_1[1]-eta_2[1]
  psi[2*dim_cov_1]=psi[2*dim_cov_1]+eta_1[2]-eta_2[2]


  psi1=psi[1:dim_cov_1]
  psi2=psi[(dim_cov_1+1):(2*dim_cov_1)]
  prob_sam_1=matrix(0,num_p1,3)
  prob_sam_1[,1]=1/(1+exp(X_1[1:num_p1,] %*% psi1+eta_1[1])+exp(X_1[1:num_p1,] %*% psi2+eta_1[2]))
  prob_sam_1[,2]=exp(X_1[1:num_p1,] %*% psi1+eta_1[1])/(1+exp(X_1[1:num_p1,] %*% psi1+eta_1[1])+exp(X_1[1:num_p1,] %*% psi2+eta_1[2]))
  prob_sam_1[,3]=exp(X_1[1:num_p1,] %*% psi2+eta_1[2])/(1+exp(X_1[1:num_p1,] %*% psi1+eta_1[1])+exp(X_1[1:num_p1,] %*% psi2+eta_1[2]))
  prob_sam_2=matrix(0,num_p2,3)
  prob_sam_2[,1]=1/(1+exp(X_1[(num_p1+1):num_sub,] %*% psi1+eta_2[1])+exp(X_1[(num_p1+1):num_sub,] %*% psi2+eta_2[2]))
  prob_sam_2[,2]=exp(X_1[(num_p1+1):num_sub,] %*% psi1+eta_2[1])/(1+exp(X_1[(num_p1+1):num_sub,] %*% psi1+eta_2[1])+exp(X_1[(num_p1+1):num_sub,] %*% psi2+eta_2[2]))
  prob_sam_2[,3]=exp(X_1[(num_p1+1):num_sub,] %*% psi2+eta_2[2])/(1+exp(X_1[(num_p1+1):num_sub,] %*% psi1+eta_2[1])+exp(X_1[(num_p1+1):num_sub,] %*% psi2+eta_2[2]))
  prob_sam=rbind(prob_sam_1,prob_sam_2)
  prob_real=matrix(0,num_sub,3)
  prob_real[,1]=1/(1+exp(X_1[1:num_sub,] %*% psi1)+exp(X_1[1:num_sub,] %*% psi2))
  prob_real[,2]=exp(X_1[1:num_sub,] %*% psi1)/(1+exp(X_1[1:num_sub,] %*% psi1)+exp(X_1[1:num_sub,] %*% psi2))
  prob_real[,3]=exp(X_1[1:num_sub,] %*% psi2)/(1+exp(X_1[1:num_sub,] %*% psi1)+exp(X_1[1:num_sub,] %*% psi2))

  #Construct the design matrix Z#
  dim_cov <- dim(X)[2]
  dim_cov_1 <- dim(X_1)[2]
  Z <- matrix(0,num_sub,2*dim_cov)
  Z[,1:dim_cov] <- (trait[,2]-prob_real[,2])*X
  Z[,(dim_cov+1):(2*dim_cov)] <- (trait[,3]-prob_real[,3])*X
  fit_1 <- lm(y~X[,-1]+Z)
  out <- summary(fit_1)
  coef_g <- out$coefficients[2,1]
  theta_1 <- t(t(out$coefficients[,1]))
  err <- out$residuals

  ## Asyptotical variance and pvalue
  theta_2 <- rbind(theta_1,t(t(psi)))
  beta <- theta_2[1:dim_cov]
  gamma1 <- theta_2[(dim_cov+1):(2*dim_cov)]
  gamma2 <- theta_2[(2*dim_cov+1):(3*dim_cov)]

  num_p1 <- length(phase_1)
  num_p2 <- length(phase_2)
  sd_beta <- MGLREG_var2(X,X_1,prob_real,prob_sam,gamma1,gamma2,trait,err,dim_cov,dim_cov_1,num_sub)
  p_value_G <- 2*(1-pnorm(abs(coef_g/sd_beta)))

  return(list(coef <- coef_g, pvalue <- p_value_G))
}

