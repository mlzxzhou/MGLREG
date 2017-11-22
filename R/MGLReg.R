#' MGLReg on multiple-group sample (single-phase).
#'
#' \code{STACPreg3} This function is to deal with the basic single-phase three-group case discussed in the paper by running the unbiased linear regression model MGLReg. The output includes the estimation, asymptotical variance and the corresponding p-value of the genetic effect.
#'
#' @param y Response. Currently, we only consider continuous variables. A vector with length n ( n is number of observations).
#'
#' @param G Genetic covariate, valued 0,1,2 according to the allele numbers for a certain SNP. A vector with length n (n is number of observations).
#'
#' @param X Environmental Covariates such as age, sex and first few PCs to adjust for strafication. An n by q matrix (n: number of observations, q: number of covariates).
#'
#' @param D Primary status. Valued 0,1,2 in the three-subgroup case, representing the belongs of obervation to a certain subgroup.
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
#' @examples
#'
#' # suppose the size of the whole population is 100000 with three subgruops and the proportions of the three groups are (10%, 15%, 75%)
#' # Generate genetic covariate G with MAF=0.2
#' MAF <- 0.3
#' G <- rbinom(100000,2,MAF)
#' s <- 800
#' beta <- c(5, 2, 1)
#' gamma1 <- c(-2,-1,-1)
#' gamma2 <- c(1,1,1)
#' # Generate X, which is mixture of normal variables
#' x_p <- c(0.2,0.8)
#' x_mean <- c(0,2)
#' x_1 <- rnorm(20000, mean = 2, sd = 2)
#' x_2 <- rnorm(80000, mean=0, sd=2)
#' X <- c(x_1,x_2)
#' data <- generate_data3(beta, gamma1, gamma2, X, G, s)
#' reg <- STACPreg3(data[,2],data[,3],data[,4],data[,1],c(0.1,0.15,0.75))
#' reg$coef; reg$pvalue
#'
#' @export

MGLReg <- function(y, G, X, D, p_til){
  num_sub <- length(y)
  cov <- cbind(rep(1,num_sub),G,X)
  dim_cov <- dim(cov)[2]
  d_0 <- which(D==0)
  d_1 <- which(D==1)
  d_2 <- which(D==2)
  trait <- matrix(0,num_sub,3)
  trait[d_0,1] <- 1
  trait[d_1,2] <- 1
  trait[d_2,3] <- 1
  pi_til <- as.numeric(table(D)/length(D))

  ##First step: calculate P(D|X)##
  #construct the vector of response D#
  trait_1 <- factor(D,levels=c("0","1","2"))
  prob_fit_DXG <- multinom(trait_1~cov[,-1],trace=FALSE,Hess = TRUE)
  out <- summary(prob_fit_DXG)
  psi1 <- out$coefficients[1,]
  psi2 <- out$coefficients[2,]
  prob_sam <- prob_fit_DXG$fitted.values
  prob_real <- matrix(0,num_sub,3)
  prob_adjust <- rep(1,3)  #The weights to adjust the prob into the prob in whole population
  prob_adjust[2] <- (pi_til[1]/pi_til[2])*(p_til[2]/p_til[1])
  prob_adjust[3] <- (pi_til[1]/pi_til[3])*(p_til[3]/p_til[1])
  prob_real[,1] <- 1/(prob_sam[,2]/prob_sam[,1]*prob_adjust[2]+prob_sam[,3]/prob_sam[,1]*prob_adjust[3]
                      +1)
  prob_real[,2] <- prob_sam[,2]/prob_sam[,1]*prob_adjust[2]/(prob_sam[,2]/prob_sam[,1]*prob_adjust[2]
                                                             +prob_sam[,3]/prob_sam[,1]*prob_adjust[3]+1)
  prob_real[,3] <- prob_sam[,3]/prob_sam[,1]*prob_adjust[3]/(prob_sam[,2]/prob_sam[,1]*prob_adjust[2]
                                                             +prob_sam[,3]/prob_sam[,1]*prob_adjust[3]+1)

  #Construct the design matrix Z#
  Z <- matrix(0,num_sub,2*dim_cov)
  Z[,1:dim_cov] <- (trait[,2]-prob_real[,2])*cov
  Z[,(dim_cov+1):(2*dim_cov)] <- (trait[,3]-prob_real[,3])*cov
  fit_1 <- lm(y~cov[,-1]+Z)
  out <- summary(fit_1)
  coef_G <- out$coefficients[2,1]
  theta_1 <- t(t(out$coefficients[,1]))
  err <- out$residuals

  ## Asyptotical variance and pvalue
  theta_2 <- rbind(theta_1,t(t(psi1)),t(t(psi2)))
  beta <- theta_2[1:dim_cov]
  gamma1 <- theta_2[(dim_cov+1):(2*dim_cov)]
  gamma2 <- theta_2[(2*dim_cov+1):(3*dim_cov)]

  sd_beta <- MGLREG_var(cov,prob_real,prob_sam,gamma1,gamma2,trait,err,dim_cov,num_sub)
  p_value_G <- 2*(1-pnorm(abs(coef_G/sd_beta)))

  return(list(coef <- coef_G, pvalue <- p_value_G))
}

