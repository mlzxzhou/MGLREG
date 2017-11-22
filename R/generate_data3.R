#' Simulate three-group sample.
#'
#' \code{generate_data3} Simulate the whole population data containing 100000 observations with certain subgroup prevalence and generate a three-group sample drawn from the whole set with group proportions (25%, 50%, 25%). The output will be a three-grouop sample for the following analysis. This simulation procedure is used in our paper.
#' @param beta coefficients of E(Y|G,X).
#'
#' @param gamma1 coefficients of E(Y|G,X,D=1)-E(Y|G,X,D=0)
#'
#' @param gamma2 coefficients of E(Y|G,X,D=2)-E(Y|G,X,D=0)
#'
#' @param G Genetic covariate, valued 0,1,2 according to the allele numbers for a certain SNP. A vector with length n (n is number of observations).
#'
#' @param X Environmental Covariates. An n by q matrix (n: number of observations, q: number of covariates).
#'
#' @param s number of observation in the sample
#'
#' @return data: The generated sample, containing secondary trait y, primary trait D, G and X
#'
#' @author Fan Zhou
#'
#' @references
#'
#' Fan Zhou, Haibo Zhou, Tengfei Li, Hongtu Zhu. "Analysis of Secondary Phenotypes in Multiple-group Association Studies"
#'
#' @examples
#'
#' #suppose the size of the whole population is 100000 with three subgruops and the proportions of the three groups are (10%, 15%, 75%)
#' # Generate genetic covariate G with MAF=0.3
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
#'
#' @export
generate_data3 <- function(beta, gamma1, gamma2, X, G, s){
  num_sub <- length(G)
  cov <- cbind(rep(1,num_sub),G,X)
  dim_cov <- dim(cov)[2]
  data1 <- matrix(0,nrow=100000,ncol=3*dim_cov+4)
  gamma <- c(gamma1,gamma2)
  data1[,3:(dim_cov+1)] <- cov[,-1]
  prob <- matrix(0,100000,3)
  group_1 <- sample(100000,15000)
  data1[group_1,1] <- 1
  com <- which(data1[,1]!=1)
  group_2 <- sample(com,75000)
  data1[group_2,1]<- 2
  gen_fit <- multinom(data1[,1]~cov[,-1],trace=FALSE)
  psi_1 <- summary(gen_fit)$coefficients[1,]
  psi_2 <- summary(gen_fit)$coefficients[2,]
  ee_1 <- exp(cov%*%t(t(psi_1)))
  ee_2 <- exp(cov%*%t(t(psi_2)))
  prob[,1] <- 1/(1+ee_1+ee_2)
  prob[,2] <- ee_1/(1+ee_1+ee_2)
  prob[,3] <- ee_2/(1+ee_1+ee_2)
  Z <- matrix(0,100000,2*dim_cov)
  for(i in 1:100000){
    rn <- rmultinom(1, 1, prob[i,])
    data1[i,1] <- which(rn!=0)-1
    trait <- rep(0,3)
    status <- which(rn!=0)-1
    trait[status+1] <- 1
    data1[i,(3*dim_cov+2):(3*dim_cov+4)] <- trait
    Z[i,1:dim_cov] <- (trait[2]-prob[i,2])*cov[i,]
    Z[i,(dim_cov+1):(2*dim_cov)] <- (trait[3]-prob[i,3])*cov[i,]
    data1[i,2] <- beta[1]+sum(data1[i,3:(dim_cov+1)]*beta[-1])+sum(Z[i,]*gamma)+rnorm(1, mean = 0, sd = 2)
  }
  data1[,(dim_cov+2):(3*dim_cov+1)] <- Z
  group_1 <- data1[data1[,1]==0,]
  group_2 <- data1[data1[,1]==1,]
  group_3 <- data1[data1[,1]==2,]
  n_1 <- dim(group_1)[1]
  n_2 <- dim(group_2)[1]
  n_3 <- dim(group_3)[1]
  nT <- n_1+n_2+n_3
  p_1 <- n_1/nT
  p_2 <- n_2/nT
  p_3 <- 1-p_1-p_2
  nT1 <- s

  Ind_1 <- sample(n_1, s*0.25, replace = FALSE)
  Ind_2 <- sample(n_2, s*0.5, replace = FALSE)
  Ind_3 <- sample(n_3, s*0.25, replace = FALSE)
  data1_train <- group_1[Ind_1,]
  data2_train <- group_2[Ind_2,]
  data3_train <- group_3[Ind_3,]
  data <- rbind(data1_train,data2_train,data3_train)

  data
}
