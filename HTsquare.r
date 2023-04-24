



HTsquared <- function(data, test){
  if(test == "one-sample"){
    munull <- 0 
    n <- length(data[,1])
    p <- length(data[1,])-1
    ybar <- apply(data[,-1], 2, mean)
    sigma <- var(data[,-1])
    v <- n - 1
    Tsq <- n * t(ybar - munull) %*% solve(sigma) %*% (ybar - munull)
    fstat = (v - p + 1)/(v * p) * Tsq
    pvalue <- 1 - pf(fstat, p, v-p+1)
    my_list <- list("Test" = "One-sample Hotelling T-squared Test", 
                    "Tsq" = Tsq, "F" = fstat, "P-value" = pvalue)
  }
  if(test == "two-sample"){
    n <- length(data[,1])/2
    p <- length(data[1,])-1
    ybar1 <- matrix(apply(data[data[,1]==1,-1], 2, mean), ncol = 1)
    ybar2 <- matrix(apply(data[data[,1]==2,-1], 2, mean), ncol = 1)
    sigma1 <- var(data[data[,1]==1,-1])
    sigma2 <- var(data[data[,1]==2,-1])
    W1 <- (n - 1) * sigma1 
    W2 <- (n - 1) * sigma2
    spool <- 1/(n + n -2) * (W1 + W2)
    Tsq <- n*n/(n+n) * t(ybar1 - ybar2) %*% solve(spool) %*% (ybar1 - ybar2)
    fstat = (n + n - p + 1)/((n + n - 2)* p) * Tsq
    pvalue <- 1 - pf(fstat, p, n + n -p)
    a <- solve(spool) %*% (ybar1 - ybar2)
      if(pvalue <= 0.05 & fstat >= df(0.1,p,n+n-2)){ Conclusion = "Reject Null mu1 = mu2"
      }else{Conclusion = "Fail to reject Null"}
    my_list <- list("Test" = "Multivariate Two-Sample Hotelling T-squared Test", 
                  "Tsq" = Tsq, "F" = fstat, "P-value" = pvalue,  
                  "Discriminant" = a, "Conclusion" = Conclusion)
  }
  if(test == "paired-sample"){
    n <- length(data[,1])/2
    p <- length(data[1,])-1
    d <- as.matrix(data[data[,1]==1,-1] - data[data[,1]==2,-1])
    dbar <- apply(d, 2, mean)
    sigma <- 1/(n-1)* t(d-dbar) %*% (d-dbar)
    Tsq <- n * t(dbar) %*% solve(sigma) %*% dbar 
    fstat <- (n - p + 1)/((n - 2)* p) * Tsq
    pvalue <- 1 - pf(fstat, p, n -p)
    a <- solve(sigma) %*% dbar
      if(pvalue <= 0.05 & fstat >= df(0.1,p,n-1)){ Conclusion = "Reject Null mudiff = 0"}
    else{ Conclusion = "Fail to reject null"}
    my_list <- list("Test" = "Paired Samples Hotelling T-squared Test", 
                  "Tsq" = Tsq, "F" = fstat, "P-value" = pvalue, 
                  "Discriminant" = a, "Conclusion" = Conclusion) }
  return(my_list)
}
