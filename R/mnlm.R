
##### Distributed Logistic Multinomial Regression  ######

mnlm <- function(counts, covars, bins=NULL, cores=1, k=log(n), ...)
{
  checked <- collapse(counts, covars, bins)
  x = checked$x
  v = checked$v
  n = nrow(x)
  mu = rowSums(x)/ncol(x)

  grun <- function(j){
    if(is.infinite(mu[j]))
     return(matrix(0,nrow=ncol(v)))
    fit <- gamlr(v, x[,j], family="poisson", fix=mu, ...)
    if(length(fit$lambda)<100) print(j)
    beta <- coef(fit,which.min(AIC(fit,k=k)))
    as.matrix(beta) }

  B <- lapply(colnames(x), grun)#, mc.cores=cores)
  B <- matrix(unlist(B),ncol=ncol(x))
  rownames(B) <- c("intercept",colnames(v))
  colnames(B) <- colnames(x)
  return(B) 
}



