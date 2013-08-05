
setClass("dmr",contains="dgCMatrix")

##### Distributed Logistic Multinomial Regression  ######
dmr <- function(counts, covars, bins=NULL, 
                cores=1, k=2, ...)
{
  checked <- collapse(counts, covars, bins)
  x <- checked$x
  v <- checked$v

  mu <- rowSums(x)/ncol(x)
  nz <- mu>0
  mu <- mu[nz]
  x <- x[nz,]
  v <- v[nz,,drop=FALSE]

  grun <- function(j){
    if(is.infinite(mu[j]))
     return(matrix(0,nrow=ncol(v)))
    fit <- gamlr(v, x[,j], family="poisson", fix=mu, ...)
    if(length(fit$lambda)<100) print(j)
    beta <- coef(fit,which.min(AIC(fit,k=k)))
    as.matrix(beta)
  }

  B <- mclapply(colnames(x), grun, mc.cores=cores)
  B <- matrix(unlist(B),ncol=ncol(x))
  rownames(B) <- c("intercept",colnames(v))
  colnames(B) <- colnames(x)
  as(as(B,"dgCMatrix"),"dmr")
}

## method predict function
predict.dmr <- function(object, newdata, 
                  type=c("link","response","reduction"), ...)
{
  if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
  if(is.data.frame(newdata)){ newdata <- as.matrix(newdata) }

  type=match.arg(type)
  if(type=="reduction"){
    m <- rowSums(newdata)
    m[m==0] <- 1
    newdata <- newdata/m
    eta <- tcrossprod(newdata,object[-1,])
    colnames(eta) <- rownames(object)[-1]
    rownames(eta) <- rownames(newdata)
  }
  else{
    eta <- t(tcrossprod(t(object[-1,]),newdata) + object[1,])
    if(type=="response"){
      expeta <- exp(eta)
      eta <- expeta/rowSums(expeta) }
    rownames(eta) <- rownames(newdata)
    colnames(eta) <- colnames(object)
  }
  as.matrix(eta)
}

setGeneric("predict")
setMethod("predict","dmr",predict.dmr)
