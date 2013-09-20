
setClass("dmrcoef",
  representation(lambda="numeric"), 
  contains="dgCMatrix")

##### Distributed Logistic Multinomial Regression  ######
dmr <- function(counts, covars, bins=NULL, 
                lambda.start=NULL, cores=1, ...)
{
  chk <- collapse(counts, covars, bins)
  x <- chk$x
  v <- chk$v

  ## calculate fixed shift
  u <- rowMeans(x)
  mu <- log(chk$n + u)

  ## set path 
  if(is.null(lambda.start))
  {   
    e0 = x-outer(u,colSums(x)/sum(u))
    g0 = abs(t(v)%*%e0)
    std = list(...)$standardize
    if(is.null(std)) std = TRUE
    if(std){
      vs <- sqrt(colSums(v^2)/nrow(v) - colMeans(v)^2)
      g0 <- g0/vs
    }
    lambda.start <- max(g0/nrow(v)) 
  }

  ## grab defaults
  argl <- list(...)
  nlambda <- ifelse(is.null(argl$nlambda),
                    formals(gamlr)$nlambda,
                    argl$nlambda)

  ## inner loop function
  grun <- function(xj){
    fit <- gamlr(v, xj, family="poisson", 
                fix=mu, 
                lambda.start=lambda.start, ...)
    if(length(fit$lambda)<nlambda) 
        write(colnames(xj),stderr())
    return(fit)
  }

  ## parallel computing
  xvar <- colnames(x)
  x <- lapply(xvar,function(j) x[,j,drop=FALSE])
  names(x) <- xvar
  mods <- mclapply(x,grun,mc.cores=cores)[xvar]
  class(mods) <- "dmr"
  attr(mods,"nobs") <- sum(chk$n)
  attr(mods,"cores") <- cores
  return(mods)
}

logLik.dmr <- function(object, ...){
  dev <- sapply(object, function(fit) fit$dev)
  df <- sapply(object, function(fit) fit$df)
  if(!inherits(dev,"matrix")){
      nl <- max(sapply(dev,length))
      dev <- sapply(dev,function(a) c(a,rep(NA,nl-length(a))))
      df <- sapply(df,function(a) c(a,rep(NA,nl-length(a))))
    }

  ll <- -0.5*dev
  attr(ll,"nobs") = attributes(object)$nobs
  attr(ll,"df") = df
  class(ll) <- "logLik"
  ll
}

coef.dmr <- function(object, select=NULL, 
  grouped=TRUE, k=2, cores=attributes(object)['cores'], ...){
  ## model selection
  if(is.null(select)){
    aic <- AIC(object,k=k)
    if(grouped) 
      select <- which.min(rowSums(aic, na.rm=TRUE))
    else{
      select <- apply(aic, 2, which.min)
      select <- sapply(select, 
        function(s) ifelse(length(s)==0,1,s))
    }
  }
  if(length(select)==1) select <- rep(select, ncol(aic))

  ## process, match names, and output
  B <- mcmapply(function(f,s) as.matrix(coef(f,s)), 
          object, select, mc.cores=cores)
  B <- as(as(B,"dgCMatrix"),"dmrcoef")
  B@lambda <- mapply(function(f,s) f$lambda[s], object, select)
  
  return(B)
}

## method predict functions
predict.dmr <- function(object, newdata, 
                    type=c("link","response","reduction"), ...){
  B <- coef(object, ...)
  predict(B,newdata=newdata,type=type)
}

predict.dmrcoef <- function(object, newdata, 
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
    return(as.data.frame(as.matrix(eta)))
  }
  else{
    eta <- t(tcrossprod(t(object[-1,]),newdata) + object[1,])
    if(type=="response"){
      expeta <- exp(eta)
      eta <- expeta/rowSums(expeta) }
    rownames(eta) <- rownames(newdata)
    colnames(eta) <- colnames(object)
    return(as.matrix(eta))
  }
}

setGeneric("predict")
setMethod("predict","dmrcoef",predict.dmrcoef)
