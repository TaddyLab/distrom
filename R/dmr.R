
setClass("dmr",representation(lambda="numeric"),contains="dgCMatrix")

##### Distributed Logistic Multinomial Regression  ######
dmr <- function(counts, covars, bins=NULL, 
                k=2, grouped=FALSE, 
                cores=1, ...)
{
  checked <- collapse(counts, covars, bins)
  x <- checked$x
  v <- checked$v
  u <- rowMeans(x)
  bn <- checked$n
  if(is.null(bn)) bn <- rep(1,length(u))

  argl <- list(...)
  if(is.null(argl$lambda.start))
  { if(grouped){  
      e0 = x-outer(u,colSums(x)/sum(u))
      g0 = abs(t(v)%*%e0)
      std = list(...)$standardize
      if(is.null(std)) std = TRUE
      if(std){
        vs <- sqrt(colSums(v^2)/nrow(v) - colMeans(v)^2)
        g0 <- g0/vs
      }
      lambda.start <- max(g0/nrow(v))
    } else{ lambda.start = Inf }
  } else{ lambda.start = argl$lambda.start }
  if(is.null(argl$nlambda)) nlambda <- 100
  else nlambda <- argl$nlambda

  grun <- function(xj){
    fit <- gamlr(v, xj, family="poisson", 
                fix=log(bn+u), lambda.start=lambda.start, nlambda=nlambda, ...)
    if(length(fit$lambda)<nlambda) print(colnames(xj))
    return(fit)
  }

  ##### parallel computing
  xvar <- colnames(x)
  x <- lapply(xvar,function(j) x[,j,drop=FALSE])
  names(x) <- xvar
  mods <- mclapply(x,grun,mc.cores=cores)
  #######

  if(grouped){
    aic <- sapply(mods, function(fit) AIC(fit,k=k))
    if(!inherits(aic,"matrix")){
      nl <- max(sapply(aic,length))
      aic <- sapply(aic,function(a) c(a,rep(Inf,nl-length(a))))
    }
    seg <- rep(which.min(rowSums(aic)),length(mods))
  } else{ 
    seg <- sapply(mods, function(fit) which.min(AIC(fit,k=k))) 
  }
  
  B <- mapply(function(f,s) as.matrix(coef(f,s)), mods, seg)
  rownames(B) <- c("intercept",colnames(v))

  lambda <- mapply(function(f,s) f$lambda[s], mods, seg)

  zebra <- match(xvar,colnames(B))
  B <- as(as(B[,zebra],"dgCMatrix"),"dmr")
  B@lambda <- as.numeric(lambda[zebra])
  names(B@lambda) <- names(lambda[zebra])
  return(B)
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
