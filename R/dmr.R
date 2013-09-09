
setClass("dmr",representation(lambda="numeric"),contains="dgCMatrix")

##### Distributed Logistic Multinomial Regression  ######
dmr <- function(counts, covars, bins=NULL, 
                k=2, grouped=FALSE, 
                cores=1, type=NULL, ...)
{
  checked <- collapse(counts, covars, bins)
  x <- checked$x
  v <- checked$v

  mu <- rowMeans(x)
  nz <- mu>0
  mu <- mu[nz]
  x <- x[nz,]
  v <- v[nz,,drop=FALSE]

  if(grouped){  
    e0 = x-outer(mu,colSums(x)/sum(mu))
    g0 = abs(t(v)%*%e0)
    std = list(...)$standardize
    if(is.null(std)) std = TRUE
    if(std){
      vs <- sqrt(colSums(v^2)/nrow(v) - colMeans(v)^2)
      g0 <- g0/vs
    }
    lambda.start <- max(g0/nrow(v))
  } else{ lambda.start = Inf }

  grun <- function(xj){
    require(Matrix)
    require(gamlr)
    fit <- gamlr(v, xj, family="poisson", 
                fix=mu, lambda.start=lambda.start, ...)
    if(length(fit$lambda)<100) print(colnames(xj))
    return(fit)
  }

  xvar <- colnames(x)
  x <- lapply(xvar, function(j) x[,j,drop=FALSE])
  names(x) <- xvar

  ## parallel computing
  if(is.null(type)){
    if(.Platform$OS.type == "unix") type <- "FORK"
    else type <- "PSOCK"
  }

  cl <- makeCluster(cores,type=type) 
  #clusterExport(cl, c("v","mu","lambda.start"), envir=environment())
  mods <- parLapplyLB(cl,x,grun)
  stopCluster(cl)

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
