##### Distributed Logistic Multinomial Regression  ######

## define class
setClass("dmrcoef",
  representation(lambda="numeric"), 
  contains="dgCMatrix")

## undocumented inner loop function
porun <- function(xj, v, mu, nlambda, ...){
  fit <- gamlr(v, xj, family="poisson", fix=mu, nlambda=nlambda, ...)
  ## below only works if you've specified an outfile in makeCluster
  if(length(fit$lambda)<nlambda) print(colnames(xj))
  return(fit)
}

## main function
dmr <- function(covars, counts, mu=NULL, bins=NULL, cl=NULL, 
                nlambda=formals(gamlr)$nlambda, ...)
{
  chk <- collapse(covars, counts, mu, bins)
  rm(covars,counts,mu)

  ## parallel computing
  stopcl = FALSE
  if(is.null(cl)){
    cl <- makeCluster(detectCores(),
                    type=ifelse(
                      .Platform$OS.type=="unix",
                      "FORK","PSOCK"))
    stopcl = TRUE
  }
  mods <- parLapply(cl,chk$x,porun,
            v=chk$v,mu=chk$mu,nlambda=nlambda,...)
  if(stopcl) stopCluster(cl)

  ## align names (probably unnecessary)
  mods <- mods[names(chk$x)]

  ## classy exit
  class(mods) <- "dmr"
  attr(mods,"nobs") <- sum(chk$nbin)
  attr(mods,"nlambda") <- nlambda
  return(mods)
}

coef.dmr <- function(object, select=NULL, k=2, ...){
  ## model selection
  if(is.null(select)){
    aic <- AIC(object,k=k)
    select <- apply(aic, 2, which.min)
    select <- sapply(select, 
      function(s) ifelse(length(s)==0,1,s))
  }
  
  ## grab coef
  B <- mapply( 
        function(f,s){
          s <- min(s,length(f$alpha))
          c(f$alpha[s],f$beta[,s]) }, 
        object, select)
 
   ## set class and double check correct naming
  B <- as(as(B[,names(object)],"dgCMatrix"),"dmrcoef")
  rownames(B) <- c("intercept",rownames(object[[1]]$b))
  B@lambda <- mapply(
    function(f,s) 
      ifelse(is.na(f$lambda[s]),f$lambda[which.min(f$lambda)],f$lambda[s]),
    object, select)
  
  return(B)
}

logLik.dmr <- function(object, ...){
  dev <- sapply(object, function(fit) fit$dev)
  df <- sapply(object, function(fit) fit$df)
  if(!inherits(dev,"matrix")){
      nl <- attributes(object)$nlambda
      dev <- sapply(dev,function(a) c(a,rep(NA,nl-length(a))))
      df <- sapply(df,function(a) c(a,rep(NA,nl-length(a))))
    }

  ll <- -0.5*dev
  attr(ll,"nobs") = attributes(object)$nobs
  attr(ll,"df") = df
  class(ll) <- "logLik"
  ll
}

## method predict functions
predict.dmr <- function(object, newdata, 
                    type=c("link","response","class","reduction"), ...){
  B <- coef(object, ...)
  predict(B,newdata=newdata,type=type)
}

predict.dmrcoef <- function(object, newdata, 
                  type=c("link","response","class","reduction"), ...)
{
  if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
  if(is.data.frame(newdata)){ newdata <- as.matrix(newdata) }


  type=match.arg(type)
  if(type=="reduction"){
    m <- rowSums(newdata)
    newdata <- newdata/(m + 1*(m==0))
    z <- tcrossprod(newdata,object[-1,])
    colnames(z) <- rownames(object)[-1]
    rownames(z) <- rownames(newdata)
    z <- as.data.frame(as.matrix(z))
    z$m <- m
    return(z)
  }
  else{
    eta <- t(tcrossprod(t(object[-1,,drop=FALSE]),newdata) + object[1,])
    if(type=="response"){
      expeta <- exp(eta)
      eta <- expeta/rowSums(expeta) }
    rownames(eta) <- rownames(newdata)
    colnames(eta) <- colnames(object)
  }
  if(type=="class"){
    c <- apply(eta,1,function(e) colnames(eta)[which.max(e)])
    return(c)
  }
  else return(as.matrix(eta))

}

setGeneric("predict")
setMethod("predict","dmrcoef",predict.dmrcoef)
