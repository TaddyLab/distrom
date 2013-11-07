##### Distributed Logistic Multinomial Regression  ######

## define class
setClass("dmrcoef",
  representation(lambda="numeric"), 
  contains="dgCMatrix")

## inner loop function
onerun <- function(xj, argl){
  argl$y <- xj
  fit <- do.call(gamlr,argl)
  ## below only works if you've specified an outfile in makeCluster
  if(length(fit$lambda)<argl$nlambda) print(colnames(xj))
  return(fit)
}

## convert Matrix to list of columns
listcol <- function(x,cl){
  rownames(x) <- NULL
  if(length(cl) > 3*ncol(x)){
    N <- length(cl)
    chunks <- round(seq.int(0,ncol(x),length.out=N+1))
    xblock <- lapply(1:N, 
      function(i) x[,(chunks[i]+1):chunks[i+1]])
    xlist <- parLapply(cl, xblock, 
      function(x) sapply(colnames(x), function(j) x[,j,drop=FALSE]))
    xlist <- unlist(xlist, recursive=FALSE)
  }
  else{ xlist <- sapply(colnames(x), function(j) x[,j,drop=FALSE]) }
 return(xlist) 
}


## main function
dmr <- function(covars, counts, mu=NULL, bins=NULL, cl=NULL, ...)
{
  #build the default argument list
  argl <- list(...)
  if(is.null(argl$family))
    argl$family="poisson"
  if(is.null(argl$nlambda))
    argl$nlambda <- formals(gamlr)$nlambda
  if(is.null(argl$verb))
    argl$verb <- FALSE

  ## start cluster
  stopcl = FALSE
  if(is.null(cl)){
    cl <- makeCluster(detectCores(), 
                    type=ifelse(
                      .Platform$OS.type=="unix",
                      "FORK","PSOCK"))
    stopcl = TRUE
  }
  if(argl$verb) print(cl)

  ## collapse and convert counts to list
  chk <- collapse(covars, counts, mu, bins)
  argl$x <- chk$v
  argl$fix <- chk$mu
  nobs <- sum(chk$nbin)
  counts <- listcol(chk$counts,cl)
  rm(chk,covars,mu)

  mods <- parLapply(cl,counts,onerun,argl=argl)
  if(stopcl) stopCluster(cl)

  ## align names (probably unnecessary)
  mods <- mods[names(counts)]

  ## classy exit
  class(mods) <- "dmr"
  attr(mods,"nobs") <- nobs
  attr(mods,"nlambda") <- argl$nlambda
  attr(mods,"mu") <- argl$fix
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
    z <- as.matrix(z)
    return(cbind(z,m=m))
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
