##### Distributed Logistic Multinomial Regression  ######

## define class
setClass("dmrcoef", contains="dgCMatrix")

## inner loop function
onerun <- function(xj, argl){
  argl$y <- xj
  if(argl$cv) fit <- do.call(cv.gamlr,argl)
  else fit <- do.call(gamlr,argl)
  ## print works only if you've specified an outfile in makeCluster
  if(length(fit$lambda)<argl$nlambda) print(colnames(xj))
  return(fit)
}

## main function
dmr <- function(cl, covars, counts, mu=NULL, bins=NULL, verb=0, cv=FALSE, ...)
{
  if(!is.null(cl)){
    if(!inherits(cl,"cluster")) stop("first argument `cl' must be NULL or a socket cluster.")
  }
  #build the default argument list
  argl <- list(...)
  argl$family="poisson"
  if(is.null(argl$nlambda))
    argl$nlambda <- formals(gamlr)$nlambda
  argl$verb <- max(verb-1,0)
  argl$cv <- cv

  ## collapse and clean
  chk <- collapse(covars, counts, mu, bins)
  if(verb)
    cat(sprintf("fitting %d observations on %d categories, %d covariates.\n",
        nrow(chk$v), ncol(chk$counts), ncol(chk$v)))
  argl$x <- chk$v
  argl$fix <- chk$mu
  nobs <- sum(chk$nbin)
  p <- ncol(chk$counts)
  vars <- colnames(chk$counts)
  ## cleanup
  rownames(argl$x) <- rownames(chk$counts) <- NULL
  counts <- chk$counts
  rm(covars,mu,chk)

  ## convert X to list
  if(verb) cat("converting counts matrix to column list...\n")
  C <- ifelse(is.null(cl),Inf,length(cl))
  if(C < p/4){
    chunks <- round(seq(0,p,length.out=C+1))
    counts <- lapply(1:C, 
      function(i) counts[,(chunks[i]+1):chunks[i+1]])
    counts <- parLapply(cl,
                counts, 
                function(x) 
                  sapply(colnames(x), 
                  function(j) x[,j,drop=FALSE]))
    counts <- unlist(counts,recursive=FALSE)
  } else{
    counts <- sapply(vars,
      function(j) counts[,j,drop=FALSE]) }

  ## lapply somehow, depending on cl and p
  if(is.null(cl)){
    if(verb) cat("running in serial.\n ")
    mods <- lapply(counts,onerun,argl=argl) 
  }
  else{
    if(verb){ 
     cat("distributed run.\n") 
     print(cl) }
    mods <- parLapply(cl,counts,onerun,argl=argl) 
  }
    
  ## classy exit
  class(mods) <- "dmr"
  attr(mods,"nobs") <- nobs
  attr(mods,"nlambda") <- argl$nlambda
  attr(mods,"mu") <- argl$fix
  return(mods)
}

coef.dmr <- function(object, ...){
  B <- do.call(cBind, lapply(object,coef, ...))
  colnames(B) <- names(object)
  B <- as(as(B,"dgCMatrix"),"dmrcoef")
  return(B)
}

## method predict functions
predict.dmr <- function(object, newdata, 
                  type=c("link","response","class"), ...){
  B <- coef(object, ...)
  predict(B,newdata=newdata,type=type)
}

predict.dmrcoef <- function(object, newdata, 
                  type=c("link","response","class"), ...)
{
  if(inherits(newdata,"simple_triplet_matrix"))
    newdata <- sparseMatrix(i=newdata$i,j=newdata$j,x=newdata$v,
      dims=dim(newdata),dimnames=dimnames(newdata))
  if(is.vector(newdata)){ newdata <- matrix(newdata, nrow=1) }
  if(is.data.frame(newdata)){ newdata <- as.matrix(newdata) }

  type=match.arg(type)
  if(type=="reduction")
    stop("type `reduction' has been replaced by the `srproj' function in the textir library.")

  eta <- t(tcrossprod(t(object[-1,,drop=FALSE]),newdata) + object[1,])
  if(type=="response"){
    expeta <- exp(eta)
    eta <- expeta/rowSums(expeta) }
  rownames(eta) <- rownames(newdata)
  colnames(eta) <- colnames(object)
  
  if(type=="class"){
    c <- apply(eta,1,function(e) colnames(eta)[which.max(e)])
    return(c)
  }
  else return(as.matrix(eta))

}

setGeneric("predict")
setMethod("predict","dmrcoef",predict.dmrcoef)
