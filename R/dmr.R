##### Distributed Logistic Multinomial Regression  ######

## define class
setClass("dmrcoef",
  representation(lambda="numeric"), 
  contains="dgCMatrix")

## inner loop function(s)
linapprox <- function(xj, argl){
  if(!is.null(argl$zeta))
    zeta <- argl$zeta 
  else 
    zeta <- 0.05
  cat(sprintf("linear approx with zeta = %g\n",zeta))

  if(zeta>0){
    argl$y <- log(xj[,1]+zeta) - argl$fix
    argl$obsweight <- xj[,1]+zeta
  } 
  else{
    argl$y <- log(xj@x) - argl$fix[xj@i+1]
    argl$obsweight <- xj@x
    argl$x <- argl$x[xj@i+1,,drop=FALSE]
  }

  argl$fix <- NULL

  fit <- do.call(gamlr,argl)
  if(length(fit$lambda)<argl$nlambda) 
    print(colnames(xj))

  return(fit)
}

onerun <- function(xj, argl){
  if(argl$family=="gaussian")
    return(linapprox(xj,argl))
  argl$y <- xj
  fit <- do.call(gamlr,argl)
  ## print works only if you've specified an outfile in makeCluster
  if(length(fit$lambda)<argl$nlambda) print(colnames(xj))
  return(fit)
}

## main function
dmr <- function(cl, covars, counts, mu=NULL, bins=NULL, verb=0, ...)
{
  #build the default argument list
  argl <- list(...)
  if(is.null(argl$family))
    argl$family="poisson"
  if(is.null(argl$nlambda))
    argl$nlambda <- formals(gamlr)$nlambda
  argl$verb <- max(verb-1,0)

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

coef.dmr <- function(object, select=NULL, k=2, ...){
  ## model selection
  if(is.null(select)){
    if(exists("AICc")) aic <- AICc(object,k=k)
    else aic <- AIC(object,k=k)
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

  if(type=="reduction")
    stop("type `reduction' has been replaced by the `srproj' function in the textir library.")
  type=match.arg(type)
  
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
