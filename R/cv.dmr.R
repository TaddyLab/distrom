##### OOS experimentation with distributed MN regression #####

## multinomial deviance utility
mndev <- function(s, x, m, v, f){
    E <- predict(f,newdata=v,select=s)
    L <- rowSums(E*x) - m*log(rowSums(exp(E)))
    -2*mean(L)
}

## outer R loop that calls dmr
cv.dmr <- function(counts, covars, 
                  lambda.start=NULL,
                  nfold=5, foldid=NULL, 
                  verb=TRUE, cl=NULL, ...){
  
  ## basic input checking
  chk <- collapse(counts,covars,listx=FALSE)
  x <- chk$x
  v <- chk$v
  m <- rowSums(x)
  nobs <- sum(chk$nbin)

  ## set shared lambda.start
  if(is.null(lambda.start))
  {  
    u <- rowMeans(x)
    e0 = as.matrix(x-outer(u,colSums(x)/sum(u)))
    g0 = abs(t(v)%*%e0)/nrow(v)
    std = list(...)$standardize
    if(is.null(std)) std = TRUE
    if(std){
      vs <- sqrt(colSums(v^2)/nrow(v) - colMeans(v)^2)
      vs[vs==0] <- 1
      g0 <- g0/vs
    }
    lambda.start <- max(g0) 
  }

  ## parallel setup
  stopcl = FALSE
  if(is.null(cl)){
    cl <- makeCluster(detectCores(),
                    type=ifelse(
                      .Platform$OS.type=="unix",
                      "FORK","PSOCK"))
    stopcl = TRUE
  }

  ## full fit and properties
  full <- dmr(x, v, lambda.start=lambda.start, cl=cl, ...)

  ## get lambda
  lambda <- sort(
    unique(unlist(lapply(full,function(f) f$lambda))),
    decreasing=TRUE)
  nlambda <- length(lambda)

  #### get full model MN deviance and df
  dev <- parSapply(cl,1:nlambda,mndev,x,m,v,full)
  names(dev) <- paste("seg",1:nlambda,sep=".")
  df <- sapply(full, 
          function(f) c(f$df,rep(NA,nlambda-length(f$df))))

  ## set the folds
  if(is.null(foldid)){
    nfold <- min(nfold,nobs)
    rando <- sample.int(nobs)
    chunks <- round(seq.int(0,nobs,length.out=nfold+1))
    foldid <- rep.int(1:nfold,times=diff(chunks))[rando]
  } else  stopifnot(length(foldid)==nobs)
  foldid <- factor(foldid)
  nfold <- nlevels(foldid)

  oos <- matrix(Inf, nrow=nfold, 
            ncol=nlambda, 
            dimnames=list(levels(foldid),
              paste("seg",1:nlambda,sep=".")))

  if(verb) cat("fold ")
  for(k in levels(foldid)){
    train <- which(foldid!=k)
    fit <- dmr(x[train,],v[train,], ..., 
              lambda.start=lambda.start)
    oos[k,] <- parSapply(cl,
                1:nlambda,mndev,
                x[-train,],m[-train],v[-train,],fit)
    if(verb) cat(sprintf("%s,",k))
  }
  
  cvm <- apply(oos,2,mean)
  cvs <- apply(oos,2,sd)/sqrt(nfold-1)
  seg.min <- which.min(cvm)
  lambda.min = lambda[seg.min]
  cv1se <- (cvm[seg.min]+cvs[seg.min])-cvm
  seg.1se <- min((1:length(cvm))[cv1se>=0])
  lambda.1se = lambda[seg.1se]

  if(stopcl) stopCluster(cl)

  if(verb) cat("done.\n")
  out <- list(dmr=full,
          nobs=nobs,
          dev=dev,
          df=df,
          lambda=lambda,
          nfold=nfold,
          foldid=foldid,
          cvm=cvm,
          cvs=cvs,
          seg.min=seg.min,
          seg.1se=seg.1se,
          lambda.min=lambda.min,
          lambda.1se=lambda.1se)

  class(out) <- "cv.dmr"
  invisible(out)
}

## S3 method functions
plot.cv.dmr <- function(x, ...){
  x$gamlr <- list(lambda=x$lambda, df=rowMeans(x$df))
  x$family <- "multinomial"
  gamlr:::plot.cv.gamlr(x, ...)
}

coef.cv.dmr <- function(object, 
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  coef(object$dmr, select=object[[seg]])
}

logLik.cv.dmr <- function(object, ...){
  ll <- -0.5*object$dev*object$nobs
  attr(ll,"nobs") = object$nobs
  attr(ll,"df") = rowSums(object$df)
  class(ll) <- "logLik"
  ll
}


