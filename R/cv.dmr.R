##### OOS experimentation with distributed MN regression #####

## multinomial deviance utility
mndev <- function(s, counts, m, v, f){
    B <- coef(f,select=s)
    naz <- which(is.finite(B[1,])) ## not all zero in training
    E <- t(tcrossprod(t(B[-1,,drop=FALSE]),v) + B[1,])
    L <- rowSums((counts*E)[,naz]) - m*log(rowSums(exp(E[,naz])))
    -2*mean(L)
}

## outer R loop that calls dmr
cv.dmr <- function(cl, covars, counts, 
                  mu=NULL, lambda.start=NULL,
                  nfold=5, foldid=NULL, 
                  verb=1, savek=FALSE, ...){
  
  ## basic input checking
  chk <- collapse(covars,counts,mu=mu)
  counts <- chk$counts
  v <- chk$v
  m <- rowSums(counts)
  nobs <- sum(chk$nbin)
  if(is.null(list(...)$bins))
    mu <- chk$mu
  else mu <- NULL

  ## set shared lambda.start
  if(verb) cat("calculating lambda.start...\n")
  if(is.null(lambda.start))
  {  
    u <- exp(chk$mu)
    e0 = as.matrix(counts-outer(u,colSums(counts)/sum(u)))
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

  ## full fit and properties
  if(verb) cat("full model fit, ")
  full <- dmr(cl, v, counts, mu=mu, 
            lambda.start=lambda.start, 
            verb=max(verb-1,0), ...)

  ## get lambda
  lambda <- sort(
    unique(unlist(lapply(full,function(f) f$lambda))),
    decreasing=TRUE)
  nlambda <- length(lambda)

  #### get full model MN deviance and df
  if(!is.null(cl))
    dev <- parSapply(cl,1:nlambda,mndev,counts,m,v,full)  
  else 
    dev <- sapply(1:nlambda,mndev,counts,m,v,full)
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
  
  if(savek) kfit <- vector(nfold, mode="list")
  else kfit <- NULL

  if(verb) cat("fold ")
  mut <- mu
  for(k in levels(foldid)){
    train <- which(foldid!=k)
    if(!is.null(mut)) mut <- mu[train]
    fit <- dmr(cl, v[train,,drop=FALSE], counts[train,], mu=mut,
            lambda.start=lambda.start, verb=max(verb-1,0), ...)
    if(savek) kfit[[k]] <- fit
    if(!is.null(cl))
      oos[k,] <- parSapply(cl,
                1:nlambda,mndev,
                counts=counts[-train,],
                m=m[-train],v=v[-train,,drop=FALSE],f=fit)
    else 
      oos[k,] <- sapply(
                1:nlambda,mndev,
                counts=counts[-train,],
                m=m[-train],v=v[-train,,drop=FALSE],f=fit)

    if(verb) cat(sprintf("%s,",k))
  }
  
  cvm <- apply(oos,2,mean)
  cvs <- apply(oos,2,sd)/sqrt(nfold-1)
  seg.min <- which.min(cvm)
  lambda.min = lambda[seg.min]
  cv1se <- (cvm[seg.min]+cvs[seg.min])-cvm
  seg.1se <- min((1:length(cvm))[cv1se>=0])
  lambda.1se = lambda[seg.1se]

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
          lambda.1se=lambda.1se,
          kfit=kfit)

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


