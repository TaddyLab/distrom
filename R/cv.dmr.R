##### OOS experimentation with distributed MN regression #####

## outer R loop that calls dmr
cv.dmr <- function(counts, covars, 
                  lambda.start=NULL,
                  nfold=5, foldid=NULL, 
                  verb=FALSE, ...){
  
  ## basic input checking
  if(inherits(counts,"factor")){
    counts <- sparseMatrix(i=1:length(counts),
      j=as.numeric(counts),
      x=rep(1,length(counts)),
      dimnames=list(names(counts),levels(counts)))
  }
  covars <- as.data.frame(as.matrix(covars))
  nobs <- nrow(covars)

  ## set lambda.start
  if(is.null(lambda.start))
  {  
    u <- rowMeans(counts)
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

  full <- dmr(counts, covars, lambda.start=lambda.start, ...)

  ## get lambda
  lambda.min.ratio <- ifelse(is.null(argl$lambda.min.ratio),
                        formals(gamlr)$lambda.min.ratio,
                        argl$lambda.min.ratio)
  nlambda <- attributes(full$nlambda)
  lambda <- exp(seq(log(lambda.start),
                log(lambda.start*lambda.min.ratio),
                length=nlambda))

  ## set the folds
  if(is.null(foldid)){
    nfold <- min(nfold,nobs)
    rando <- sample.int(nobs)
    chunks <- round(seq.int(0,nobs,length.out=nfold+1))
    foldid <- rep.int(1:nfold,times=diff(chunks))[rando]
  } else  stopifnot(length(foldid)==nobs)
  foldid <- factor(foldid)
  nfold <- nlevels(foldid)

  nlambda <- attributes(full)$nlambda

  oos <- matrix(Inf, nrow=nfold, 
            ncol=nlambda, 
            dimnames=list(levels(foldid),
              paste("seg",1:nlambda,sep=".")))

  if(verb) cat("fold ")
  for(k in levels(foldid)){
    train <- which(foldid!=k)
    fit <- dmr(counts[train,],covars[train], 
              ..., lambda.start=lambda.start)

    for(s in 1:nlambda){
      eta <- predict(fit, covars[-train,], select=s)
      ylo <- cbind(1:nrow(eta),counts[-train,])
      d <- log(rowSums(exp(eta)))*rowSums(counts[-train,])
        D <- D - rowSums(E*Yspend[-ind,])
      eta[ylo] - log(rowSums(exp(eta))) 
      oos[k,s] <- mean(-2*d)
    }
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
          nfold=nfold,
          foldid=foldid,
          lambda=lambda,
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

  argl = list(...)

  argl$x <- log(x$lambda)
  argl$y <- x$cvm
  argl$type <- "n"

  if(is.null(argl$xlab)) argl$xlab="log lambda"
  if(is.null(argl$ylab)) argl$ylab="multinomial deviance"
  if(is.null(argl$pch)) argl$pch=20
  if(is.null(argl$col)) argl$col=4

  cvlo <- x$cvm-x$cvs
  cvhi <- x$cvm+x$cvs

  if(is.null(argl$ylim)) 
    argl$ylim=range(c(cvlo,cvhi),finite=TRUE)
  if(is.null(argl$xlim))
    argl$xlim=range(argl$x[is.finite(argl$y)])


  suppressWarnings(do.call(plot, argl))
  segments(x0=argl$x, y0=cvlo, y1=cvhi, col="grey70")
  argl$type <- NULL
  suppressWarnings(do.call(points, argl))

  abline(v=log(x$lambda.min), lty=3, col="grey20")
  abline(v=log(x$lambda.1se), lty=3, col="grey20")

  xdf <- rowMeans(sapply(x,
    function(f) c(f$df,rep(tail(f$df,1),nlambda-length(f$df)))))

  dfi <- unique(round(seq(1,nlambda,length=ceiling(length(axTicks(1))))))
  axis(3,at=log(lambda[dfi]), labels=round(xdf[dfi])-1,tick=FALSE)
}

coef.cv.dmr <- function(object, 
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  coef(object$dmr, select=object[[seg]])
}

predict.cv.gamlr <- function(object, newdata,
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  predict.dmr(object$dmr, newdata, select=object[[seg]], ...)
}



## internal lhd adjustment
mnadjust <- function(object){
  chk <- attributes(object)$data
  if(is.null(chk)) 
    stop("You need to run dmr with store=TRUE to get grouped deviance.")

  ## undo saturated poisson adjustment
  satd <- chk$x@x*log(chk$x@x) - chk$x@x
  satd[is.nan(satd)] <- 0
  satd <- sum(satd)

  ## calculate difference in normalizing constants
  m <- rowSums(chk$x)
  nlambda <- attributes(object)$nlambda
  dshift <- sapply(1:nlambda, 
    function(s){
        ee <- exp(predict(object,chk$v,select=s))
        sum(ee) - sum(m*log(rowSums(ee))) })
  return(-2*(dshift + satd))

  aic <- aic + mnadjust(object)
      select <- which.min(rowSums(aic, na.rm=TRUE))
}


