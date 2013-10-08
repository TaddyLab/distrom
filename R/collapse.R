  
#####  argument checking and binning #####

collapse <- function(x,v,bins=NULL){

	if(inherits(v,c("Matrix","simple_triplet_matrix")))
		v <- as.matrix(v)
	v <- as.data.frame(v)

  if(inherits(x,"data.frame")){
    if(ncol(x)>1) x <- as.matrix(x)
    else x <- factor(x[,1])
  }
  if(inherits(x,"factor")){
    x <- sparseMatrix(i=1:length(x),
      j=as.numeric(x),
      x=rep(1,length(x)),
      dimnames=list(names(x),levels(x)))
  }
  if(inherits(x,"simple_triplet_matrix")){
    x <- sparseMatrix(i=x$i,j=x$j,x=x$v,
    dims=dim(x),dimnames=dimnames(x))
  }
  x=as(x,"dgCMatrix") 
  if(is.null(colnames(x))) colnames(x) <- 1:ncol(x)

  if(nrow(x) != nrow(v)) 
    stop("counts and covars have a different number of observations")
  if(is.null(bins)) 
    return(list(x=x,v=v,n=rep.int(1,nrow(x))))

  qs <- (0:bins)/bins
  cutit <- function(vj){
    if(length(unique(vj))<=bins) return(factor(vj))
    return(cut(vj, 
            breaks=unique(quantile(vj, qs)), 
            include.lowest=TRUE))
  }
  B <- apply(v,2,cutit)
  I <- interaction(as.data.frame(B), drop=TRUE)
  vbin <- apply(v,2,function(vj) tapply(as.numeric(vj), I, mean))
  nbin <- table(I)

  xstm <- summary(x)
  xbin <- sparseMatrix(i=as.numeric(I)[xstm$i],
    j = xstm$j, x=xstm$x,
    dims=c(nlevels(I),ncol(x)),
    dimnames=list(levels(I),colnames(x)))

  return(list(x=xbin,v=vbin,n=nbin))

}



