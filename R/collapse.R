  
#####  argument checking and binning #####

collapse <- function(v,x,mu=NULL,bins=NULL,listx=TRUE){

	if(inherits(v,c("simple_triplet_matrix"))){
    v <- sparseMatrix(i=v$i,j=v$j,x=v$v,
            dims=dim(v),dimnames=dimnames(v))
  }
	if(is.null(dim(v))) v <- as.matrix(v)
  rownames(v) <- NULL

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
  rownames(x) <- NULL
  if(is.null(colnames(x))) colnames(x) <- 1:ncol(x)

  if(nrow(x) != nrow(v)) 
    stop("counts and covars have a different number of observations")

  ## uncollapsed exit
  if(is.null(bins)){
    if(is.null(mu)) mu <- suppressWarnings(log(rowMeans(x)+1))
    if(length(mu)==1) mu <- rep(mu,nrow(x))
    mu[is.infinite(mu)] <- -1e6 # a zero row total
    nbin <- rep(1,nrow(x))
    if(listx) x <- sapply(colnames(x),
                function(j) x[,j,drop=FALSE],simplify=FALSE)
    return(list(v=v,x=x,nbin=nbin,mu=mu))
  }

  ## binning
  qs <- (0:bins)/bins
  cutit <- function(vj){
    if(length(unique(vj))<=bins) return(factor(vj)) 
    return(cut(vj, 
            breaks=unique(quantile(vj, qs)), 
            include.lowest=TRUE))
  }
  B <- apply(v,2,cutit) 
  I <- interaction(as.data.frame(B), drop=TRUE)
  v <- apply(v,2,function(vj) tapply(as.numeric(vj), I, mean))
  nbin <- table(I)

  xstm <- summary(x)
  x <- sparseMatrix(i=as.numeric(I)[xstm$i],
    j = xstm$j, x=xstm$x,
    dims=c(nlevels(I),ncol(x)),
    dimnames=list(levels(I),colnames(x)))
  if(!is.null(mu)) warning("pre-specified mu is ignored after binning.")
  mu <- suppressWarnings(log(rowMeans(x)+1))
  if(listx) x <- sapply(colnames(x),
              function(j) x[,j,drop=FALSE],simplify=FALSE)

  return(list(v=v,x=x,nbin=nbin,mu=mu))

}



