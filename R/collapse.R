  
#####  argument checking and binning #####

collapse <- function(v,counts,mu=NULL,bins=NULL){

	if(inherits(v,c("simple_triplet_matrix"))){
    v <- sparseMatrix(i=v$i,j=v$j,x=v$v,
            dims=dim(v),dimnames=dimnames(v))
  }
	if(is.null(dim(v))) v <- as.matrix(v)
  rownames(v) <- NULL

  if(inherits(counts,"data.frame")){
    if(ncol(counts)>1) counts <- as.matrix(counts)
    else counts <- factor(counts[,1])
  }
  if(inherits(counts,"factor")){
    counts <- sparseMatrix(i=1:length(counts),
      j=as.numeric(counts),
      x=rep(1,length(counts)),
      dimnames=list(names(counts),levels(counts)))
  }
  if(inherits(counts,"simple_triplet_matrix")){
    counts <- sparseMatrix(i=counts$i,j=counts$j,x=counts$v,
      dims=dim(counts),dimnames=dimnames(counts))
  }
  counts=as(counts,"dgCMatrix") 
  p <- ncol(counts)
  if(is.null(colnames(counts))) colnames(counts) <- 1:p
  n <- nrow(counts)
  if(n != nrow(v)) 
    stop("counts and covars have a different number of observations")

  ## uncollapsed exit
  if(is.null(bins)){
    if(is.null(mu)) mu <- getmu(rowSums(counts),p)
    if(length(mu)==1) mu <- rep(mu,n)
    return(list(v=v,counts=counts,nbin=rep(1,n),mu=mu))
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

  cstm <- summary(counts)
  counts <- sparseMatrix(i=as.numeric(I)[cstm$i],
    j = cstm$j, x=cstm$x,
    dims=c(nlevels(I),p),
    dimnames=list(levels(I),colnames(counts)))
  if(!is.null(mu)) warning("pre-specified mu is ignored after binning.")
  mu <- getmu(rowSums(counts),p)

  return(list(v=v,counts=counts,nbin=nbin,mu=mu))

}

getmu <- function(m,p){ suppressWarnings(log(m+1)) }

