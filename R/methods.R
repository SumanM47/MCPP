## PRINT METHOD

#' @title Print function for MCPP class objects
#'
#' @param x object of class MCPP
#' @param digits number of digits to print
#' @param ci level of confidence interval, between 0 and 1
#' @param ... additional parameters supplied to the function
#'
#' @return prints outputs from the object
#'
#' @exportS3Method

print.MCPP <- function(x,digits=3,ci=0.95,...){

  la <- (1 - ci)/2; ua <- 1-la

  lambdaCest <- colMeans(x$lambdaC)
  lambdaCmed <- apply(x$lambdaC,2,"median")
  lambdaCsd <- apply(x$lambdaC,2,"sd")
  lambdaCll <- apply(x$lambdaC,2,quantile,probs=la)
  lambdaCul <- apply(x$lambdaC,2,quantile,probs=ua)

  lambdaCtab <- cbind(lambdaCest,lambdaCmed,lambdaCsd,lambdaCll,lambdaCul)
  rownames(lambdaCtab) <- unique(x$setup$papir[,1])
  colnames(lambdaCtab) <- c("EST","MED","SD","LL","UL")


  if(!is.null(x$lambdaO)){
  lambdaOest <- colMeans(x$lambdaO)
  lambdaOmed <- apply(x$lambdaO,2,"median")
  lambdaOsd <- apply(x$lambdaO,2,"sd")
  lambdaOll <- apply(x$lambdaO,2,quantile,probs=la)
  lambdaOul <- apply(x$lambdaO,2,quantile,probs=ua)

  lambdaOtab <- cbind(lambdaOest,lambdaOmed,lambdaOsd,lambdaOll,lambdaOul)
  rownames(lambdaOtab) <- x$setup$other
  colnames(lambdaOtab) <- c("EST","MED","SD","LL","UL")
  }

  mu0est <- colMeans(x$mu0)
  mu0med <- apply(x$mu0,2,"median")
  mu0sd <- apply(x$mu0,2,"sd")
  mu0ll <- apply(x$mu0,2,quantile,probs=la)
  mu0ul <- apply(x$mu0,2,quantile,probs=ua)

  mu0tab <- cbind(mu0est,mu0med,mu0sd,mu0ll,mu0ul)
  rownames(mu0tab) <- x$setup$popair[,2]
  colnames(mu0tab) <- c("EST","MED","SD","LL","UL")

  hest <- colMeans(x$h)
  hmed <- apply(x$h,2,"median")
  hsd <- apply(x$h,2,"sd")
  hll <- apply(x$h,2,quantile,probs=la)
  hul <- apply(x$h,2,quantile,probs=ua)

  htab <- cbind(hest,hmed,hsd,hll,hul)
  rownames(htab) <- x$setup$popair[,2]
  colnames(htab) <- c("EST","MED","SD","LL","UL")

  cat("\nAnalysis of Multivariate Cluster Point Process \n")
  cat("Number of scans = ",x$setup$nreps,"\n",sep="")
  cat("Initial Burn-in samples = ",x$setup$burn,"\n",sep="")
  cat("Thinning interval =",x$setup$thin,"\n",sep="")
  cat("Number of posterior samples = ",nrow(x$lambdaC),"\n",sep="")

  cat("\nParent-Offspring Pairs in analysis: \n")
  print(x$setup$popair,quote=F)

  cat("\nParent Intensity \n")
  print(round(lambdaCtab,digits))

  cat("\nOffspring Density \n")
  print(round(mu0tab,digits))

  cat("\nBandwidth \n")
  print(round(htab,digits))

  cat("\n Other Taxa (if present) \n")
  if(!is.null(x$lambdaO)){
    print(round(lambdaOtab,digits))
  }

  invisible()
}

## SUMMARY METHOD ##

#' @title Summary function for MCPP class objects
#'
#' @param x object of class MCPP
#' @param digits number of digits to print
#' @param ci level of confidence interval, between 0 and 1
#' @param ... additional parameters supplied to the function
#'
#' @return Summary outputs from the object
#'
#' @exportS3Method

summary.MCPP <- function(x,digits=3,ci=0.95,...){
  la <- (1-ci)/2; ua <- 1-la

  lambdaCest <- colMeans(x$lambdaC)
  lambdaCmed <- apply(x$lambdaC,2,"median")
  lambdaCsd <- apply(x$lambdaC,2,"sd")
  lambdaCll <- apply(x$lambdaC,2,quantile,probs=la)
  lambdaCul <- apply(x$lambdaC,2,quantile,probs=ua)

  lambdaCtab <- cbind(lambdaCest,lambdaCmed,lambdaCsd,lambdaCll,lambdaCul)
  rownames(lambdaCtab) <- unique(x$setup$papir[,1])
  colnames(lambdaCtab) <- c("EST","MED","SD","LL","UL")

  lambdaOtab <- NULL
  if(!is.null(x$lambdaO)){
    lambdaOest <- colMeans(x$lambdaO)
    lambdaOmed <- apply(x$lambdaO,2,"median")
    lambdaOsd <- apply(x$lambdaO,2,"sd")
    lambdaOll <- apply(x$lambdaO,2,quantile,probs=la)
    lambdaOul <- apply(x$lambdaO,2,quantile,probs=ua)

    lambdaOtab <- cbind(lambdaOest,lambdaOmed,lambdaOsd,lambdaOll,lambdaOul)
    rownames(lambdaOtab) <- x$setup$other
    colnames(lambdaOtab) <- c("EST","MED","SD","LL","UL")
  }

  mu0est <- colMeans(x$mu0)
  mu0med <- apply(x$mu0,2,"median")
  mu0sd <- apply(x$mu0,2,"sd")
  mu0ll <- apply(x$mu0,2,quantile,probs=la)
  mu0ul <- apply(x$mu0,2,quantile,probs=ua)

  mu0tab <- cbind(mu0est,mu0med,mu0sd,mu0ll,mu0ul)
  rownames(mu0tab) <- x$setup$popair[,2]
  colnames(mu0tab) <- c("EST","MED","SD","LL","UL")

  hest <- colMeans(x$h)
  hmed <- apply(x$h,2,"median")
  hsd <- apply(x$h,2,"sd")
  hll <- apply(x$h,2,quantile,probs=la)
  hul <- apply(x$h,2,quantile,probs=ua)

  htab <- cbind(hest,hmed,hsd,hll,hul)
  rownames(htab) <- x$setup$popair[,2]
  colnames(htab) <- c("EST","MED","SD","LL","UL")

  value <- list("lambdaC"=lambdaCtab,"alpha"=mu0tab,"h"=htab,"lambdaO"=lambdaOtab,"class"=x$class,"setup"=x$setup)
  class(value) <- "summ.MCPP"
  return(value)
}

## PRINT SUMMARY METHOD ###

#' @title Print function for summary of MCPP class objects
#'
#' @param x summary object of class MCPP
#' @param digits number of digits to print
#' @param ... additional parameters supplied to the function
#'
#' @return prints outputs from the summary object
#'
#' @exportS3Method

print.summ.MCPP <- function(x,digits=3,...){
  cat("\nAnalysis of Multivariate Cluster Point Process \n")
  cat("Number of scans = ",x$setup$nreps,"\n",sep="")
  cat("Initial Burn-in samples = ",x$setup$burn,"\n",sep="")
  cat("Thinning interval =",x$setup$thin,"\n",sep="")
  cat("Number of posterior samples = ",nrow(x$lambdaC),"\n",sep="")

  cat("\nParent-Offspring Pairs in analysis: \n")
  print(x$setup$popair,quote=F)

  cat("\nParent Intensity \n")
  print(round(x$lambdaC,digits))

  cat("\nOffspring Density \n")
  print(round(x$mu0,digits))

  cat("\nBandwidth \n")
  print(round(x$h,digits))

  cat("\n Other Taxa (if present) \n")
  if(!is.null(x$lambdaO)){
    print(round(x$lambdaO,digits))
  }

  invisible()
}
