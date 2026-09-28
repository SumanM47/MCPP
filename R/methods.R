## PRINT METHOD

#' @title Print function for MCPP class objects
#'
#' @param x object of class MCPP
#' @param digits number of digits to print
#' @param ci level of confidence interval, between 0 and 1
#' @param ... additional parameters supplied to the function
#'
#' @import stats
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

  lambdaCtab <- cbind(lambdaCest,lambdaCsd,lambdaCll,lambdaCmed,lambdaCul)
  rownames(lambdaCtab) <- x$setup$uparent
  colnames(lambdaCtab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))


  if(!is.null(x$lambdaO)){
  lambdaOest <- colMeans(x$lambdaO)
  lambdaOmed <- apply(x$lambdaO,2,"median")
  lambdaOsd <- apply(x$lambdaO,2,"sd")
  lambdaOll <- apply(x$lambdaO,2,quantile,probs=la)
  lambdaOul <- apply(x$lambdaO,2,quantile,probs=ua)

  lambdaOtab <- cbind(lambdaOest,lambdaOsd,lambdaOll,lambdaOmed,lambdaOul)
  rownames(lambdaOtab) <- x$setup$other
  colnames(lambdaOtab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))
  }

  mu0est <- colMeans(x$mu0)
  mu0med <- apply(x$mu0,2,"median")
  mu0sd <- apply(x$mu0,2,"sd")
  mu0ll <- apply(x$mu0,2,quantile,probs=la)
  mu0ul <- apply(x$mu0,2,quantile,probs=ua)

  mu0tab <- cbind(mu0est,mu0med,mu0sd,mu0ll,mu0med,mu0ul)
  rownames(mu0tab) <- x$setup$popair[,2]
  colnames(mu0tab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))

  hest <- colMeans(x$h)
  hmed <- apply(x$h,2,"median")
  hsd <- apply(x$h,2,"sd")
  hll <- apply(x$h,2,quantile,probs=la)
  hul <- apply(x$h,2,quantile,probs=ua)

  htab <- cbind(hest,hsd,hll,hmed,hul)
  rownames(htab) <- x$setup$popair[,2]
  colnames(htab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))

  cat("\nAnalysis of Multivariate Cluster Point Process \n")
  cat("Number of scans = ",x$setup$nreps,"\n",sep="")
  cat("Initial Burn-in samples = ",x$setup$burn,"\n",sep="")
  cat("Thinning interval =",x$setup$thin,"\n",sep="")
  cat("Number of posterior samples = ",nrow(x$lambdaC),"\n",sep="")

  cat("Parent-Offspring Pairs in analysis: \n")
  print(x$setup$popair,quote=F)

  cat("Parent Intensity \n")
  print(round(lambdaCtab,digits))

  cat("Offspring Density \n")
  print(round(mu0tab,digits))

  cat("Bandwidth \n")
  print(round(htab,digits))

  cat("Other Taxa (if present) \n")
  if(!is.null(x$lambdaO)){
    print(round(lambdaOtab,digits))
  }

  cat("DIC \n")
  print(round(x$DIC,digits))

  invisible()
}

## SUMMARY METHOD ##

#' @title Summary function for MCPP class objects
#'
#' @param object object of class MCPP
#' @param digits number of digits to print
#' @param ci level of confidence interval, between 0 and 1
#' @param ... additional parameters supplied to the function
#'
#' @import stats
#'
#' @return Summary outputs from the object
#'
#' @exportS3Method

summary.MCPP <- function(object,digits=3,ci=0.95,...){
  la <- (1-ci)/2; ua <- 1-la

  lambdaCest <- colMeans(object$lambdaC)
  lambdaCmed <- apply(object$lambdaC,2,"median")
  lambdaCsd <- apply(object$lambdaC,2,"sd")
  lambdaCll <- apply(object$lambdaC,2,quantile,probs=la)
  lambdaCul <- apply(object$lambdaC,2,quantile,probs=ua)

  lambdaCtab <- cbind(lambdaCest,lambdaCsd,lambdaCll,lambdaCmed,lambdaCul)
  rownames(lambdaCtab) <- object$setup$uparent
  colnames(lambdaCtab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))

  lambdaOtab <- NULL
  if(!is.null(object$lambdaO)){
    lambdaOest <- colMeans(object$lambdaO)
    lambdaOmed <- apply(object$lambdaO,2,"median")
    lambdaOsd <- apply(object$lambdaO,2,"sd")
    lambdaOll <- apply(object$lambdaO,2,quantile,probs=la)
    lambdaOul <- apply(object$lambdaO,2,quantile,probs=ua)

    lambdaOtab <- cbind(lambdaOest,lambdaOmed,lambdaOsd,lambdaOll,lambdaOmed,lambdaOul)
    rownames(lambdaOtab) <- object$setup$other
    colnames(lambdaOtab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))
  }

  mu0est <- colMeans(object$mu0)
  mu0med <- apply(object$mu0,2,"median")
  mu0sd <- apply(object$mu0,2,"sd")
  mu0ll <- apply(object$mu0,2,quantile,probs=la)
  mu0ul <- apply(object$mu0,2,quantile,probs=ua)

  mu0tab <- cbind(mu0est,mu0med,mu0sd,mu0ll,mu0med,mu0ul)
  rownames(mu0tab) <- object$setup$popair[,2]
  colnames(mu0tab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))

  hest <- colMeans(object$h)
  hmed <- apply(object$h,2,"median")
  hsd <- apply(object$h,2,"sd")
  hll <- apply(object$h,2,quantile,probs=la)
  hul <- apply(object$h,2,quantile,probs=ua)

  htab <- cbind(hest,hsd,hll,hmed,hul)
  rownames(htab) <- object$setup$popair[,2]
  colnames(htab) <- c("Mean","SD",paste(100*la,"%",sep=""),"50%",paste(100*ua,"%",sep=""))

  value <- list("lambdaC"=lambdaCtab,"alpha"=mu0tab,"h"=htab,"lambdaO"=lambdaOtab, "DIC"=object$DIC,"class"=object$class,"setup"=object$setup)
  class(value) <- "summ.MCPP"
  return(value)
}

### PRINT SUMMARY METHOD ###

#' @title Print function for summary of MCPP class objects
#'
#' @param x summary object of class MCPP
#' @param digits number of digits to print
#' @param ... additional parameters supplied to the function
#'
#' @import stats
#'
#' @return prints outputs from the summary object
#'
#' @exportS3Method

print.summ.MCPP <- function(x,digits=3,...){
  cat("Summary of MCPP object",x,"\n",sep=" ")
  cat("Number of scans = ",x$setup$nreps,"\n",sep="")
  cat("Initial Burn-in samples = ",x$setup$burn,"\n",sep="")
  cat("Thinning interval =",x$setup$thin,"\n",sep="")
  cat("Number of posterior samples = ",nrow(x$lambdaC),"\n",sep="")

  cat("DIC = ",x$DIC,"\n",sep="")

  cat("Parent-Offspring Pairs in analysis: \n")
  print(x$setup$popair,quote=F)

  cat("Parent Intensity \n")
  print(round(x$lambdaC,digits))

  cat("Offspring Density \n")
  print(round(x$alpha,digits))

  cat("Bandwidth \n")
  print(round(x$h,digits))

  cat("Other Taxa (if present) \n")
  if(!is.null(x$lambdaO)){
    print(round(x$lambdaO,digits))
  }

  invisible()
}


### PLOT METHOD ###

plot.MCPP(x, trace=TRUE, density=TRUE){
  oldpar <- par(no.readonly = TRUE)

  on.exit(par(old.par))

  nsamp <- nrow(x$lambdaC)
  nparC <- ncol(x$lambdaC)
  nparmu0 <- ncol(x$mu0)
  nparh <- ncol(x$h)
  nparO <- 0
  if(!is.null(x$lambdaO)){nparO <- ncol(x$lambdaO)}

  totpar <- nparC + nparmu0 + nparh + nparO

  ylim <- c(0,0)

  par(mfrow=c(nparC,2),ask=TRUE)
  for(i in 1:nparC){
    obj2use <- x$lambdaC[,i]
    hyp <- x$setup$hyeprParams$lambda
    ddpost <- density(obj2use)
    xvals <- sort(obj2use)
    ddpr <- dgamma(xvals,shape=hyp[1],rate=hyp[2])
    miny <- min(min(ddpost$y),min(ddpr))
    maxy <- max(max(ddpost$y),max(ddpr))
    plot(obj2use,xlab="index",ylab="",main=paste("lambdaC_",i,sep=""))
    plot(ddpost,ylim=c(miny-0.05,maxy + 0.05),xlab="",ylab="density",main=paste("Density of lambdaC_",i,sep=""))
    lines(xvals,ddpr,col="grey",lty=2)
    legend("topright",c("prior","posterior"),col=c("grey","black"),lty=1:2,cex=0.6)
  }

  par(mfrow=c(nparmu0,2),ask=TRUE)
  for(i in 1:nparmu0){
    obj2use <- x$lambdamu0[,i]
    hyp <- x$setup$hyeprParams$alpha
    ddpost <- density(obj2use)
    xvals <- sort(obj2use)
    ddpr <- dgamma(xvals,shape=hyp[1],rate=hyp[2])
    miny <- min(min(ddpost$y),min(ddpr))
    maxy <- max(max(ddpost$y),max(ddpr))
    plot(obj2use,xlab="index",ylab="",main=paste("alpha_",i,sep=""))
    plot(ddpost,ylim=c(miny-0.05,maxy + 0.05),xlab="",ylab="density",main=paste("Density of alpha_",i,sep=""))
    lines(xvals,ddpr,col="grey",lty=2)
    legend("topright",c("prior","posterior"),col=c("grey","black"),lty=1:2,cex=0.6)
  }

  par(mfrow=c(nparh,2),ask=TRUE)
  for(i in 1:nparh){
    obj2use <- x$h[,i]
    hyp <- x$setup$hyeprParams$h
    ddpost <- density(obj2use)
    xvals <- sort(obj2use)
    ddpr <- 2*dnorm(xvals,mean=0,sd=hyp)
    miny <- min(min(ddpost$y),min(ddpr))
    maxy <- max(max(ddpost$y),max(ddpr))
    plot(obj2use,xlab="index",ylab="",main=paste("h_",i,sep=""))
    plot(ddpost,ylim=c(miny-0.05,maxy + 0.05),xlab="",ylab="density",main=paste("Density of h_",i,sep=""))
    lines(xvals,ddpr,col="grey",lty=2)
    legend("topright",c("prior","posterior"),col=c("grey","black"),lty=1:2,cex=0.6)
  }

  if(nparO > 0){
    par(mfrow=c(nparO,2),ask=TRUE)
    for(i in 1:nparO){
      obj2use <- x$lambdaO[,i]
      hyp <- x$setup$hyeprParams$lambdaO
      ddpost <- density(obj2use)
      xvals <- sort(obj2use)
      ddpr <- dgamma(xvals,shape=hyp[1],rate=hyp[2])
      miny <- min(min(ddpost$y),min(ddpr))
      maxy <- max(max(ddpost$y),max(ddpr))
      plot(obj2use,xlab="index",ylab="",main=paste("lambdaO_",i,sep=""))
      plot(ddpost,ylim=c(miny-0.05,maxy + 0.05),xlab="",ylab="density",main=paste("Density of lambdaO_",i,sep=""))
      lines(xvals,ddpr,col="grey",lty=2)
      legend("topright",c("prior","posterior"),col=c("grey","black"),lty=1:2,cex=0.6)
    }
  }


  invisible()
}
