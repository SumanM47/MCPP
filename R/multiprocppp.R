#' @name multiprocppp
#' @title R function to generate from multiprocess poisson point process
#'
#' @param win owin class object, the observation window
#' @param npar number of parent taxa which are not offsprings themselves
#' @param noff number of offspring taxa. This includes taxa which may be parents to other offspring taxa
#' @param nxtra number of unrelated taxa. Defaults to 0
#' @param lambdapar intensity parameter for the parent processes. Must be a positive scalar or a positive vector of size npar
#' @param mu0vec numeric vector of positive numbers, vector of offspring densities
#' @param hvec numeric vector of positive numbers, vector of bandwidth parameters
#' @param kern character, choice of kernel function. Permitted values are "Gaussian", "Uniform" or "Cauchy"
#' @param lambdaxtra scalar or vector of positive numbers, intensities of the unrelated entities processes
#' @param paroffrel a character matrix of size noff x 2. The second column denotes the offspring and the first column denotes its parent.
#' @param name_par optional character, name of the parent
#' @param names_off character vector, names of the offspring entities
#' @param  names_xtra scalar or vector of characters, names of the unrelated entities
#'
#' @import spatstat.geom
#' @import spatstat.random
#' @import stats
#' @import grDevices
#'
#' @return object of class ppp
#' @export
multiprocppp <- function(win,
                         npar, noff, nxtra=0,
                         lambdapar, mu0vec, hvec, kern="Gaussian", lambdaxtra,
                         paroffrel,
                         names_par=NULL, names_off=NULL,names_xtra=NULL){
  ## Initial Checks
  if(!spatstat.geom::is.owin(win)) stop("win must be object of class owin")
  if(npar <= 0) stop("npar must be a positive integer")
  if(any(lambdapar <= 0) | length(lambdapar) != npar) stop("lambdapar must be a vector of positive reals of size npar")
  if(noff <= 0) stop("noff must be a positive integer")
  if(length(mu0vec)!=1 & length(mu0vec)!=noff) stop("mu0vec must be a positive scalar or vector of length noff")
  if(length(hvec)!=1 & length(hvec)!=noff) stop("hvec must be a positive scalar or vector of length noff")
  if(any(mu0vec<=0)) stop("mu0vec must be positive")
  if(any(hvec<=0)) stop("hvec must be positive")
  if(nxtra < 0) stop("nxtra must be a non-negative integer")
  if(nxtra > 0 & missing(lambdaxtra)) stop("lambdaxtra is missing. Provide a positive scalar or vector of length nxtra")
  if(nxtra > 0 & (length(lambdaxtra)!=1 & length(lambdaxtra)!=nxtra))  stop("Provide a positive scalar or vector of length nxtra")
  if(nxtra > 0 & (any(lambdaxtra <=0))) stop("lambdaxtra must be positive")
  if(!is.null(names_par) & (length(name_par)!=npar | !is.character(name_par))) stop("names_par must a character string of size npar")
  if(!is.null(names_off) & (length(names_off)!=noff & !is.character(names_off))) stop("names_off must be a character vector of size noff")
  if(nxtra > 0 & !is.null(names_xtra) & (length(names_xtra)!=nxtra | !is.character(names_xtra))) stop("names_xtra must be a character vector of size nxtra")
  if(!is.matrix(paroffrel) | nrow(paroffrel)!= noff | ncol(paroffrel)!= 2 | !is.character(paroffrel))stop("paroffrel must be a noff x 2 matrix of characters")
  if(unique(paroffrel[,1]) != npar | unique(paroffrel[,2]) != noff)stop("all parent-offspring pairs must be listed in paroffrel")


  ## Setting labels if not already given
  if(is.null(name_par)){name_par <- paste0("parent",1:npar,sep=""); paroffrel[,1] <- paste0("parent",paroffrel[,1],sep="")}
  if(is.null(names_off)){
    names_off <- paste0("offspring",1:noff,sep=""); paroffrel[,2] <- paste0("offspring",paroffrel[,2],sep="")
  }
  if(nxtra > 0 & is.null(names_xtra)){
    names_xtra <- paste0("xtra",1:nxtra,sep="")
  }

  ## Checking to have a consistent parent-offspring ordering
  flag <- rep(0,length(noff))
  for(i in 1:nrow(paroffrel)){
    par <- paroffrel[i,1]
    offlistpos <- which(paroffrel[,2]==par)
    if(length(offlistpos)>1)stop("cannot have a taxon be offspring to multiple parents")
    if(length(offlistpos)==1){if(offlistpos == i)stop("A taxon canont be its own offspring")}
    if(length(offlistpos)==1){if(offlistpos > i) flag[i] <- offlistpos}
  }

  if(any(flag!=0)){
    inds <- which(flag!=0)
    for(ii in inds){
      paroffrel[c(ii,flag[ii]),] <- paroffrel[c(flag[ii],ii),]
    }
  }

  ## Managing inputs
  if(length(lambdapar) == 1) lambdapar <- rep(lambdapar,npar)
  names(lambdapar) <- names_par
  if(length(mu0vec)==1) mu0vec <- rep(mu0vec,noff)
  names(mu0vec) <- names_off
  if(length(hvec)==1) hvec <- rep(hvec,noff)
  names(hvec) <- names_off
  if(nxtra > 0 & length(lambdaxtra)==1) lambdaxtra <- rep(lambdaxtra,nxtra)
  names(lambdaxtra) <- names_xtra

  if(!(kern %in% c("Gaussian","Cauchy","Uniform"))){stop("Only Gaussian, Uniform and Cauchy kernels are supported")}
  if(kern == "Gaussian"){
    rkern <- rkern.norm
  }
  if(kern == "Cauchy"){
    rkern <- rkern.cauchy
  }
  if(kern == "Uniform"){
    rkern <- rkern.unif
  }

  ntot <- npar+noff+nxtra

  alllist <- vector(mode="list",length=ntot)

  #require('spatstat')

  ## Creating a factor of labels for marking
  types <- c(name_par,names_off,names_xtra)
  factortype <- factor(types, levels = types)

  names(alllist) <- types

  ## First create the parent processes

  for(i in 1:npar){
    parname <- names_par[i]
    parinds <- which(paroffrel[,1] == parname)
    offnames <- paroffrel[parind,2]
    maxmu0vec <- max(mu0vec[offnames])
    lambdapar2use <- lambdapar[i]
    X <- spatstat.random::rpoispp(lambdapar2use*(1-exp(-maxmu0vec)),win=win)
    alllist[[parname]] <- list("x"=X$x,"Y"=X$y, "marks"=rep(factortype[i],X$n))
  }

  for(i in 1:noff){
    offname <- paroffrel[i,2]
    parname <- paroffrel[i,1]
    parentx <- alllist[[parname]]$x
    parenty <- alllist[[parname]]$y
    np <- length(parentx)
    mu02use <- mu0vec[offname]
    h2use <- hvec[offname]

    find <- which(factortype == offname)

    csize <- rpois(np,mu02use)
    numoff <- sum(csize)
    x0 <- rep.int(parentx,csize)
    y0 <- rep.int(parenty,csize)

    lll <- rkern(1,x0,y0,h2use)
    xoff <- c(lll$Bx)
    yoff <- c(lll$By)

    retain <- spatstat.geom::inside.owin(xoff,yoff,win)

    xoff <- xoff[retain]
    yoff <- yoff[retain]
    nn <- length(xoff)
    alllist[[offname]] <- list("x"=xoff,"y"=yoff,"marks"=rep(factortype[find],nn))
  }

  if(nxtra > 0){
    for(i in 1:nxtra){
      xname <- names_xtra[i]
      lambdax2use <- lambdaxtra[i]
      Y <- spatstat.random::rpoispp(lambda2use,lamx=NULL,
                                    win=win,
                                    nsim=1,drop=TRUE,ex=NULL,warnwin=TRUE)
      find <- which(factortype == xname)
      alllist[[xname]] <- list("x"=Y$x,"Y"=Y$y, "marks"=rep(factortype[find],Y$n))
    }
  }

  allx <- unlist(lapply(alllist,"[[","x"),use.names=FALSE)
  ally <- unlist(lapply(alllist,"[[","y"),use.names=FALSE)
  allmarks <- unlist(lapply(alllist,"[[","marks"),use.names=FALSE)

  Xf <- spatstat.geom::ppp(x=allx,y=ally,marks=allmarks, window=win)

  ## Permuting to make it look randomize; not needed
  permu <- sample(Xf$n)
  Xf <- Xf[permu]


  # X <- spatstat.random::rpoispp(lambdapar*(1-exp(-max(mu0vec))),win=win)
  # parentx <- X$x; parenty <- X$y; np <- X$n
  # X <- spatstat.geom::setmarks(X,factortype[1])

  ## Create the offspring process
  ## Specify parent locations in input for the lambda function
  ## Superimpose on the already created marked point process
  # for (i in 1:noff) {
  #
  #   csize <- rpois(np,mu0vec[i])
  #   numoff <- sum(csize)
  #   x0 <- rep.int(parentx,csize)
  #   y0 <- rep.int(parenty,csize)
  #
  #   lll <- rkern(1,x0,y0,hvec[i])
  #   xoff <- c(lll$Bx)
  #   yoff <- c(lll$By)
  #
  #   ## dd <- matrix(rnorm(2*numoff,0,hvec[i]),ncol=2)
  #   ## xy <- xy.coords(dd)
  #   ## dx <- xy$x
  #   ## dy <- xy$y
  #   ## xoff <- x0 + dx
  #   ## yoff <- y0 + dy
  #
  #   retain <- spatstat.geom::inside.owin(xoff,yoff,win)
  #
  #   xoff <- xoff[retain]
  #   yoff <- yoff[retain]
  #
  #   Y <- spatstat.geom::ppp(xoff,yoff,window=win)
  #
  #   Y <- spatstat.geom::setmarks(Y,factortype[1+i])
  #   X <- spatstat.geom::superimpose(X, Y, W = X$window, check = FALSE)
  # }
  ## Generate from other processes, if any
  # if(nxtra > 0){
  #   for(i in 1:nxtra){
  #     Y <- spatstat.random::rpoispp(lambdaxtra[i],lamx=NULL,
  #                                   win=win,
  #                                   nsim=1,drop=TRUE,ex=NULL,warnwin=TRUE)
  #     Y <- spatstat.geom::setmarks(Y,factortype[1+noff+i])
  #     X <- spatstat.geom::superimpose(X, Y, W = X$window, check = FALSE)
  #   }
  # }

  ## Permuting to make it look randomize; not needed
  # permu <- sample(X$n)
  # X <- X[permu]

  return(Xf)
}
