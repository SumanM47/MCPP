#' @name multiprocppp
#' @title R function to generate from multiprocess poisson point process
#'
#' @param win owin class object, the observation window
#' @param lambdapar positive scalar, the intensity parameter for the parent process
#' @param name_par optional character, name of the parent
#' @param noff positive scalar, number of offsprings
#' @param mu0vec numeric vector of positive numbers, vector of offspring densities
#' @param hvec numeric vector of positive numbers, vector of bandwidth parameters
#' @param kern character, choice of kernel function. Permitted values are "Gaussian", "Uniform" or "Cauchy"
#' @param names_off character vector, names of the offspring entities
#' @param nxtra non-negative integer, number of unrelated entities
#' @param lambdaxtra scalar or vector of positive numbers, intensities of the unrelated entities processes
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
                         lambdapar, name_par=NULL,
                         noff,mu0vec,hvec,kern="Gaussian",names_off=NULL,
                         nxtra=0,lambdaxtra,names_xtra=NULL){
  ## Initial Checks
  if(lambdapar <= 0 | length(lambdapar) != 1) stop("lambdapar must be a positive scalar")
  if(noff <= 0) stop("noff must be a positive integer")
  if(length(mu0vec)!=1 & length(mu0vec)!=noff) stop("mu0vec must be a positive scalar or vector of length noff")
  if(length(hvec)!=1 & length(hvec)!=noff) stop("hvec must be a positive scalar or vector of length noff")
  if(sum(mu0vec<=0)!=0) stop("mu0vec must be positive")
  if(sum(hvec<=0)!=0) stop("hvec must be positive")
  if(!spatstat.geom::is.owin(win)) stop("win must be object of class owin")
  if(nxtra < 0) stop("nxtra must be a non-negative integer")
  if(nxtra > 0 & missing(lambdaxtra)) stop("lambdaxtra is missing. Provide a positive scalar or vector of length nxtra")
  if(nxtra > 0 & (length(lambdaxtra)!=1 & length(lambdaxtra)!=nxtra))  stop("Provide a positive scalar or vector of length nxtra")
  if(nxtra > 0 & (sum(lambdaxtra <=0)!=0)) stop("lambdaxtra must be positive")
  if(!is.null(name_par) & (length(name_par)!=1 | !is.character(name_par))) stop("name_par must a character string of size 1")
  if(!is.null(names_off) & (length(names_off)!=noff & !is.character(names_off))) stop("names_off must be a character vector of size noff")
  if(nxtra > 0 & !is.null(names_xtra) & (length(names_xtra)!=nxtra | !is.character(names_xtra))) stop("names_xtra must be a character vector of size nxtra")

  ## Setting labels if not already given
  if(is.null(name_par)) name_par <- "parent"
  if(is.null(names_off)){
    names_off <- paste0("offspring",1:noff,sep="")
  }
  if(nxtra > 0 & is.null(names_xtra)){
    names_xtra <- paste0("xtra",1:nxtra,sep="")
  }

  ## Managing inputs
  if(length(mu0vec)==1) mu0vec <- rep(mu0vec,noff)
  if(length(hvec)==1) hvec <- rep(hvec,noff)
  if(nxtra > 0 & length(lambdaxtra)==1) lambdaxtra <- rep(lambdaxtra,nxtra)

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

  #require('spatstat')

  ## Creating a factor of labels for marking
  types <- c(name_par,names_off,names_xtra)
  factortype <- factor(types, levels = types)

  ## First create the parent process
  X <- spatstat.random::rpoispp(lambdapar*(1-exp(-max(mu0vec))),win=win)
  parentx <- X$x; parenty <- X$y; np <- X$n
  X <- spatstat.geom::setmarks(X,factortype[1])

  ## Create the offspring process
  ## Specify parent locations in input for the lambda function
  ## Superimpose on the already created marked point process
  for (i in 1:noff) {

    csize <- rpois(np,mu0vec[i])
    numoff <- sum(csize)
    x0 <- rep.int(parentx,csize)
    y0 <- rep.int(parenty,csize)

    lll <- rkern(1,x0,y0,hvec[i])
    xoff <- c(lll$Bx)
    yoff <- c(lll$By)

    ## dd <- matrix(rnorm(2*numoff,0,hvec[i]),ncol=2)
    ## xy <- xy.coords(dd)
    ## dx <- xy$x
    ## dy <- xy$y
    ## xoff <- x0 + dx
    ## yoff <- y0 + dy

    retain <- spatstat.geom::inside.owin(xoff,yoff,win)

    xoff <- xoff[retain]
    yoff <- yoff[retain]

    Y <- spatstat.geom::ppp(xoff,yoff,window=win)

    Y <- spatstat.geom::setmarks(Y,factortype[1+i])
    X <- spatstat.geom::superimpose(X, Y, W = X$window, check = FALSE)
  }
  ## Generate from other processes, if any
  if(nxtra > 0){
    for(i in 1:nxtra){
      Y <- spatstat.random::rpoispp(lambdaxtra[i],lamx=NULL,
                   win=win,
                   nsim=1,drop=TRUE,ex=NULL,warnwin=TRUE)
      Y <- spatstat.geom::setmarks(Y,factortype[1+noff+i])
      X <- spatstat.geom::superimpose(X, Y, W = X$window, check = FALSE)
    }
  }

  ## Permuting to make it look randomize; not needed
  permu <- sample(X$n)
  X <- X[permu]

  return(X)
}
