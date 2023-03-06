#' @title Internal function to check if a set of points is within a given window
#'
#' @param x numeric vector, x coordinates of the points
#' @param y numeric vector, y coordinates of the points
#' @param w owin class object, observation window
#'
#' @import sp
#' @import spatstat.geom
#'
#' @return numeric vector of ones and zeros with
#'  one meaning the point is within the window

check_within_Window <- function(x,y,w){
  if(w$type == "rectangle"){
    xb <- w$xrange; xlb <- min(xb); xub <- max(xb)
    yb <- w$yrange; ylb <- min(yb); yub <- max(yb)
    out <- as.vector((x <= xub & x >= xlb)*(y <= yub & y >= ylb))
  }
  if(w$type=="polygonal"){
    #require(sp)
    polyx <- w$bdry[[1]]$x
    polyy <- w$bdry[[1]]$y
    #dimx <- dim(x)
    ind <- sp::point.in.polygon(x,y,polyx,polyy)
    out <- as.numeric(ind>0)
  }
  return(out)
}


#' @title Internal function to compute the sum of the log terms
#'
#' @param Y_x numeric vector, x coordinates of the offspring process
#' @param Y_y numeric vector, y coordinates of the offspring process
#' @param C_x numeric vector, x coordinates of the parent process
#' @param C_y numeric vector, y coordinates of the parent process
#' @param h positive scalar, bandwidth
#'
#' @useDynLib MCPP bignum2func
#'
#' @return real scalar

getbignum2 <- function(Y_x,Y_y,C_x,C_y,h){
  out <-  .Call("bignum2func",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y),as.numeric(h))
  return(out)
}

#' @title Internal function to get initial values for the bandwidth parameters
#'
#' @param Y_x numeric vector, x coordinates of the offspring process
#' @param Y_y numeric vector, y coordinates of the offspring process
#' @param C_x numeric vector, x coordinates of the parent process
#' @param C_y numeric vector, y coordinates of the parent process
#'
#' @useDynLib MCPP gethsamp
#'
#' @return positive scalar for initial value of h

gethpars <- function(Y_x,Y_y,C_x,C_y){
  hsamp <- .Call("gethsamp",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y))
  hbar <- mean(hsamp)

  return(hbar)
}


#' @title Internal function to subset the data based on genus name
#'
#' @param dat ppp class object, the entire dataset
#' @param genus character, name of the genus to be subsetted
#'
#' @return ppp class object, subset of the data with only the specified genus

subdata <- function(dat, genus){
  dd <- dat
  ind <- which(dat$marks==as.character(genus))
  dd$n <- length(ind)
  dd$x <- dat$x[ind]
  dd$y <- dat$y[ind]
  dd$markformat <- "none"
  dd$marks <- NULL

  return(dd)
}


#' @title R function to draw posterior samples from MCPP
#'
#' @param obj ppp class object, the dataset
#' @param parent_genus character vector, scalar if all parent types are same or a vector of name of the parent genuses
#' @param offspring_genus character vector, name of the offspring genuses
#' @param bdtype numeric scalar, can take values 1 or 2. 1 keeps the window as is; 2 creates a convex hull and uses that as the window
#' @param B positive integer, number of Monte Carlo samples to use to approximate the integral term
#' @param jitter logical, if TRUE, the initial value for the bandwidth parameter is jittered. Useful for running multiple chains
#' @param hyperParams named list of hyperparameters needed for the MCPP model, components are, lambda (two-vector of shape and rate hyperparamaters for the Gamma prior to the parent intensities), alpha (two-vector of shape and rate hyperparamaters for the Gamma prior to the alpha parameter of the offspring densities), lambdaO (two-vector of shape and rate hyperparamaters for the Gamma prior to the unrelated taxa intensities), h (scalar proportion of the maximum window distance to be set equal to the 99th percentile of the Half-normal prior to the bandwidth parameter)
#' @param startValues named list of starting values for the alpha and h parameters, components are, alpha (positive vector of initial values for the alpha parameter, same size as the number of offspring taxa) and h (positive vector of initial values for the bandwidth parameter, same size as the number of offspring taxa). Default is NULL in which case initial values are computed
#' @param mcmcParams named list of parameters required for MCMC sampling, components are, run (named list with entries numReps: total number of scans, thin: the thinning interval length and burninPerc: percentage of scans to be considered as burn-in), storage (ignored for MCPP) and tuning (ignored for MCPP)
#' @param store_res logical, indicator to store results periodically
#' @param path path to the storage file where the periodical results are to be saved
#'
#' @import spatstat.geom
#' @import sp
#'
#' @useDynLib MCPP gethsamp
#' @useDynLib MCPP bignum2func
#'
#' @return list of posterior samples for the corresponding parameters, object of class MCPP
#'
#' @export

MCPP_run <- function(obj,parent_genus,offspring_genus,
                             bdtype=1,B=100,
                             jitter=FALSE,
                             hyperParams=list("lambda"=c(0.01,0.01),
                                              "alpha"=c(0.01,0.01),
                                              "lambdaO"=c(0.01,0.01),
                                              "h"=0.05),
                             startValues=NULL,
                             mcmcParams=list("run"=list(numReps=20000,
                                                        thin=1,
                                                        burninPerc=0.5),
                                             "storage"=NULL,
                                             "tuning"=NULL),
                             store_res=FALSE,path=NULL){
  if(store_res & is.null(path)){stop("Please specify intermediary output storage file")}
  # require(fields)
  # require(spatstat)

  ## Create window according to boundary type
  ## bdtype=1 (default) is to use as is; 2 and 3 are convex and concave boundaries

  if(bdtype==2){
    cvxhull <- rev(chull(obj$x,obj$y))
    wcvx <- owin(poly=list(x=obj$x[cvxhull],y=obj$y[cvxhull]))
    obj$window <- wcvx
  }

  ## Parse hyperParam Inputs
  al <- hyperParams$lambda[1]
  bl <- hyperParams$lambda[2]
  am <- hyperParams$alpha[1]
  bm <- hyperParams$alpha[2]
  ao <- hyperParams$lambdaO[1]
  bo <- hyperParams$lambdaO[2]
  hclimp <- hyperParams$h

  ## Parse mcmcParams Inputs
  iters <- mcmcParams$run$numReps
  burn <- round(iters*mcmcParams$run$burninPerc)
  thin <- mcmcParams$run$thin

  ## Get Window Info
  W <- obj$window
  xmax <- max(W$xrange)
  xmin <- min(W$xrange)

  ymax <- max(W$yrange)
  ymin <- min(W$yrange)
  maxd <- sqrt((xmax-xmin)^2 + (ymax-ymin)^2)
  hclim = hclimp*maxd;
  hsd <- hclim/qnorm(0.995)
  W_area <- spatstat.geom::area(W)

  ##Sorting out obj Input

  alltaxa <- unique(obj$marks)
  atn <- length(alltaxa)

  not <- length(offspring_genus)
  np <- length(parent_genus)
  if(np!=1 & np!=not){stop("Number of parent taxa must be one or same as the number of offspring taxa")}
  notp <- length(unique(c(parent_genus,offspring_genus)))
  nxt <- atn - notp

  up <- unique(as.character(parent_genus))
  upno <- setdiff(up,as.character(offspring_genus))
  nupno <- length(upno)
  upnoid <- sapply(upno,function(v){min(which(parent_genus==v))})

  if(length(parent_genus)==1){parent_genus <- rep(parent_genus,not)}

  Cpart_x <- Cpart_y <- list(not)
  Ypart_x <- Ypart_y <- list(not)
  n_Ypart <- rep(0,not); n_Cpart <- rep(0,not)

  ## Subsetting the data to parent, offspring and other taxa info
  for(ind in 1:not){
    parent_dat <- subdata(obj,as.character(parent_genus[ind]))
    Cpart_x[[ind]] <- parent_dat$x; Cpart_y[[ind]] <- parent_dat$y; n_Cpart[ind] <- parent_dat$n
  }


  ## List of offspring subdata
  for(ind in 1:not){
    od <- subdata(obj,as.character(offspring_genus[ind]))
    # offspring_datall[[ind]] <- od
    Ypart_x[[ind]] <- od$x; Ypart_y[[ind]] <- od$y; n_Ypart[ind] <- od$n
  }


  ## List of Other subdata
  if(nxt > 0){
    xtra_genus <- setdiff(alltaxa,union(parent_genus,offspring_genus))
    Opart_x <- Opart_y <- list(nxt)
    n_Opart <- rep(0,nxt)

    for (ind in 1:nxt){
      xd <- subdata(obj,as.character(xtra_genus[ind]))
      # offspring_datall[[ind]] <- od
      Opart_x[[ind]] <- xd$x; Opart_y[[ind]] <- xd$y; n_Opart[ind] <- xd$n
    }
  }


  hvec <- rep(0,not)
  if(!is.null(startValues$h)){hvec <- startValues$h}
  Bignum1vec <- Bignum2vec <- rep(0,not)
  for(ind in 1:not){
    Y_x <- Ypart_x[[ind]]; Y_y <- Ypart_y[[ind]]; n_Ys <- n_Ypart[ind]
    C_x <- Cpart_x[[ind]]; C_y <- Cpart_y[[ind]]; n_Cs <- n_Cpart[ind]
    if(is.null(startValues$h)){hs <- gethpars(Y_x,Y_y,C_x,C_y)}else{hs <- hvec[ind]}
    if(jitter){
      hs <- jitter(hs)
    }
    hvec[ind] <- hs
    Bigmat_x <- C_x + hs*matrix(rnorm(n_Cs*B),n_Cs,B); Bigmat_y <- C_y + hs*matrix(rnorm(n_Cs*B),n_Cs,B)
    Bignum1vec[ind] <- sum(check_within_Window(Bigmat_x,Bigmat_y,W))/B

    Bignum2vec[ind] <- getbignum2(Y_x,Y_y,C_x,C_y,hs)
  }

  if(is.null(startValues$alpha)){mu0vec <- n_Ypart[1:not]/Bignum1vec}else{mu0vec <- startValues$alpha}

  ll.lh <- -mu0vec*Bignum1vec + Bignum2vec

  ## Bookkeping
  v.h <- 0.0075
  v.beta <- 0.5

  ## Storage
  keep.lambdaC <- matrix(0,(iters-burn)/thin,nupno)
  keep.mu0 <- matrix(0,(iters-burn)/thin,not)
  keep.h <- matrix(0,(iters-burn)/thin,not)
  if(nxt > 0){
    keep.lambdaO <- matrix(0,(iters-burn)/thin,nxt)
  }
  if(store_res){
    temp.h <- temp.mu0 <- matrix(0,10000,not); temp.lambdaC <- matrix(0,10000,nupno)
    if(nxt > 0){
      temp.lambdaO <- matrix(0,10000,nxt)
    }
    store.h <- store.mu0 <- store.lambdaC <- store.lambdaO <- NULL
  }

  pratio <- xtra_ep <- 0

  ## GO!
  for(i in 1:iters){



    ## Update lambdaC
    lambdaC <- rgamma(nupno,al + n_Cpart[upnoid], bl + W_area)

    if(nxt > 0){
      ## Update lambdaO
      lambdaO <- rgamma(nxt,ao + n_Opart, bo + W_area)
    }


    ## Update alpha_j
    mu0vec <- rgamma(not,am + n_Ypart, bm + Bignum1vec)

    ## Update h_j
    for(j in 1:not){
      Y_x <- Ypart_x[[j]]; Y_y <- Ypart_y[[j]]; n_Ys <- n_Ypart[j]
      C_x <- Cpart_x[[j]]; C_y <- Cpart_y[[j]]; n_Cs <- n_Cpart[j]
      Bignum1s <- Bignum1vec[j]; Bignum2s <- Bignum2vec[j]
      hs <- hvec[j]

      ll.lhs <- -mu0vec[j]*Bignum1s + Bignum2s
      ll.lh[j] <- ll.lhs

      can_hs <- hs + v.h*rnorm(1)

      if(can_hs > 0){
        pratio <- -0.5*((can_hs/hsd)^2 - (hs/hsd)^2)
        Bigmat_x <- C_x + can_hs*matrix(rnorm(n_Cs*B),n_Cs,B); Bigmat_y <- C_y + can_hs*matrix(rnorm(n_Cs*B),n_Cs,B)
        can_Bignum1 <- sum(check_within_Window(Bigmat_x,Bigmat_y,W))/B

        can_Bignum2 <- getbignum2(Y_x,Y_y,C_x,C_y,can_hs)

        can_ll.lhs <- -mu0vec[j]*can_Bignum1 + can_Bignum2

        a_lhs <- can_ll.lhs - ll.lhs + pratio
        if(log(runif(1)) < a_lhs){
          Bignum1vec[j] <- can_Bignum1
          Bignum2vec[j] <- can_Bignum2
          hvec[j] <- can_hs
          ll.lh[j] <- can_ll.lhs
        }
      }
    }


    if(i > burn & (i-burn)%%thin == 0){
      index <- (i-burn)/thin
      keep.lambdaC[index,] <- lambdaC
      keep.mu0[index,] <- mu0vec
      keep.h[index,] <- hvec
      if(nxt > 0){
        keep.lambdaO[index,] <- lambdaO
      }
    }
    if(store_res){
      ind2 <- 1 + ((i-1)%%10000)
      temp.lambdaC[ind2,] <- lambdaC
      temp.mu0[ind2,] <- mu0vec
      temp.h[ind2,] <- hvec
      if(nxt > 0){temp.lambdaO[ind2,] <- lambdaO}
      if(i%%10000==0){
        store.lambdaC <- rbind(store.lambdaC,temp.lambdaC)
        store.mu0 <- rbind(store.mu0,temp.mu0)
        store.h <- rbind(store.h,temp.h)
        if(nxt > 0){
          store.lambdaO <- rbind(store.lambdaO,temp.lambdaO)
        }
        temp.h <- temp.mu0 <- matrix(0,10000,not); temp.lambdaC <- matrix(0,10000,nupno)
        if(nxt > 0){
          temp.lambdaO <- matrix(0,10000,nxt)
        }
        lst <- list("alllambdaC" = store.lambdaC,"allmu0"=store.mu0,"allh"=store.h,"alllambdaO"=store.lambdaO)
        save(lst,file=path)
      }
    }

  }
  if(nxt==0){keep.lambdaO <- NULL}
  setup <- list()
  setup$nreps <- iters
  setup$burn <- burn
  setup$thin <- thin
  setup$popair <- cbind(parent_genus,offspring_genus)
  setup$others <- NULL
  if(nxt>0){setup$others <- xtra_genus}
  lst <- list("lambdaC"=keep.lambdaC,"mu0"=keep.mu0,"h"=keep.h,"lambdaO"=keep.lambdaO,"setup" <- setup)
  class(lst) <- "MCPP"
  return(lst)
}
