#' @name check_within_window
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
#' @noRd

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


#' @name getbignum2.norm
#' @title Internal function to compute the sum of the log terms for Gaussian kernel
#'
#' @param Y_x numeric vector, x coordinates of the offspring process
#' @param Y_y numeric vector, y coordinates of the offspring process
#' @param C_x numeric vector, x coordinates of the parent process
#' @param C_y numeric vector, y coordinates of the parent process
#' @param h positive scalar, bandwidth
#'
#' @useDynLib MCPP, .registration=TRUE
#'
#' @return real scalar
#' @noRd

getbignum2.norm <- function(Y_x,Y_y,C_x,C_y,h){
  out <-  .Call("bignum2func_norm",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y),as.numeric(h))
  return(out)
}

#' @name getbignum2.cauchy
#' @title Internal function to compute the sum of the log terms for Cauchy kernel
#'
#' @param Y_x numeric vector, x coordinates of the offspring process
#' @param Y_y numeric vector, y coordinates of the offspring process
#' @param C_x numeric vector, x coordinates of the parent process
#' @param C_y numeric vector, y coordinates of the parent process
#' @param h positive scalar, bandwidth
#'
#' @useDynLib MCPP, .registration=TRUE
#'
#' @return real scalar
#' @noRd

getbignum2.cauchy <- function(Y_x,Y_y,C_x,C_y,h){
  out <-  .Call("bignum2func_cauchy",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y),as.numeric(h))
  return(out)
}

#' @name getbignum2.unif
#' @title Internal function to compute the sum of the log terms for Uniform Kernel
#'
#' @param Y_x numeric vector, x coordinates of the offspring process
#' @param Y_y numeric vector, y coordinates of the offspring process
#' @param C_x numeric vector, x coordinates of the parent process
#' @param C_y numeric vector, y coordinates of the parent process
#' @param h positive scalar, bandwidth
#'
#' @useDynLib MCPP, .registration=TRUE
#'
#' @return real scalar
#' @noRd

getbignum2.unif <- function(Y_x,Y_y,C_x,C_y,h){
  out <-  .Call("bignum2func_unif",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y),as.numeric(h))
  return(out)
}

#' @name gethpars
#' @title Internal function to get initial values for the bandwidth parameters
#'
#' @param Y_x numeric vector, x coordinates of the offspring process
#' @param Y_y numeric vector, y coordinates of the offspring process
#' @param C_x numeric vector, x coordinates of the parent process
#' @param C_y numeric vector, y coordinates of the parent process
#'
#' @useDynLib MCPP, .registration=TRUE
#'
#' @return positive scalar for initial value of h
#' @noRd

gethpars <- function(Y_x,Y_y,C_x,C_y){
  hsamp <- .Call("gethsamp",as.numeric(Y_x),as.numeric(Y_y),as.numeric(C_x),as.numeric(C_y))
  hbar <- mean(hsamp)

  return(hbar)
}


#' @name subdata
#' @title Internal function to subset the data based on genus name
#'
#' @param dat ppp class object, the entire dataset
#' @param genus character, name of the genus to be subsetted
#'
#' @return ppp class object, subset of the data with only the specified genus
#' @noRd

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

#' @name rkern.norm
#' @title Internal function for generating from a Gaussian kernel
#'
#' @param B positive integer, number of samples to draw for each parent location
#' @param Cx real vector, x coordinates of the centers
#' @param Cy real vector, y coordinates of the centers
#' @param h positive real scalar, bandwidth parameter
#'
#' @import stats
#' @return list with two matrices for the x and y coordinates of the generated points
#' @noRd

rkern.norm <- function(B,Cx,Cy,h){
  nCs <- length(Cx)
  Bx <- Cx + h*matrix(rnorm(nCs*B),nCs,B)
  By <- Cy + h*matrix(rnorm(nCs*B),nCs,B)
  out <- list("Bx"=Bx,"By"=By)
  return(out)
}

#' @name rkern.cauchy
#' @title Internal function for generating from a Cauchy kernel
#'
#' @param B positive integer, number of samples to draw for each parent location
#' @param Cx real vector, x coordinates of the centers
#' @param Cy real vector, y coordinates of the centers
#' @param h positive real scalar, bandwidth parameter
#'
#' @import stats
#' @return list with two matrices for the x and y coordinates of the generated points
#' @noRd

rkern.cauchy <- function(B,Cx,Cy,h){
  nCs <- length(Cx)
  Z <- sqrt(matrix(rchisq(nCs*B),nCs,B))
  Bx <- Cx + (h/Z)*matrix(rnorm(nCs*B),nCs,B)
  By <- Cy + (h/Z)*matrix(rnorm(nCs*B),nCs,B)
  out <- list("Bx"=Bx,"By"=By)
  return(out)
}

#' @name rkern.unif
#' @title Internal function for generating from a Uniform kernel
#'
#' @param B positive integer, number of samples to draw for each parent location
#' @param Cx real vector, x coordinates of the centers
#' @param Cy real vector, y coordinates of the centers
#' @param h positive real scalar, bandwidth parameter
#'
#' @import stats
#' @return list with two matrices for the x and y coordinates of the generated points
#' @noRd

rkern.unif <- function(B,Cx,Cy,h){
  nCs <- length(Cx)
  R <- h*matrix(runif(nCs*B),nCs,B)
  theta <- 2*pi*matrix(runif(nCs*B),nCs,B)
  Bx <- Cx + R*cos(theta)
  By <- Cy + R*sin(theta)
  out <- list("Bx"=Bx,"By"=By)
  return(out)
}
