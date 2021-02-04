#' Calculate winding numbers.
#'
#' The winding numbers (u,v) for a pair (m,n) satisfy mv-nu=Delta for Delta=+/- 1, and if 1<m,n they are in the range
#' 0<u<m, 0<v<n. windingNumbers returns 0,1 when called with 1,n and m-1,0 when called with m,1.
#' Winding numbers can only be found if m and n are coprime and the function throws an error if not.
#' If
#' @param m,n Positive integer
#' @param Delta Plus or minus 1.
#' @export
windingNumbers <- function(m,n,Delta=1) {

  if (m<1 || n < 1) {stop("arguments must be strictly positive")}
  if (hcf(m,n)!=1) {stop(sprintf("%g and %g are not coprime",m,n))}

  if (n<m) {
    vu <- windingNumbers(n,m,-Delta)
    uv <- c(vu[2],vu[1])
    return(uv)
  } else {
    uv <- if (m==1 && n==1 ) {
      c(1/2,1/2)
    } else if (m==1) {
      if(Delta==1) {c(0,1)} else {c(1,n-1)}
    } else if (Delta==1) {
      vu <- jCongruence(m,n) # returns the wrong way round...
      c(vu[2],vu[1])
    } else if (Delta==-1) {
      uv <- jCongruence(n,m)
    }
    return(uv)
  }
}

#' The sign of the area of the generating parallelogram
#'
#' Interpreting this as an area assumes that the rise is equal to one
#' @param m,n Integer. Orders of the generating parastichy
#' @param d Divergence
#' @export
Deltamnd <- function(m,n,d) {
  if (m==0) return(n)
  if (m!=n) return( round(n*d)*m - round(m*d)*n)
  if (m==1 & n==1 & d==1/2) return(1)
}



#' The range of possible divergences for a parastichy pair of a given type
#'
#' If a parastichy pair is generating or opposed there is a constraint on the ange of possible divergences
#' which this function returns
#' @param m,n Integer parastichy orders
#' @param Delta +/- 1
#' @param scale If "unity", the divergence (between 0 and 1); if "twomn" the divergence multiplied by 2*m*n
#' @param type One of "generating", "pseudogenerating", "opposed", or "generatingopposed"
#' @export
divergenceInterval <- function(m,n,Delta,scale="unity",type="generating") {
  interval <- switch(type,
                "generating" = {
                  interval <- dGeneratingInterval(m,n,Delta)$Interval
                  if (scale=="twomn") interval <- interval *  (2 * m * n);
                  interval
                },
                "opposed" = {
                  interval <- dOpposedInterval(m,n,"twomn") # returns on 2mn scale
                  if (scale=="unity") interval <- interval/ (2 * m * n);;
                  interval
                },
                "pseudogenerating" = {
                  interval <- dPseudogeneratingInterval(m,n)
                  if (scale=="unity") interval <- interval/ (2 * m * n);
                  interval
                },
                "generatingopposed" = {
                  interval <- dGeneratingOpposedInterval(m,n,Delta)
                  if (scale=="twomn") interval <- interval *  (2 * m * n);
                  interval
                }

  )
  interval
}

###### probably better reimplemented as computing the interval for each type directly then taking intersections

dGeneratingInterval <- function(m,n,Delta) {
  uv <- windingNumbers(m,n,Delta)
  u <- uv[1];v <- uv[2]
  Lm <- (u-1/2)/m; Rm <- (u+1/2)/m
  Ln <- (v-1/2)/n; Rn <- (v+1/2)/n
  if (n<= m-2) {
    Itype <- "LmRm"
  } else if (n==m-1 || n==m+1) {
    if (Delta==1) {
      Itype <- "LnRm"
    } else {
      Itype <- "LmRn"
    }
  } else if (n==m) {
    stopifnot (n==1)
    Itype <- "half"
  } else if (n>=m+2) {
    Itype <- "LnRn"
  } else{ stop(sprintf("%d %d \n",m,n))}
  intervalL <- NA;intervalR <- NA
  if (Itype == "half") { intervalL <- 1/2; intervalR <- 1/2 }
  if (Itype =="LnRm" || Itype =="LnRn") intervalL <- Ln
  if (Itype =="LmRm" || Itype =="LmRn") intervalL <- Lm
  if (Itype =="LnRm" || Itype =="LmRm") intervalR <- Rm
  if (Itype =="LnRn" || Itype =="LmRn") intervalR <- Rn
  return(list(Itype=Itype,Interval=c(intervalL,intervalR)))
}
# GeneratingIntervalDelta <- function(m,n,Delta) {





dGeneratingOpposedInterval <- function(m,n,Delta) {
  if (m > n) {
    interval <- dGeneratingOpposedInterval (n,m,Delta)
  } else if (m==1 & n==1) {
    interval <- c(1/2,1/2)
  } else if (m==1 & n>1) {
    if (Delta==1) {
      interval <- c(1/(2*n),1/n)
    } else {
      interval <- 1 - c(1/(2*n),1/n)
    }
  } else if (1 < m & m < n) {
    uv <- windingNumbers(m,n,Delta)
    u <- uv[1];v <- uv[2]
    interval <- c(u/m,v/n)
  } else {
    stop(sprintf("%d %d \n",m,n))
  }
  if(interval[2]<interval[1]) interval <- rev(interval)
  return(interval)
}

dOpposedInterval <- function(m,n,type="twomn") {
  if (m==1 & n==1) {
    intervalBoundaries <- c(c(-1/2,-1/2),c(1/2,1/2))
  } else {
    if (hcf(m,n) != 1) {stop("Can't compute opposed intervals for non coprime m,n")}
    twomn <- 2 * m * n
    # we work with a scaled d'=2mnd,  0< d' < 2mn
    # first construct for d<1/2
    mchanges <- n*seq(from=1,to= m )
    nchanges <- m*seq(from=1,to= n-1 ) # because m and n coprime mchanges and nprimes distinct until mn

    bothchanges <- sort(c(mchanges,nchanges)) # there are m+n-1 of these
    if ( (m+n -1) %% 2 ==1) {
      # the final boundary is the start of a nonopposed
      bothchanges <- bothchanges[-length(bothchanges)]
    }
    # now replicate into d>1/2
    intervalBoundaries <- c(bothchanges,rev(twomn-bothchanges))
}
  #     intervalBoundaries <- integer()
#     while( length(c(mchanges,nchanges))>0 ) {
#       nextmchange <- mchanges[1]
#       nextnchange <- nchanges[1]
#       if (nextmchange==nextnchange) {
#         mchanges <- mchanges[-1]
#         nchanges <- nchanges[-1] # if they both change sign together, no change of opposition
#       } else if (nextmchange<nextnchange) {
#         mchanges <- mchanges[-1]
#         intervalBoundaries <- c(intervalBoundaries,nextmchange)
#       } else {
#         nchanges <- nchanges[-1]
#         intervalBoundaries <- c(intervalBoundaries,nextnchange)
#       }
#     }
#  }
  # back to 0,1
  stopifnot(length(intervalBoundaries)%%2==0)
  #intervalBoundaries <- split(intervalBoundaries,rep(1:(length(intervalBoundaries)%/%2),each=2))
  interval <- matrix(intervalBoundaries,ncol=2,byrow=TRUE)

  return(interval)

}


dPseudogeneratingInterval <- function(m,n,type="twomn") {
  interval <- c(0,min(m,n))
  return(interval)
}


#' Return the Delta for which the generating interval (m,n,Delta) is in 0<d<1/2
#' @param m,n Integer
DeltaForLHSd <- function(m,n) { # force the choice of 0<d<1/2
  Delta <- 1
  gi <- dGeneratingInterval(m,n,1)
  if (gi$Interval[1] >= 0.5 ) {
    Delta <- -1
  }
  Delta
}

#' Return the interval for which m,n is generating in 0<d<1/2
#'
#' @param m,n Integer
dGeneratingIntervalLHS <- function(m,n) {
  if (hcf(m,n)!=1) return(NA)
  gi <- dGeneratingInterval(m,n,1)
  if (gi$Interval[1] >= 0.5 ) {
    gi <- dGeneratingInterval(m,n,-1)
  }
  return(gi$Interval)
}

#' Return the interval for which m,n is generating and opposed in 0<d<1/2
#'
#' @param m,n Integer
dGeneratingOpposedIntervalLHS <- function(m,n) {
  if (hcf(m,n)!=1) return(NA)
  gi <- dGeneratingOpposedInterval(m,n,1)
  if (gi[1] >= 0.5 ) {
    gi <- dGeneratingOpposedInterval(m,n,-1)
  }
  return(gi)
}


#' Check if a parastichy pair is opposed
#'
#' Used by test suite
#'  @param m,n Integer. Orders of a generating parastichy pair
#' @param d Divergence of the lattice. 0<d<1 if scale=unity
#' @param scale twomn or unity
isOpposed <- function(m,n,d,scale="twomn") {
  if (scale=="twomn") { d <- d/(2*m*n)}
  xm <- m * d - round(m*d)
  xn <- n * d - round(n*d)
  return(sign(xm * xn)<0)
}

#' Cylindrical coordinate
#'  @param m Integer.
#' @param d Divergence of the lattice. 0<d<1 if scale=unity
#' @param scale twomn or unity
cylindricalX <- function(m,d,scale="unity") {
  if (scale=="twomn") { d <- d/(2*m*n)}
  xm <- m * d - round(m*d)
  return(xm)
}


