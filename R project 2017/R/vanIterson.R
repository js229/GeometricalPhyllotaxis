# # what is this?
# PP <- function(h,d) {
#   PGhd <- newPhyllotaxisGenetic(Rise=h,Divergence=2*pi*d,Rho=(1/(2*pi)),Jugacy=1,L=1,lhsOrigin=FALSE)
#   GPPrincipalParastichyNumbers(PGhd)
# }

#' The square of the rise h for a van Iterson diagram
#' @param m,n Integer. Orders of a generating parastichy pair
#' @param d Divergence of the lattice. 0<d<1.
vanItersonhSquared <- function(m,n,d) {
	dm <- m*d - round(m*d)
	dn <- n*d - round(n*d)
	h2 <- (dm^2 - dn^2)/ (n^2-m^2)
	h2[h2<0] <- NA
	h2
}


#' The rise for a given divergence in a van Iterson curve for m,n
#' @param m,n Integer. Orders of a generating parastichy pair
#' @param d Divergence of the lattice. 0<d<1.
#' @export
hvanIterson <- function(m,n,d) {
  h2 <- vanItersonhSquared(m,n,d)
  h <- sqrt(h2)
  h
}

#' The divergence-rise curve in a van Iterson diagram for m,n
#' @param m,n Integer. Orders of a generating parastichy pair
#' @param d Divergence of the lattice. 0<d<1.
#' @export
dhvanIterson <- function(m,n,ontype,lprtype) {
  range <- getdrange(m,n,ontype,lprtype)
  h2 <- vanItersonhSquared(m,n,range)
  h <- sqrt(h2)
  res <- cbind(range,h)
  res <- res[!is.na(res[,2]),]
  res
}

#' The square of the disk radius D for a van Iterson diagram
#' @param m,n Integer. Orders of a generating parastichy pair
#' @param d Divergence of the lattice. 0<d<1.
#' @export
vanItersonDSquared <- function(m,n,d){
	dm <- m*d - round(m*d)
	dn <- n*d - round(n*d)
	D2  <- (n^2* dm^2 - m^2 *dn^2)/ (n^2-m^2)
	D2[D2<0] <- NA
	D2
}

#' The (h,d) coordinates of the (m,n,m+n) triple point in a van Iterson diagram
#' @param m,n Integers
#' @export
vanItersonTriplePoint <- function(m,n) {
	h2 <- (3/4)*1/(m*m+m*n+n*n)^2
	h <- sqrt(h2)
	uv <- windingNumbers (m,n,Delta=1)
	wm <- uv[1]
	dm <- (1/2)*(m+2*n)/(m*m+m*n+n*n)
	d <- (dm + wm)/m
	if (d>1/2) { d<- 1-d}
	c(h,d)
}



#' @export
getdrange <- function(m,n,ontype,lprtype) {

  grange <- dGeneratingIntervalLHS(m,n)
  orange <- dGeneratingOpposedIntervalLHS(m,n)

  # on the primary interval (unique in 0<d<1/2) m,n is the primary pair
  # left and right are the portions to the left or right of this where still generating
  # first classify by whether we are left, primary or right
  if (lprtype=="left") {
    rrange <- c(grange[1],vanItersonTriplePointL(m,n))
  } else if (lprtype=="primary") {
    rrange <- c(vanItersonTriplePointL(m,n),vanItersonTriplePointR(m,n))
  } else if (lprtype=="right") {
    rrange <- c(vanItersonTriplePointR(m,n),grange[2])
  } else stop("lprtype")
  # now which portion is generating
  if (ontype=="opposed") {
    orrange <- intervalIntersect(rrange,orange)
  } else  if (ontype=="nonopposed"){
    orrange <- intervalDelete(rrange,orange)
  } else stop("ontype")
  if (!any(is.na(orrange))) {
    resrange <- seq(orrange[1],orrange[2],length=1001)
  } else {
    resrange <- NA
  }
  resrange
}

#' The van Iterson triple point for m,n with smaller 0<d<1/2
#'
#' @param m,n Integer
#' @export
vanItersonTriplePointL <- function(m,n) {
  if (m==1 & n==2) return(vanItersonTriplePoint(1,2)[2])
  stopifnot(m<=n)
  pm <- sort(c(m,n,n-m))
  dpm  <- vanItersonTriplePoint(pm[1],pm[2])[2]
  dpp <- vanItersonTriplePoint(m,n)[2]
  min(dpm,dpp)
}
#' The van Iterson triple point for m,n with larger 0<d<1/2
#'
#' @param m,n Integer
#' @export
vanItersonTriplePointR <- function(m,n) {
		if (m==1 & n==2) return(1/2)
		stopifnot(m<=n)
		pm <- sort(c(m,n,n-m))
		dpm  <- vanItersonTriplePoint(pm[1],pm[2])[2]
		dpp <- vanItersonTriplePoint(m,n)[2]
		max(dpm,dpp)
}


#' The intersection of two intervals
#'
#' Throws an error if the resulting interval is disjoint
#' @param a,b Vectors of length two representing real intervals
#' @export
intervalIntersect <- function(a,b) {
  a <- sort(a)
  b <- sort(b)
  if (b[1]>a[2] || b[2]<a[1]) {return(NA) }
  c(max(a[1],b[1]),min(a[2],b[2]))
}

#' The difference of two intervals
#'
#' Throws an error if the resulting interval is disjoint
#' @param a,b Vectors of length two representing real intervals
#' @export
intervalDelete <- function(a,b) {
  if (a[1]< b[1] && b[2] < a[2]) stop("Disjoint intervals")
  if (b[2]<=a[1] || b[1]>= a[2]) {return (a)}
  if (b[1]<= a[1] && b[2]>=a[2]) {return (c(NA,NA))}
  if (b[2]< a[2]) { return (c(b[2],a[2])) } else { return (c(a[1],b[1])) }
}



vanItersonTriplePointR <- function(m,n) {
  if (m==1 & n==2) return(1/2)
  stopifnot(m<=n)
  pm <- sort(c(m,n,n-m))
  dpm  <- vanItersonTriplePoint(pm[1],pm[2])[2]
  dpp <- vanItersonTriplePoint(m,n)[2]
  max(dpm,dpp)
}



#' The point in d,h space where a touching circles lattice becomes pseudogenerating
#'
#' @param m Integer
#' @export
dhvanItersonPseudoGenerating <- function(m) {
  # if n=m+1 the equidistant line for m,n terminates when m ceases to be visible
  # when m=1/(2d)
  d <- 1/(2*m)
  h <- hvanIterson(m,m+1,d)
  c(d,h)
}
