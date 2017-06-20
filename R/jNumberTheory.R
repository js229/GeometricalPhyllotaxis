congruence <- function (a, l = 1)  { # from Hankin's elliptic library
    l <- as.integer(l)
    m <- as.integer(a[1])
    n <- as.integer(a[2])
    zero <- as.integer(0)
    one <- as.integer(1)
    if (m == zero & n == one) {
        return(NULL)
    }
    if (m == one & n == zero) {
        return(c(NA, 1))
    }
    if (m == 1) {
        return(rbind(c(a), c(1, n + l), c(0, l)))
    }
    if (n == 1) {
        return(rbind(c(a), c(m - l, 1), c(l, 0)))
    }
    q1 <- which((+l + (1:m) * n)%%m == zero)
    if (!any(q1)) {
        q1 <- NA
    }
    q2 <- which((-l + (1:n) * m)%%n == zero)
    if (!any(q2)) {
        q2 <- NA
    }
    out <- rbind(a, cbind(q1, q2))
    rownames(out) <- NULL
    colnames(out) <- NULL
    return(out)
}

jFactorize <- function(n) # renamed from factorize
{
  ## Purpose:  Prime factorization of integer(s) 'n'
  ## -------------------------------------------------------------------------
  ## Arguments: n vector of integers to factorize (into prime numbers)
  ##	--> needs 'prime.sieve'
  ## >> Better would be: Define class 'primefactors' and "multiply" method
  ##			 then use this function recursively only "small" factors
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 26--30 Jan 96
  if(all(n < .Machine$integer.max))
      n <- as.integer(n)
  else {
      warning("factorizing large int ( > maximal integer )")
      n <- round(n)
  }
  N <- length(n)
  M <- trunc(sqrt(max(n))) #-- check up to this prime number
  ##-- for M > 100 to 200: should DIVIDE first and then go on ..
  ##-- Here, I am just (too) lazy:
  k <- length(pr <- prime.sieve(maxP = M))
  nDp <- outer(pr, n, FUN = function(p,n) n %% p == 0) ## which are divisible?
  ## dim(nDp) = (k,N) ;
  ## Divide those that are divisible :
  ## quot <- matrix(n,k,N,byrow=T)[nDp] %/% matrix(pr,k,N)[nDp]
  ## quot <- rep(n,rep(k,N))[nDp] %/% rep(pr,N)[nDp]
  res <- vector("list",length = N)
  names(res) <- n
  for(i in 1:N) { ## factorize	n[i]
    nn <- n[i]
    if(any(Dp <- nDp[,i])) { #- Dp: which primes are factors
      nP <- length(pfac <- pr[Dp]) # all the small prime factors
      if(exists("DEBUG")&& DEBUG) cat(nn," ")
    } else { # nn is a prime
      res[[i]] <- cbind(p = nn, m = 1)
  #    print("direct prime", nn)
      next # i
    }
    m.pr <- rep(1,nP)# multiplicities
    Ppf <- prod(pfac)
    while(1 < (nn <- nn %/% Ppf)) { #-- have multiple or only bigger factors
      Dp <- nn %% pfac == 0
      if(any(Dp)) { # have more smaller factors
	m.pr[Dp] <- m.pr[Dp] + 1
	Ppf <- prod(pfac[Dp])
      } else { #-- the remainder is a bigger prime
	pfac <- c(pfac,nn)
	m.pr <- c(m.pr, 1)
	break # out of while(.)
      }
    }
    res[[i]] <- cbind(p = pfac,m = m.pr)

  } # end for(i ..)

  res
}

#' Highest common factor of two integers
#' @param n1, n2 Integers
#' @export
hcf <- function(n1,n2) {
	nvec <- c(n1,n2)
	isInt <- all(sapply(nvec,function(n) isTRUE(all.equal(n-round(n),0))))
	if (!isInt) {
		warning(sprintf("Nonintegers (%g,%g) in hcf",n1,n2))
		return(1)
	}
	if (n1*n2==0) {
#		warning(sprintf("Zero (%g,%g) in hcf",n1,n2))
		return(0)
	}
	nvec <- abs(nvec)


	if ( min(nvec)==1 ) {
		hcf <- 1
	} else {
		eta1.fac <- data.frame(jFactorize(nvec[1])[[1]])
		eta2.fac <- data.frame(jFactorize(nvec[2])[[1]])
		eta.pow <- merge(eta1.fac,eta2.fac,by="p")
		if (nrow(eta.pow)==0) {
			hcf <- 1
		} else {
			eta.pow$minp <- min(eta.pow$m.x,eta.pow$m.y)
			eta.pow$val <- eta.pow$p ^ eta.pow$minp
			hcf <- prod(eta.pow$val)
		}
	}
	hcf
}


jCongruence <- function(m,n) { #x,y s.t. x * m - y * n = 1
	h <- hcf(m,n);
	if (h>1) {
		warning(sprintf("%g and %g are not coprime: hcf %d",m,n,h))
	}
	mh <- m/h; nh <- n/h
	res <- congruence(c(mh,nh))[2,]
	x <- res[2]; y <- res[1]
	c(x,y)
}

testCongruence <- function(m,n) {
	res <- jCongruence(m,n)
	cat(sprintf("%d * %d - %d * %d = %d \n", m ,res[1],n,res[2],m*res[1] - n*res[2]))
}
#testCongruence(4,11)
#testCongruence(4,8)
do.test.Congruence <- function() {
for (i in 1:10) {
	for (j in 1:10) {
		GeometricalPhyllotaxis:::testCongruence(i,j)
	}
}
}
#do.test.Congruence

prime.sieve <- function(p2et = c(2,3,5), maxP = pM^2)
{
  ## Purpose: Produce ALL prime numbers from 2, 3.., using 2,3,5,7,...
  ## -------------------------------------------------------------------------
  ## Arguments: p2et: primes c(2,3,5,..., pM);
  ##		maxP : want primes up to maxP
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date: 26 Jan 96, 15:08
  if(any(p2et[1:2] != 2:3) || is.unsorted(p2et) || !is.numeric(p2et))
	stop("argument 'p2et' must contain SORTED primes 2,3,..")
  k <- length(p2et)
  pM <- p2et[k]
  if(maxP <= pM+1) p2et #- no need to compute more
  else if(maxP > pM^2)	prime.sieve(prime.sieve(p2et), maxP = maxP)
  else { #-- pM < maxP <= pM^2
    r <- seq(from = pM+2, to = maxP, by = 2)
    for(j in 1:k)
      if(0 == length(r <- r[r%% p2et[j] != 0])) break
    c(p2et,r)
  }
}

#'
#' Give the continued fraction representation
#' @param x Real
#' @param nmax Integer. Maximum number of terms to return
continuedFraction <- function(x, nmax=10) {
  tt <- rep(NA,nmax); uu <- rep(NA,nmax);
  tt[1] <- floor(x)
  uu[1] <- 1/x- tt[1]
  for (k in 2:nmax) {
    if (uu[k-1]==0) break;

    tt[k] <- floor(1/uu[k-1])
    uu[k] <- uu[k-1]-tt[k-1]
  }
}
