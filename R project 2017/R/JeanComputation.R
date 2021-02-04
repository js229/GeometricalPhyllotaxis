#' Perform the computation of Jean(1994) p40
#' @param Divergence
#' @param kmax
#' @param doOpposed
#' @export
JeanComputation <- function(divergence, kmax=21,doOpposed=TRUE) {
	# this is an implementation of the computation on p40 of Jean,
	# who negelects to add that the dividing line is defined by the sign
	# of the $x$ coordinate

	# it doesn't work correctly for rational divergences
	A <- data.frame(k=1:kmax)
	A$kdivergence <- A$k*divergence
	A$kd <- round(A$kdivergence)
	A$dk <- A$kdivergence - A$kd
	A$hcf <- NA
	for (k in 1:kmax) A$hcf[k] <- hcf(A$kd[k],k)
	A <- A[A$hcf<=1,]
	A$fracval <- A$kd/A$k
	A$ksequence <- 1:nrow(A)
	A <- A[order(A$fracval),]
	showit <- function(Fstart) {

		 cat(	"\n",sapply(Fstart[[1]],function(k){sprintf("%3g",A$kd[A$k==k])}),
			"|",
			sapply(Fstart[[2]],function(k){sprintf("%3g",A$kd[A$k==k])}),
			"\n",
			sapply(Fstart[[1]],function(k){sprintf("%3g",k)}),
			"|",
			sapply(Fstart[[2]],function(k){sprintf("%3g",k)}),"\n"
		)
	}

	Fset <- list()
	for(k in 1:(nrow(A)-1)) {
		Aset <- A[A$ksequence<=k+1,]
		Fset[[k]] <- split(Aset$k,Aset$dk<0)
		showit(Fset[[k]])
	}
	visiblePairs <- do.call(rbind,sapply(Fset, function(k){
		kd <- c(k$`FALSE`,k$`TRUE`);
		kpairs <- cbind(kd[-length(kd)],kd[-1]);
	}))
	visiblePairs <- unique(visiblePairs)
	opposedvisiblePairs <- do.call(rbind,lapply(Fset, function(k){
		kd <- c(k$`FALSE`[length(k$`FALSE`)],k$`TRUE`[1]);
	}))
	opposedvisiblePairs <- unique(opposedvisiblePairs )
	if (doOpposed) {
		return(opposedvisiblePairs)
	} else {
		return(visiblePairs)
	}

}


