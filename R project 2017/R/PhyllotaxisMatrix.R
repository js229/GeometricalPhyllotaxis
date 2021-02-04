

#######################################################################
# The rest of this file is about a matrix representation of the lattice



setClass("PhyllotaxisMatrix",representation(
	nScale="numeric",dScale="numeric", Coord="matrix"),  # Coord will be integer,
	contains="Cylinder"
	)


getRise.PM <- function(object) { # for a PhyllotaxisMatrix
  mn <- getParastichyNumbers(object)
  rise <- hcf(mn[1],mn[2])
  rise <- rise * object@nScale/object@dScale
  rise
}

getParastichyNumbers<- function(object,order=c(1,2)) {
	object@Coord[order,2]
}

getParastichyVector <- function(object,order) {
	mat <- as(object,"matrix")
	mat[order,]
}

getRho.PM <- function(object,Jugacy) {# for a PhyllotaxisMatrix
	mn <- getParastichyNumbers(object)
	m <- mn[1]; n <- mn[2]
	if (missing(Jugacy)) {J <- getJugacy(object)} else { J<-Jugacy}
	P1 <- getParastichyVector(object,1)
	P2 <- getParastichyVector(object,2)
	rho <- J* abs(( n * P1[1] - m * P2[1]))/(2 * pi)
	stopifnot(rho != 0)
	rho
}

setMethod("getRho","PhyllotaxisMatrix",getRho.PM)


extractHCFIntoRise <- function(object) {
	mn <- getParastichyNumbers(object)
	if ( all(mn==floor(mn))) {
		h <- hcf(mn[1],mn[2])
		object@Coord <- object@Coord/h
		object@nScale <- object@nScale*h
	} else {
		stop ("noninteger values")
	}
	object
}


setAs("matrix","PhyllotaxisMatrix", function(from) {
	object <- new("PhyllotaxisMatrix", nScale=1,dScale=1, Coord=from  )
	rho <- getRho.PM(object,Jugacy=1)
	object@Rho <- rho
	object <- extractHCFIntoRise(object)
	object
})

setAs("PhyllotaxisMatrix","matrix", function(from) 	from@Coord*from@nScale/from@dScale)


###################################################
### chunk number 39: defpmtest
###################################################
#line 1146 "GeometricPhyllotaxisTestVignette.Rnw"
#matrix.test <- matrix(c(2,1,-1,2),ncol=2,byrow=TRUE)
#PM.test <- as(matrix.test,"PhyllotaxisMatrix")



summary.PhyllotaxisMatrix.matrix <- function(object,...) {
	cat(	sprintf("|a\tb|  = %-3g/%-3g *  |%3g\t%3g|\n",object@nScale,object@dScale,object@Coord[1,1],object@Coord[1,2]))
	cat(	sprintf("|c\td|               |%3g\t%3g|\n",object@Coord[2,1],object@Coord[2,2]))
}

#setMethod("print","PhyllotaxisMatrix",function(x,...)summary.PhyllotaxisMatrix.matrix(x))
#summary.PhyllotaxisMatrix.matrix(PM.test)

setAs("PhyllotaxisGenetic","PhyllotaxisMatrix",function(from) {
	alpha <- getDivergence(from)
	rise <- getRise(from)
	a <-  alpha / rise
	b <- 1
	nover <- floor(1/alpha)+1
	c <- ((nover * alpha) %% (1)) * 1 / rise
	d <- nover
	new("PhyllotaxisMatrix",as(from,"Cylinder"),
		 nScale=rise,dScale=1,Coord=matrix(c(a,b,c,d),byrow=TRUE,ncol=2)
		)
})
#PGM.test <- as(PG.test,"PhyllotaxisMatrix")

#PGM.test.R <- as(PG.test.R,"PhyllotaxisMatrix")

getCongruence <- function(object) {
	mn <- (getParastichyNumbers(object))
	absmn <- abs(mn)
	kl <- jCongruence(absmn[1],absmn[2]) # assumes all +ve
	kl <- kl * sign(mn)
	stopifnot(kl[1] * mn[1] - kl[2] * mn[2] ==1)
	kl
}


getDivergence.PM <- function(object) {# for a PhyllotaxisMatrix
	# for this, insist that P numbers are positive and exactly integer
	# 	find a point at k P1 - l P2 with a rise of 1
	kl <- getCongruence (object)
	if (length(kl)!=2) {return(NA)}
	k <- kl[1]; l <- kl[2]
	# use positive parastichy vectors
	P1 <- getParastichyVector(object,1) ; if (P1[2]<0) { P1 <- -P1 }
	P2 <- getParastichyVector(object,2) ; if (P2[2]<0) { P2 <- -P2 }


# what's its x-ccord ?
	rho <- getRho(object)
	alpha <-  (k * P1[1] - l * P2[1]) / rho
	alpha <- alpha %% (2 * pi)

	alpha
}

setMethod("getDivergence","PhyllotaxisMatrix",getDivergence.PM )


setMethod("getJugacy","PhyllotaxisMatrix",function(object) {1} )


summary.PhyllotaxisMatrix.genetic <- function(object,...) {
	mn <- getParastichyNumbers(object)
	kl <-  getCongruence(object)
	cat(	sprintf("(m,n)=(%g,%g); (%g*%g)-(%g*%g) = %g\n",mn[1],mn[2],kl[1],mn[1],kl[2],mn[2],mn[1]*kl[1]-mn[2]*kl[2]))
	cat(sprintf("Rise %g\n",getRise(object)))
	cat(sprintf("Minimal cylinder radius %g circumference %g \n",getRho(object),2*pi*getRho(object)))
	alpha <- getDivergence(object)
	d <- alpha/(2*pi)
	if (d>1/2) {d<- d-(1/2)}
	cat(sprintf("Divergence %g (%g degrees; Jean d %g)\n",
		alpha , alpha*360/(2*pi),d)
	)
	d1 <- abs(kl[2]/mn[1])
	d2 <- abs(kl[1]/mn[2])
	cat(sprintf("Visible opposed divergence interval (%g,%g)\n", d1 * 2 *pi,d2* 2 * pi))
	}

#setMethod("print","PhyllotaxisMatrix",function(x,...) {summary.PhyllotaxisMatrix.matrix(x);summary.PhyllotaxisMatrix.genetic(x)}		)

makeGenetic <- function(from) {
	rise <- getRise(from)
	rho <- getRho(from)
	alpha <-getDivergence(from)
	J <- getJugacy(from)
	Cylinder <- as(from,"Cylinder")
	PG <- newPhyllotaxisGenetic(Rise=rise,Divergence=alpha,Cylinder=Cylinder,Jugacy=J)
	PG }

setAs("PhyllotaxisMatrix","PhyllotaxisGenetic",function(from)makeGenetic(from))


setMethod("plot",signature(x="PhyllotaxisMatrix",y="numeric"),
	function(x,y,plotPrincipals=TRUE,...) {
		PG <- as(x,"PhyllotaxisGenetic")
		plotPhyllotaxis(PG,y,plotPrincipals=plotPrincipals,...)
	})




###################################################
### chunk number 56: doreduction
###################################################


dot  <- function(u,v) { sum(u*v)}
mod <- function(u) { sqrt(dot(u,u))}

isPrincipalMatrix <- function(object) {
	u = getParastichyVector(object,1)
	v = getParastichyVector(object,2)
	dot(u,v)>= 0 & mod(u) <= mod(v) & mod(v)<=mod(u-v)
}


isOpposedMatrix <- function(object) {
	u = getParastichyVector(object,1)
	v = getParastichyVector(object,2)
	if (sign(u[1]) < 0 ) { u <- -u }
	if (sign(v[1]) < 0 ) { v <- -v }
	sign(u[2])*sign(v[2]) <= 0
}

swapParastichyVectors <- function(from) {
	object <- from
	object@Coord <- object@Coord[ c(2,1),]
	object
}

makeDotPositive <- function(object,trace=FALSE) {
	u = getParastichyVector(object,1)
	v = getParastichyVector(object,2)
	if (dot(u,v)<0) {
		if (trace) {cat("M") }
		object <- addkUtolV (object,0,-1)
	}
	object
}

makeUTheShorter <- function(object,trace=FALSE) {
	u = getParastichyVector(object,1)
	v = getParastichyVector(object,2)
	if ( sum(u*u) > sum(v * v)) {
		if (trace) { cat("S") }
		object <- swapParastichyVectors(object)
	}
	object
}

addkUtolV <- function(object,k,l) { # u,v -> u, k u + l v
	Coord <- object@Coord
	Coord[2,] <- k * Coord[1,] + l * Coord[2,]
	object@Coord <- Coord
	object
}

reduceV <- function(object,trace=FALSE) { #u,v -> u,   u - v
	if(trace) { cat("R")}
	addkUtolV(object,1,-1)
}

reverseXY <- function(object) {
	object@Coord <- - object@Coord
	object
}

makePrincipalPhyllotaxisMatrix <- function(object,trace=FALSE) {
	object <- makeUTheShorter(object,trace)
	object <- makeDotPositive(object,trace)
	while(!isPrincipalMatrix(object)) {
		if(trace) {		cat("\n=====\n");print(object);cat("=====\n")}
		object <- reduceV(object,trace)
		object <- makeUTheShorter(object,trace)
		object <- makeDotPositive(object,trace)
	}
	if(trace){cat("\n")}
	# finally, adopt convention that first principal vector has positive y component
	pysign <- getParastichyVector(object,1)[2] <  0
	if (pysign) {
		object <- reverseXY(object)
	}
	object
}

summary.PhyllotaxisMatrix.reduction <- function(object,...) {
	if (isPrincipalMatrix(object)) {
		cat("Principal vector matrix\n")
	} else {
		cat("Not a principal vector matrix\n")
	}
	if (isOpposedMatrix(object)) {
		cat("Opposed\n")
	} else {
		cat("Not opposed\n")
	}



	u <- getParastichyVector(object,1);v=getParastichyVector(object,2)
	cat(sprintf("u.v: %g; u^2 %g; v^2: %g; (u-v)^2: %g; (u+v)^2: %g\n",sum(u*v),sum(u*u),sum(v*v),sum((u-v)*(u-v)),sum((u+v)*(u+v))))
	cat(sprintf("Delta^2 %g \n",sum(u*u)*sum(v*v)-sum(u*v)^2))
}

summary.PhyllotaxisMatrix <- function(object,...) {
	summary.PhyllotaxisMatrix.matrix(object)
	summary.PhyllotaxisMatrix.reduction(object,...)
	summary.PhyllotaxisMatrix.genetic (object,...)
}
setMethod("print","PhyllotaxisMatrix",function(x,...)summary.PhyllotaxisMatrix(x))


GPPrincipalPhyllotaxisMatrix <- function(object) {
	PM  <- as(object,"PhyllotaxisMatrix")
	PM <- makePrincipalPhyllotaxisMatrix(PM)
	PM
}
GPPrincipalParastichyNumbers <- function(object) {
	PM  <- as(object,"PhyllotaxisMatrix")
	PM <- makePrincipalPhyllotaxisMatrix(PM)
	getParastichyNumbers(PM)
}


