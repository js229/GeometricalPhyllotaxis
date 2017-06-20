# todo - change default representation of Rho to be cylindrical circumference


###########################################
# for graphic display, but we embed one of these
# in our geometric objects so we know how to plot them if asked
# also stores the cylinder radius which is potentially of
# geometric importance but we fix its circumference at 1
#' An S4 class representing a cylinder and how a lattice should be displayed on it.
#' Node zero of the lattice has coordinates (0,0).  We display a window of size -1/2<x<1/2, 0<z<L, and by default node zero is at the centre of
#' this window, but this can be changed with the origin slot.
#' @slot Rho Radius of the cylinder unused
#' @slot L display height
#' @slot origin eg (-1/2,-L/2) will place node 0 at the bottom left of the display window
##' @slot nRepeats number of time to repeat periodically left and right NOT RELIABLE
#' @export
setClass("Cylinder",representation(Rho="numeric",origin="numeric",L="numeric"),prototype(Rho=1/(2*pi),origin=c(0,0),L=50))

#' Range of horizontal coordinate in the plotted cylindrical lattice
#' @param Cylinder an object inheriting from class Cylinder
#' @export
getXlim <- function(Cylinder) {
	xlim <- c(-1/2,1/2)
	xlim
}

#' Range of vertical coordinate in the plotted cylindrical lattice
#' @param Cylinder an object inheriting from class Cylinder
#' @export
getYlim <- function(Cylinder) {
	ylim <- c(-Cylinder@L/2,Cylinder@L/2)
	ylim
}

#' The unscaled radius of the cylinder
#' @param Cylinder an object inheriting from class Cylinder
getRho <- function(object,Jugacy=1) { cat(sprintf("getRho default generic called with object of class %s\n",class(object)))}
#' @describeIn getRho
setMethod("getRho","Cylinder",function(object) object@Rho)



.setup.grid <- function(Cylinder,doRepeat=0,doAxes=TRUE,doNewPage=TRUE,...) {
	# doRepeat is ingnored for now
	if (doNewPage)  { grid::grid.newpage() }
	pushViewport(plotViewport(c(0,0,0,0))) # better be square

	# we need equiscaled coordinates and want height L
	L <- Cylinder@L
	if (L<1) {
		warning("Cylinder too wide to equiscale plot; increase y")
	}

	ylim <- getYlim(Cylinder)
	xlim <- getXlim(Cylinder)

# use ylim as viewport xrange to ensure equiscaling
	pushViewport(dataViewport(
		xData = ylim,extension=0,
		yData = ylim,
		name="plotRegion"))
	if (doAxes) {
		doAxes <- as.numeric(doAxes)
		if (doAxes==1) {
			grid.xaxis(at=c(-1/2,0,1/2));
		}
	}
}

#setup.grid(Cylinder,doAxes=2)

.plot.Cylinder.grid <- function(Cylinder,
	doRepeat=0,doNewPage=TRUE,gpCylinder,...) {
	.setup.grid(Cylinder,doRepeat=doRepeat,doNewPage=doNewPage,...)
	L <- 	Cylinder@L

	if (missing(gpCylinder)) { gpCylinder <- gpar(fill="lightgreen",lty=0) }
	grid.rect(x=0,width=1,y=0,height=L,default.units="native",just="centre",gp=gpCylinder)

}

#' A class representing lattices defined by their divergence and rise
#' @slot Rise Numeric
#' @slot Divergence Numeric
#' @slot Jugacy Numeric (untested)
#' @slot Cylinder object of class Cylinder
#' @export
setClass("PhyllotaxisGenetic",
	representation(Rise="numeric",Divergence="numeric" , Jugacy="numeric"),
	contains="Cylinder"
	)

#' Create a new PhyllotaxisGenetic object
#' @param Rise Numeric
#' @param Divergence Numeric
#' @param Rho cylinder radius - Untested. Used if Cylinder is missing.
#' @param Jugacy untested
#' @param L Display height of cylinder. Used if Cylinder is missing.
#' @param A Cylinder object.
#' @export
newPhyllotaxisGenetic <- function(Rise=1,Divergence,Rho=1/(2*pi),Jugacy=1,L,Cylinder,...) {
	if (missing(L)) { L<- 50*Rise }
	if (missing(Cylinder)) { Cylinder <- new("Cylinder",Rho=Rho,L=L,...) }
	PG <- new("PhyllotaxisGenetic",Cylinder,Rise=Rise,Divergence=Divergence,Jugacy=Jugacy)
	}





#' @describein PhyllotaxisGenetic Get the divergence of a lattice object
#' @export
getRise <- function(object) { cat(sprintf("getRise default generic called with object of class %s\n",class(object)))}
getJugacy <- function(object) { cat(sprintf("getJugacy default generic called with object of class %s\n",class(object)))}
getDivergence <- function(object) { cat(sprintf("getJugacy default generic called with object of class %s\n",class(object)))}

#' @describein PhyllotaxisGenetic Get the divergence of a lattice object
#' @export
setMethod("getRise","PhyllotaxisGenetic",function(object) object@Rise)


#' @describein PhyllotaxisGenetic Get the divergence of a lattice object
#' @export
setMethod("getDivergence","PhyllotaxisGenetic",function(object) object@Divergence)
setMethod("getRho","PhyllotaxisGenetic",function(object) object@Rho)
setMethod("getJugacy","PhyllotaxisGenetic",function(object) object@Jugacy)

setGeneric("setRise<-",function(object,value)standardGeneric("setRise<-"))
setMethod("setRise<-","PhyllotaxisGenetic",function(object,value){
	object@Rise<-value;object})



.tox <- function(x){dx<- x%%1 ; ifelse(dx>1/2,dx-1,dx)}

# setup for plotPhyllotaxis
# should really abstract out the idea of drawing points, lines, circles, polgyons on a cylinder
# and silently handle both periodicity in primary strip and repeat strips if wanted...
plot.Cylinder.Lattice.grid <- function(PG,doRepeat=0,doNumbers=FALSE,pointSize=0.2,
	doCircles=FALSE,
	doFacets=FALSE,...) {
	# plot.Cylinder sets us up to have a square of with coords -L/2,L/2, with our cylinder origin at (0,0)


	# ie from the genetic spiral
	rise <- getRise(PG)
	divergence <- getDivergence(PG)
	J <- getJugacy(PG)

	Cylinder <- as(PG,"Cylinder")
	L <- Cylinder@L
	origin <- Cylinder@origin

	nodeZeroFromBottom <- origin[2]- (-L/2) # In [0,L] if node 0 is to be visible, but not required
	nodeAtBottom <- - nodeZeroFromBottom %/% rise
	sequenceLength <- L %/% rise
	bodge <- 5 # to avoid calculating how many extra we need
	sequenceNumbers <- seq(from=nodeAtBottom-bodge,length=sequenceLength+2*bodge)
	# that tells us all the nodes visible on the cylinder, but there
	# might be visible circles around nodes just above and below the edges
	# the bodge above is just because the calc is wrong
	# on the other hand the bodge below is to include circles which are visible but whose centres aren't
	circleSequenceNumbers <- seq(from=nodeAtBottom-bodge,length=sequenceLength+2*bodge)


	# need to modify latticePoint for J>1
	# and vectorise while we're at it
	xypts <- latticePoint(PG,sequenceNumbers)
	xypts.circle <- latticePoint(PG,circleSequenceNumbers)


	# that gives us x coordinates in the range [-1/2,1/2]
	xvalues <- xypts[,1]
	yvalues <- xypts[,2]

	xvalues.circle <- xypts.circle[,1]
	yvalues.circle <- xypts.circle[,2]
	if (doNumbers & J >1) {stop ("Multijugate unimplemented")}


	# for circles and facets we clip tightly to the cylinder
	# plot these first
	grid.clip(x=0,y=0,width=1,height=L,just="centre",default.units="native")

	PM <- GPPrincipalPhyllotaxisMatrix(PG)
	if (doCircles) {
		P1 <- getParastichyVector(PM,1)
		radius <- sqrt(sum(P1^2))/2
		outers <- c(-1,0,1)
		for (outer in outers) {
			oxpts <- xvalues.circle+  outer
			grid.circle(x=oxpts,y=yvalues.circle ,r=radius,default.units="native")
		}
	}
	if (doFacets) {
		P1 <- getParastichyVector(PM,1)
		P2 <- getParastichyVector(PM,2)
		outers <- c(-1,0,1)
		for (outer in outers) {
			oxpts <- xvalues.circle +  outer
			lxpts <- rep(oxpts,each=4)
			lypts <- rep(yvalues.circle,each=4)
			lpts <- cbind(lxpts,lypts)

			poly <-  matrix(
				rep(c( (P1+P2)/2,(-P1+P2)/2,(-P1-P2)/2,(P1-P2)/2),times=length(oxpts))
				,ncol=2,byrow=TRUE)
			polypts <- lpts + poly

			grid.polyline(x=polypts[,1],y=polypts[,2],id.lengths=rep(4,length(oxpts)),default.units="native")
		}
	}


	# if we want to show the periodic repeat, then we clip to the vertical section of cylinder, allowing nodes to appear on the repeats to the side
	# but not above and below the green edges
	# also allows nodes at top and bottom of cylinder to be fully visible

	grid.clip(x=0,y=0,width=1,height=L,just="centre",default.units="native")

	grid.points(x=xvalues,y=yvalues ,default.units="native",pch=16,gp=gpar(cex=pointSize))
	if (doRepeat>0) {
		outers <- c(-doRepeat:-1,1:doRepeat)
		for (outer in outers) {
			oxpts <- xvalues
			oxpts <- (oxpts +  outer)
			grid.points(x=oxpts,y=yvalues ,default.units="native",pch=1,gp=gpar(cex=pointSize))
		}

	}
	# we leave that clipping in place for the label numbers as they may be displaced to the side a little

	if (doNumbers) {
	  grid.clip(x=0,y=0,width=1,height=1,just=c("left","bottom"))
	  grid.text(label=paste("",sequenceNumbers),x=xvalues,y=yvalues ,default.units="native",
				just=c("left","center"))
	}

}

#' Plot a phyllotactic lattice
#' @param PG A lattice object which can be coerced to class PhyllotaxisGenetic
#' @param doRepeat =0,
#' @param doNumbers =FALSE,
#' @param pointSize =0.2,
#' @param doCircles =FALSE,
#' @param doFacets =FALSE
#' @param plotPrincipals Default FALSE. If TRUE, plot the principal parastichy vectors.
#' @export
plotPhyllotaxis <- function(PG,plotPrincipals=FALSE,...) {
		PC <- as(PG,"Cylinder")
		.plot.Cylinder.grid(PC,...)
		plot.Cylinder.Lattice.grid(PG,...)
		if (plotPrincipals) {
			plot.Principal.Vectors(PG,plotPrincipals)
			# not defined until later!!
		}
}



#' Coordinates of a point in the cylindrical lattice within the display window
#' @param PG A lattice object
#' @param n Sequence number of the point to return
#' @export
latticePoint <- function(PG,n) {
	Divergence <- PG@Divergence
	rise <- PG@Rise
	origin <- PG@origin

	x <- n * Divergence
	x <- x + origin[1]
	d <- .tox(x)
	h <- n * rise
	h <- h + origin[2]
	cbind(d,h)
}

getParastichyVectorOfPoint <- function(PG,n) {
  vec <-  latticePoint(PG,n)- latticePoint(PG,0)
  vec[1] <- .tox(vec[1])
  vec
}


clipToCylinder <- function(PG) {
	L <- PG@L
	grid.clip(x=0,y=0,width=1,height=L,just="centre",default.units="native")
}

clipToPlane <- function(PG) {
  L <- PG@L
  grid.clip(x=0,y=0,width=L,height=L,just="centre",default.units="npc")
}



plotParastichy <- function(PG,parastichyCount,through,type="cylinder",gp) {
	rise <- PG@Rise
	L <- PG@L
	if (missing(through)) {
		for (through in seq(0,parastichyCount-1)) {
			plotParastichy(PG,parastichyCount,through=through,type=type,gp)
		}
	}

	throughNode <- latticePoint(PG,through)
	parastichyVector <- getParastichyVectorOfPoint(PG,parastichyCount)
	slope <- parastichyVector[2]/parastichyVector[1]
	cylinderRise <- slope * 1 # nb may be negative

	lhs <- -1/2
	z.lhs <-throughNode[2] - slope * (throughNode[1] - (lhs)) # height of the parastichy through that node at lhs of cylinder


	if(type=="cylinder") {
	  .plot.parastichylines(PG,lhs=-1/2,rhs=1/2,z.lhs,cylinderRise,gp)
	} else { #"plane"
	  .plot.parastichylines(PG,lhs=-3/2,rhs=-1/2,z.lhs,cylinderRise,gp)
	  .plot.parastichylines(PG,lhs= 1/2,rhs= 3/2,z.lhs,cylinderRise,gp)
	}

}

	.plot.parastichylines <- function(PG, lhs, rhs, z.lhs, cylinderRise,gp) {
	  # z.lhs may not be in -L/2,L/2!
	  L <- PG@L
	  nabove <- 1 + ((L / 2 - z.lhs) %/% abs(cylinderRise))
	  nbelow <- max(0,1 + ((z.lhs - (-L / 2)) %/% abs(cylinderRise)))

	  ntodo <- c(-seq(length = nbelow), 0, seq(nabove))
	  for (n in ntodo) {
	    grid.lines(
	      x = c(lhs, rhs),
	      y = c(z.lhs + n * cylinderRise, z.lhs + (n + 1) * cylinderRise),
	      default.units = "native",
	      gp = gp
	    )
	  }

}

	#' Plot a segement of a parastichy line without wrapping
	#' @param PG A lattice object
	#' @param parastichyCount The order of the parastichy
	#' @param through The sequence number of the point to pass through
	#' @param length The length of the parastichy line (in cylinder units, ugh)
	#' @param gp A gpar object
	#' @export
	plotParastichySegment <- function(PG,parastichyCount,through,length,gp) {
	  rise <- PG@Rise
	  L <- PG@L
	  # don't worry about wrapping round cylinder
	  throughNode <- latticePoint(PG,through)
	  parastichyVector <- getParastichyVectorOfPoint(PG,parastichyCount)

	  points.x <- throughNode[1]+length*c(0,parastichyVector[1])
	  points.y <- throughNode[2]+length*c(0,parastichyVector[2])

	  grid.lines(x=points.x,y=points.y, default.units = "native",gp = gp)


	}


#' Plot a lattice and a parastichy pair
#' @param PG A lattice object
#' @param m Order of the parastichy to plot
#' @param n Order of the parastichy to plot
#' @param doParallelogram UNUSED
#' @param dodash UNUSED
#' @export
	plotParastichyPairs <- function(PG,y,m,n,doParallelogram=0,dodash=TRUE) {
	  plotPhyllotaxis(PG,doNumbers=TRUE,doRepeat=0,doAxes=FALSE)
	  plotParastichy (PG,parastichyCount=m,gp=gpar(col="blue"))
	  if (dodash) {
	    plotParastichy (PG,parastichyCount=m,gp=gpar(col="blue",lty=2),type="plane")
	    plotParastichy (PG,parastichyCount=n,gp=gpar(col="red",lty=2),type="plane")
	  }
	  plotParastichy (PG,parastichyCount=n,gp=gpar(col="red"))

	  if (!is.na(doParallelogram)) {
	    plotParallelogram(PG,m,n,through=doParallelogram)
	  }
	}

plotParallelogram <- function (PG,m,n,through) {

	mvec <- getParastichyVectorOfPoint(PG,n=m)
	nvec <- getParastichyVectorOfPoint(PG,n=n)
	ogin <- latticePoint(PG,n=through)

	pgram.x <- ogin[1]+ c( 0,mvec[1],mvec[1]+nvec[1],nvec[1],0)
	pgram.y <- ogin[2]+ c( 0,mvec[2],mvec[2]+nvec[2],nvec[2],0)

	grid.lines(x=pgram.x,y=pgram.y,default.units="native",gp=gpar(col="black",lwd=2))
}



plot.Parastichy.Vectors.grid <- function(Phyllotaxis,plotPrincipals) { #
  origin <- Phyllotaxis@origin
	P1 <- getParastichyVector(Phyllotaxis,1)+origin
	P2 <- getParastichyVector(Phyllotaxis,2)+origin
	L <- Phyllotaxis@L
	if (is.numeric(plotPrincipals)) {
	  warning("Code with numeric plotprinciapls does not work")
		P1scaled <- P1 * ceiling((L/2)/abs(P1[2]))
		P2scaled <- P2 * ceiling((L/2)/abs(P2[2]))
		rho <- getRho(Phyllotaxis)
		grid.clip(x=-pi*rho,y=-L/2,width=2*pi*rho,height=L,just=c("left","bottom"),default.units="native")
		for (offset in -plotPrincipals : plotPrincipals) {
			v0 <- -P1scaled + P2 * offset
			v1 <- P1scaled + P2*offset
			grid.lines(x=c(v0[1],v1[1]),y=c(v0[2],v1[2]),default.units="native",gp=gpar(col="red"))
			v0 <- -P2scaled + P1 * offset
			v1 <- P2scaled + P1*offset
			grid.lines(x=c(v0[1],v1[1]),y=c(v0[2],v1[2]),default.units="native",gp=gpar(col="blue"))
			}
	} else { # just the vectors
		grid.lines(x=c(origin[1],P1[1]),y=c(origin[2],P1[2]),default.units="native",gp=gpar(col="red"))
		grid.lines(x=c(origin[1],P2[1]),y=c(origin[2],P2[2]),default.units="native",gp=gpar(col="blue"))
	}
}

plot.Principal.Vectors<- function(object,plotPrincipals=TRUE,...) {
		object.PM <- as(object,"PhyllotaxisMatrix")
		object.PP <- makePrincipalPhyllotaxisMatrix (object.PM )
		plot.Parastichy.Vectors.grid(object.PP,plotPrincipals)
}



