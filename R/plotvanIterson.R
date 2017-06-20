
mnshow2 <- function(m,n,labellines,onlyPrimary,primaryColour="black") {
  v1 <- dhvanIterson(m,n,ontype="opposed",lprtype="primary")
  lines(v1,lwd=2,col=primaryColour)
  v2 <- dhvanIterson(m,n,ontype="nonopposed",lprtype="primary")
  lines(v2,lwd=2,col=primaryColour,lty=3)
  v1 <-v1[order(v1[,2]),]
  lpos <- min(which(v1[,2]>median(v1[,2])))

  colval <- "darkgray"
  if (!onlyPrimary) {
    v3 <- dhvanIterson(m,n,ontype="opposed",lprtype="left")
    lines(v3,lwd=1,col=colval ,lty=1)
    v4 <- dhvanIterson(m,n,ontype="nonopposed",lprtype="left")
    lines(v4,lwd=1,col=colval ,lty=3)
    v5 <- dhvanIterson(m,n,ontype="opposed",lprtype="right")
    lines(v5,lwd=1,col=colval ,lty=1)
    v6 <- dhvanIterson(m,n,ontype="nonopposed",lprtype="right")
    lines(v6,lwd=1,col=colval ,lty=3)
    #	v <- rbind(v1,v2,v3,v4,v5,v6)
  }
  #	vmax <- v[which.max(v[,2]),] + c(0, .005)

  if (labellines) {
    text(v1[lpos,1],v1[lpos,2],label=sprintf(" %d=%d",m,n),cex=0.5,adj=0)
    #		text(vmax[1],vmax[2],label=sprintf("%d=%d",m,n),cex=0.5)
  }
}

#' Plot elements of lattice space including the van Iterson tree
#'
#' @param vpairs A data frame in which the first two columns are pair of parastichy numbers. If the third column is FALSE, do not label the pair
#' @param labellines Label the branches of the van Iterson tree. Default FALSE
#' @param labelregions Label the regions by their principal parastichy values
#' @param vpairlabels A data frame in which the first two columns are the parastichy numbers for labelling the region and the second two the (d,h) coordinates where the labels sohuld be placed
#' @param xlim,ylim A pair of reals defining the plotting area to use
#' @param onlyPrimary If true, do not plot lines on which a pair is van Iterson but not primary
#' @param doCircles If true, plot entire circles on which van Iterson branches lie
#' @param primaryColour The colour for the primary branch
#' @param doTriplePoints If true, plot a marker at each triplepoint
#' @export
plotvanIterson <- function(vpairs,labellines=FALSE,labelregions=TRUE,vpairlabels,xlim,ylim,onlyPrimary=FALSE,doCircles=TRUE,
                   primaryColour="black",doTriplePoints=FALSE) {
  ylab <- "Rise h"
  ymax <- 0.3;

  if (missing(ylim)) {ylim <- c(0,ymax)}
  if (missing(xlim)) {xlim <- c(0,.5)}
  plot(0,xlim=xlim,ylim=ylim,ylab=ylab,type="n",xlab="Divergence d")
  palette(rainbow(5))

  if (doCircles) {
    for (ix in 1:nrow(vpairs)){
      m <- vpairs[ix,1]
      n <- vpairs[ix,2]
      uv <- windingNumbers(m,n,DeltaForLHSd(m,n))
      u <- uv[1];v<-uv[2]
      R <- 1/(n^2-m^2)
      d0 <- (n*v - m*u)/(n^2-m^2)
      theta <- seq(0,pi,length=101)
      dh <- cbind(R*cos(theta)+d0, R*sin(theta))
      lines(dh[,1],dh[,2],col="lightgrey")
    }
  }



  if (doTriplePoints){ # plot triple points
    for (ix in 1:nrow(vpairs)){
      m <- vpairs[ix,1]
      n <- vpairs[ix,2]
      hd <- vanItersonTriplePoint(m,n);
      h <- hd[1];d <- hd[2]
      #		text(d,h,label=sprintf("%d%d",m,n))
      points(d,h,pch=1)
    }
  }

  for (ix in 1:nrow(vpairs)) {
    vpr <- vpairs; if (!labellines) vpr$showLabel <- FALSE
    mnshow2(vpr[ix,1],vpr[ix,2],labellines=vpr$showLabel[ix],onlyPrimary=onlyPrimary,
            primaryColour=primaryColour)
  }

  if (!onlyPrimary) {
    colval <- "darkgray"
    mrange <- setdiff(unique(vpairs[,1]),1)
    for (m in mrange) {
      dh <- dhvanItersonPseudoGenerating(m)
      lines(c(dh[1],dh[1]),c(dh[2],ymax),col=colval,lty=3)
    }
  }

  if (labelregions) {
    text(vpairlabels[,3],vpairlabels[,4],label=vpairlabels[,5],cex=.7)
  }

}
