#' Plot the regions on which d is generating/opposed for a range of m,n pairs
#' @param mnlist A list of integer (m,n) pairs
#' @export

plotGeneratingIntervals <- function(mnlist) {
  gix <- mnlist
  gixnames <- sapply(gix,function(x)sprintf("(%d,%d)",x[1],x[2]))
  ng <- length(gix)
  yrange <- c(0,ng+1)
  xrange <- c(0,1)
  #pushViewport(plotViewport(c(2,2,2,2))) # better be square
  grid.newpage()
  pushViewport(dataViewport(
    xData = xrange ,extension=0.01,
    yData = yrange ,
    name="plotRegion"))
  for (ix in 1:ng) {
    grid.lines(x=c(0,1),y=c(ix,ix),default.units="native")
    grid.text(x=.5,y=ix+.5,label=gixnames[[ng-ix+1]],default.units="native",
              just="center")
    grid.lines(x=c(0,0),y=ix+c(-.1,.1),default.units="native")
    grid.lines(x=c(0.5,0.5),y=ix+c(-.1,.1),default.units="native")
    grid.lines(x=c(1,1),y=ix+c(-.1,.1),default.units="native")
    mn <- gix[[ng-ix+1]]
    m <- mn[1]; n <- mn[2]
    uv <- windingNumbers(m,n,Delta=1)
    gi <- divergenceInterval (m,n,Delta=1,type="generating")
    grid.rect(x=gi[1],y=ix,width=gi[2]-gi[1],height=.2,default.units="native",
              gp=gpar(fill=grey(.6)),just=c("left"))
    gi <- divergenceInterval (m,n,Delta=-1,type="generating")
    grid.rect(x=gi[1],y=ix,width=gi[2]-gi[1],height=.2,default.units="native",
              gp=gpar(fill=grey(.4)),just=c("left"))
    oi <- divergenceInterval(m,n,Delta=1,type="generatingopposed")
    grid.rect(x=oi[1],y=ix,width=oi[2]-oi[1],height=.2,default.units="native",
              gp=gpar(fill="red"),just=c("left"))
    oi <- divergenceInterval(m,n,Delta=-1,type="generatingopposed")
    grid.rect(x=oi[1],y=ix,width=oi[2]-oi[1],height=.2,default.units="native",
              gp=gpar(fill="red"),just=c("left"))

  }
}
