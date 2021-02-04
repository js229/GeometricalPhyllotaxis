library(GeometricalPhyllotaxis)
context("Tests of congruence and interval creation")

test_that("Winding numbers work", {
  expect_equal(windingNumbers(1,1,1),c(.5,.5))
  expect_equal(windingNumbers(4,7,1),c(1,2))
  expect_equal(windingNumbers(1,8,1),c(0,1))
  expect_equal(windingNumbers(9,1,1),c(8,1))
  expect_equal(windingNumbers(19,55,1),c(10,29))
  expect_equal(windingNumbers(1,1,-1),c(.5,.5))
  expect_equal(windingNumbers(4,7,-1),c(3,5))
  expect_equal(windingNumbers(1,8,-1),c(1,7))
  expect_equal(windingNumbers(9,1,-1),c(1,0))

  # used in book
  expect_equal(windingNumbers(1,10,1),c(0,1))
  expect_equal(windingNumbers(13,1,1),c(12,1))
  expect_equal(windingNumbers(6,5,1),c(1,1))
  expect_equal(windingNumbers(5,6,1),c(4,5))
  expect_equal(windingNumbers(34,55,1),c(21,34))
  expect_equal(windingNumbers(55,89,1),c(21,34))


  expect_equal(windingNumbers(19,55,1),rev(windingNumbers(55,19,-1)))
})

test_that("generating intervals have special case properties", {


  expect_equal(divergenceInterval (m=1,n=2,Delta=1,scale="twomn",type="generating"),c(1,2))
  expect_equal(divergenceInterval (m=1,n=11,Delta=1,scale="twomn",type="generating"),c(1,3))
  m<- 11; expect_equal(divergenceInterval (m=m,n=1,Delta=1,scale="twomn",type="generating"),c(2*m-3,2*m-1))
  m<- 11; expect_equal(divergenceInterval (m=m,n=m+1,Delta=1,scale="twomn",type="generating"),(2*m-1)*c(m,m+1))
  n<- 7; expect_equal(divergenceInterval (m=2,n=n,Delta=1,scale="twomn",type="generating"),c(2*n,2*n+4))
  n<- 10; expect_equal(divergenceInterval (m=3,n=n,Delta=1,scale="twomn",type="generating"),c(4*n-1,4*n+5))
  # 34 = F(9)
  m<- 34; n<- 55; expect_equal(divergenceInterval (m=m,n=n,Delta=1,scale="twomn",type="generating"),m*c(2*m-1,2*m+1))
  # different for odd or even n in F(n)
  m<- 55; n<- 89; k<- 34;
  expect_equal(divergenceInterval (m=m,n=n,Delta=1,scale="twomn",type="generating"),m*c(2*k-1,2*k+1))

  pdi <- divergenceInterval (m=1,n=4,scale="unity",type="pseudogenerating")
  dinpi <- pdi[1]+ 0.3*pdi[2]
  expect_equal(Deltamnd(1,4,1),0)

  dint <- divergenceInterval (m=4,n=11,Delta=1,scale="unity",type="generating"); dindi <- dint[1]+.63*(dint[2]-dint[1])
  expect_equal(windingNumbers(4,11,1) , c(round(4*dindi), round(11*dindi)))



})


test_that("generating intervals are actually generating", {
  gi <- divergenceInterval (m=4,n=5,Delta=1,scale="unity",type="generating")
  gi2 <- divergenceInterval (m=4,n=5,Delta=1,scale="twomn",type="generating")
  expect_equal(gi2[2],2*4*5*gi[2])
  dingi <- gi[1]+ 0.25*(gi[2]-gi[1])
  expect_equal(Deltamnd(4,5,dingi),1)
  expect_equal(Deltamnd(4,5,.6),2)



  pdi <- divergenceInterval (m=1,n=4,scale="unity",type="pseudogenerating")
  dinpi <- pdi[1]+ 0.3*pdi[2]
  expect_equal(Deltamnd(1,4,1),0)

  dint <- divergenceInterval (m=4,n=11,Delta=1,scale="unity",type="generating"); dindi <- dint[1]+.63*(dint[2]-dint[1])
  expect_equal(windingNumbers(4,11,1) , c(round(4*dindi), round(11*dindi)))



  })


test_that("opposed intervals work", {

  # pick a d in the predicted interval to test:
  dinoi <- function(m,n) {
    oiList <- divergenceInterval (m=m,n=n,Delta=1,scale="unity",type="opposed")
    oi <- oiList[3,]
    dinoi <- oi[1]+ 0.6*(oi[2]-oi[1])
    dinoi
  }

  m<-4;n<-5; expect_equal(GeometricalPhyllotaxis:::isOpposed(m=m,n=n,dinoi(m,n),scale="unity"),TRUE)
  m<-34;n<-55; expect_equal(GeometricalPhyllotaxis:::isOpposed(m=m,n=n,dinoi(m,n),scale="unity"),TRUE)
  m<-34;n<-19; expect_equal(GeometricalPhyllotaxis:::isOpposed(m=m,n=n,dinoi(m,n),scale="unity"),TRUE)


})

test_that("Older generating tests",{


  checkgenerating <- function(m,n,d) {
    if (m==n && m==1  & isTRUE(all.equal(d,1/2))) { return(TRUE) }
    return(abs(abs(Deltamnd(m,n,d))-1)<1e-8 ) # should be an integer so not much point to the 1e-8
  }


  checkInterval<- function (m,n,Delta){
    res <- divergenceInterval(m,n,Delta,type="generating")
    Interval <- res
    intervalL <- Interval[1]
    intervalR <- Interval[2]
    # if (Itype=="half") {
    # 	check <- checkgenerating(m,n,0.5)
    # 	if (!check) {stop(sprintf("m=%d n=%d I=(%g,%g)  fails\n",m,n,intervalL,intervalR))}
    # 	return(TRUE)
    # }
    stopifnot(intervalR>=intervalL)
    stopifnot(intervalL<=1)
    stopifnot(intervalR>=0)
    dmin <- max(0,intervalL)+1e-8
    dmax <- min(1,intervalR)-1e-8
    dnumber <- 5
    drange <- seq(dmin,dmax,length=dnumber)
    checks <-sapply(drange,function(dtest)checkgenerating(m,n,dtest))
    if(any(!checks)) {
      stop(sprintf("m=%d n=%d I=(%g,%g) d=%g fails\n",m,n,intervalL,intervalR,drange[!checks]))
    }
  }

  #checkInterval(1,1,1)
  checkInterval(1,2,1)
  checkInterval(1,2,-1)
  checkInterval(1,3,1)
  checkInterval(1,3,-1)
  checkInterval(1,4,1)
  checkInterval(1,4,-1)
  checkInterval(1,5,1)
  checkInterval(1,5,-1)
  checkInterval(2,3,1)
  checkInterval(2,3,-1)
  checkInterval(2,5,1)
  checkInterval(2,5,-1)

})

test_that("Old generating tests",{


  checkerInterval<- function (m,n,Delta){
    res <- divergenceInterval (m,n,Delta,type="generating")
    Interval <- res
    #	Itype <- res$Itype
    mntwo <- 2 * m *n
    Interval <- Interval * mntwo
    intervalL <- round(Interval[1] - 1 )
    intervalR <- round(Interval[2]+ 1)
    from <- seq(intervalL,intervalR)
    df <- data.frame(to=from[-length(from)],from=from[-1])
    df$mid <- (df$to+df$from)/2
    df$d <- df$mid/mntwo
    df$Delta <- sapply(df$d, function(d)Deltamnd(m,n,d))
    df$Generating <-     abs(df$Delta) ==1
    df$dm <- sapply(df$d, function(d) (d*m)-round(d*m))
    df$dn <- sapply(df$d, function(d) (d*n)-round(d*n))
    df$ShouldBeGenerating <- df$d > Interval[1]/mntwo & df$d< Interval[2]/mntwo
    df$OK <- df$Generating == df$ShouldBeGenerating
    all(df$OK)
  }


  expect_equal(checkerInterval(4,5,1),TRUE)
  expect_equal(divergenceInterval (m=4,n=5,scale="twomn",type="generatingopposed",Delta=1),c(30,32))
})

test_that("van Iterson",{
  expect_equal(hvanIterson(1,2,c(.4,.4)),c(.2,.2))


  expect_equal(intervalDelete(c(0,1),c(-2,-1)),c(0,1))
  expect_equal(intervalIntersect(c(0,3),c(2,1)),c(1,2))

})
