findAGeneratingPair <- function (d) {
  p0 <- c(1,0)
  p1 <- c(d,1)
  return( matrix(c(p0,p1),byrow=TRUE,ncol=2))
}

DeltaGP <- function(L,d) {
  m <- L[1,2];n <-L[2,2]
  Deltamnd(m,n,d)
}

CylindricalCoordinates <- function(L,h) {
  L %*% matrix(c(1,0,0,h),ncol=2) #sigh
}

ParastichyVectors <- function(L,h) {
  mnvec <- CylindricalCoordinates(L,h)
  list(mnvec[1,], mnvec[2,])
}


ParastichyVectorLengths <- function(L,h) {
  mnvec <- ParastichyVectors(L,h)
  sapply(mnvec,function(x)sum(x*x))
}

UShorterThanV <- function(L,h){
  uvLengths <- ParastichyVectorLengths(L,h)
  return(uvLengths[1]<uvLengths[2])
}


UdotV <- function(L,h) {
  mnvec <- ParastichyVectors(L,h)
  mvec <- mnvec[[1]]; nvec <- mnvec[[2]]
  res <- (mvec %*% nvec);dim(res)<-NULL
  res
}


mirrorV <- function(L){
  L[2,]<- -L[2,]
  L
}

dot  <- function(u,v) { sum(u*v)}
mod <- function(u) { sqrt(dot(u,u))}

isPrincipalPair <- function(L,h) {
  mnvec <- ParastichyVectors(L,h)
  u <- mnvec[[1]]; v <- mnvec[[2]]
  dot(u,v)>= 0 & mod(u) <= mod(v) & mod(v)<=mod(u-v)
}

swapUV <- function(L){
  L[c(2,1),]
}
reduceV <- function(L){
  cat("R");
  L[2,] <- L[2,]-L[1,]
  L
}


# u <- c(-4/72, 4/26); u2 <- u%*%u
# v <- c(13/72, 5/26); v2 <- v%*%v
# s <- v-u; s2 <- s%*% s

h <- 1/31
d <- 17/72

makeUShorter <- function(L,h) {
  if (!UShorterThanV(L,h)) {
    cat("S");
    L <- swapUV(L)
  }
  L
}

makeDotPositive <- function(L,h) {
  if (UdotV(L,h) >0) {
    cat("M");
    L <- mirrorV(L)
  }
  L
}



getPrincipalPair  <- function(d,h) {
  L <- findAGeneratingPair(d)
  L <- makeUShorter(L,h)
  L <- makeDotPositive(L,h)
  while(!isPrincipalPair(L,h)) {
    print(L)
    L <- reduceV(L)
    L <- makeUShorter(L,h)
    L <- makeDotPositive(L,h)

    if(abs(L[2,2])>40) { break()}
  }
  L
}
getPrincipalPair(d=17/72,h=1/30)




