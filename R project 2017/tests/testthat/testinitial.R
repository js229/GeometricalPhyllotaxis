library(GeometricalPhyllotaxis)
context("Tests from package creation")

test_that("Rise bug fixed", {
  # there was a bug when Rise>1
  PG.test.R <- newPhyllotaxisGenetic(Rise=2.5,Divergence=1.1,Rho=0.795775,Jugacy=1)
  expect_equal(getRise(PG.test.R),2.5)
  expect_equal(getJugacy(PG.test.R),1)
})

test_that("Create objects", {
  expect_is(new("Cylinder",Rho=0.795775),"Cylinder")
})


test_that("Parastichy points appear in the right place", {
  library(grid)
  PG.test <- newPhyllotaxisGenetic(Rise=1,Divergence=2.4,Rho=0.795775,Jugacy=1,L=10,origin=c(0,0))
  plotPhyllotaxis (PG.test,doRepeat=0,doNumbers=FALSE,doAxes=FALSE,pointSize=.3)
})
