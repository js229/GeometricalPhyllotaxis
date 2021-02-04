library(GeometricalPhyllotaxis)
context("Tests of inverse matrix creation")

#
#
# test_that("Make matrix from genetic",{
#   PG.test.R <- newPhyllotaxisGenetic(Rise=2.5,Divergence=1.1,Rho=0.795775,Jugacy=1)
#   PGM.test.R <- as(PG.test.R,"PhyllotaxisMatrix")
#   expect_is(PGM.test.R,"PhyllotaxisMatrix")
#
#   PG.test <- newPhyllotaxisGenetic(Rise=1,Divergence=2.4,Rho=0.795775,Jugacy=1,L=10,lhsOrigin=FALSE,bottomOrigin=FALSE)
#   #Does it roundtrip? Yes.
#   PMG.test <- as(PG.test,"PhyllotaxisMatrix")
#   PGPG.test <- (as(PMG.test,"PhyllotaxisGenetic"))
#   PMG.test.R <- as(PG.test.R,"PhyllotaxisMatrix")
#   PGPG.test.R <- as(PMG.test.R,"PhyllotaxisGenetic")
#   expect_equal(PGPG.test.R ,PG.test.R)
#   expect_equal(PGPG.test ,PG.test)
#
# })
#
# test_that("Parastichy points appear in the right place", {
#   library(grid)
#   PG.test <- newPhyllotaxisGenetic(Rise=1,Divergence=2.4,Rho=0.795775,Jugacy=1,L=10,lhsOrigin=FALSE,bottomOrigin=TRUE)
#   plotPhyllotaxis (PG.test,doRepeat=0,doNumbers=FALSE,doAxes=FALSE,pointSize=.3)
# })
