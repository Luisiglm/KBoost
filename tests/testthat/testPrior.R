library(KBoostv)

test_that("Prior works using a Input list of Symbols", {

  Inp = c("TP53","CTCF", "MYBL2","FOXA1","MDM2","DKK1","DGAT2","ZFP64","PNPO")
  weight = 0.6
  P =  get_prior_Gerstein(Inp,1:5,weight,0.5 )
  expect_equal(dim(P)[2], 5)
  # Check the number of rows in the Prior should be equal to input.
  expect_equal(dim(P)[1], length(Inp))
})

test_that("get_tfs_human should not work if there aren't any transcription factors in input", {
  Inp = c("MDM2","DKK1","DGAT2","PNPO")
  Nomen = "Symbol"
  expect_error(get_tfs_human(Inp))
})

