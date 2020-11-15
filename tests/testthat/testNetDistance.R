library(KBoostv)

# Check that the net distance function works.

# We will check Outdegree
test_that("The Distance between nodes is calculated correctly", {

  # Use the IRMA Net as an example
  irma_net = KBoostv::IRMA_Gold
  d = net_dist_bin(irma_net,1:5,0.1)
  # Write down manually calculated distance
  d_ = matrix(0,5,5)
  d_[,1] = c(0,1,2,3,3)
  d_[,2] = c(2,0,1,2,2)
  d_[,3] = c(1,2,0,1,1)
  d_[,4] = c(1,2,3,0,4)
  d_[,5] = c(Inf, Inf, Inf, Inf, 0)
  # Check that the distance between the same nodes is the same.
  expect_equal(sum(d_==d), 25)

})

#Let's use another example and check.

test_that("The Distance between nodes is calculated correctly", {

  # Use the IRMA Net as an example
  ex_net = matrix(0,3,3)
  ex_net[,1] = c(0,1,0)
  ex_net[,2] = c(1,0,1)
  ex_net[,3] = c(0,1,0)
  d = net_dist_bin(ex_net,1:3,0.1)
  # Write down manually calculated distance
  d_ = matrix(0,3,3)
  d_[,1] = c(0,1,2)
  d_[,2] = c(1,0,1)
  d_[,3] = c(2,1,0)
  # Check that the distance between the same nodes is the same.
  expect_equal(sum(d_==d), 9)

})
