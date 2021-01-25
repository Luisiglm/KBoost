library(KBoost)

# Check that the net summarzing function works.

# We will check Outdegree and Indegree
test_that("The Outdegree and Indegree is calculated correctly in IRMA Net", {

  # Use the IRMA Net as an example
  irma_net = KBoost::IRMA_Gold
  s_n = net_summary_bin(irma_net)
  # Write down manually calculated Outdegree and Indegree
  o_ = c(1,1,3,1,0)
  i_ = c(2,1,1,1,1)
  o_ = sort(o_, decreasing = TRUE)
  # Check that the distance between the same nodes is the same.
  expect_equal(sum(o_==s_n$Outdegree), 5)
  expect_equal(sum(i_==s_n$Indegree), 5)

})


