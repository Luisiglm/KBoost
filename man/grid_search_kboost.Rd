% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grid_search_kboost.R
\name{grid_search_kboost}
\alias{grid_search_kboost}
\title{Function to perform a grid search and find the hyperparameters.}
\usage{
grid_search_kboost(dataset, vs, gs, ite)
}
\arguments{
\item{dataset}{1 for IRMA or 2 for DREAM4 multifactorial.}

\item{vs}{The range of values of v. All values need to be between 0 and 1.}

\item{gs}{The range of values of g. All values need to be larger than 0.}

\item{ite}{An integer that is the number of iterations, fixed in this case.}
}
\value{
list with auprs and aurocs of different values of vs and gs and ite.
}
\description{
Function to perform a grid search and find the hyperparameters.
}
\examples{
res = grid_search_kboost(1,c(0.1,0.5,1),c(1,10,60,100),3)

}
