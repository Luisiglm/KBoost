% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/d4_mfac.R
\name{d4_mfac}
\alias{d4_mfac}
\title{Function to obtain the AUPR and AUROC in the DREAM4 Multifactorial Challenge.}
\usage{
d4_mfac(v, g, ite, write_res)
}
\arguments{
\item{v}{a number between 0 and 1 that is the shrinkage parameter}

\item{g}{a number larger than 0, width parameter for the RBF Kernel}

\item{ite}{an integer with number of iterations.}

\item{write_res}{a logical to indicate if the tables should be written.}
}
\value{
list with auroc and auprs of the DREAM4 multifactorial challenge.
}
\description{
Function to obtain the AUPR and AUROC in the DREAM4 Multifactorial Challenge.
}
\examples{
res = d4_mfac()

}
