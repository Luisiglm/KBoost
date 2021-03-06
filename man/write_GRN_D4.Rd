% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_GRN_D4.R
\name{write_GRN_D4}
\alias{write_GRN_D4}
\title{Function to write output in DREAM4 Challenge Format.}
\usage{
write_GRN_D4(GRN, TFs, filename)
}
\arguments{
\item{GRN}{a GxK gene regulatory network.}

\item{TFs}{a K set of indixes of G that are TFs.}

\item{filename}{a string with the filename.}
}
\value{
a file with the network written as a file.
}
\description{
Function to write output in DREAM4 Challenge Format.
}
\examples{
data(D4_multi_1)
Net = kboost(D4_multi_1)
write_GRN_D4(Net$GRN, seq_len(100), "D4_multi_1_network.txt")
}
