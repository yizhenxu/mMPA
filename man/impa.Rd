% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/impa.R
\name{impa}
\alias{impa}
\title{Mini Pooling Test Number}
\usage{
impa(v, s, K = 5, vf_cut = 1000, lod = 0)
}
\arguments{
\item{v}{Vector of viral load.}

\item{s}{Vector of risk score with the same length of viral load.}

\item{K}{Pool size, default is 5.}

\item{vf_cut}{Cutoff for individual viral failure, default is 1000.}

\item{lod}{Vector of true VL of those undetectable, default is 0.}
}
\value{
A vector with length equals the number of pools. Each element indicates the estimated number of tests for the corresponding pool.
}
\description{
This function allows you to compute the average number of assays needed per pool within each provided pooled cluster.
}
\examples{
d = Simdata
V = d$VL # Viral Load
S = d$S # Risk Score
impa(V, S, K = 5, vf_cut = 1000, lod = 0)
}
\keyword{Pooling,}
\keyword{assays.}
\keyword{average}
\keyword{number}
\keyword{of}

