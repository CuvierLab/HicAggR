% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MeanScale.R
\name{MeanScale}
\alias{MeanScale}
\title{Scale values by mean.}
\usage{
MeanScale(x)
}
\arguments{
\item{x}{<numeric>: Numerical vector.}
}
\value{
Scaled numeric vector.
}
\description{
Scale values with mean.
}
\details{
MeanScale
}
\examples{
set.seed(655213)
x <- rnorm(500, 500)
set.seed(522613)
y.num <- rnorm(500, 100)
plot(density(x), col = "red", xlim = c(min(y.num), max(x)))
lines(density(y.num), col = "green")
plot(density(MeanScale(x)), col = "red", xlim = c(min(MeanScale(y.num)), max(MeanScale(x))))
lines(density(MeanScale(y.num)), col = "green")

}
\keyword{internal}
