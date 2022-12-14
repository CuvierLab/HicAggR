% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AddAttr.R
\name{AddAttr}
\alias{AddAttr}
\title{Add list as attributes.}
\usage{
AddAttr(x = NULL, attrs = NULL, overwrite = FALSE)
}
\arguments{
\item{x}{<any>: An object to which attributes are to be added.}

\item{attrs}{<list>: A named list of new attributes.}

\item{overwrite}{<logical>: Whether an overwrite is required on attributes with the same name.(Default FALSE)}
}
\value{
The object with new attributes.
}
\description{
Add list as attributes to any object with or without overwrite.
}
\details{
AddAttr
}
\examples{
x <- seq_len(10)
x <- AddAttr(x, list(dim = c(2, 5)))
x
x <- AddAttr(x, list(dim = c(5, 2)))
x
x <- AddAttr(x, list(dim = c(5, 2)), overwrite = TRUE)
x

}
\keyword{internal}
