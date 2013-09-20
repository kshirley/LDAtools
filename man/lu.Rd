\name{lu}
\alias{lu}
\title{Compute the number of unique elements in a vector}
\usage{
  lu(x, ...)
}
\arguments{
  \item{x}{a vector or a data frame or an array or NULL.}

  \item{...}{expressions evaluated in the context of
  \code{x} and then fed to \code{\link{unique}}}
}
\description{
  This function is an alias for length(unique()).
}
\examples{
x <- c(1, 1, 2, 2, 3)
lu(x)
}

