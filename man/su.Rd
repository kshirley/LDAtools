\name{su}
\alias{su}
\title{Sort the unique elements in a vector}
\usage{
  su(x, ...)
}
\arguments{
  \item{x}{a vector or a data frame or an array or NULL.}

  \item{...}{expressions evaluated in the context of
  \code{x} and then fed to \code{\link{unique}}}
}
\description{
  This function is an alias for sort(unique()).
}
\examples{
x <- c(1, 1, 2, 2, 3)
su(x)
}

