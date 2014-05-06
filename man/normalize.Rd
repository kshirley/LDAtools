\name{normalize}
\alias{normalize}
\title{Normalize}
\usage{
normalize(x)
}
\arguments{
  \item{x}{a vector or a data frame or an array or NULL.}
}
\description{
This function is an alias for x/sum(x).
}
\examples{
x <- c(1, 1, 2, 2, 3)
normalize(x)
}

