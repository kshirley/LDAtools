\name{ent}
\alias{ent}
\title{Entropy}
\usage{
ent(x)
}
\arguments{
  \item{x}{a vector that contains either a normalized or
  un-normalized vector of probabilities of a discrete
  (categorical) distribution}
}
\description{
This function calculates the entropy of a discrete
(categorical) probability distribution
}
\examples{
x <- c(1, 1, 2, 2, 3)
ent(x)
y <- c(1/9, 1/9, 2/9, 2/9, 3/9)
ent(y) # should be same as ent(x)
z <- c(1/3, 1/3, 1/6, 1/6)
ent(z)
z2 <- c(1/3, 1/3, 1/6, 1/6, 0)
ent(z2) # should be the same as ent(z), since this version of entropy assumes 0*log(0) = 0
}

