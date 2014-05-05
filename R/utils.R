#' @title Compute the number of unique elements in a vector
#'
#' @description This function is an alias for length(unique()).
#'
#' @param x a vector or a data frame or an array or NULL.
#' @param ... expressions evaluated in the context of \code{x} and then fed to \code{\link{unique}}
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' lu(x)
lu <- function(x, ...) {
  length(unique(x, ...))
}

#' @title Sort the unique elements in a vector
#'
#' @description This function is an alias for sort(unique()).
#'
#' @param x a vector or a data frame or an array or NULL.
#' @param ... expressions evaluated in the context of \code{x} and then fed to \code{\link{unique}}
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' su(x)

su <- function(x, ...) {
  sort(unique(x, ...))
}

#' @title Normalize
#'
#' @description This function is an alias for x/sum(x).
#'
#' @param x a vector or a data frame or an array or NULL.
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' normalize(x)

normalize <- function(x) x/sum(x)

#' @title Entropy
#'
#' @description This function calculates the entropy of a discrete (categorical) probability distribution
#'
#' @param x a vector that contains either a normalized or un-normalized vector of probabilities of a 
#' discrete (categorical) distribution
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' ent(x)
#' y <- c(1/9, 1/9, 2/9, 2/9, 3/9)
#' ent(y) # should be same as ent(x)
#' z <- c(1/3, 1/3, 1/6, 1/6)
#' ent(z)
#' z2 <- c(1/3, 1/3, 1/6, 1/6, 0)
#' ent(z2) # should be the same as ent(z), since this version of entropy assumes 0*log(0) = 0

ent <- function(x) -sum(x[x != 0]/sum(x[x != 0])*log(x[x != 0]/sum(x[x != 0])))


