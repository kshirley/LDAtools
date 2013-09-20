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

#' @title Sum
#'
#' @description This function is an alias for sum(x, na.rm=TRUE).
#'
#' @param x a vector or a data frame or an array or NULL.
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' sum.na(x)

sum.na <- function(x) sum(x, na.rm=TRUE)

#' @title Norm
#'
#' @description This function is an alias for x/sum(x).
#'
#' @param x a vector or a data frame or an array or NULL.
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' norm(x)

norm <- function(x) x/sum(x)

#' @title Entropy
#'
#' @description This function calculates entropy
#'
#' @param x a vector or a data frame or an array or NULL.
#' @export
#' @examples
#' x <- c(1, 1, 2, 2, 3)
#' ent(x)

ent <- function(x) -sum(x[x!=0]/sum(x[x!=0])*log(x[x!=0]/sum(x[x!=0])))
