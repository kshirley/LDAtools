\name{flag.exact}
\alias{flag.exact}
\title{Flag the documents that exactly match a pre-specified list of strings}
\usage{
flag.exact(data, exact, verbose = FALSE, quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus.
  Each element should correspond to a 'document'.}

  \item{exact}{a character vector in which each element is
  a string, phrase, or longer snippet of text that you wish
  to discard, if the element matches the entire content of
  a document.}

  \item{verbose}{logical. Track the categories of exact
  matches. For instance, if a document exactly matches the
  third element of \code{exact}, then the corresponding
  value returned will be 3.}

  \item{quiet}{logical. Should a summary of the
  preprocessing steps be printed to the screen?}
}
\value{
category an integer vector of the same length as
\code{data}, where, if verbose=TRUE, 0 indicates that the
document did not match any of the strings in \code{exact},
and an integer j = 1, ..., K indicates that a document was
an exact match to the jth element of \code{exact}, and if
verbose=FALSE, an indicator vector of whether the document
exactly matched any of the elements of \code{exact}
(without indicating which element it matched).
}
\description{
If there are certain (typically very short) documents that
occur frequently in your data, and you wish to remove them
from the data before you fit the LDA model, this function
can be used to flag those documents. It's a trivial
operation, but it's a useful reminder that users should
visually inspect their data before running LDA (so as to
throw out documents that don't require topic modeling in
the first place).
}
\examples{
data <- c("bla bla bla", "foo", "bar", "text")
match.exact <- c("foo", "junk")
flag.exact(data, match.exact, verbose=FALSE, quiet=FALSE) # c(0, 1, 0, 0)
flag.exact(data, match.exact, verbose=TRUE, quiet=FALSE) # c(0, 2, 0, 0)
}
\seealso{
flag.partial
}

