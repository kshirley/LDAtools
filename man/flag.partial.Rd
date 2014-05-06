\name{flag.partial}
\alias{flag.partial}
\title{Flag the documents that contain an occurrence of one or more strings
 from a pre-specified list of strings}
\usage{
flag.partial(data, partial, verbose, quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus.
  Each element should correspond to a 'document'.}

  \item{partial}{a character vector in which each element
  is a string, phrase, or longer snippet of text that you
  wish to discard, if the element matches a subset of a
  document.}

  \item{verbose}{logical. Track the categories of partial
  matches. For instance, if a document partially matches
  the third element of \code{partial}, then the
  corresponding value returned will be 3.}

  \item{quiet}{logical. Should a summary of the
  preprocessing steps be printed to the screen?}
}
\value{
category an integer vector of the same length as
\code{data}, where, if verbose=TRUE, 0 indicates that the
document did not match any of the strings in
\code{partial}, and an integer j = 1, ..., K indicates that
a document was a partial match to the jth element of
\code{partial}, and if verbose=FALSE, an indicator vector
of whether the document partially matched any of the
elements of \code{partial} (without indicating which
element it matched).
}
\description{
Often a data set contains documents that you wish to remove
before you fit the LDA model, and these documents share a
common "boilerplate" string or phrase (along with
potentially unique information). This function can be used
to flag those documents. Similar to the function
\code{flag.exact}, this is a very simple operation that may
be more useful as a signal to the user that he or she
should visually inspect the data before running LDA (so as
to remove documents that don't require topic modeling in
the first place).
}
\examples{
data <- c("Automatic Message: Account 12 ...",
"Automatic Message: Account 314 ...",
"A document with unknown content",
"Boilerplate text: Customer 1532 ...")
match.exact <- c("Automatic Message:", "Auto Text:", "Boilerplate text")
flag.partial(data, match.partial, verbose=FALSE, quiet=FALSE) # c(1, 1, 0, 1)
flag.partial(data, match.partial, verbose=TRUE, quiet=FALSE) # c(1, 1, 0, 3)
}
\seealso{
flag.exact
}

