\name{flag.exact}
\alias{flag.exact}
\title{Flag the documents that exactly match a pre-specified list of strings}
\usage{
  flag.exact(data, exact, D, verbose = FALSE,
    quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus.
  Each element should correspond to a 'document'.}

  \item{exact}{a character vector in which each element is
  a string, phrase, or longer snippet of text that you wish
  to discard (if the element(s) match the entire content of
  a document).}

  \item{D}{the original number of documents}

  \item{verbose}{logical vector. Should a summary of the
  documents being trimmed be reported?}
}
\value{
  category an integer vector of the same length as
  \code{data}, where 0 indicates that the document did not
  match any of the strings in \code{match.exact}, and an
  integer j = 1, ..., K that indicates that a document was
  an exact match to the jth element of \code{match.exact}.
}
\description{
  If there are certain (typically very short) documents
  that occur frequently in your data, and you wish to
  remove them from the data before you fit the LDA model,
  this function can be used to flag those documents. It's a
  trivial operation, but it's a useful reminder that users
  should visually inspect their data before running LDA (so
  as to throw out documents that don't require topic
  modeling in the first place). An example in customer care
  rep notes are the phrases "Not offered" (in reference to
  deals that were specifically not offered to the customer
  during the phone call) and "Too expensive" (in reference
  to an offer that was declined because if was too
  expensive).
}
\examples{
data <- c("bla bla bla", "foo", "bar", "text")
match.exact <- c("foo", "junk")
flag.exact(data, match.exact) # c(0, 2, 0, 0)
}
\seealso{
  flag.partial
}

