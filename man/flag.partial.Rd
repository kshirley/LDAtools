\name{flag.partial}
\alias{flag.partial}
\title{Flag the documents that contain an occurrence of one or more strings from a pre-specified list of strings}
\usage{
  flag.partial(data, partial, D, verbose, quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus.
  Each element should correspond to a 'document'.}

  \item{partial}{a character vector in which each element
  is a string, phrase, or longer snippet of text that you
  wish to discard (if the element(s) match a subset of a
  document).}

  \item{D}{the original number of documents}

  \item{verbose}{logical vector. Track the specific
  categories of 'flagged' documents?}
}
\value{
  category an integer vector of the same length as
  \code{data}, where 0 indicates that the document did not
  contain any of the strings in \code{partial}, and an
  integer j = 1, ..., K that indicates that a document
  contained the jth element of \code{partial}.
}
\description{
  Often a data set contains documents that you wish to
  remove before you fit the LDA model, and these documents
  share a common "boilerplate" string or phrase (along with
  potentially unique information). This function can be
  used to flag those documents. Similar to the function
  \code{flag.exact}, this is a very simple operation that
  may be more useful as a signal to the user that he or she
  should visually inspect the data before running LDA (so
  as to remove documents that don't require topic modeling
  in the first place). An example of a phrase that would
  indicate a document should be discarded is "iPhone Unlock
  Approved", which is a common phrase followed by a unique
  transaction ID. Each document containing this string is
  unique (so we can't use \code{match.exact} to flag these
  documents), but we don't want to include them in a topic
  model.
}
\examples{
data <- c("Automatic Message: Account 12 ...", "Automatic Message: Account 314 ...", "A document with unknown content",
"Boilerplate text: Customer 1532 ...")
match.exact <- c("Automatic Message:", "Auto Text:", "Boilerplate text")
flag.partial(data, match.partial) # c(1, 1, 0, 3)
}
\seealso{
  flag.exact
}

