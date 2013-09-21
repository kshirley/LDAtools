\name{preprocess}
\alias{preprocess}
\title{Filter raw documents according to various options}
\usage{
  preprocess(data, exact = NULL, partial = NULL,
    subs = NULL, stopwords = NULL, cutoff = 2,
    verbose = FALSE, quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus.
  Each element should correspond to a 'document'.}

  \item{exact}{a character vector in which each element is
  a string, phrase, or longer snippet of text that you wish
  to discard (if the element(s) match the entire content of
  a document).}

  \item{partial}{a character vector in which each element
  is a string, phrase, or longer snippet of text that you
  wish to discard (if the element(s) match a subset of a
  document).}

  \item{subs}{character vector where the odd-numbered
  element(s) are removed from the corpus and the subsequent
  even-numbered element are inserted in their place.}

  \item{stopwords}{character vector of words that should be
  ignored.}

  \item{cutoff}{The number of times a token should appear
  in order to be considered for prediction.}

  \item{verbose}{logical. Track the categories of partial
  and exact matches. For instance, if a document exactly
  matches the third element of \code{exact}, then the
  corresponding value of \code{category} is 3.}

  \item{quiet}{logical. Should a summary of the documents
  being trimmed be reported?}
}
\value{
  Returns a list of length four. The first element,
  \code{word.id}, is an integer vector with a unique value
  for each token in the corpus. The numbering corresponds
  to the ordering of third element \code{vocab}. The second
  element, \code{doc.id}, is an integer vector which
  corresponds to the document each token belongs to. The
  fourth element, \code{category} has length equal to the
  number of documents. If the value of an element in this
  vector is 0, then the corresponding document was
  retained. Otherwise, it was discarded. If the value is
  positive, it was an exact or partial match and if
  \code{verbose == TRUE} then the value points to the
  relevant element of \code{exact} or \code{partial}. If
  the value is -1, then the document had no relevant words
  after removing \code{stopwords} and applying the
  \code{cutoff}.
}
\description{
  Conduct a series of preprocessing techniques. By default,
  a very limited amount of preprocessing will occur (just
  basic stuff like remove documents that are blank or NA).
  The user must enter character vector arguments such as
  \code{exact}, \code{partial}, \code{stopwords} and/or
  \code{subs} in order to conduct more advanced techniques.
}
\examples{
data(APcorpus) #Dollars and dates already removed
data(stopwords)
input <- filter(data=text, stopwords=stopwords)
input <- filter(data=text, stopwords=stopwords, verbose=TRUE)
}

