\name{preprocess.newdocs}
\alias{preprocess.newdocs}
\title{Filter raw version of new documents based on previously fit model}
\usage{
  preprocess.newdocs(data = character(),
    vocab = character(), exact = NULL, partial = NULL,
    subs = NULL, verbose = FALSE, quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus.
  Each element should correspond to a 'document'.}

  \item{vocab}{a character vector containing the vocabulary
  of the fitted topic model from which we will estimate the
  topic proportions for the new documents entered in
  \code{data}.}

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

  \item{verbose}{logical. Track the categories of partial
  and exact matches. For instance, if a document exactly
  matches the third element of \code{exact}, then the
  corresponding value of \code{category} is 3.}

  \item{quiet}{logical. Should a summary of the documents
  being trimmed be reported?}
}
\value{
  Returns a list of length three. The first element,
  \code{word.id}, is an integer vector with a unique value
  for each token in the corpus. The numbering corresponds
  to the ordering of \code{vocab}. The second element,
  \code{doc.id}, is an integer vector which corresponds to
  the document each token belongs to. The third element,
  \code{category} has length equal to the number of
  documents. If the value of an element in this vector is
  0, then the corresponding document was retained.
  Otherwise, it was discarded. If the value is positive, it
  was an exact or partial match and if \code{verbose ==
  TRUE} then the value points to the relevant element of
  \code{exact} or \code{partial}. If the value is -1, then
  the document had no relevant words after removing
  \code{stopwords} and applying the \code{cutoff}.
}
\description{
  This function performs the same preprocessing steps as
  the function \code{filter()}, except this time for a set
  of new documents whose topic proportions we wish to
  estimate given the topics from a previously fit model.
  The key difference is that the vocabulary won't be
  constructed based on the words that occur in the new
  documents, but rather it will be entered as input from
  the previously fit model.
}
\examples{
data(APcorpus) #Dollars and dates already removed
data(stopwords)
input <- filter(data=text, stopwords=stopwords)
input <- filter(data=text, stopwords=stopwords, verbose=TRUE)
}

