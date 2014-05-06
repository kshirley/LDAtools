\name{remap.terms}
\alias{remap.terms}
\title{This function replaces instances of specified terms with existing
or new terms in a corpus of documents}
\usage{
remap.terms(term.map = data.frame(), term.id = integer(),
  vocab = character())
}
\arguments{
  \item{term.map}{data.frame with two columns, the first of
  which contains the terms to be replaced, and the second
  of which contains their replacements. If a replacement is
  not in the current vocabulary, it will be added to the
  vocabulary.}

  \item{term.id}{an integer vector containing the term ID
  number of every token in the corpus. Should take values
  between 1 and W, where W is the number of terms in the
  vocabulary.}

  \item{vocab}{a character vector of length W, containing
  the terms in the vocabulary. This vector must align with
  \code{term.id}, such that a term.id of 1 indicates the
  first element of \code{vocab}, a term.id of 2 indicates
  the second element of \code{vocab}, etc.}
}
\value{
Returns a list of length two. The first element,
\code{new.vocab}, is a character vector containing the new
vocabulary. The second element, \code{new.term.id} is the
new vector of term ID numbers for all tokens in the data,
taking integer values from 1 to the length of the new
vocabulary.
}
\description{
After tokenization, use this function to replace all
occurrences of a given term with a new term, whether the
new term is in the existing vocabulary or not.
}

