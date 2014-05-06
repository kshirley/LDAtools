\name{preprocess.newdocs}
\alias{preprocess.newdocs}
\title{Preprocess raw version of new documents based on previously fit model}
\usage{
preprocess.newdocs(data = character(), vocab = character(), exact = NULL,
  partial = NULL, subs = NULL, verbose = FALSE, quiet = FALSE)
}
\arguments{
  \item{data}{a character vector containing the raw corpus,
  where each element is a document.}

  \item{vocab}{a character vector containing the vocabulary
  of the fitted topic model from which we will estimate the
  topic proportions for the new documents entered in
  \code{data}.}

  \item{exact}{a (case-sensitive) character vector in which
  each element is a string, phrase, or longer snippet of
  text that results in a document being discarded from the
  data if the entire document matches an element of
  \code{exact}.}

  \item{partial}{a (case-sensitive) character vector in
  which each element is a string, phrase, or longer snippet
  of text that results in a document being discarded from
  the data if any part of the document matches an element
  of \code{partial}.}

  \item{subs}{character vector of regular expressions where
  the odd-numbered element(s) are removed from the corpus
  and the subsequent even-numbered element are inserted in
  their place. These substitutions are performed using the
  \code{gsub()} function after forcing the raw text to
  lowercase.}

  \item{verbose}{logical. If set to TRUE the function will
  retain the indices of the elements of \code{exact} and
  \code{partial} that were matched. For instance, if a
  document exactly matches the third element of
  \code{exact}, then the corresponding value of
  \code{category} will be 3, if \code{verbose = TRUE}}

  \item{quiet}{logical. Should a summary of the
  preprocessing steps be printed to the screen?}

  \item{stem}{logical. Should the porter stemmer be used to
  stem the tokens in the vocabulary?}

  \item{hash}{a length-1 character vector indicating the
  prefix of substitution replacements that should be
  replaced with a '#' symbol after tokenizing. Set to "ent"
  by default, where "ent" stands for "entity", and is often
  used as a prefix to a substitution replacement for a
  class of terms, like dollar amounts ("entdollaramount")
  and timestamps ("entdatestamp", "enttimeofday"), etc.}
}
\value{
Returns a list of length three. The first element,
\code{term.id}, is an integer vector containing the index
in the vocabulary of each token in the corpus. If the 4th
token in the corpus is "tree" and "tree" is the 50th
element of the vocabulary, then the 4th element of term.id
will be 50, for example. The second element, \code{doc.id},
is an integer vector which corresponds to the document each
token belongs to. The third element, \code{category} has
length equal to the number of documents. If the value of an
element in this vector is 0, then the corresponding
document was retained. Otherwise, it was discarded. If the
value is positive, it was an exact or partial match and if
\code{verbose == TRUE} then the value points to the
relevant element of \code{exact} or \code{partial}. If the
value is -1, then the document contained no tokens in the
vocabulary.
}
\description{
This function performs the same preprocessing steps as the
function \code{preprocess()}, except this time for a set of
new documents whose topic proportions we wish to estimate
given the topics from a previously fit model. The key
difference is that the vocabulary won't be constructed
based on the words that occur in the new documents, but
rather it will be entered as input from the previously fit
model.
}

