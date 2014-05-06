\name{preprocess}
\alias{preprocess}
\title{Preprocess raw documents according to various options}
\usage{
preprocess(data, exact = NULL, partial = NULL, subs = NULL,
  stopwords = NULL, cutoff = 2, verbose = FALSE, quiet = FALSE,
  stem = FALSE, hash = "ent")
}
\arguments{
  \item{data}{a character vector containing the raw corpus,
  where each element is a document.}

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

  \item{stopwords}{character vector of tokens that should
  be excluded from the vocabulary.}

  \item{cutoff}{The minimum number of times a token must
  appear in the corpus in order to be included in the
  vocabulary.}

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
Returns a list of length five. The first element,
\code{term.id}, is an integer vector containing the index
in the vocabulary of each token in the corpus. If the 4th
token in the corpus is "tree" and "tree" is the 50th
element of the vocabulary, then the 4th element of term.id
will be 50. The second element, \code{doc.id}, is an
integer vector which corresponds to the document each token
belongs to. The third element, \code{vocab}, is the
vocabulary of the corpus, which contains all the terms
(i.e. unique tokens) in the data. It is sorted in
decreasing order of term frequency by default. The fourth
element, \code{category} has length equal to the number of
input documents in \code{data}. If the value of an element
in this vector is 0, then the corresponding document was
retained. Otherwise, it was discarded. If the value is
positive, it was an exact or partial match and if
\code{verbose == TRUE} then the value points to the
relevant element of \code{exact} or \code{partial}. If the
value is -1, then the document contained no tokens in the
vocabulary after removing \code{stopwords} and applying the
\code{cutoff}. The fifth element, \code{call}, is a named
list returning the arguments supplied to the
\code{preprocess} function: exact, partial, subs,
stopwords, stem, and cutoff.
}
\description{
Conduct a series of preprocessing steps on raw documents.
By default, a very limited amount of preprocessing will
occur (just basic stuff like remove documents that are
blank or NA and then tokenize documents by separating by
whitespace). The user can optionally filter documents,
perform global substitutions using regular expressions,
remove stopwords, and perform stemming.
}
\examples{
data(APcorpus)
data(stopwords)
input <- preprocess(data=APcorpus, exact=NULL, partial=NULL, subs=NULL,
                    stopwords=stopwords, cutoff=5, verbose=FALSE,
                    quiet=FALSE, stem=FALSE, hash="ent")
}

