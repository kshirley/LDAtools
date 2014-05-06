\name{bigram.table}
\alias{bigram.table}
\title{Compute table of bigrams}
\usage{
bigram.table(term.id = integer(), doc.id = integer(), vocab = character(),
  n = integer())
}
\arguments{
  \item{term.id}{an integer vector containing the term ID
  number of every token in the corpus. Should take values
  between 1 and W, where W is the number of terms in the
  vocabulary.}

  \item{doc.id}{an interger vector containing the document
  ID number of every token in the corpus. Should take
  values between 1 and D, where D is the total number of
  documents in the corpus.}

  \item{vocab}{a character vector of length W, containing
  the terms in the vocabulary. This vector must align with
  \code{term.id}, such that a term.id of 1 indicates the
  first element of \code{vocab}, a term.id of 2 indicates
  the second element of \code{vocab}, etc.}

  \item{n}{an integer specifying how large the bigram table
  should be. The function will return the top n most
  frequent bigrams. This argument is here because the
  number of bigrams can be as large as W^2.}
}
\value{
a dataframe with three columns and \code{n} rows,
containing the bigrams (column 2), their frequencies
(column 3), and their rank in decreasing order of frequency
(column 1). The table is sorted by default in decreasing
order of frequency.
}
\description{
This function counts the bigrams in the data. It's based on
the vector of term IDs and document IDs -- that is, the
vocabulary has already been established, and this function
simply counts occurrences of consecutive terms in the data.
}

