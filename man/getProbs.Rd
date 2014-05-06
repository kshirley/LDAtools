\name{getProbs}
\alias{getProbs}
\title{Compute topic-word and document-topic probability distribution matrices, and re-label topic indices}
\usage{
getProbs(word.id = numeric(), doc.id = numeric(), topic.id = numeric(),
  vocab = character(), alpha = 0.01, beta = 0.01,
  sort.topics = c("None", "byDocs", "byTerms"), K = integer())
}
\arguments{
  \item{word.id}{a numeric vector with the token id of each
  token occurrence in the data.}

  \item{doc.id}{a numeric vector containing the document id
  number of each token occurrence in the data.}

  \item{topic.id}{a numeric vector with a unique value for
  each topic.}

  \item{vocab}{a character vector of the unique words
  included in the corpus. The length of this vector should
  match the max value of \code{word.id}.}

  \item{alpha}{Dirichlet hyperparameter. See
  \link{fitLDA}.}

  \item{beta}{Dirichlet hyperparameter. See \link{fitLDA}.}

  \item{sort.topics}{Sorting criterion for topics.
  Supported methods include: "byDocs" to sort topics by the
  number of documents for which they are the most probable
  or "byTerms" to sort topics by the number of terms within
  topic.}
}
\value{
A list of two matrices and one vector. The first matrix is,
\code{phi.hat}, contains the distribution over tokens for
each topic, where the rows correspond to topics. The second
matrix, \code{theta.hat}, contains the distribution over
topics for each document, where the rows correspond to
documents. The vector returned by the function,
\code{topic.id}, is the vector of sampled topics from the
LDA fit, with topic indices re-labeled in decreasing order
of frequency by the \code{sort.topics} argument.
}
\description{
This function assumes the ordering of \code{word.id},
\code{doc.id}, \code{topic.id} matters! That is, the first
element of \code{word.id} corresponds to the first element
of \code{doc.id} which corresponds to the first element of
\code{topic.id}. Similarly, the second element of tokens
corresponds to the second element of \code{doc.id} which
corresponds to the second element of \code{topic.id} (and
so on). Also, the ordering of the elements of \code{vocab}
are assumed to correspond to the elements of
\code{word.id}, so that the first element of \code{vocab}
is the token with \code{word.id} equal to 1, the second
element of \code{vocab} is the token with \code{word.id}
equal to 2, etc.
}
\examples{
data(APinput)
#takes a while
\dontrun{o <- fitLDA(APinput$word.id, APinput$doc.id)}
data(APtopics) #load output instead for demonstration purposes
probs <- getProbs(word.id=APinput$word.id, doc.id=APinput$doc.id, topic.id=APtopics$topics,
                   vocab=APinput$vocab)
head(probs$phi.hat[,1:5])
head(probs$theta.hat)
}

