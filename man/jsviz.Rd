\name{jsviz}
\alias{jsviz}
\title{Create JSON object from fitted topic model to visualize using d3}
\usage{
  jsviz(text = character(), doc.id = integer(),
    word.id = integer(), topic.id = integer(),
    vocab = character(), K = integer(), k.clusters = 1,
    lambda = 0.5, n.terms = 30, file = character())
}
\arguments{
  \item{text}{A character vector of all training documents
  used to fit the LDA model}

  \item{doc.id}{An integer vector of the document ID
  numbers for each token occurrence in the data}

  \item{word.id}{An integer vector of the token ID numbers
  for each token occurrence in the data}

  \item{topic.id}{An integer vector of the topic ID numbers
  for each token occurrence in the data (from the fitted
  topic model)}

  \item{vocab}{A character vector of the tokens in the
  vocabulary}

  \item{K}{The number of topics in the fitted topic model}

  \item{k.clusters}{The number of clusters into which to
  group the topics}

  \item{lambda}{A number in [0,1] to govern the weighted
  average that defines a given token's relevance for a
  given topic. Default to 0.5.}

  \item{n.terms}{The number of terms to display on the
  right panel of the interactive visualization for each
  topic}

  \item{file}{The name of the file to which to write the
  JSON object containing all the relevant information for
  the interactive visualization}
}
\description{
  Create JSON object from fitted topic model to visualize
  using d3
}
\examples{
data(APinput, package="ldatools")
names(APinput)  # [1] "word.id"  "doc.id"   "vocab"    "category"
data(APtopics, package="ldatools")
names(APtopics)  # [1] "topics" "loglik"
data(APcorpus)
length(APcorpus)  # [1] 2250

# select only the documents that were used in training the model (i.e. remove a few very short, or empty docs)
text <- APcorpus[APinput$category == 0]

# Run the function and write the JSON file:
z <- jsviz(text=text, doc.id=APinput$doc.id, word.id=APinput$word.id, topic.id=APtopics$topics, vocab=APinput$vocab,
          K=30, k.clusters=1, lambda=0.5, n.terms=30, file="path/lda.json")
# Now serve index.html from path/
}

