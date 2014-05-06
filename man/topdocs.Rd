\name{topdocs}
\alias{topdocs}
\title{Find representative documents for each topic}
\usage{
topdocs(theta, docs, n = 30)
}
\arguments{
  \item{theta}{Matrix of document-topic probabilities.
  Could be taken from the output of \link{getProbs}}

  \item{docs}{A character vetor where each element is a
  document from the corpus. The length should be equal to
  the first dimension of \code{theta}.}

  \item{n}{The number of documents to return within each
  topic.}
}
\description{
Find representative documents for each topic
}
\examples{
\dontrun{
 data(APinput)
 data(APtopics)
 data(APcorpus)
 probs <- getProbs(APinput$word.id, APinput$doc.id, APtopics$topic, APinput$vocab, sort.topics="byTerms")
 top.docs <- topdocs(probs$theta.hat, APcorpus[APinput$category == 0], n=5)
 #write the file for uploading to shiny...
 write.table(top.docs, file=paste0(getwd(), "/top5docs.txt"), sep="\\t", row.names=FALSE)
 #save(top.docs, file="~/LDAtool/data")
 #sanity check (peaks on plot below should line up on same topic)
 corpus <- APcorpus[APinput$category == 0]
 idx <- which(corpus \%in\% top.docs[,"Topic1"])
 plot(probs$theta.hat[idx[1],], type="l")
 for (i in 2:length(idx)) {
   lines(probs$theta.hat[idx[i],], col=i)
 }
}
}

