#' @title Traceplot of log-likelihood.
#'
#' @description Trace plot of the log-likelihood (the x-axis is MCMC iteration, and the y-axis is the log-likelihood). 
#' The MCMC algorithm has converged when this plot levels out and the log-likelihood time series plot looks like white noise. 
#' If the log-likelihood is still increasing at the last iteration, you should go back and run the algorithm longer (using \code{fitLDA}.)
#'
#' @param data a matrix with log-likelihood values. Each column of the matrix refers to a different MCMC chain.
#' @param zoom should a second plot of the latter iterations be output.
#' @param prop If \code{zoom == TRUE}, this number (between 0 and 1) sets the proportion of final iterations to plot.
#' @param start number of iterations (if any) performed before these chains.
#' @export
#' 
#' @examples
#' one.chain <- matrix(rnorm(1000, mean=-500, sd=100), ncol=1) #hopefully our logliklihood is essentially random noise!
#' plotLoglik(one.chain)
#' par(mfrow=c(1,2))
#' plotLoglik(one.chain, zoom=TRUE)
#' three.chains <- matrix(rnorm(3000, mean=-500, sd=100), ncol=3)
#' plotLoglik(three.chains)
#' plotLoglik(three.chains, zoom=TRUE, prop=.1, start=5000)
#' 

plotLoglik <- function(data=matrix(), zoom=FALSE, prop=0.5, start=0) {
  #par(mfrow=c(1,2))
  #on.exit(par(mfrow=c(1, 1)))
  # First, a plot of all K series, from first iteration to last:
  n <- dim(data)[1]
  base <- floor(log(-data[n, 1], base=10))
  y.vec <- data[,1]*10^(-base)
  rg <- range(y.vec)
  plot(1:n, y.vec, las=1, xlab="Iteration", ylab=paste("Log-likelihood (x 10^-", base, ")", sep=""), 
       ylim=rg, xaxt="n", type="l", col=1, yaxt="n")
  ticks <- seq(0, n, length=7)
  axis(1, at=ticks, labels=round(ticks)+start)
  ticks2 <- seq(rg[1], rg[2], length=6)
  axis(2, at=ticks2, labels=round(ticks2), las=1)
  K <- dim(data)[2]
  if (K > 1) {
    for (k in 2:K) {
      y.vec <- data[,k]*10^(-base)
      lines(1:n, y.vec, col=k)
    }
  }
  # Second, a plot of all K series, from second half of iterations only (to zoom in on this section to assess convergence):
  if (zoom) {
    second <- floor(n*(1-prop)):n
    rg <- range(data[second, ]*10^(-base))
    plot(second, data[second, 1]*10^(-base), las=1, xlab="Iteration", 
         ylab=paste("Log-likelihood (x 10^-", base, ")", sep=""), ylim=rg, xaxt="n", type="l", col=1, yaxt="n")
    ticks <- round(seq(n/2, n, length=7))
    axis(1, at=ticks, labels=ticks+start)
    ticks2 <- seq(rg[1], rg[2], length=6)
    axis(2, at=ticks2, labels=round(ticks2), las=2)
    if (K > 1){
      for (k in 2:K){
        lines(second, data[second, k]*10^(-base), col=k)
      }
    }
  }
}


#' @title Compute topic-word and document-topic probability distribution matrices, and re-label topic indices
#'
#' @description This function assumes the ordering of \code{word.id}, \code{doc.id}, \code{topic.id} matters! 
#' That is, the first element of \code{word.id} corresponds to the first element of \code{doc.id} which corresponds to the first 
#' element of \code{topic.id}. Similarly, the second element of tokens corresponds to the second element of \code{doc.id} 
#' which corresponds to the second element of \code{topic.id} (and so on). Also, the ordering of the elements of \code{vocab}
#' are assumed to correspond to the elements of \code{word.id}, so that the first element of \code{vocab} is the token with \code{word.id}
#' equal to 1, the second element of \code{vocab} is the token with \code{word.id} equal to 2, etc.
#'
#' @param word.id a numeric vector with the token id of each token occurrence in the data.
#' @param doc.id a numeric vector containing the document id number of each token occurrence in the data.
#' @param topic.id a numeric vector with a unique value for each topic.
#' @param vocab a character vector of the unique words included in the corpus. The length of this vector should match the max value of \code{word.id}.
#' @param alpha Dirichlet hyperparameter. See \link{fitLDA}.
#' @param beta Dirichlet hyperparameter. See \link{fitLDA}.
#' @param sort.topics Sorting criterion for topics. Supported methods include: "byDocs" to sort topics by the 
#' number of documents for which they are the most probable or "byTerms" to sort topics by the number of terms within topic.
#' 
#' @return A list of two matrices and one vector. The first matrix is, \code{phi.hat}, contains the distribution over tokens for each topic,
#' where the rows correspond to topics. The second matrix, \code{theta.hat}, contains the distribution over topics for each document, where
#' the rows correspond to documents. The vector returned by the function, \code{topic.id}, is the vector of sampled topics from the LDA fit, 
#' with topic indices re-labeled in decreasing order of frequency by the \code{sort.topics} argument.
#' @export
#' @examples
#' data(APinput)
#' #takes a while
#' \dontrun{o <- fitLDA(APinput$word.id, APinput$doc.id)}
#' data(APtopics) #load output instead for demonstration purposes
#' probs <- getProbs(word.id=APinput$word.id, doc.id=APinput$doc.id, topic.id=APtopics$topics,
#'                    vocab=APinput$vocab)
#' head(probs$phi.hat[,1:5])
#' head(probs$theta.hat)
#' 


getProbs <- function(word.id=numeric(), doc.id=numeric(), topic.id=numeric(), vocab=character(), 
                     alpha=0.01, beta=0.01, sort.topics=c("None", "byDocs", "byTerms"), K=integer()) {
  # K is number of topics
  stopifnot(sort.topics[1] %in% c("None", "byDocs", "byTerms"))
  #stopifnot(sort.terms[1] %in% c("none", "freq", "distinct", "saliency"))
  if (!all(sort.topics == c("None", "byDocs", "byTerms")) & length(sort.topics) > 1) stop("Please enter only one topic sorting choice")
  #if (!all(sort.terms == c("none", "freq", "distinct", "saliency")) & length(sort.terms) > 1) stop("Please enter only one term sorting choice")
  N <- length(word.id)
  stopifnot(N == length(doc.id), N == length(topic.id))
  # compute phi, the matrix of topic-word probability distributions
  df <- table(topic.id, word.id)
  W <- max(word.id)
  D <- max(doc.id)
  stopifnot(W == length(vocab))
  
  CTW <- matrix(0, K, W)
  CTW[, as.numeric(colnames(df))] <- df
  
  # compute theta, the matrix of document-topic probability distributions
  df <- table(doc.id, topic.id)
  CDT <- matrix(0, D, K)
  CDT[as.numeric(rownames(df)),] <- df

  # compute the topic to which the most tokens from each doc are assigned to
  main.topic <- max.col(CDT)
  main.topic.table <- table(main.topic)
  
  # compute posterior point estimates of phi.hat and theta.hat:
  CTW.b <- CTW + beta
  phi.hat <- CTW.b/apply(CTW.b, 1, sum)
  CDT.a <- CDT + alpha
  theta.hat <- CDT.a/apply(CDT.a, 1, sum)
  
  #set relevant names for the two matrices
  #rownames(phi.hat) <- rownames(phi.hat, do.NULL=FALSE, prefix= "Topic")
  colnames(phi.hat) <- vocab
  #colnames(theta.hat) <- colnames(theta.hat, do.NULL=FALSE, prefix= "Topic")
  
  #sort topics (if necessary)
  topic.o <- 1:K
  if (sort.topics[1] == "byDocs") {
    # order the topics by the number of documents for which they are the main topic:
    topic.o <- order(main.topic.table, decreasing=TRUE)
    main.topic <- match(main.topic, topic.o)
  }
  if (sort.topics[1] == "byTerms") {
    topic.o <- order(apply(CTW, 1, sum), decreasing=TRUE)
    main.topic <- match(main.topic, topic.o)
  }
  phi.hat <- phi.hat[topic.o, ]
  theta.hat <- theta.hat[, topic.o]
  topic.id <- match(topic.id, topic.o)

  #sort terms (if necessary)
  if (FALSE) {  # put this here to skip this block of code (Kenny - 9/18/13)
  term.o <- NULL
  if (sort.terms[1] != "none") {
    word.tab <- table(word.id)
    if (sort.terms[1] == "freq") {
      term.o <- order(word.tab, decreasing=TRUE)
    } else {
      if (sort.terms[1] %in% c("distinct", "saliency")) {
        topic.tab <- table(topic.id)
        pt <- topic.tab/sum(topic.tab)
        t.w <- t(t(phi.hat)/apply(phi.hat, 2, sum)) #P(T|w)
        kernel <- t.w*log(t.w/as.vector(pt))
        distinct <- apply(kernel, 2, sum)
      }
      if (sort.terms == "distinct") term.o <- order(distinct, decreasing=TRUE)
      if (sort.terms == "saliency") {
        pw <- word.tab/sum(word.tab)
        saliency <- pw*distinct
        term.o <- order(saliency, decreasing=TRUE)
      }
    }
  }
  } # End sort terms section
  #if (!is.null(term.o)) phi.hat <- phi.hat[,term.o]
  #if (sort.topics[1] != "byDocs") main.topic=NULL  # now we have main.topic no matter what the sort.topic option
  rownames(phi.hat) <- 1:K
  colnames(theta.hat) <- 1:K
  colnames(phi.hat) <- vocab
  return(list(phi.hat=phi.hat, theta.hat=theta.hat, topic.id=topic.id, main.topic=main.topic, topic.order=topic.o))
}


#' @title Compute distinctiveness and saliency of the words in the vocabulary for a given topic model
#'
#' @description This function computes 'distinctiveness' and 'saliency' as defined in Chuang, et al (2012), 'Termite'.
#'
#' @param word.id an integer vector with the token id of each token occurrence in the data.
#' @param topic.id an integer vector with the topic id for each token in the data.
#' @param phi.hat a numeric (K x W) matrix with topic distributions on each of the K rows
#' 
#' @return A list of two length-W vectors, the distinctiveness and saliency of each of the W tokens in the vocab
#'
#' @export
#' 

token.rank <- function(word.id=integer(), topic.id=integer(), phi.hat=numeric()) {
  word.tab <- table(word.id)
  topic.tab <- table(topic.id)
  pt <- topic.tab/sum(topic.tab)
  t.w <- t(t(phi.hat)/apply(phi.hat, 2, sum)) #P(T|w)
  kernel <- t.w*log(t.w/as.vector(pt))
  distinct <- as.numeric(apply(kernel, 2, sum))
  pw <- word.tab/sum(word.tab)
  saliency <- as.numeric(pw*distinct)
  return(list(distinct=distinct, saliency=saliency))
}



#' @title Plot probable tokens for a given topic
#'
#' @param phi numeric vector with the probability of each token for a given topic.
#' @param vocab character vector of the vocabulary for the corpus
#' @param n.token the number of tokens to plot, where the default is \code{n.token = 20}.
#' @param lambda a parameter between 0 and 1 to control how tokens are ranked within topics
#' @param p the marginal probabilities of the tokens in the vocabulary
#' @param ... additional arguments to the plot() function
#'
#' @details The ranking of tokens within topics is based on a weighted average of the probability of a 
#' token (given the topic) and the lift, where the lift of a token is defined as the probability of the token 
#' (given the topic) divided by the marginal probability of the token (i.e. across all topics). 
#' The ranking that determines the top \code{n.token} tokens to plot is simply 
#' \code{lambda * log(p(token)) + (1 - lambda) * log(p(token | topic)/p(token))}.
#' 
#' Note: the ordering of \code{phi}, \code{vocab}, and \code{p} must be the same 
#' (i.e. the nth element of each vector must correspond to the same token)
#' 
#' @export
#'
#' @examples
#' data(APinput)
#' data(APtopics) #load output instead for demonstration
#' probs <- getProbs(word.id=APinput$word.id, doc.id=APinput$doc.id, topic.id=APtopics$topics, 
#'                vocab=APinput$vocab)
#'  #THE ORDERING OF phi, vocab and p MUST MATCH!
#' tokens <- factor(APinput$vocab[APinput$word.id], levels=colnames(probs$phi.hat))
#' token.tab <- table(tokens)
#' p <- token.tab/sum(token.tab)
#' plotTokens(phi=probs$phi.hat[1,], vocab=names(p), n.tokens=30, lambda=1/3, p)
#' # plot all the topics!
#' \dontrun{
#'  for (i in seq_along(probs$phi.hat[,1])) {
#'    plotWords(probs$phi.hat[i,], tokens)
#'  }
#' }

plotTokens <- function(phi=vector(), vocab=character(), n.tokens=20, lambda=0.5, p=vector(), ...) {  
  # draw figure of the top words for each topic:
  if (length(dim(phi)) > 1) {
    warning("phi should be a vector. Only the first topic will be used for plotting.")
    if (which.max(dim(phi)) == 2) { #guess that max dimension holds the terms
      phi <- phi[1,]
    } else {
      phi <- phi[,1]
    }
  }
  if (length(vocab) == 0) warning("vocab is missing.")
  if (length(p) == 0) warning("marginal probability vector for tokens is missing.")
  if (length(phi) != length(vocab)) warning("mismatch between length of probability vector for this topic and length of vocab")
  if (length(phi) != length(p)) warning("mismatch between length of probability vector for this topic and length of marginal probability vector for all tokens")
  #if (length(vocab) != length(p)) warning("mismatch between length of probability vector for this topic and length of vocab")
  #vocab <- unique(tokens)
  #if (length(vocab) == 0) {
  #  vocab <- names(phi)
  #  if (length(vocab) == 0) {
  #    warning("Vocabulary (names of phi) are missing. Please supply them.")
  #  }
  #}
  #N <- length(tokens)
  #k <- length(phi)
  W <- length(vocab)
  
  # Get overall proportion for each word
  #p <- table(factor(tokens, levels=vocab))/N
  #p <- p[names(phi)]
  #stopifnot(all(names(p) == names(phi))) #order matters!
  w <- lambda*log(phi) + (1 - lambda)*log(phi/p)
  o <- order(w, decreasing=TRUE)
  top <- phi[o][1:n.tokens]
  topwords <- vocab[o][1:n.tokens]
  topprobs <- as.numeric(top)
  #po <- p[o]
  #p.vec <- po[names(po) %in% topwords]
  p.vec <- p[o][1:n.tokens]
  
  # set up basic plotting variables:
  x.vec <- c(0, cumsum(nchar(topwords))[-n.tokens]) + 1:n.tokens - 1
  y.vec <- topprobs
  yl <- "Probability"
  xl <- "Word Rank"
  xr <- sum(nchar(topwords)) + n.tokens
  midwords <- x.vec + nchar(topwords)/2 - 2
  # start the plot:
  plot(x.vec, y.vec, type="n", las=1, ylab=yl, xlab=xl, ylim=c(0, max(y.vec, p.vec)), xlim=c(0, xr), xaxt="n", ...)
  abline(h=0, col=gray(0.7))
  axis(1, at=midwords, labels=1:n.tokens)
  # plot the word probabilities within this topic:
  for (j in 1:n.tokens) {
    for (k in 1:nchar(topwords[j])) {
      text(x.vec[j] + k - 1, y.vec[j], substr(topwords[j], k, k), family="mono", pos=2)
    }
  }
  # plot the word probabilities across topics in light gray:
  for (j in 1:n.tokens){
    for (k in 1:nchar(topwords[j])){
      text(x.vec[j] + k - 1, p.vec[j], substr(topwords[j], k, k), family="mono", col=gray(0.7), pos=2)
    }
  }
  legend("topright", inset=0.01, col=c(1, gray(0.7)), legend=c("P(Word | Topic)", "P(Word) across all topics"), lty=1, lwd=2)
}

#' Find representative documents for each topic
#' 
#' @param theta Matrix of document-topic probabilities. Could be taken from the output of \link{getProbs}
#' @param docs A character vetor where each element is a document from the corpus. 
#' The length should be equal to the first dimension of \code{theta}.
#' @param n The number of documents to return within each topic.
#' @export
#' @examples
#' \dontrun{
#'  data(APinput)
#'  data(APtopics)
#'  data(APcorpus)
#'  probs <- getProbs(APinput$word.id, APinput$doc.id, APtopics$topic, APinput$vocab, sort.topics="byTerms")
#'  top.docs <- topdocs(probs$theta.hat, APcorpus[APinput$category == 0], n=5)
#'  #write the file for uploading to shiny...
#'  write.table(top.docs, file=paste0(getwd(), "/top5docs.txt"), sep="\t", row.names=FALSE)
#'  #save(top.docs, file="~/LDAtool/data")
#'  #sanity check (peaks on plot below should line up on same topic)
#'  corpus <- APcorpus[APinput$category == 0]
#'  idx <- which(corpus %in% top.docs[,"Topic1"])
#'  plot(probs$theta.hat[idx[1],], type="l")
#'  for (i in 2:length(idx)) {
#'    lines(probs$theta.hat[idx[i],], col=i)
#'  }
#' }

topdocs <- function(theta, docs, n=30){
  stopifnot(dim(theta)[1] == length(docs))
#   indices <- apply(theta, 2, function(x) which(order(x, decreasing=TRUE) <= n))
#   df <- data.frame(array(docs[indices], dim=c(n, dim(theta)[2])), stringsAsFactors=FALSE)
#   names(df) <- colnames(theta)
#   return(df)
  k <- dim(theta)[2]
  tops <- data.frame(matrix(NA, nrow=n, ncol=k))
  for (i in 1:k){
    o <- order(theta[,i], decreasing=TRUE)
    top <- docs[o][1:n]
    tops[,i] <- top
  }
  names(tops) <- colnames(theta)
  return(tops)
}



#' Compute symmetric version of Kullback-Leibler (KL) divergence between two categorical distributions
#' 
#' @param x The vector of probabilities in the first distribution
#' @param y The vector of probabilities in the second distribution
#'
#' @export

# symmetric version of kullback-leibler divergence:
KL <- function(x, y) {
  #if (length(x) != length(y)) break # must have same support
  #if (sum(x) !=1 | sum(y) != 1) break # probabilities must sum to one
  #if (sum(x <= 0) > 0 | sum(y <= 0) > 0) break # probabilities must be greater than zero
  #if (sum(x >= 1) > 0 | sum(y >= 1) > 0) break # probabilities must be less than one
  0.5*sum(x*log(x/y)) + 0.5*sum(y*log(y/x))
}


#' Create a list of required objects from fitted topic model to visualize using d3
#' 
#' @param text A character vector of all training documents used to fit the LDA model
#' @param doc.id An integer vector of the document ID numbers for each token occurrence in the data
#' @param word.id An integer vector of the token ID numbers for each token occurrence in the data
#' @param topic.id An integer vector of the topic ID numbers for each token occurrence in the data 
#' (from the fitted topic model)
#' @param vocab A character vector of the tokens in the vocabulary
#' @param K The number of topics in the fitted topic model
#' @param k.clusters The number of clusters into which to group the topics
#' @param lambda A number in [0,1] to govern the weighted average that defines a given token's relevance
#' for a given topic. Default to 0.5.
#' @param n.terms The number of terms to display on the right panel of the interactive visualization for each topic
#' @param n.docs The number of example documents to display for each topic
#'
#' @export
#' @examples
#'
#' data(APinput, package="ldatools")
#' names(APinput)  # [1] "word.id"  "doc.id"   "vocab"    "category"
#' data(APtopics, package="ldatools")
#' names(APtopics)  # [1] "topics" "loglik"
#' data(APcorpus)
#' length(APcorpus)  # [1] 2250
#'
#' # select only the documents that were used in training the model (i.e. remove a few very short, or empty docs)
#' text <- APcorpus[APinput$category == 0]
#'
#' # Run the function and write the JSON file:
#' z <- jsviz(text=text, doc.id=APinput$doc.id, word.id=APinput$word.id, topic.id=APtopics$topics, vocab=APinput$vocab, 
#'           K=30, k.clusters=1, lambda=0.5, n.terms=30, n.docs=10)
#' 
#' # Write the list to a JSON object and place in a directory from which to serve the d3 webpage:
#' library(RJSONIO)
#' z.out <- toJSON(z)
#' cat(z.out, file="path/lda.json")
#' # Now serve index.html from path/

jsviz <- function(text=character(), doc.id=integer(), word.id=integer(), topic.id=integer(), 
                  vocab=character(), K=integer(), k.clusters=1, lambda=0.5, n.terms=30, n.docs=10) {

  # Set some relevant local variables and run a few basic checks:
  N <- length(word.id)
  W <- length(vocab)
  D <- length(text)
  if (D != max(doc.id)) print ("Number of documents not equal to max Document ID")
  if (N != length(doc.id)) print ("Number of doc.id elements not equal to number of word.id elements")
  if (N != length(topic.id)) print ("Number of topic.id elements not equal to number of word.id elements")

  # Get estimated probability matrices theta and phi, and main.topic for each document:
  dat <- getProbs(word.id, doc.id, topic.id, vocab, alpha=0.01, beta=0.01, sort.topics="byTerms", K=K)

  # re-set topic.id and name it 'topics':
  topics <- dat$topic.id

  # compute the distinctiveness and saliency of the tokens:
  tr <- token.rank(word.id, dat$topic.id, dat$phi.hat)
  distinct <- tr$distinct
  saliency <- tr$saliency

  # set some more variables for the rest of this function:
  tokens <- vocab[word.id] # the list of tokens (printed out as actual words)
  docs <- doc.id
  words.per.doc <- as.numeric(table(doc.id))

  # Create a data.frame with the top M documents for each topic (with at least w.min words per document)
  # Data.frame will have two columns, first column is labels of the form "Topick" with k=topic id number:
  w.min <- 5 # minimum number of words per document for showing exmaple documents
  top.docs <- rep("", n.docs*K)
  for (k in 1:K){
    sel <- dat$main.topic == k & words.per.doc > w.min
    o <- order(dat$theta.hat[sel, k], decreasing=TRUE)
    top.docs[(k - 1)*n.docs + 1:min(sum(sel), n.docs)] <- text[sel][o][1:min(sum(sel), n.docs)]
  }

  ### 1 ### The first main object to be returned by this function:
  doc.df <- data.frame(Category=paste("Topic", rep(1:K, each=n.docs), sep=""), Document=as.character(top.docs), stringsAsFactors=FALSE)

  #calculate 'distance' between topics given the top 450 most frequent tokens
  avg <- as.numeric(table(word.id))/length(word.id)
  phi.hat <- dat$phi.hat

  # set distance computation cutoff at minimum number of tokens whose cumulative probability exceeds 80%
  o <- order(table(word.id), decreasing=TRUE) # order tokens by overall frequency
  dist.cutoff <- min(which(cumsum(avg[o]) > 0.80))
  phi <- phi.hat[, o[1:dist.cutoff]]
  d <- dist(phi, KL)
  fit.cmd <- cmdscale(d, k=2)
  x <- fit.cmd[,1]
  y <- fit.cmd[,2]
  lab <- gsub("Topic", "", names(x))
  loc.df <- data.frame(x, y, topics=lab, stringsAsFactors=FALSE)

  # compute % of tokens that come from each topic:
  p.topics <- table(topics)/length(topics)
  topics.df <- data.frame(topics=1:K, Freq=as.numeric(p.topics*100))

  # join the MDS location data.frame with the Topic Frequency data.frame
  mds.df <- data.frame(loc.df, Freq=topics.df[, "Freq"])  

  # workaround errors if no clustering is done (ie, input$kmeans == 1)
  mds.df$cluster <- 1
  centers <- data.frame(x=0, y=0)
  if (k.clusters > 1) {  # and clustering info (if relevant)
    cl <- kmeans(cbind(x, y), centers=k.clusters)
    mds.df$cluster <- factor(cl$cluster)
    centers <- data.frame(cl$centers)
  }

  # map tokens to a cluster
  frame <- data.frame(tokens, topics, stringsAsFactors=FALSE)
  topic.cl.ID <- mds.df[c("topics", "cluster")]  #why does subsetting here make the topic -> cluster assignment go haywire???
  #framed <- plyr::join(frame, topic.cl.ID, by="topics")
  framed <- data.frame(frame, cluster=topic.cl.ID[frame[, "topics"], "cluster"])
  # I think the second way is faster.


  ##############################################################################
  ### Create a df with the info neccessary to make the default OR new bar chart when selecting a topic or cluster.
  ### This functionality requires that we keep track of the top input$nTerms within each cluster and topic (as well as overall).
  ### The ranking of words under a topic is done via a weighted average of the lift and probability of the word given the topic.
  ### The ranking of words under a cluster is done via a similar weighted average (summing over the relevant topics)

  phi.t <- t(dat$phi.hat)
  weight <- lambda*log(phi.t) + (1 - lambda)*log(phi.t/avg) 
    
  #get the most relevant terms for each topic:
  top.terms <- NULL
  for (i in 1:K) {
    weights <- weight[, i]
    o <- order(weights, decreasing=TRUE)
    terms <- rownames(phi.t)[o][1:n.terms]
    top.terms <- c(top.terms, terms)
  }
  term.labs <- rep(paste0("Topic", 1:K), each=n.terms)
  topic.df <- data.frame(Term=top.terms, Category=term.labs, stringsAsFactors=FALSE)
    
  # get the top terms for each cluster:
  t.weights <- as.numeric(table(topics))/N
  clust.terms <- NULL
  if (k.clusters == 1) {
    # if no clustering is done, we don't want to change the 'most informative words' upon hover
    clust.terms <- rownames(phi.t)[1:n.terms]
  } else {
    for (i in 1:k.clusters) {
      # grab topics that belong to the current cluster
      topicz <- which(topic.cl.ID[, "cluster"] == i)
      sub.phi <- phi.t[, topicz]
      sub.theta <- dat$theta.hat[, topicz]
      #only sum if multiple columns exist
      if (!is.null(dim(sub.phi))) {
        sub.phi <- apply(t(sub.phi)*t.weights[topicz], 2, sum)  # weighted by topic term frequency
        sub.theta <- apply(t(sub.theta)*t.weights[topicz], 2, sum)  # weighted by topic term frequency
      }
      weight <- lambda*log(sub.phi) + (1 - lambda)*log(sub.phi/avg)
      o <- order(weight, decreasing=TRUE)
      terms <- rownames(phi.t)[o][1:n.terms]
      clust.terms <- c(clust.terms, terms)
    }
  }
  term.labs <- rep(paste0("Cluster", 1:k.clusters), each=n.terms)
  clust.df <- data.frame(Term=clust.terms, Category=term.labs, stringsAsFactors=FALSE)
    
  # Order the terms for the "default" view by decreasing saliency:
  top.df <- data.frame(Term=vocab[order(saliency, decreasing=TRUE)][1:n.terms], Category="Default")
    
  all.df <- rbind(topic.df, clust.df, top.df)
  all.df$Freq <- 0
    
  # now we have all the top ranking words within each cluster and topic
  # next, we find the frequency of those words in each category
  all.words <- unique(all.df$Term)
  all.frame <- subset(framed, tokens %in% all.words)
  counts <- table(as.character(all.frame$tokens), all.frame$topics)
  counts2 <- table(as.character(all.frame$tokens), all.frame$cluster)
    
  for (i in 1:K) {
    idx <- which(all.df$Category == paste0("Topic", i))
    all.df$Freq[idx] <- counts[all.df$Term[idx], i]
  }
    
  for (i in 1:k.clusters) {
    idx <- which(all.df$Category == paste0("Cluster", i))
    all.df$Freq[idx] <- counts2[all.df$Term[idx], i]
  }
    
  totals <- table(as.character(all.frame$tokens))
  idx <- which(all.df$Category == "Default")
  all.df$Freq[idx] <- totals[all.df$Term[idx]]
  all.df$Total <- as.integer(totals[all.df$Term])
    
  # relative frequency (in percentages) over topics for each possible term
  probs <- t(apply(counts, 1, function(x) round(100*x/sum(x)))) # round() gets closer to 100, although sometimes over
  topic.probs <- data.frame(probs, stringsAsFactors=FALSE)
  topic.probs$Term <- rownames(probs)
  topic.table <- data.frame(Term = rep(rownames(probs), K), Topic=rep(1:K, each=length(all.words)),
                            Freq = as.numeric(as.matrix(topic.probs[, 1:K])))
  return(list(mdsDat=mds.df, mdsDat2=topic.table, barDat=all.df, docDat=doc.df,
            centers=centers, nClust=k.clusters))
}





