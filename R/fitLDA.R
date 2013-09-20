#' @title Fit LDA model via Gibbs sampler
#'
#' @description This function implements the Gibbs sampling method described by Griffiths and Steyvers (2004). The Gibbs sampler
#' portion of the function is a call to C code. Note that we only return the latent topic assignments (for each token) from the last iteration.
#' Thus, memory limitations aren't really an issue. However, the run time is O(num.chains*n.iter*N*k) where \code{n.chains} is number of MCMC chains,
#' \code{n.iter} is the number of iterations, N is the total number of tokens in teh data, and k is the number of topics.
#' It is possible to resume a Gibbs sampler from a previous fit
#' by using the topics from that fit to initiate the next set of iterations using \code{topics.init}.
#'
#' @param word.id Unique token ID. Can be taken directly from the output of \code{filter}.
#' @param doc.id Unique document ID. Can be taken directly from the output of \code{filter}.
#' @param k number of topics.
#' @param n.chains number of MCMC chains.
#' @param n.iter number of iterations.
#' @param topics.init A vector of topics to initially assign. The Markov property of MCMC allows one to input the topic assignments from the last iteration of a previous model fit. 
#' Note that this vector should be the same length of the \code{word.id} vector times the number of chains. 
#' @param alpha Dirichlet hyperparameter
#' @param beta Dirichlet hyperparameter
#' 
#' @return A list of length two. The first element is the sampled latent topic value from the last iteration (for each token).
#' The second element is a vector with the log-likelihood values for every iteration of the gibbs sampler. 
#' 
#' @references Griffiths and Steyvers (2004). Finding Scientific Topics. Proceedings of the National Academy of Sciences. 
#'
#' @export
#' @useDynLib ldatools
#'
#' @examples
#' data(APinput)
#' #takes a while
#' \dontrun{o <- fitLDA(APinput$word.id, APinput$doc.id, k=20)}


fitLDA <- function(word.id=integer(), doc.id=integer(), k=10, n.chains=1, n.iter=1000, topics.init=NULL, alpha=0.01, beta=0.01) {
  # insert a bunch of checks here to make sure the inputs are valid:
  stopifnot(!is.null(word.id))
  stopifnot(!is.null(doc.id))
  stopifnot(length(word.id) == length(doc.id))

  # force inputs to their proper types so C doesn't crash:
  word.id <- as.integer(word.id)
  doc.id <- as.integer(doc.id)
  N <- as.integer(length(word.id))
  W <- as.integer(max(word.id))
  D <- as.integer(max(doc.id))
  # number of topics (note that T is a special character in R - using it as an object name can lead to unintended results)
  k <- as.integer(k) 
  n.chains <- as.integer(n.chains)  
  n.iter <- as.integer(n.iter)
  if (is.null(topics.init)) {
  	topics.init <- as.integer(sample(1:k, N*n.chains, replace=TRUE))
  } else {
  	topics.init <- as.integer(topics.init)
  	stopifnot(length(topics.init) == N*n.chains)
  }
  # number of log-likelihood evaluations that will be returned (depend on n.iter and thin)
  loglik <- as.double(numeric(n.chains*n.iter))
  # call the C function from R:
  f <- .C("lda", word_id_R=word.id, doc_id_R=doc.id, T_R=k, n_chains_R=n.chains, n_iter_R=n.iter, 
            topics_init_R=topics.init, alpha_R=as.double(alpha), beta_R=as.double(beta), 
            loglik_R=loglik, N_R=N, W_R=W, D_R=D, package = 'ldatools')
  topics.out <- matrix(f$topics_init_R, nrow=N, ncol=n.chains, byrow=TRUE)
  loglik.out <- matrix(f$loglik_R, nrow=n.iter, ncol=n.chains, byrow=TRUE)
  #Throw away the first element since this is always 0 and not informative
  loglik.out <- loglik.out[-1,]
  return(list(topics = topics.out, loglik = loglik.out))
}




#' @title Estimate topics for new documents using a Gibbs sampler
#'
#' @description This function estimates topic proportions for a new corpus of documents, using the
#' the vocabulary and the topic-token probability distributions from a previously fit LDA topic model. The function
#' samples the latent topics for each token in the new corpus using a Gibbs sampler, and returns the latent topics
#' from the last iteration.
#'
#' @param word.id Unique token ID. Can be taken directly from the output of \code{filter}.
#' @param doc.id Unique document ID. Can be taken directly from the output of \code{filter}.
#' @param k number of topics.
#' @param n.chains number of MCMC chains.
#' @param n.iter number of iterations.
#' @param topics.init A vector of topics to initially assign. The Markov property of MCMC allows one to input the topic assignments from the last iteration of a previous model fit. 
#' Note that this vector should be the same length of the \code{word.id} vector times the number of chains. 
#' @param alpha Dirichlet hyperparameter
#' @param phi The \code{T} x \code{W} matrix containing the topic-token probability distributions for each of the \code{T} topics in the 
#' previously fit topic model.
#' 
#' @return A list of length two. The first element is the sampled latent topic value from the last iteration (for each token).
#' The second element is a vector with the log-likelihood values for every iteration of the Gibbs sampler. 
#' 
#' @export
#' @useDynLib ldatools
#'
#' @examples
#' data(APinput)
#' #takes a while
#' \dontrun{o <- fitLDA(APinput$word.id, APinput$doc.id, k=20)}


predictLDA <- function(word.id=integer(), doc.id=integer(), k=10, n.chains=1, n.iter=1000, topics.init=NULL, alpha=0.01, phi) {
  # insert a bunch of checks here to make sure the inputs are valid:
  stopifnot(!is.null(word.id))
  stopifnot(!is.null(doc.id))
  stopifnot(length(word.id) == length(doc.id))
  #stopifnot(length(dim(phi)) == 2)
  if (length(dim(phi)) != 2) stop("phi must be a matrix (with one row per topic and one column per token in the vocabulary)")

  # force inputs to their proper types so C doesn't crash:
  word.id <- as.integer(word.id)
  doc.id <- as.integer(doc.id)
  N <- as.integer(length(word.id))
  W <- dim(phi)[2]
  D <- as.integer(max(doc.id))
  k <- as.integer(k) 
  n.chains <- as.integer(n.chains)  
  n.iter <- as.integer(n.iter)
  if (is.null(topics.init)) {
  	topics.init <- as.integer(sample(1:k, N*n.chains, replace=TRUE))
  } else {
  	topics.init <- as.integer(topics.init)
  	stopifnot(length(topics.init) == N*n.chains)
  }
  # number of log-likelihood evaluations that will be returned (depend on n.iter and thin)
  loglik <- as.double(numeric(n.chains*n.iter))
  # call the C function from R:
  f <- .C("ldapredict", word_id_R=word.id, doc_id_R=doc.id, T_R=k, n_chains_R=n.chains, n_iter_R=n.iter, 
            topics_init_R=topics.init, alpha_R=as.double(alpha), loglik_R=loglik, N_R=N, W_R=W, D_R=D, 
            phi_R=as.double(phi), package = 'ldatools')
  topics.out <- matrix(f$topics_init_R, nrow=N, ncol=n.chains, byrow=TRUE)
  loglik.out <- matrix(f$loglik_R, nrow=n.iter, ncol=n.chains, byrow=TRUE)
  #Throw away the first element since this is always 0 and not informative
  loglik.out <- loglik.out[-1,]
  return(list(topics = topics.out, loglik = loglik.out))
}








