#' @title Preprocess raw documents according to various options
#'
#' @description Conduct a series of preprocessing steps on raw documents.
#' By default, a very limited amount of preprocessing will occur
#' (just basic stuff like remove documents that are blank or NA and then 
#' tokenize documents by separating by whitespace). The user can optionally
#' filter documents, perform global substitutions
#' using regular expressions, remove stopwords, and perform stemming.
#'
#' @param data a character vector containing the raw corpus, where each 
#' element is a document.
#'
#' @param exact a (case-sensitive) character vector in which 
#' each element is a string, phrase, or longer snippet of text
#' that results in a document being discarded from the data if the entire 
#' document matches an element of \code{exact}.
#'
#' @param partial a (case-sensitive) character vector in which 
#' each element is a string, phrase, or longer snippet of text
#' that results in a document being discarded from the data if any part 
#' of the document matches an element of \code{partial}.
#'
#' @param subs character vector of regular expressions where the odd-numbered 
#' element(s) are removed from the corpus 
#' and the subsequent even-numbered element are inserted in their place. 
#' These substitutions are performed
#' using the \code{gsub()} function after forcing the raw text to lowercase.
#'
#' @param stopwords character vector of tokens that should be excluded from 
#' the vocabulary.
#'
#' @param cutoff The minimum number of times a token must appear in the corpus 
#' in order to be included in the vocabulary.
#'
#' @param quiet logical. Should a summary of the preprocessing steps be 
#' printed to the screen?
#'
#' @param stem logical. Should the porter stemmer be used to stem the tokens 
#' in the vocabulary?
#'
#' @param hash a length-1 character vector indicating the prefix of substitution 
#' replacements that should be replaced with a '#' symbol after tokenizing. 
#' Set to "ent" by default, where "ent" stands for "entity",
#' and is often used as a prefix to a substitution replacement for a class 
#' of terms, like dollar amounts ("entdollaramount") and timestamps
#' ("entdatestamp", "enttimeofday"), etc. 
#' 
#' @return Returns a list of length five.
#' The first element, \code{term.id}, is an integer vector containing the
#' index in the vocabulary of each token in the corpus. If the 4th token in 
#' the corpus is "tree" and "tree" is the 50th element of the vocabulary, then
#' the 4th element of term.id will be 50.
#' The second element, \code{doc.id}, is an integer vector which corresponds to 
#' the document each token belongs to.
#' The third element, \code{vocab}, is the vocabulary of the corpus, which
#' contains all the terms (i.e. unique tokens) in the data. It is sorted in 
#' decreasing order of term frequency by default.
#' The fourth element, \code{category} has length equal to the number of input 
#' documents in \code{data}. If the value of an element in this vector is 0, 
#' then the corresponding document was retained. Otherwise, it was discarded. 
#' If the value is positive, it was an exact or partial match 
#' and if \code{verbose == TRUE} then the value points to the relevant element 
#' of \code{exact} or \code{partial}. If the value is -1, then the document
#' contained no tokens in the vocabulary after removing \code{stopwords} and 
#' applying the \code{cutoff}.
#' The fifth element, \code{call}, is a named list returning the 
#' arguments supplied to the \code{preprocess} function: exact, partial, 
#' subs, stopwords, stem, and cutoff.
#' 
#' @export
#'
#' @examples
#' data(APcorpus)
#' data(stopwords)
#' input <- preprocess(data=APcorpus, exact=NULL, partial=NULL, subs=NULL, 
#'                     stopwords=stopwords, cutoff=5,
#'                     quiet=FALSE, stem=FALSE, hash="ent")

preprocess <- function(data, exact=NULL, partial=NULL, subs=NULL, 
                       stopwords=NULL, cutoff=2, verbose=FALSE, quiet=FALSE,
                       stem=FALSE, hash="ent") {
  D.orig <- length(data)
  
  # Print part 1:
  if (!quiet) {
  	cat("Filtering documents:\n\n")
  }

  # track and get rid of NA/blank documents:
  sel.na <- as.numeric(is.na(data))
  sel.blank <- as.numeric(data == "" & !is.na(data))
  if (!quiet) {
    count.na <- sum(sel.na)
    count.blank <- sum(sel.blank, na.rm=TRUE)
    cat(paste0(sprintf("%.1f", round(count.na/D.orig, 3)*100), "% (",
               count.na,"/", D.orig, ") of documents are NA"), "\n")
    cat(paste0(sprintf("%.1f", round(count.blank/D.orig, 3)*100), "% (",
               count.blank,"/", D.orig, ") of documents are blank"), "\n")
  }
  
  # Discard documents that exactly match elements of 'exact' or contain 
  # strings that are elements of 'partial':
  sel.exact <- flag.exact(data, exact, verbose=FALSE, quiet)$category
  sel.partial <- flag.partial(data, partial, verbose=FALSE, quiet)$category
  if (!quiet) {
    cat(paste0(sprintf("%.1f", round(sum(sel.exact)/D.orig, 3)*100), "% (",
               sum(sel.exact), "/", D.orig, ") of documents discarded
               as exact matches"), "\n")
    cat(paste0(sprintf("%.1f", round(sum(sel.partial)/D.orig, 3)*100), "% (",
               sum(sel.partial), "/", D.orig, ") of documents discarded
               as partial matches"), "\n")
  }
  
  # Discard the filtered documents and track the category that each one 
  # belonged to:
  tmp.mat <- cbind(sel.exact, sel.partial, sel.na, sel.blank)
  keep.docs <- apply(tmp.mat, 1, sum) == 0
  dat <- data[keep.docs]
  # stable sorting of ties using max.col(..., ties.method='first')
  # for when a document matches multiple elements of 'partial':
  category <- max.col(tmp.mat, ties.method="first")
  category[keep.docs] <- 0
  
  # Re-set the total number of documents:
  D <- length(dat)

  # Print out progress to console:
  if (!quiet) {
    cat(paste0(sprintf("%.1f", 100 - round(D/D.orig, 3)*100), 
               "% of documents discarded; that is, ", D, " out of ", 
               D.orig, " documents remaining"), "\n")  
  }
  
  # coerce documents to lowercase:
  if (!quiet) {
  	cat("\nForcing documents to lowercase\n\n")
  }
  temp <- tolower(dat)

  # make substitutions:
  if (!is.null(subs)) {
    if (!quiet) {
  	  cat("Performing regular expression substitutions\n\n")
    }
    n.subs <- length(subs)/2
    if (n.subs != floor(n.subs)){ 
      warning("The length of the subs object should be divisible by 2.")
    }
    for (i in 1:n.subs) {
      temp <- gsub(subs[i*2-1], subs[i*2], temp)
    }
  }

  # handle the punctuation:
  if (!quiet) {
  	cat("Tokenizing documents by separating on whitespace and punctuation\n\n")
  }

  # remove apostrophes
  temp <- gsub("\'", "", temp)
  # change any punctuation to single space
  temp <- gsub("[[:punct:]]+", " ", temp)
  # change any multiple space to single space
  temp <- gsub("[[:space:]]+", " ", temp)
  # strip space at beginning of document
  temp <- gsub("^+[[:space:]]", "", temp)
  # strip space at end of document
  temp <- gsub("[[:space:]]+$", "", temp)
  
  # tokenize documents:
  wl <- strsplit(temp," ", useBytes=TRUE)
  token.vec <- unlist(wl)
  
  # throw out 1-letter terms (but keep single-digit terms)
  if (!quiet) cat("Automatically discard 1-letter terms\n\n")
  token.vec <- token.vec[nchar(token.vec) > 1 | token.vec %in% 0:9]
  
  # Form a table of terms:
  term.table <- table(token.vec)
  term.table <- sort(term.table, decreasing=TRUE)
  
  if (!quiet) cat(paste("Total # of Terms: ", length(term.table), sep=""), "\n")
  if (!quiet) cat(paste("Total # of Tokens: ", sum(term.table), sep=""), "\n")

  # remove stopwords:
  if (!is.null(stopwords)) {
    if (!quiet) {
      cat("Removing stop words from vocabulary\n")
    }
    stops <- names(term.table) %in% stopwords
    n.tokens.removed <- sum(term.table[stops])
    term.table <- term.table[!stops]
    if (!quiet) {
      cat(paste0("Removed ", sum(stops), " stop words from provided list of ", 
                 length(stopwords), " stop words"), "\n")
      cat(paste0("Removed ", n.tokens.removed, 
                 " stop word occurrences from data"), "\n")
      cat(paste0("Total Remaining Terms: ", length(term.table)), "\n")
      cat(paste0("Total Remaining Tokens: ", sum(term.table)), "\n")
    }
  }
  
  # Remove terms that occur fewer than cutoff times from the vocabulary:
  pred.tab <- term.table[term.table >= cutoff]
  vocab <- names(pred.tab)
  if (!quiet) {
    cat(paste0("Removing ", sum(term.table < cutoff), 
               " terms that appear less than ", cutoff, 
               " times in the data"), "\n")
    cat(paste0("Total Remaining Terms: ", length(pred.tab)), "\n")
    cat(paste0("Total Remaining Tokens: ", sum(pred.tab)), "\n")
  }
  # now the vocabulary has been defined

  # compute number of tokens per document (some may be 0)
  n <- sapply(wl, length)
  
  # create document ID vector:
  doc.id <- rep(1:D, n)

  # Now delete tokens that are not in the vocabulary:
  keep <- token.vec %in% vocab
  doc.id <- doc.id[keep]

  # set 'category' value to -1 if all tokens from a document were removed 
  # because of various preprocessing steps
  # i.e. substitutions, punctuation, stopwords, or the cutoff setting 
  # caused the document to contain zero tokens in the vocabulary:
  category[category == 0][(1:D) %in% unique(doc.id) == FALSE] <- -1
  if (!quiet) {
  	cat(paste0(sum(category == -1), "additional documents removed because
  	           they consisted entirely of punctuation or rare terms
  	           that are not in the vocabulary."))
  }
  
  # every document must have at least one token in the vocabulary:
  doc.id <- match(doc.id, unique(doc.id))
  term.id <- match(token.vec[keep], vocab)
  
  # Insert hash symbol '#' for tokens that originated as replacements 
  # in the substitution table:
  if (!is.null(subs)) {
    hash.subs <- subs[seq(2, by=2, length=length(subs))]
    hash.subs <- gsub("^[[:space:]]+", "", hash.subs)
    hash.subs <- gsub("[[:space:]]+$", "", hash.subs)
    vocab[vocab %in% hash.subs] <- gsub("^ent", "#", 
                                        vocab[vocab %in% hash.subs])
  }
  
  # Do stemming:
  if (stem) {
    V.old <- length(vocab)
  	vocab.stemmed <- vocab

    # Only perform stemming on non-substitution tokens:
    g <- grep("^#", vocab)
    if (length(g) > 0) {
  	  vocab.stemmed[-g] <- wordStem(vocab[-g])
  	} else {
  	  vocab.stemmed <- wordStem(vocab)
  	}

    # Re-map token ids to new vocabulary:
    vocab <- unique(vocab.stemmed)
    term.id <- match(vocab.stemmed[term.id], vocab)
    V.new <- length(vocab)
    if (!quiet) cat(paste0("Vocab of length ", V.old, 
                    " reduced to length ", V.new, 
                    " after stemming."), "\n")
  }
  if (!quiet) {
    cat(paste0("\n"))
    cat(paste0("Final Summary:\n"))
    cat(paste0("Total number of tokens in data, N = ", 
               length(term.id), ".\n"))
    cat(paste0("Total number of documents remaining, D = ", 
               max(doc.id), ".\n"))
    cat(paste0("Total number of documents discarded is ", 
               D.orig - max(doc.id), ".\n"))
    cat(paste0("Total number of terms in vocabulary, W = ", 
               length(vocab), ".\n"))
  }
  return(list(term.id=term.id, doc.id=doc.id, vocab=vocab, category=category, 
              call=list(exact=exact, partial=partial, subs=subs, 
                        stopwords=stopwords, stem=stem, cutoff=cutoff)))
}





#' @title Flag the documents that exactly match a pre-specified list of strings
#'
#' @description If there are certain (typically very short) documents that 
#' occur frequently in your data, and 
#' you wish to remove them from the data before you fit the LDA model, this
#' function can be used to flag those
#' documents. It's a trivial operation, but it's a useful reminder that 
#' users should visually inspect their
#' data before running LDA (so as to throw out documents that don't require 
#' topic modeling in the first place).
#'
#' @param data a character vector containing the raw corpus. Each element 
#' should correspond to a 'document'.
#' 
#' @param exact a character vector in which each element is a string, phrase, 
#' or longer snippet of text that you wish to discard, if the element matches
#' the entire content of a document.
#'
#' @param verbose logical. Track the categories of exact matches. For instance,
#' if a document exactly matches the third element of \code{exact}, 
#' then the corresponding value returned will be 3.
#' 
#' @param quiet logical. Should a summary of the preprocessing steps be printed
#' to the screen?
#'
#' @return category an integer vector of the same length as \code{data}, where,
#' if verbose=TRUE, 0 indicates that the document did not
#' match any of the strings in \code{exact}, and an integer j = 1, ..., K 
#' indicates that a document was an exact match to the jth element of 
#' \code{exact}, and if verbose=FALSE, an indicator vector
#' of whether the document exactly matched any of the elements of 
#' \code{exact} (without indicating which element
#' it matched).
#'
#' @seealso flag.partial
#'
#' @export
#'
#' @examples
#' data <- c("bla bla bla", "foo", "bar", "text")
#' match.exact <- c("foo", "junk")
#' flag.exact(data, match.exact, verbose=FALSE, quiet=FALSE) # c(0, 1, 0, 0)
#' flag.exact(data, match.exact, verbose=TRUE, quiet=FALSE) # c(0, 2, 0, 0)

flag.exact <- function(data, exact, verbose=FALSE, quiet=FALSE) {
  D <- length(data)
  l <- length(exact)
  if (l == 0) {
    return(list(category=numeric(D)))
  } else {
    if (!quiet) cat("Looking for exact matches", "\n")
    if (verbose) {
      flagged <- matrix(0, D, l)
      for (i in 1:l) {
        flagged[, i] <- data == exact[i]
        if (!quiet) {
          rmd <- sum(flagged[,i])
          cat(paste0(sprintf("%.1f", round(rmd/D, 3)*100), "% (", 
                     rmd,"/", D, ") of notes are \'", exact[i], "\'"), 
                     "\n")
        }
      }
    } else {
      flagged <- as.numeric(data %in% exact)
    }  
  }
  return(list(category=flagged))
}

#' @title Flag the documents that contain an occurrence of one or more strings
#'  from a pre-specified list of strings
#'
#' @description Often a data set contains documents that you wish to remove 
#' before you fit the LDA model, and these
#' documents share a common "boilerplate" string or phrase (along with 
#' potentially unique information).
#' This function can be used to flag those
#' documents. Similar to the function \code{flag.exact}, this is a very simple 
#' operation that may be more useful as a
#' signal to the user that he or she should visually inspect the
#' data before running LDA (so as to remove documents that don't require topic 
#' modeling in the first place).
#'
#' @param data a character vector containing the raw corpus. Each element 
#' should correspond to a 'document'.
#'
#' @param partial a character vector in which each element is a string, phrase,
#' or longer snippet of text
#' that you wish to discard, if the element matches a subset of a document.
#' 
#' @param verbose logical. Track the categories of partial matches. 
#' For instance, if a document partially matches the third element of 
#' \code{partial}, then the corresponding value returned will be 3.
#'
#' @param quiet logical. Should a summary of the preprocessing steps be 
#' printed to the screen?
#'
#' @return category an integer vector of the same length as \code{data}, 
#' where, if verbose=TRUE, 0 indicates that the document did not
#' match any of the strings in \code{partial}, and an integer j = 1, ..., K 
#' indicates that a document was a partial match to the jth element of 
#' \code{partial}, and if verbose=FALSE, an indicator vector
#' of whether the document partially matched any of the elements of 
#' \code{partial} (without indicating which element it matched).
#'
#' @seealso flag.exact
#'
#' @export
#'
#' @examples
#' data <- c("Automatic Message: Account 12 ...", 
#' "Automatic Message: Account 314 ...", 
#' "A document with unknown content", 
#' "Boilerplate text: Customer 1532 ...")
#' match.exact <- c("Automatic Message:", "Auto Text:", "Boilerplate text")
#' flag.partial(data, match.partial, verbose=FALSE, quiet=FALSE) # c(1, 1, 0, 1)
#' flag.partial(data, match.partial, verbose=TRUE, quiet=FALSE) # c(1, 1, 0, 3)

flag.partial <- function(data, partial, verbose, quiet=FALSE) {
  D <- length(data)
  l <- length(partial)
  if (l == 0) {
    return(list(category=numeric(D)))
  } else {
    if (!quiet) cat("Looking for partial matches", "\n")
    if (verbose) {
      flagged <- matrix(0, D, l)
      for (i in 1:l) {
        flagged[, i][grep(partial[i], data)] <- 1
        if (!quiet) {
          cat(paste0(sprintf("%.1f", round(sum(flagged[,i])/D, 3)*100), 
                     "% (", sum(flagged[, i]), "/", D, ") of notes are \'", 
                     partial[i], "\'"), "\n")
        }
      }
    } else {
      flagged <- as.numeric(grepl(paste(partial, collapse="|"), data))
    }  
  }
  return(list(category=flagged))
}


#' @title Preprocess raw version of new documents based on previously fit model
#'
#' @description This function performs the same preprocessing steps as the 
#' function \code{preprocess()}, except this time for a set of new
#' documents whose topic proportions we wish to estimate given the topics 
#' from a previously fit model. The key difference is that the vocabulary 
#' won't be constructed based on the words that occur in the new documents, but 
#' rather it will be entered as input from the previously fit model.
#'
#' @param data a character vector containing the raw corpus, where each 
#' element is a document.
#'
#' @param vocab a character vector containing the vocabulary of the fitted 
#' topic model from which we will estimate the topic proportions for the new 
#' documents entered in \code{data}.
#'
#' @param exact a (case-sensitive) character vector in which 
#' each element is a string, phrase, or longer snippet of text
#' that results in a document being discarded from the data if the entire 
#' document matches an element of \code{exact}.
#'
#' @param partial a (case-sensitive) character vector in which 
#' each element is a string, phrase, or longer snippet of text
#' that results in a document being discarded from the data if any part 
#' of the document matches an element of \code{partial}.
#'
#' @param subs character vector of regular expressions where the odd-numbered 
#' element(s) are removed from the corpus 
#' and the subsequent even-numbered element are inserted in their place. 
#' These substitutions are performed
#' using the \code{gsub()} function after forcing the raw text to lowercase.
#'
#' @param verbose logical. If set to TRUE the function will retain the indices 
#' of the elements of \code{exact} and \code{partial} that were matched. 
#' For instance, if a document exactly matches the third element of 
#' \code{exact}, then the corresponding value of \code{category} will be 3, 
#' if \code{verbose = TRUE}
#'
#' @param quiet logical. Should a summary of the preprocessing steps be 
#' printed to the screen?
#'
#' @param stem logical. Should the porter stemmer be used to stem the tokens 
#' in the vocabulary?
#'
#' @param hash a length-1 character vector indicating the prefix of substitution 
#' replacements that should be replaced with a '#' symbol after tokenizing. 
#' Set to "ent" by default, where "ent" stands for "entity",
#' and is often used as a prefix to a substitution replacement for a class 
#' of terms, like dollar amounts ("entdollaramount") and timestamps
#' ("entdatestamp", "enttimeofday"), etc. 
#' 
#' @return Returns a list of length three.
#' The first element, \code{term.id}, is an integer vector containing the
#' index in the vocabulary of each token in the corpus. If the 4th token in 
#' the corpus is "tree" and "tree" is the 50th element of the vocabulary, then
#' the 4th element of term.id will be 50, for example.
#' The second element, \code{doc.id}, is an integer vector which corresponds 
#' to the document each token belongs to. 
#' The third element, \code{category} has length equal to the number of 
#' documents. If the value of an element in this vector is 0, 
#' then the corresponding document was retained. Otherwise, it was discarded. 
#' If the value is positive, it was an exact or partial match 
#' and if \code{verbose == TRUE} then the value points to the relevant element 
#' of \code{exact} or \code{partial}. If the value is -1, then the document 
#' contained no tokens in the vocabulary.
#' 
#' @export

preprocess.newdocs <- function(data=character(), vocab=character(), exact=NULL,
                               partial=NULL, subs=NULL, verbose=FALSE, 
                               quiet=FALSE) {
  D.orig <- length(data)
  
  # track and get rid of NA/blank documents:
  sel.na <- as.numeric(is.na(data))
  sel.blank <- as.numeric(data=="")
  if (!quiet) {
    count.na <- sum(sel.na)
    count.blank <- sum(sel.blank)
    cat(paste0(sprintf("%.1f", round(count.na/D.orig, 3)*100), "% (",
               count.na,"/", D.orig, ") of notes are NA"), "\n")
    cat(paste0(sprintf("%.1f", round(count.blank/D.orig, 3)*100), "% (",
               count.blank,"/", D.orig, ") of notes are blank"), "\n")
  }
  
  sel.exact <- flag.exact(data, exact, verbose, quiet)$category
  sel.partial <- flag.partial(data, partial, verbose, quiet)$category
  
  # Discard the 'bad' notes and track the category each document belonged to:
  tmp.mat <- cbind(sel.exact, sel.partial, sel.na, sel.blank)
  keep.docs <- apply(tmp.mat, 1, sum) == 0
  dat <- data[keep.docs]

  # stable sorting of ties using max.col(..., ties.method='first')
  category <- max.col(tmp.mat, ties.method="first")
  category[keep.docs] <- 0
  
  D <- length(dat)
  
  if (!quiet) {
    cat(paste0(sprintf("%.1f", 100 - round(D/D.orig, 3)*100), 
               "% of notes removed; that is, ", D, 
               " out of ", D.orig, " docs remaining"), "\n")
  }
  
  # coerce documents to lowercase
  temp <- tolower(dat)

  # make substitutions
  if (!is.null(subs)) {
    n.subs <- length(subs)/2
    if (n.subs != floor(n.subs)) {
      warning("The length of the subs object should be divisible by 2.")
    }
    for (i in 1:n.subs) temp <- gsub(subs[i*2-1], subs[i*2], temp)
  }
  # handle the punctuation:
  # remove apostrophes
  temp <- gsub("\'", "", temp)
  # change any punctuation to single space
  temp <- gsub("[[:punct:]]+", " ", temp)
  # change any multiple space to single space
  temp <- gsub("[[:space:]]+", " ", temp)
  # strip space at beginning of note
  temp <- gsub("^+[[:space:]]", "", temp)
  # strip space at end of note
  temp <- gsub("[[:space:]]+$", "", temp)
  
  # get word vector:
  wl <- strsplit(temp," ", useBytes=TRUE)
  word.vec <- unlist(wl)

  # number of words per document (some may be 0)
  n <- sapply(wl, length)
  doc.id <- rep(1:D, n)
  keep <- word.vec %in% vocab
  doc.id <- doc.id[keep]
  category[category == 0][(1:D) %in% unique(doc.id) == FALSE] <- -1
  
  # every document must have at least one token in the vocabulary:
  doc.id <- match(doc.id, unique(doc.id))
  term.id <- match(word.vec[keep], vocab)

  return(list(term.id=term.id, doc.id=doc.id, category=category))
}




#' @title Compute the lower bound of the perplexity of a topic model fit 
#' using LDA
#'
#' @description This is an 'experimental' function that computes the lower 
#' bound of the perplexity of the training data in an LDA topic model. We 
#' claim that the perplexity of the training data is minimized when alpha and 
#' beta, the priors for the document-topic distributions and the topic-term
#' distributions, respectively, approach zero, and the number of topics is 
#' equal to min(D, W), where D is the number of training documents and W is 
#' the size of the vocabulary. I'll have to write up the idea some day.
#'
#' @param term.id an integer vector containing the term ID number of every 
#' token in the corpus. Should take values between 1 and W, where W is the 
#' number of terms in the vocabulary.
#' 
#' @param alpha the dirichlet prior parameter for the document-topic
#' multionomial distributions
#'
#' @param beta the dirichlet prior parameter for the topic-term multionomial 
#' distributions
#'
#' @param term.frequency an integer vector containing the counts of the 
#' number of occurences of each term in the vocabulary
#'
#' @param doc.frequency an integer vector containing the number of tokens per 
#' document, whose length is equal to the total number of documents in the 
#' corpus.
#'
#' @return bounds a numeric vector of length two containing the upper and 
#' lower bound of perplexity. The lower bound is computed by the C function. 
#' The upper bound is simple -- it is just the perplexity of a 1-topic model.
#' Included as output just for convenience.
#'
#' @export

perplexity.bounds <- function(term.id=integer(), alpha=double(), beta=double(),
                              term.frequency=integer(), doc.frequency=integer(),
                              print=50) {
                   	
  # Compute a few global parameters from the inputs:
  N <- as.integer(length(term.id))
  W <- as.integer(max(term.id))
  D <- as.integer(length(doc.frequency))

  # Set the placeholder for the log-probability vector that the C function 
  # will compute:
  logp <- as.double(numeric(N))

  # Call the C function:
  f <- .C("minperplex", token_id=as.integer(term.id), N_R=N, W_R=W, 
          D_R=D, alpha_R=as.double(alpha), 
          beta_R=as.double(beta), logp=logp, 
          tokentable=as.integer(term.frequency), 
          ndocs=as.integer(doc.frequency), print_R=as.integer(print), 
          package='ldatools')
  lower.bound <- exp(-sum(f$logp)/N)
  upper.bound <- exp(-sum(log(normalize(term.frequency + beta)[term.id]))/N)
  return(c(upper.bound, lower.bound))
}



#' @title Compute table of bigrams
#'
#' @description This function counts the bigrams in the data. It's based on 
#' the vector of term IDs and document IDs -- that is, the vocabulary has 
#' already been established, and this function simply counts occurrences of
#' consecutive terms in the data.
#'
#' @param term.id an integer vector containing the term ID number of every
#' token in the corpus. Should 
#' take values between 1 and W, where W is the number of terms in the 
#' vocabulary.
#'
#' @param doc.id an interger vector containing the document ID number of 
#' every token in the corpus. Should
#' take values between 1 and D, where D is the total number of documents in 
#' the corpus.
#'
#' @param vocab a character vector of length W, containing the terms 
#' in the vocabulary. This vector must align with \code{term.id}, such that a 
#' term.id of 1 indicates the first element of \code{vocab}, a term.id
#' of 2 indicates the second element of \code{vocab}, etc.
#'
#' @param n an integer specifying how large the bigram table should be. The 
#' function will return the top n most frequent bigrams. This argument is here 
#' because the number of bigrams can be as large as W^2.
#'
#' @return a dataframe with three columns and \code{n} rows, containing the 
#' bigrams (column 2), their frequencies (column 3), and their rank in 
#' decreasing order of frequency (column 1). The table is sorted
#' by default in decreasing order of frequency.
#'
#' @export

# Function to compute table of bigrams:
bigram.table <- function(term.id=integer(), doc.id=integer(), 
                         vocab=character(), n=integer()) {
  stopifnot(n > 0)
  N <- length(term.id)
  stopifnot(N > 0)
  bg.vec <- paste(vocab[term.id][1:(N - 1)], vocab[term.id][2:N], sep="-")
  n.doc <- as.numeric(table(doc.id))
  del <- cumsum(n.doc)
  tb.bigram <- data.frame(table(bg.vec[-del]))
  tb.bigram <- tb.bigram[order(tb.bigram[, 2], decreasing=TRUE), ][1:min(n, 
               dim(tb.bigram)[1]), ]
  colnames(tb.bigram) <- c("Bigram", "Frequency")
  rownames(tb.bigram) <- 1:dim(tb.bigram)[1]
  tb.bigram <- data.frame(Rank=1:dim(tb.bigram)[1], 
                          Bigram=tb.bigram[, "Bigram"], 
                          Frequency=tb.bigram[, "Frequency"])
  return(tb.bigram)
}




#' @title This function replaces instances of specified terms with existing 
#' or new terms in a corpus of documents
#'
#' @description After tokenization, use this function to replace all 
#' occurrences of a given term with a new term, whether the new term is 
#' in the existing vocabulary or not.
#'
#' @param term.map data.frame with two columns, the first of which contains 
#' the terms to be replaced, and the second of which contains their 
#' replacements. If a replacement is not in the current vocabulary, it will
#' be added to the vocabulary.
#'
#' @param term.id an integer vector containing the term ID number of every 
#' token in the corpus. Should take values between 1 and W, where W is the
#' number of terms in the vocabulary.
#'
#' @param vocab a character vector of length W, containing the terms in the
#' vocabulary. This vector must align with \code{term.id}, such that a 
#' term.id of 1 indicates the first element of \code{vocab}, a term.id
#' of 2 indicates the second element of \code{vocab}, etc.
#'
#' @return Returns a list of length two. The first element, \code{new.vocab}, 
#' is a character vector containing the new vocabulary. The second element, 
#' \code{new.term.id} is the new vector of term ID numbers
#' for all tokens in the data, taking integer values from 1 to the length 
#' of the new vocabulary.
#' 
#' @export

remap.terms <- function(term.map=data.frame(), term.id=integer(), 
                        vocab=character()) {

  # Set up some warnings:
  if(class(term.map) != "data.frame") {
  	stop("Error: term.map must be a data.frame")
  }
  if(dim(term.map)[2] != 2) {
  	stop("Error: term.map must have two columns, one for the terms to 
  	     be replaced, and one for the replacement terms")
  }
  
  # create a vector of all tokens in data:
  tokens <- vocab[term.id]
  
  # Make the substitutions:
  J <- dim(term.map)[1]
  for (i in 1:J){
    if (term.map[i, 1] %in% vocab) {
      print(paste0("Replacing ", term.map[i, 1], " with ", term.map[i, 2]))
      tokens[tokens == term.map[i, 1]] <- term.map[i, 2]
    } else {
      print(paste0(term.map[i, 1], " not in vocabulary; no replacement 
                   to perform."))
    }
  }
  token.tab <- table(tokens)
  token.tab <- sort(token.tab, decreasing=TRUE)
  new.vocab <- names(token.tab)
  new.term.id <- match(tokens, new.vocab)

  cat(paste0(length(vocab), " terms in old vocabulary"), "\n")
  cat(paste0(length(new.vocab), " terms in new vocabulary"), "\n")
  return(list(new.vocab=new.vocab, new.term.id=new.term.id))
}



#' @title Replace specified bigrams with terms representing the bigrams
#'
#' @description After tokenization, use this function to replace all 
#' occurrences of a given bigram with a single token representing the bigram,
#' and 'delete' the occurrences of the two individual tokens that comprised 
#' the bigram (so that it is still a generative model for text).
#'
#' @param bigrams A character vector, each element of which is a bigram 
#' represented by two terms separated by a hyphen, such as 'term1-term2'.
#' Every consecutive occurrence of 'term1' and 'term2' in the data
#' will be replaced by a single token representing this bigram.
#'
#' @param doc.id an interger vector containing the document ID number of 
#' every token in the corpus. Should
#' take values between 1 and D, where D is the total number of documents in 
#' the corpus.
#'
#' @param term.id an integer vector containing the term ID number of every
#' token in the corpus. Should 
#' take values between 1 and W, where W is the number of terms in the 
#' vocabulary.
#'
#' @param vocab a character vector of length W, containing the terms 
#' in the vocabulary. This vector must align with \code{term.id}, such that a 
#' term.id of 1 indicates the first element of \code{vocab}, a term.id
#' of 2 indicates the second element of \code{vocab}, etc.
#'
#' @return Returns a list of length three. The first element, 
#' \code{new.vocab}, is a character vector
#' containing the new vocabulary. The second element, \code{new.term.id} 
#' is the new vector of term ID numbers
#' for all tokens in the data, taking integer values from 1 to the length 
#' of the new vocabulary. The third
#' element is \code{new.doc.id}, which is the new version of the document 
#' id vector. If any of the specified bigrams
#' were present in the data, then \code{new.term.id} and \code{new.doc.id} 
#' will be shorter vectors than the original \code{term.id} and 
#' \code{doc.id} vectors.
#' 
#' @export

collapse.bigrams <- function(bigrams=character(), doc.id=integer(), 
                             term.id=integer(), vocab=character()) {

  # Set a few parameters:
  D <- max(doc.id)
  N <- length(doc.id)
  W <- length(vocab)
  W.orig <- W

  # Just grab unique values of bigrams:
  bigrams <- sort(unique(bigrams))

  # Set up some warnings:
  g <- grep("-", bigrams)
  if (length(g) != length(bigrams)) {
  	stop("Error: each bigram in 'bigrams' must contain a hyphen")
  }
  bg <- strsplit(bigrams, "-", fixed=TRUE)
  if (sum(unlist(lapply(bg, length)) == 2) != length(bigrams)) {
    cat(paste0("  Warning: the following bigrams did not contain two 
               non-empty character strings 
               separated by a single hyphen, and were deleted:\n\n  ", 
               bigrams[unlist(lapply(bg, length)) != 2]))
  	bigrams <- bigrams[unlist(lapply(bg, length)) == 2]
  	bg <- bg[unlist(lapply(bg, length)) == 2]
  }

  # Get the first and second term in each bigram:
  term1 <- unlist(bg)[seq(1, by=2, length=length(bg))]
  term2 <- unlist(bg)[seq(2, by=2, length=length(bg))]

  # concatenate these two vectors:
  all.terms <- c(term1, term2)

  # Send out a warning message if some of the tokens that comprise the 
  # bigrams don't occur in the vocabulary:
  if (sum(all.terms %in% vocab) != length(all.terms)) {
    w <- which(all.terms %in% vocab == FALSE)
    cat(paste0("  Note: the following term(s) were not found in the 
               vocabulary:\n\n  ", paste(token1[w], collapse=", "), 
               "\n"))
  }

  # get indices of vocab for each token:
  m1 <- match(term1, vocab)
  m2 <- match(term2, vocab)

  # Start replacing bigrams:
  nb <- length(bigrams)
  for (i in 1:nb) {
    if (!is.na(m1[i]) & !is.na(m2[i])) {
      # only bother is both terms are in the vocab

      # Indices of last token for each document:
      num.docs <- cumsum(table(doc.id))
      N <- length(term.id)
      W <- length(vocab)    
    
      # Find occurrences of bigrams:
      sel.bigram <- term.id[1:(N - 1)] == m1[i] & term.id[2:N] == m2[i] & 
                    1:(N - 1) %in% num.docs == FALSE
      if (sum(sel.bigram) > 0) {  # this bigram occurred at least once
        term.id[sel.bigram] <- W + 1
        term.id <- term.id[-(which(sel.bigram) + 1)]
        doc.id <- doc.id[-(which(sel.bigram) + 1)]
        vocab <- c(vocab, bigrams[i])
      }
    }
  }

  # Almost done, but there is the possibility of zero-occurrence terms.
  # If all the occurrences of a term were within a bigram, then they were 
  # all converted, and the original count
  # for that term will now be zero. Let's delete such terms from the vocabulary 
  # (and re-sort according to frequency):
  
  # create a vector of all tokens in data:
  tokens <- vocab[term.id]
  term.tab <- table(tokens)
  term.tab <- sort(term.tab, decreasing=TRUE)
  new.vocab <- names(term.tab)
  new.term.id <- match(tokens, new.vocab)

  cat(paste0("  ", W.orig, " terms in old vocabulary"), "\n")
  cat(paste0("  ", length(new.vocab), " terms in new vocabulary"), "\n")
  return(list(new.vocab=new.vocab, new.term.id=new.term.id, 
              new.doc.id=doc.id))
}





