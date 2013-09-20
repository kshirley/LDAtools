#' @title Filter raw documents according to various options
#'
#' @description Conduct a series of preprocessing techniques. By default, a very limited amount of preprocessing will occur
#' (just basic stuff like remove documents that are blank or NA). The user must enter character vector arguments such as
#' \code{exact}, \code{partial}, \code{stopwords} and/or \code{subs} in order to conduct more advanced techniques.
#'
#' @param data a character vector containing the raw corpus. Each element should correspond to a 'document'.
#' @param exact a character vector in which each element is a string, phrase, or longer snippet of text
#' that you wish to discard (if the element(s) match the entire content of a document).
#' @param partial a character vector in which each element is a string, phrase, or longer snippet of text
#' that you wish to discard (if the element(s) match a subset of a document).
#' @param subs character vector where the odd-numbered element(s) are removed from the corpus 
#' and the subsequent even-numbered element are inserted in their place.
#' @param stopwords character vector of words that should be ignored.
#' @param cutoff The number of times a token should appear in order to be considered for prediction.
#' @param verbose logical. Track the categories of partial and exact matches. 
#' For instance, if a document exactly matches the third element of \code{exact}, 
#' then the corresponding value of \code{category} is 3.
#' @param quiet logical. Should a summary of the documents being trimmed be reported?
#' 
#' @return Returns a list of length four. The first element, \code{word.id}, is an integer vector with a unique value 
#' for each token in the corpus. The numbering corresponds to the ordering of third element \code{vocab}. The second
#' element, \code{doc.id}, is an integer vector which corresponds to the document each token belongs to. The fourth element,
#' \code{category} has length equal to the number of documents. If the value of an element in this vector is 0, 
#' then the corresponding document was retained. Otherwise, it was discarded. If the value is positive, it was an exact or partial match 
#' and if \code{verbose == TRUE} then the value points to the relevant element of \code{exact} or \code{partial}. If the value is -1, then
#' the document had no relevant words after removing \code{stopwords} and applying the \code{cutoff}.
#' 
#' @export
#' @examples
#' data(APcorpus) #Dollars and dates already removed
#' data(stopwords)
#' input <- filter(data=text, stopwords=stopwords)
#' input <- filter(data=text, stopwords=stopwords, verbose=TRUE)

preprocess <- function(data, exact=NULL, partial=NULL, subs=NULL, stopwords=NULL, cutoff=2, verbose=FALSE, quiet=FALSE) {
  D.orig <- length(data)
  
  # track and get rid of NA/blank documents:
  sel.na <- as.numeric(is.na(data))
  sel.blank <- as.numeric(data == "" & !is.na(data))
  if (!quiet) {
    count.na <- sum(sel.na)
    count.blank <- sum(sel.blank, na.rm=TRUE)
    cat(paste0(sprintf("%.1f", round(count.na/D.orig, 3)*100), "% (",
               count.na,"/", D.orig, ") of notes are NA"), "\n")
    cat(paste0(sprintf("%.1f", round(count.blank/D.orig, 3)*100), "% (",
               count.blank,"/", D.orig, ") of notes are blank"), "\n")
  }
  
  sel.match <- flag.exact(data, exact, D=D.orig, verbose, quiet)
  sel.grep <- flag.partial(data, partial, D=D.orig, verbose, quiet)
  
  # Discard the 'bad' notes and track the category each document belonged to:
  tmp.mat <- cbind(sel.match, sel.grep, sel.na, sel.blank)
  good <- apply(tmp.mat, 1, sum) == 0
  dat <- data[good]
  vec <- apply(tmp.mat, 1, which.max)
  vec[good] <- 0
  
  D <- length(dat)
  
  if (!quiet) {
    cat(paste0(sprintf("%.1f", 100 - round(D/D.orig, 3)*100),"% of notes removed; that is, ", D, " out of ", D.orig, " docs remaining"), "\n")  
  }
  
  # coerce documents to lowercase
  temp <- tolower(dat)
  # make substitutions
  if (!is.null(subs)) {
    n.subs <- length(subs)/2
    if (n.subs != floor(n.subs)) warning("The length of the subs object should be divisible by 2.")
    for (i in 1:n.subs) temp <- gsub(subs[i*2-1], subs[i*2], temp)
  }
  # handle the punctuation:
  temp <- gsub("\'", "", temp) # remove apostrophes
  temp <- gsub("[[:punct:]]+", " ", temp) # change any punctuation to single space
  temp <- gsub("[[:space:]]+", " ", temp) # change any multiple space to single space
  temp <- gsub("^+[[:space:]]", "", temp) # strip space at beginning of note
  temp <- gsub("[[:space:]]+$", "", temp) # strip space at end of note
  
  # get word vector:
  wl <- strsplit(temp," ", useBytes=TRUE)
  word.vec <- unlist(wl)
  
  # throw out 1-letter words (but keep single-digit words)
  word.vec <- word.vec[nchar(word.vec) > 1 | word.vec %in% 0:9]
  
  # Form a table of words
  word.tab <- table(word.vec)
  w.tab <- sort(word.tab, decreasing=TRUE)
  
  if (!quiet) cat(paste("Total Unique Tokens: ", length(w.tab), sep=""), "\n")
  if (!quiet) cat(paste("Total Tokens: ", sum(w.tab), sep=""), "\n")

  # remove the first set of stopwords, sw1, from w.tab:
  if (!is.null(stopwords)) {
    stops <- names(w.tab) %in% stopwords
    n.tokens.removed <- sum(w.tab[stops])
    w.tab <- w.tab[!stops]
    if (!quiet) {
      cat(paste0("Removed ", sum(stops), " stopwords from provided list of ", length(stopwords), " stopwords"), "\n")
      cat(paste0("Removed ", n.tokens.removed, " stopword occurrences from data"), "\n")
      cat(paste0("Total Remaining Unique Tokens: ", length(w.tab)), "\n")
      cat(paste0("Total Remaining Tokens: ", sum(w.tab)), "\n")
    }
  }
  
  # Use just words that occur more than cutoff times for prediction:
  pred.tab <- w.tab[w.tab > cutoff]
  vocab <- names(pred.tab)
  if (!quiet) {
    cat(paste0("Removing ", sum(w.tab <= cutoff), " tokens that appear ", cutoff, " times or less in the data"), "\n")
    cat(paste0("Total Remaining Unique Tokens: ", length(pred.tab)), "\n")
    cat(paste0("Total Remaining Tokens: ", sum(pred.tab)), "\n")
  }
  n <- sapply(wl, length)  # number of words per document (some may be 0)
  doc.id <- rep(1:D, n)
  words <- unlist(wl)
  idx <- words %in% vocab
  doc.id <- doc.id[idx]
  vec[vec==0][(1:D) %in% unique(doc.id) == FALSE] <- -1
  
  # every document must have at least one token in the vocabulary:
  doc.id <- match(doc.id, unique(doc.id))
  word.id <- match(words[idx], vocab)

  return(list(word.id=word.id, doc.id=doc.id, vocab=vocab, category=vec))
}





#' @title Flag the documents that exactly match a pre-specified list of strings
#'
#' @description If there are certain (typically very short) documents that occur frequently in your data, and 
#' you wish to remove them from the data before you fit the LDA model, this function can be used to flag those
#' documents. It's a trivial operation, but it's a useful reminder that users should visually inspect their
#' data before running LDA (so as to throw out documents that don't require topic modeling in the first place).
#' An example in customer care rep notes are the phrases "Not offered" (in reference to deals that were specifically
#' not offered to the customer during the phone call) and "Too expensive" (in reference to an offer that was declined
#' because if was too expensive).
#'
#' @param data a character vector containing the raw corpus. Each element should correspond to a 'document'.
#' @param exact a character vector in which each element is a string, phrase, or longer snippet of text
#' that you wish to discard (if the element(s) match the entire content of a document).
#' @param D the original number of documents
#' @param verbose logical vector. Should a summary of the documents being trimmed be reported?
#'
#' @return category an integer vector of the same length as \code{data}, where 0 indicates that the document did not
#' match any of the strings in \code{match.exact}, and an integer j = 1, ..., K that indicates that a document was an
#' exact match to the jth element of \code{match.exact}.
#'
#' @seealso flag.partial
#'
#' @export
#'
#' @examples
#' data <- c("bla bla bla", "foo", "bar", "text")
#' match.exact <- c("foo", "junk")
#' flag.exact(data, match.exact) # c(0, 2, 0, 0)
flag.exact <- function(data, exact, D, verbose=FALSE, quiet=FALSE) {
  l <- length(exact)
  if (l == 0) {
    return(numeric(D))
  } else {
    if (!quiet) cat("Looking for exact matches", "\n")
    if (verbose) {
      flagged <- matrix(0, D, l)
      for (i in 1:l) {
        flagged[,i] <- data == exact[i]
        if (!quiet) {
          rmd <- sum(flagged[,i])
          cat(paste0(sprintf("%.1f", round(rmd/D, 3)*100), "% (", rmd,"/", D, ") of notes are \'", exact[i], "\'"), "\n")
        }
      }
    } else {
      flagged <- as.numeric(data %in% exact)
    }  
  }
  return(flagged)
}

#' @title Flag the documents that contain an occurrence of one or more strings from a pre-specified list of strings
#'
#' @description Often a data set contains documents that you wish to remove before you fit the LDA model, and these
#' documents share a common "boilerplate" string or phrase (along with potentially unique information).
#' This function can be used to flag those
#' documents. Similar to the function \code{flag.exact}, this is a very simple operation that may be more useful as a
#' signal to the user that he or she should visually inspect the
#' data before running LDA (so as to remove documents that don't require topic modeling in the first place). An example
#' of a phrase that would indicate a document should be discarded is "iPhone Unlock Approved", which is a common phrase
#' followed by a unique transaction ID. Each document containing this string is unique (so we can't use \code{match.exact}
#' to flag these documents), but we don't want to include them in a topic model.
#'
#' @param data a character vector containing the raw corpus. Each element should correspond to a 'document'.
#' @param partial a character vector in which each element is a string, phrase, or longer snippet of text
#' that you wish to discard (if the element(s) match a subset of a document).
#' @param D the original number of documents
#' @param verbose logical vector. Track the specific categories of 'flagged' documents?
#' @return category an integer vector of the same length as \code{data}, where 0 indicates that the document did not
#' contain any of the strings in \code{partial}, and an integer j = 1, ..., K that indicates that a document contained
#' the jth element of \code{partial}.
#'
#' @seealso flag.exact
#'
#' @export
#'
#' @examples
#' data <- c("Automatic Message: Account 12 ...", "Automatic Message: Account 314 ...", "A document with unknown content", 
#' "Boilerplate text: Customer 1532 ...")
#' match.exact <- c("Automatic Message:", "Auto Text:", "Boilerplate text")
#' flag.partial(data, match.partial) # c(1, 1, 0, 3)

flag.partial <- function(data, partial, D, verbose, quiet=FALSE) {
  l <- length(partial)
  if (l == 0) {
    return(numeric(D))
  } else {
    if (!quiet) cat("Looking for partial matches", "\n")
    if (verbose) {
      flagged <- matrix(0, D, l)
      for (i in 1:l) {
        flagged[,i][grep(partial[i], data)] <- 1
        if (!quiet) {
          cat(paste0(sprintf("%.1f",round(sum(flagged[,i])/D, 3)*100),"% (",
                    sum(flagged[,i]),"/", D,") of notes are \'", partial[i], "\'"), "\n")
        }
      }
    } else {
      flagged <- as.numeric(grepl(paste(partial, collapse="|"), data))
    }  
  }
  return(flagged)
}


#' @title Filter raw version of new documents based on previously fit model
#'
#' @description This function performs the same preprocessing steps as the function \code{filter()}, except this time for a set of new
#' documents whose topic proportions we wish to estimate given the topics from a previously fit model. The key
#' difference is that the vocabulary won't be constructed based on the words that occur in the new documents, but rather it will be entered
#' as input from the previously fit model.
#'
#' @param data a character vector containing the raw corpus. Each element should correspond to a 'document'.
#' @param vocab a character vector containing the vocabulary of the fitted topic model from which we will estimate
#' the topic proportions for the new documents entered in \code{data}.
#' @param exact a character vector in which each element is a string, phrase, or longer snippet of text
#' that you wish to discard (if the element(s) match the entire content of a document).
#' @param partial a character vector in which each element is a string, phrase, or longer snippet of text
#' that you wish to discard (if the element(s) match a subset of a document).
#' @param subs character vector where the odd-numbered element(s) are removed from the corpus 
#' and the subsequent even-numbered element are inserted in their place.
#' @param verbose logical. Track the categories of partial and exact matches. 
#' For instance, if a document exactly matches the third element of \code{exact}, 
#' then the corresponding value of \code{category} is 3.
#' @param quiet logical. Should a summary of the documents being trimmed be reported?
#' 
#' @return Returns a list of length three. The first element, \code{word.id}, is an integer vector with a unique value 
#' for each token in the corpus. The numbering corresponds to the ordering of \code{vocab}. The second
#' element, \code{doc.id}, is an integer vector which corresponds to the document each token belongs to. The third element,
#' \code{category} has length equal to the number of documents. If the value of an element in this vector is 0, 
#' then the corresponding document was retained. Otherwise, it was discarded. If the value is positive, it was an exact or partial match 
#' and if \code{verbose == TRUE} then the value points to the relevant element of \code{exact} or \code{partial}. If the value is -1, then
#' the document had no relevant words after removing \code{stopwords} and applying the \code{cutoff}.
#' 
#' @export
#' @examples
#' data(APcorpus) #Dollars and dates already removed
#' data(stopwords)
#' input <- filter(data=text, stopwords=stopwords)
#' input <- filter(data=text, stopwords=stopwords, verbose=TRUE)

preprocess.newdocs <- function(data=character(), vocab=character(), exact=NULL, partial=NULL, subs=NULL, verbose=FALSE, quiet=FALSE) {
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
  
  sel.match <- flag.exact(data, exact, D=D.orig, verbose, quiet)
  sel.grep <- flag.partial(data, partial, D=D.orig, verbose, quiet)
  
  # Discard the 'bad' notes and track the category each document belonged to:
  tmp.mat <- cbind(sel.match, sel.grep, sel.na, sel.blank)
  good <- apply(tmp.mat, 1, sum)==0
  dat <- data[good]
  vec <- apply(tmp.mat, 1, which.max)
  vec[good] <- 0
  
  D <- length(dat)
  
  if (!quiet) {
    cat(paste0(sprintf("%.1f", 100 - round(D/D.orig, 3)*100),"% of notes removed; that is, ", D, " out of ", D.orig, " docs remaining"), "\n")  
  }
  
  # coerce documents to lowercase
  temp <- tolower(dat)
  # make substitutions
  if (!is.null(subs)) {
    n.subs <- length(subs)/2
    if (n.subs != floor(n.subs)) warning("The length of the subs object should be divisible by 2.")
    for (i in 1:n.subs) temp <- gsub(subs[i*2-1], subs[i*2], temp)
  }
  # handle the punctuation:
  temp <- gsub("\'", "", temp) # remove apostrophes
  temp <- gsub("[[:punct:]]+", " ", temp) # change any punctuation to single space
  temp <- gsub("[[:space:]]+", " ", temp) # change any multiple space to single space
  temp <- gsub("^+[[:space:]]", "", temp) # strip space at beginning of note
  temp <- gsub("[[:space:]]+$", "", temp) # strip space at end of note
  
  # get word vector:
  wl <- strsplit(temp," ", useBytes=TRUE)
  word.vec <- unlist(wl)

  n <- sapply(wl, length)  # number of words per document (some may be 0)
  doc.id <- rep(1:D, n)
  idx <- word.vec %in% vocab
  doc.id <- doc.id[idx]
  vec[vec==0][(1:D) %in% unique(doc.id) == FALSE] <- -1
  
  # every document must have at least one token in the vocabulary:
  doc.id <- match(doc.id, unique(doc.id))
  word.id <- match(word.vec[idx], vocab)

  return(list(word.id=word.id, doc.id=doc.id, category=vec))
}











