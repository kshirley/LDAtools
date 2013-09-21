library(shiny)
library(LDAviz)
library(proxy)
library(reshape)
library(plyr)
source("getProbsOld.R") #Kenny made some changes to getProbs -- use the old version -- getProbsOld

options(shiny.maxRequestSize=100*1024^2) #change default file upload size from 5MB to 100MB

KL <- function(x, y) { #compute Kullback-Leibler divergence
  .5*sum(x*log(x/y)) + .5*sum(y*log(y/x))
}

shinyServer(function(input, output) {
  
  getSample <- reactive( {
    data(APtopdocs, package="ldatools")
    data(APinput, package="ldatools")
    data(APtopics, package="ldatools")
    #inp <- get("input", env=globalenv())  #input is an unfortunate name...
    return(list(docs=APinput$doc.id, topics=APtopics$topics[,1], words=APinput$word.id, 
                vocab=APinput$vocab, topdocs=top.docs))
  })
  
  getLocal <- reactive( {
    if (!is.null(input$file1)) {
      path <- input$file1$datapath 
      dat <- read.table(path, header=TRUE, sep="\t")
      #dat$tokens should already be a factor with alphabetical levels (otherwise, this wont work)
      vocab <- levels(dat$tokens)
      word.id <- as.integer(dat$tokens) 
      dat <- list(docs=dat$docs, topics=dat$topics, words=word.id, vocab=vocab)
      if (!is.null(input$file2)) {
        path <- input$file2$datapath 
        topdocs <- read.table(path, as.is=TRUE, header=FALSE) # one column with n.docs * n.topics rows
        dat$topdocs <- topdocs
      }
      return(dat)
    } else NULL
  })
  
  getData <- reactive( {
    if (input$newDat) {
      data <- getLocal()
    } else {
      data <- getSample()
    }
  })
  
  filter <- reactive({
    data <- getData()
    topdocs <- data$topdocs
    docs <- data$docs
    topics <- data$topics
    words <- data$words
    vocab <- data$vocab
    tokens <- vocab[words]
    if (isTRUE(input$tool2 %in% vocab)){
      idx <- !tokens %in% input$tool2
      tokens <- tokens[idx]
      vocab <- vocab[!vocab %in% input$tool2]
      words <- match(tokens, vocab)
      topics <- topics[idx]
      docs <- docs[idx]
      #filter the actual documents too?
    } 
    return(list(docs=docs, topics=topics, words=words, vocab=vocab, topdocs=topdocs))
  })
  
  getMatrices <- reactive({ 
    stuff <- filter()
    getProbsOld(word.id=as.integer(stuff$words), doc.id=as.integer(stuff$docs), topic.id=as.integer(stuff$topics), 
            vocab=as.character(stuff$vocab), sort.topics="byDocs", sort.terms="saliency")
#     getProbsOld(word.id=as.integer(stuff$words), doc.id=as.integer(stuff$docs), topic.id=as.integer(stuff$topics), 
#              vocab=as.character(stuff$vocab), sort.topics="byDocs", K=30)
  })
  
  output$mdsDat <- reactive({
    stuff <- filter()
    dat <- getMatrices()
    k <- min(dim(dat$phi.hat)) 
    if (!is.null(stuff)){
      tokens <- as.character(stuff$vocab[stuff$words])
      docs <- as.integer(stuff$docs)
      topics <- as.integer(stuff$topics)
      if (is.null(stuff$topdocs)) {
        doc.df <- data.frame(Category=paste0("Topic", rep(1:k)), Document="No Documents Provided!")
      } else {
        #the ordering of topics are different!
        names(stuff$topdocs) <- colnames(dat$theta.hat)
        top.docz <- unlist(stuff$topdocs)
        #the row index is appended the end of the name - so remove it!
        names(top.docz) <- substr(names(top.docz), 1, nchar(names(top.docz))-1)
        doc.df <- data.frame(Category=names(top.docz), Document=as.character(top.docz), stringsAsFactors=FALSE)
      }
    }
  
    topic.order <- as.numeric(gsub("Topic", "", colnames(dat$theta.hat)))
    #Convert topic numbering to match the ordering determined by the sort.topics option of getProbs 
    topics <- match(topics, topic.order)
    #Overwrite topic assignments so that Topic1 is now 'most represented' and the last topic is 'least represented'.
    colnames(dat$theta.hat) <- paste0("Topic", 1:k)
    rownames(dat$phi.hat) <- paste0("Topic", 1:k)
    
    #calculate 'distance' between topics given the top 450 most frequent tokens
    tab <- table(tokens)
    pw <- tab/sum(tab)
    top <- names(tab)[order(tab, decreasing=FALSE) > 450] #add option to change this arbitrary cutoff?
    a <- apply(dat$phi.hat, 2, sum)
    phi.hat <- rbind(dat$phi.hat, avg=a/sum(a))
    phi <- phi.hat[, colnames(phi.hat) %in% top]
    d <- dist(phi, KL)
    fit <- cmdscale(d, k=2)
    x <- fit[,1]
    y <- fit[,2]
    lab <- gsub("Topic", "", names(x))
    loc.df <- data.frame(x, y, topics=lab, stringsAsFactors=FALSE)
    tab.topics <- table(topics)
    p.topics <- tab.topics/sum(tab.topics)
    topics.df <- data.frame(topics=as.numeric(names(p.topics)), Freq=as.numeric(p.topics*100))
    mds.df <- plyr::join(loc.df, topics.df, by="topics")
    #workaround errors if no clustering is done (ie, input$kmeans == 1)
    mds.df$cluster <- 1
    centers <- data.frame(x=0, y=0)
    if (input$kmeans > 1) { #and clustering info (if relevant)
      cl <- kmeans(cbind(x, y), input$kmeans)
      mds.df$cluster <- factor(cl$cluster)
      centers <- data.frame(cl$centers)
    }
    
    # map tokens to a cluster
    frame <- data.frame(tokens, topics, stringsAsFactors=FALSE)
    topic.cl.ID <- mds.df[c("topics", "cluster")]  #why does subsetting here make the topic -> cluster assignment go haywire???
    framed <- plyr::join(frame, topic.cl.ID, by="topics")
    
    ##############################################################################
    ### Create a df with the info neccessary to make the default OR new bar chart when selecting a topic or cluster.
    ### This functionality requires that we keep track of the top input$nTerms within each cluster and topic (as well as overall).
    ### The ranking of words under a topic is done via a weighted average of the lift and probability of the word given the topic.
    ### The ranking of words under a cluster is done via a similar weighted average (summing over the relevant topics)
    
    phi.t <- t(dat$phi.hat)
    pws <- as.numeric(pw[rownames(phi.t)]) #reorder frequencies to match the ordering of phi
    weight <- input$lambda*log(phi.t) + (1-input$lambda)*log(phi.t/pws) 
    
    #get the top terms and top documents for each topic
    top.terms <- NULL
    for (i in 1:k) {
      weights <- weight[,i]
      o <- order(weights, decreasing=TRUE)
      terms <- rownames(phi.t)[o][1:input$nTerms]
      top.terms <- c(top.terms, terms)
    }
    term.labs <- rep(paste0("Topic", 1:k), each=input$nTerms)
    topic.df <- data.frame(Term=top.terms, Category=term.labs, stringsAsFactors=FALSE)
    
    #get the top terms and top documents for each cluster
    clust.terms <- NULL
    if (input$kmeans == 1) {
      #if no clustering is done, we don't want to change the 'most informative words' upon hover
      clust.terms <- rownames(phi.t)[1:input$nTerms]
    } else {
      for (i in 1:input$kmeans) {
        #grab topics that belong to the current cluster
        topicz <- subset(topic.cl.ID, cluster == i & topics != "avg")$topics
        topicnames <- paste0("Topic", topicz)
        sub.phi <- phi.t[,topicnames]
        sub.theta <- dat$theta.hat[,topicnames]
        #only sum if multiple columns exist
        if (!is.null(dim(sub.phi))) {
          sub.phi <- apply(sub.phi, 1, sum)
          sub.theta <- apply(sub.theta, 1, sum)
        }
        weight <- input$lambda*log(sub.phi) + (1-input$lambda)*log(sub.phi/pws)
        o <- order(weight, decreasing=TRUE)
        terms <- rownames(phi.t)[o][1:input$nTerms]
        clust.terms <- c(clust.terms, terms)
      } 
    }
    term.labs <- rep(paste0("Cluster", 1:input$kmeans), each=input$nTerms)
    clust.df <- data.frame(Term=clust.terms, Category=term.labs, stringsAsFactors=FALSE)
    
    #use the ordering of terms determined by getProbs option sort.terms for the 'default' filter
    top.df <- data.frame(Term=rownames(phi.t)[1:input$nTerms], Category="Default")
    
    all.df <- rbind(topic.df, clust.df, top.df)
    all.df$Freq <- 0
    
    #now we have all the top ranking words within each cluster and topic
    #next, we find the frequency of those words in each category
    all.words <- unique(all.df$Term)
    all.frame <- subset(framed, tokens %in% all.words)
    counts <- table(as.character(all.frame$tokens), all.frame$topics)
    counts2 <- table(as.character(all.frame$tokens), all.frame$cluster)
    
    for (i in 1:k) {
      idx <- which(all.df$Category == paste0("Topic", i))
      all.df$Freq[idx] <- counts[all.df$Term[idx], i]
    }
    
    for (i in 1:input$kmeans) {
      idx <- which(all.df$Category == paste0("Cluster", i))
      all.df$Freq[idx] <- counts2[all.df$Term[idx], i]
    }
    
    totals <- table(as.character(all.frame$tokens))
    idx <- which(all.df$Category == "Default")
    all.df$Freq[idx] <- totals[all.df$Term[idx]]
    all.df$Total <- as.integer(totals[all.df$Term])
    
    #relative frequency (in percentages) over topics for each possible term
    probs <- t(apply(counts, 1, function(x) as.integer(100*x/sum(x))))
    topic.probs <- data.frame(probs, stringsAsFactors=FALSE)
    topic.probs$Term <- rownames(probs)
    topic.table <- melt(topic.probs, id.vars="Term", variable_name = "Topic")
    topic.table$Topic <- as.numeric(gsub("X", "", topic.table$Topic))
    names(topic.table) <- gsub("value", "Freq", names(topic.table))
    return(list(mdsDat=mds.df, mdsDat2=topic.table, barDat=all.df, docDat=doc.df,
                centers=centers, nClust=input$kmeans))
  })

output$dat <- renderPrint({
#treat me like your R console!
  stuff <- filter()
  dat <- getMatrices()
  k <- min(dim(dat$phi.hat)) 
  if (!is.null(stuff)){
    tokens <- as.character(stuff$vocab[stuff$words])
    docs <- as.integer(stuff$docs)
    topics <- as.integer(stuff$topics)
    if (is.null(stuff$topdocs)) {
      doc.df <- data.frame(Category=paste0("Topic", rep(1:k)), Document="No Documents Provided!")
    } else {
      #the ordering of topics are different!
      names(stuff$topdocs) <- colnames(dat$theta.hat)
      top.docz <- unlist(stuff$topdocs)
      #the row index is appended the end of the name - so remove it!
      names(top.docz) <- substr(names(top.docz), 1, nchar(names(top.docz))-1)
      doc.df <- data.frame(Category=names(top.docz), Document=as.character(top.docz), stringsAsFactors=FALSE)
    }
  }
  doc.df
})

})
