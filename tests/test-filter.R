library(testthat)
library(LDAtools)

input1 <- filter(data=text, exact=exact, partial=partial, subs=subs, stopwords=stopwords, verbose=FALSE, quiet=TRUE)
input2 <- filter(data=text, exact=exact, partial=partial, subs=subs, stopwords=stopwords, verbose=FALSE, quiet=FALSE)
input3 <- filter(data=text, exact=exact, partial=partial, subs=subs, stopwords=stopwords, verbose=TRUE, quiet=TRUE)
input4 <- filter(data=text, exact=exact, partial=partial, subs=subs, stopwords=stopwords, verbose=TRUE, quiet=FALSE)

expect_true(identical(input1, input2))
#categories are fundamentally different between input2 and input3
expect_true(identical(input2$word.id, input3$word.id))
expect_true(identical(input2$doc.id, input3$doc.id))
expect_true(identical(input2$vocab, input3$vocab))
expect_true(identical(input3, input4))

#sanity check that sort.topics="byDocs" works
dat <- getProbs(word.id=input$word.id, doc.id=input$doc.id, topic.id=topics, 
                vocab=input$vocab, sort.topics="byDocs")
dat2 <- getProbs(word.id=input$word.id, doc.id=input$doc.id, topic.id=topics, 
                 vocab=input$vocab)
top.topic <- max.col(dat2$theta.hat)
tabs <- sort(table(top.topic), decreasing=TRUE)
expect_true(all(names(tabs) == gsub("Topic", "", colnames(dat$theta.hat))))

#sanity check that sort.topics="byTerms" works
dat <- getProbs(word.id=input$word.id, doc.id=input$doc.id, topic.id=topics, 
                vocab=input$vocab, sort.topics="byTerms")
dat2 <- getProbs(word.id=input$word.id, doc.id=input$doc.id, topic.id=topics, 
                 vocab=input$vocab, sort.topics="")
tokens <- input$vocab[input$word.id]
tab <- table(tokens, topics)
ord <- sort(apply(tab, 2, sum), decreasing=TRUE)
expect_true(all(names(ord) == gsub("Topic", "", colnames(dat$theta.hat))))



