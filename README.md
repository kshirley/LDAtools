LDAtools
========

R package for fitting topic models. 

<b>Update</b>: We changed the name of this repo from 'LDAviz' to 'LDAtools' on 5/5/2014. It will remain a repo with tools for preprocessing raw text and fitting LDA topic models in R (with C code as the back-end to run the collapsed Gibbs sampler). Two notes:

* For visualizing the output of a topic model, please check out our repo <a href='https://github.com/cpsievert/LDAvis' target='_blank'>LDAvis</a>hosted by Carson Sievert. All future work on <em>visualizing</em> topic models will be done in this repo.

* For fitting topic models, there are other software packages available, including MALLET and the R packages 'topicmodels' and 'lda', that are much more popular and better-tested (for speed and accuracy) than this package. This package was developed, more or less, (1) for practice building an R package and (2) to learn about LDA, rather than to become a widely-used package for others. So thanks for checking this out, but we'd recommend MALLET or other existing R packages for fitting topic models, and cpsievert/LDAvis for visualizing of topic models.

<b>Older README</b>:

### Online examples:

* [xkcd comics](http://glimmer.rstudio.com/cpsievert/xkcd/) with [full analysis](http://bit.ly/19Dmedr)
* [Associated Press articles](http://glimmer.rstudio.com/cpsievert/LDAviz/)

### Installation:

```library(devtools); install_github("LDAtools", "kshirley"); library(LDAtools)```

### Run AP example locally:

Make sure you have the following packages installed:

```install.packages(c("plyr","reshape", "proxy", "shiny"))```

```library(shiny); runApp(system.file('shiny', 'hover', package='LDAtools'))```

More documentation to come...
