LDAviz
======

R package for fitting and visualizing topic models. 

### Online examples:

* [xkcd comics](http://glimmer.rstudio.com/cpsievert/xkcd/) with [full analysis](http://bit.ly/19Dmedr)
* [Associated Press articles](http://glimmer.rstudio.com/cpsievert/LDAviz/)

### Installation:

```library(devtools); install_github("LDAviz", "kshirley"); library(LDAviz)```

### Run AP example locally:

Make sure you have the following packages installed:

```install.packages(c("plyr","reshape", "proxy", "shiny"))```

```library(shiny); runApp(system.file('shiny', 'hover', package='LDAviz'))```

More documentation to come...
