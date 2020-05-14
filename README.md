# condi-recur

Repository with the code used for the paper "Bayesian Joint Modelling of Recurrence and Survival: a Conditional Approach" by Willem van den Boom, Marta Tallarita and Maria De Iorio.


## Description of files

* [`main.Rmd`](main.Rmd) contains the code that runs the Gibbs sampler and creates the plots. It loads [`Gibbs_sampler.R`](Gibbs_sampler.R). When knitting, it loads computationally expensive Gibbs sampling results. Therefore, knitting it requires one to first run the relevant R chunks involving Gibbs sampling.

* [`main.html`](main.html) is the result of knitting [`main.Rmd`](main.Rmd). It does not contain any plots as those are saved as PDF files when [`main.Rmd`](main.Rmd) is run.

* [`Gibbs_sampler.R`](Gibbs_sampler.R) consists of an R function that implements the Gibbs sampler for the main method from the paper.