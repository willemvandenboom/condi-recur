# condi-recur

Repository with the code used for the paper "Bayesian inference on the number of recurrent events: A joint model of recurrence and survival" by Willem van den Boom, Maria De Iorio and Marta Tallarita ([arXiv:2005.06819](https://arxiv.org/abs/2005.06819)).


## Description of files

* [`main.Rmd`](main.Rmd) contains the code that runs the Gibbs sampler and creates the plots for the application to the atrial fibrillation data. It loads [`Gibbs_sampler.R`](Gibbs_sampler.R). When knitting, it loads computationally expensive Gibbs sampling results. Therefore, knitting it requires one to first run the relevant R chunks involving Gibbs sampling.

* [`main.html`](main.html) is the result of knitting [`main.Rmd`](main.Rmd). It does not contain any plots as those are saved as PDF files when [`main.Rmd`](main.Rmd) is run.

* [`simulation.Rmd`](simulation.Rmd) and [`simulation.html`](simulation.html) contain a simulation study in the same fashion as [`main.Rmd`](main.Rmd) and [`main.html`](main.html) contain the application to atrial fibrillation data.

* [`Gibbs_sampler.R`](Gibbs_sampler.R) consists of an R function that implements the Gibbs sampler for the main method from the paper.