# condi-recur

Repository with the code used for the paper "Bayesian inference on the number
of recurrent events: A joint model of recurrence and survival" by Willem van
den Boom, Maria De Iorio and Marta Tallarita
([arXiv:2005.06819](https://arxiv.org/abs/2005.06819) and
[doi:10.1177/09622802211048059](https://doi.org/10.1177/09622802211048059)).


## Description of files

* [`readmission.Rmd`](readmission.Rmd) contains the code that runs the Gibbs
sampler and creates the plots for the application to the colorectal cancer
data. It loads [`Gibbs_sampler.R`](Gibbs_sampler.R). When knitting, it loads
computationally expensive Gibbs sampling results. Therefore, knitting it
requires one to first run the relevant R chunks involving Gibbs sampling.

* [`simulation.Rmd`](simulation.Rmd) contains simulation studies in the same
fashion as [`main.Rmd`](main.Rmd) contains the application to the colorectal
cancer data.

* [`simulation_Gompertz.Rmd`](simulation_Gompertz.Rmd) performs a simulation
study with survival times generated from a Gompertz distribution.

* [`AF.Rmd`](AF.Rmd) contains the code for the application to the atrial
fibrillation data.

* [`Gibbs_sampler.R`](Gibbs_sampler.R) consists of an R function that
implements the Gibbs sampler for the main method from the paper.
