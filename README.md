# Code for "Efficient smoothness selection for nonparametric Markov-switching models via quasi restricted maximum likelihood"

This repository contains `R` code for reproducing the case studies and simulation experiments of the paper <a href="https://arxiv.org/abs/2411.11498" target="_blank">"Efficient smoothness selection for nonparametric Markov-switching models via quasi restricted maximum likelihood"</a>.

In the paper, we develop an efficient method for selecting the smoothness parameter in nonparametric Markov-switching models called **quasi restricted maximum likelihood** (qREML).
The method is implemented in the `qreml()` function contained in the R package <a href="https://janoleko.github.io/LaMa/" target="_blank">`LaMa`</a>.

The structure of the repository is really simple as all `.R` files are self-contained. 
The folder `case_studies` contains the three case studies from the paper and the folder `simulation_experiments` contains one `.R` file for reproducing the simulation experiments.

The figures presented in the paper are provided in the subfolders `figs` and all data necessary will either be downloaded automatically by running the case study code, or is included in the folder `data` and will be loaded automatically.

When using <a href="https://janoleko.github.io/LaMa/" target="_blank">`LaMa`</a>, please cite the package as follows:

Koslik Jan-Ole (2024). LaMa: Fast Numerical Maximum Likelihood Estimation for Latent Markov Models. R package version 2.0.1 https://CRAN.R-project.org/package=LaMa.

Or in `R` type citation(package = "LaMa").