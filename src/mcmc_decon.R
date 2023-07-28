### mcmc_decon.R --- MCMC implementation of linear regression
##
## Filename: mcmc_decon.R
## Author: Zach Maas
## Created: Wed Apr  5 12:14:11 2023 (-0600)
##
######################################################################
##
### Commentary:
##
## This file contains code to run a MCMC linear regression model on
## counts tables for deconvolution of large systems of counts data.
## This is presumptively a very slow choice of algorithm, but it is
## something that's flexible and will give us more error bounds.
##
######################################################################
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or (at
## your option) any later version.
##
## This program is distributed in the hope that it will be useful, but
## WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
## General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with GNU Emacs.  If not, see <https://www.gnu.org/licenses/>.
##
######################################################################
##
### Code:

## Load in the tidyverse
library(tidyverse)
## Load in the rstan package
library(rstan)

## Stan model recommendations
options(mc.cores=parallel::detectCores())
rstan_options(auto_write = TRUE)

## Simulate data for 6 samples
num_samples <- 6
size <- 100
X <- matrix(runif(n=num_samples*size,min=0,max=10),ncol=num_samples,nrow=size)
Y <- rowSums(X*(1/num_samples))

## Load in the stan model
count_data <- list(x=X,y=Y,N=nrow(X),K=ncol(X))
fit <- stan('src/linreg.stan', data = count_data, iter=10000, algorithm="NUTS")
posterior <- get_posterior_mean(fit)
size_fit <- ncol(posterior)
weights <- posterior[1:ncol(X)+1,size_fit]

######################################################################
### mcmc_decon.R ends here
