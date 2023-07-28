### analyze_collated_results.R --- Analyze collated results from simulation
##
## Filename: analyze_collated_results.R
## Author: Zach Maas
## Created: Tue Apr 18 11:48:26 2023 (-0600)
##
######################################################################
##
### Commentary:
##
## This file contains code to generate plots based off of collated
## simulations done on FIJI, which are relatively slow.
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

library("tidyverse")

## Load the data from /data/collated_simulations.csv
collated_results <- read_csv("./data/collated_simulations.csv",
                             col_names = c("N1", "N2", "N3",
                                           "N4", "N5", "N6",
                                           "N7", "N8", "N9",
                                           "iteration", "method", "epsilon"),
                             col_types = c("d", "d", "d",
                                           "d", "d", "d",
                                           "d", "d", "d",
                                           "i", "f", "d")) %>%
    mutate(N = N1 + N2 + N3 + N4 + N5 + N6 + N7 + N8 + N9,
           method = factor(method)) %>%
    arrange(method)

ggplot(data = collated_results %>% filter(method %in% c("mcmc", "nnls"))) +
    geom_line(aes(x = iteration, y = N1, linetype = method, color = "HCT116")) +
    geom_line(aes(x = iteration, y = N2, linetype = method, color = "HeLa"  )) + # nolint
    geom_line(aes(x = iteration, y = N3, linetype = method, color = "MCF7"  )) +
    geom_line(aes(x = iteration, y = N4, linetype = method, color = "K562"  )) +
    geom_line(aes(x = iteration, y = N5, linetype = method, color = "ESC"   )) +
    geom_line(aes(x = iteration, y = N6, linetype = method, color = "Kasumi")) +
    geom_line(aes(x = iteration, y = N7, linetype = method, color = "CD4"   )) +
    geom_line(aes(x = iteration, y = N8, linetype = method, color = "Jurkat")) +
    geom_line(aes(x = iteration, y = N9, linetype = method, color = "BJ5TA" )) +
    labs(x = "Iteration", y = "Estimated Proportion")



######################################################################
### analyze_collated_results.R ends here
