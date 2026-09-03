library(cmdstanr)

set_cmdstan_path("~/cmdstan/2.38.0")

file <- file.path("/users/project1/pt01268/cmdstanr", "mod107sim.stan")

mod107sim <- cmdstan_model(file)
