library(cmdstanr)
library(tidybayes)
library(tidyr)
library(bettermc)
library(dplyr)
setwd("/users/project1/pt01268/cmdstanr")
set_cmdstan_path("~/cmdstan/2.38.0")

file <- file.path("/users/project1/pt01268/cmdstanr", "mod107sim.stan")
mod107sim <- cmdstan_model(file)

obj <- readRDS("input_data.rds")
datastruct_sim = obj$datastruct_sim
logkobs_future = obj$logkobs_future
current_best = obj$current_best
x_candidates = obj$x_candidates

nAnalytes_sim =datastruct_sim$nAnalytes_sim 
nColumns =  datastruct_sim$nColumns

init_sim <- function(){
  list( etaStd  = t(matrix(rep(0,nAnalytes_sim*3), nrow=nAnalytes_sim)),
        dzeta1  = matrix(rep(0,nColumns*nAnalytes_sim), nrow=nAnalytes_sim),
        dzeta2  = matrix(rep(0,nColumns*nAnalytes_sim), nrow=nAnalytes_sim),
        dzeta3  = matrix(rep(0,nColumns*nAnalytes_sim), nrow=nAnalytes_sim))
}

log_sum_exp <- function(x, y) {
  m <- pmax(x, y)
  m + log(exp(x - m) + exp(y - m))
}

log1p_exp <- function(x) {
  ifelse(x > 0,
         x + log1p(exp(-x)),  # for large x
         log1p(exp(x)))       # for small or negative x
}

funlogki_fi_v102 <- function(logkw, S1,  S2, cgamma, alpha, pHmpKa, xlogkw, xS1, xS2, fi) {
  
  log10 <- 2.302585092994046  # ln(10)
  
  logfix1 <- - (1+xS2) * xS1 * fi / (1 + xS2 * fi) 
  logfix2 <- - (1+S2) * S1  * fi / (1 + S2 * fi)
  
  
  t2 <- (pHmpKa  + alpha * fi) * log10
  
  a <- (xlogkw + cgamma * fi + logfix1) * log10
  c <- (logkw + logfix2) * log10 + t2
  
  logk <- log_sum_exp(a, c) / log10 -
    log1p_exp(t2) / log10
  
  return(logk)
}

funlogki_pKa_v102 <- function(logkw, S1, logS2, cgamma, fid,
                              alpha, pHmpKa, dlogkw, dS1, dlogS2, fi) {
  
  S2 <- 10^logS2
  xlogkw <- logkw+dlogkw
  xS1 <- S1+dS1
  xS2 <- 10^(logS2+dlogS2)
  
  logkd <- funlogki_fi_v102(
    logkw, S1,  S2, cgamma,
    alpha, pHmpKa, xlogkw, xS1, xS2, fid
  )
  
  logk_raw <- funlogki_fi_v102(
    logkw, S1, S2,  cgamma,
    alpha, pHmpKa, xlogkw, xS1, xS2, fi
  )
  
  w = plogis((fid - fi) / 0.001);
  logk = w*logkd+(1-w)*logk_raw
}

utility_fn <- function(param_draws, fi_new) {
  param_draws %>%
    tidybayes::spread_draws(param[i,c,..]) %>%
    pivot_wider(
      id_cols = c(.draw,c),
      names_from = i,
      values_from = starts_with("param")) %>%
    tidyr::expand_grid(fi = fi_new) %>%
    mutate(logk1 = funlogki_pKa_v102(param.1_1, param.2_1, param.3_1,
                                     param.4_1, param.5_1, param.6_1,
                                     param.7_1, param.8_1, param.9_1,
                                     param.10_1,fi)) %>%
    mutate(logk2 = funlogki_pKa_v102(param.1_2, param.2_2, param.3_2,
                                     param.4_2, param.5_2, param.6_2,
                                     param.7_2, param.8_2, param.9_2,
                                     param.10_2,fi)) %>%
    mutate(
      d  = abs(logk1 - logk2),
      lo = 10^(pmin(logk1, logk2)),
      hi = 10^(pmax(logk1, logk2)), # midpoint for selectivity of 1.2, zero for 1.04 and 1 for 1.5
      utylity = plogis((d-0.07918125)/0.015)*plogis((lo-2)/0.25)*plogis((10-hi)/1.5)) %>%
    select(fi, c, utylity)
}

expected_utility <- function(draws, fi_new) {
  utility_fn(draws,fi_new) %>%
    group_by(c, fi) %>%
    summarize(prob = mean(utylity)) %>% ungroup()
}

safe_worker <- function(m, .stanmodel = stanmodel, 
                        .datastruct_sim  = datastruct_sim, 
                        .logkobs_sim_exp = logkobs_sim_exp,
                        .current_best = current_best) {
  
  logkobsi              = .logkobs_sim_exp %>% dplyr::filter(draw == m)
  datastruct_sim_update = .datastruct_sim
  datastruct_sim_update$nObs_sim = datastruct_sim_update$nObs_sim + length(logkobsi$n_logkobs)
  datastruct_sim_update$logkobs_sim = c(datastruct_sim_update$logkobs_sim, logkobsi$n_logkobs)
  datastruct_sim_update$fi_sim = c(datastruct_sim_update$fi_sim, logkobsi$n_fi)
  datastruct_sim_update$analyte_sim = c(datastruct_sim_update$analyte_sim, logkobsi$n_analyte)
  datastruct_sim_update$column_sim = c(datastruct_sim_update$column_sim, logkobsi$n_column)
  
  tryCatch({

    fit_updated <- .stanmodel$sample(
      data = datastruct_sim_update,
      init = init_sim,
      iter_warmup = 500,
      iter_sampling = 1000,
      chains = 1,
      parallel_chains = 1,
      refresh = 100,
      show_messages = FALSE,
      adapt_delta=0.95
    )
    
    param_updated       <- fit_updated$draws("param")
    exp_utility_updated <- expected_utility(param_updated, seq(0, 1, 0.05))
    updated_best        <- max(exp_utility_updated$prob)
    
    gains <- updated_best - .current_best
    gains
    
  }, error = function(e) {
    message(sprintf("Simulation %s failed: %s", m, e$message))
    gains = NULL
  })
}

n_future = 192

task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
logkobs_sim_exp = logkobs_future %>% filter(n_exp == x_candidates[task_id])
 
gains_list <- bettermc::mclapply(
    1:n_future,
    \(x) safe_worker(x, .stanmodel = mod107sim, 
                     .datastruct_sim  = datastruct_sim, 
                     .logkobs_sim_exp = logkobs_sim_exp,
                     .current_best = current_best),
    mc.cores = 48,
    mc.set.seed = TRUE,
    mc.preschedule = FALSE, mc.allow.fatal = NULL,
    mc.timeout.elapsed = 180, mc.force.fork = TRUE,
    mc.retry = 2
  )
 
print(unlist(Filter(Negate(is.null), gains_list)))
evsi <- mean(unlist(Filter(Negate(is.null), gains_list)),na.rm = TRUE)

saveRDS(evsi, file = sprintf("evsi_%03d-%s.rds", task_id, obj$case_names))
