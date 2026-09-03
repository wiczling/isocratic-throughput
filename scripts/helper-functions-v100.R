plot_prior_post_intervals <- function(fitp,fit,params_to_plot,transformations=list()){
  
  d_prior <- bayesplot::mcmc_intervals_data(fitp$draws(params_to_plot),transformations=transformations)
  d_post <- bayesplot::mcmc_intervals_data(fit$draws(params_to_plot),transformation=transformations)
  
  ggplot() +
    geom_linerange(
      data = d_prior,
      aes(y = parameter, xmin = ll, xmax = hh),
      linewidth = 0.6,
      colour = "grey80"
    ) +
    geom_linerange(
      data = d_prior,
      aes(y = parameter, xmin = l, xmax = h),
      linewidth = 2,
      colour = "grey60"
    ) +
    geom_point(
      data = d_prior,
      aes(y = parameter, x = m),
      colour = "grey40",
      size = 2
    ) +
    geom_linerange(
      data = d_post,
      aes(y = parameter, xmin = ll, xmax = hh),
      linewidth = 0.6,
      colour = "#1f77b4"
    ) +
    geom_linerange(
      data = d_post,
      aes(y = parameter, xmin = l, xmax = h),
      linewidth = 2,
      colour = "#1f77b4"
    ) +
    geom_point(
      data = d_post,
      aes(y = parameter, x = m),
      colour = "#1f77b4",
      size = 2
    ) +
    labs(x = NULL, y = NULL)
}

log_sum_exp <- function(x, y) {
  m <- pmax(x, y)
  m + log(exp(x - m) + exp(y - m))
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

funlogki_sil <- function(logkw, S1, logks, logS2, fid, pKasmpH, alpha, fi) {
  S2 <- 10^logS2
  delta = 0.1;
  logfixo = -(1+S2)*S1*fid/(1+S2*fid);
  logfix  = -(1+S2)*S1*fi /(1+S2*fi );
  logksfixo = logks - log1p_exp((pKasmpH-alpha*fid)*log(10))/log(10);
  logksfix  = logks - log1p_exp((pKasmpH-alpha*fi )*log(10))/log(10);
  logko = logksfixo+log1p_exp((logkw-logksfixo+logfixo)*log(10))/log(10);
  logk  = logksfix +log1p_exp((logkw-logksfix +logfix )*log(10))/log(10);
  
  return(logko - delta * log1p_exp((logko-logk) / delta))
}

log1p_exp <- function(x) {
  ifelse(x > 0,
         x + log1p(exp(-x)),  # for large x
         log1p(exp(x)))       # for small or negative x
}

`%notin%` <- Negate(`%in%`)

extend_fg_hierarchy <- function(fg_hierarchy, new_groups) {
  for (group in new_groups) {
    if (!("name" %in% names(group)) || !("smarts" %in% names(group))) {
      warning("Invalid custom functional group structure.")
      next
    }
    fg_hierarchy <- append(fg_hierarchy, list(group))
  }
  return(fg_hierarchy)
}
make_filename_safe <- function(name) {
  name <- iconv(name, from = "latin1", to = "UTF-8")
  name |>
    tolower() |>
    (\(x) gsub("[^a-z0-9]+", "_", x))() |>
    (\(x) gsub("^_|_$", "", x))() |>
    (\(x) substr(x, 1, 100))()
}

# Convert to lower triangle with indices
similarity_to_ltr_fun <- function(similarity_matrix){
  
  lower_tri_indices <- which(lower.tri(similarity_matrix, diag = TRUE), arr.ind = TRUE)
  lower_tri_df <- data.frame(
    row = lower_tri_indices[, 1],
    col = lower_tri_indices[, 2],
    similarity = similarity_matrix[lower_tri_indices] )
  
  return(lower_tri_df)}

make_similarity_matrix_fun <- function(lower_tri_df) {
  n <- max(lower_tri_df$row, lower_tri_df$col)
  similarity_matrix <- matrix(0, nrow = n, ncol = n)
  for (k in 1:nrow(lower_tri_df)) {
    i <- lower_tri_df$row[k]
    j <- lower_tri_df$col[k]
    similarity <- lower_tri_df$similarity[k]
    similarity_matrix[i, j] <- similarity
    similarity_matrix[j, i] <- similarity
  }
  diag(similarity_matrix) <- 1; 
  return(similarity_matrix)}



stan_init_from_flat <- function(x) {
  stopifnot(is.numeric(x), !is.null(names(x)))
  
  out <- list()
  nms <- names(x)
  base_names <- unique(sub("\\[.*", "", nms))
  
  for (nm in base_names) {
    sel <- grepl(paste0("^", nm, "(\\[|$)"), nms)
    vals <- x[sel]
    idx_raw <- gsub(paste0("^", nm, "\\[|\\]$"), "", nms[sel])
    
    # ---- scalar ----
    if (all(idx_raw == "")) {
      out[[nm]] <- vals
      next
    }
    
    # ---- parse indices safely ----
    idx_list <- lapply(idx_raw, function(s) {
      if (s == "") return(integer(0))
      as.integer(strsplit(s, ",", fixed = TRUE)[[1]])
    })
    
    # remove NA indices
    idx_list <- lapply(idx_list, function(ii) ii[!is.na(ii) & ii >= 1])
    
    # determine number of dimensions
    ndims <- max(lengths(idx_list))
    if (ndims == 0) {
      out[[nm]] <- vals
      next
    }
    
    # determine dimensions safely
    dims <- integer(ndims)
    for (d in seq_len(ndims)) {
      dims[d] <- max(
        sapply(idx_list, function(ii) if (length(ii) >= d) ii[d] else 1L),
        na.rm = TRUE
      )
      if (!is.finite(dims[d]) || dims[d] < 1) dims[d] <- 1L
    }
    
    # ---- create array ----
    arr <- array(NA_real_, dim = dims)
    
    # ---- fill array ----
    for (i in seq_along(vals)) {
      ii <- idx_list[[i]]
      if (length(ii) < ndims) ii <- c(ii, rep(1L, ndims - length(ii)))
      arr[matrix(ii, nrow = 1)] <- vals[i]
    }
    
    # drop to vector if 1D
    if (length(dims) == 1) arr <- as.vector(arr)
    
    out[[nm]] <- arr
  }
  
  out
}
