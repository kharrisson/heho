extract_pr_pos <- function(model, pars) {
  
  mod_tmp <- as.matrix(model)
  mat_sub <- mod_tmp[, pars]
  
  pr_pos <- apply(mat_sub, 2, function(x) sum(x > 0) / length(x))
  names(pr_pos) <- pars
  
  pr_pos
  
}
