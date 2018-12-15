line_plot <- function(model, pars, labels = NULL, bwd = 0.05, xlim = NULL, ...) {
  
  mod_summary <- summary(model, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...)
  
  plot_vals <- mod_summary[, c(1, 3:8)]
  
  xlim_set <- xlim
  if (is.null(xlim_set))
    xlim_set <- range(plot_vals)
  
  if (!missing(pars))
    plot_vals <- plot_vals[match(pars, rownames(plot_vals)), ]
  
  x_set <- rev(seq_len(nrow(plot_vals)))
  
  plot(plot_vals[, "50%"], x_set, type = "n", bty = "l", xlab = "", ylab = "",
       yaxt = "n", xlim = xlim_set, ylim = c(0.5, max(x_set) + 0.5))
  
  for (i in seq_along(x_set)) {
    
    lines(c(plot_vals[i, "2.5%"], plot_vals[i, "97.5%"]), c(x_set[i], x_set[i]), lwd = 2.1)

    lines(c(plot_vals[i, "10%"], plot_vals[i, "90%"]), c(x_set[i], x_set[i]), lwd = 4.1)
    
  }
  
  points(plot_vals[, "50%"], x_set, pch = 16, cex = 1.5)
  
  lines(c(0, 0), c(0, max(x_set) + 1), lty = 2)
  
  if (is.null(labels))
    labels <- rownames(plot_vals)

  axis(2, at = x_set, labels = labels, las = 1)
  
  mtext("Estimate", side = 1, adj = 0.5, line = 2.5)
  #mtext("Predictor", side = 2, adj = 0.5, line = 5.9)
  
}


box_plot <- function(model, pars, labels = NULL, bwd = 0.05, xlim = NULL, ...) {
  
  mod_summary <- summary(model, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...)
  
  plot_vals <- mod_summary[, c(1, 3:8)]
  
  xlim_set <- xlim
  if (is.null(xlim_set))
    xlim_set <- range(plot_vals)
  
  if (!missing(pars))
    plot_vals <- plot_vals[match(pars, rownames(plot_vals)), ]
  
  x_set <- rev(seq_len(nrow(plot_vals)))
  
  plot(plot_vals[, "50%"], x_set, type = "n", bty = "l", xlab = "", ylab = "",
       yaxt = "n", xlim = xlim_set, ylim = c(0.5, max(x_set) + 0.5),
       cex.axis = 1.2)
  
  for (i in seq_along(x_set)) {
    
    lines(c(plot_vals[i, "2.5%"], plot_vals[i, "10%"]), c(x_set[i], x_set[i]), lwd = 2)
    lines(c(plot_vals[i, "90%"], plot_vals[i, "97.5%"]), c(x_set[i], x_set[i]), lwd = 2)
    
    lines(c(plot_vals[i, "10%"], plot_vals[i, "90%"]), c(x_set[i] - bwd, x_set[i] - bwd), lwd = 1.5)
    lines(c(plot_vals[i, "10%"], plot_vals[i, "90%"]), c(x_set[i] + bwd, x_set[i] + bwd), lwd = 1.5)
    lines(c(plot_vals[i, "10%"], plot_vals[i, "10%"]), c(x_set[i] - bwd, x_set[i] + bwd), lwd = 1.5)
    lines(c(plot_vals[i, "90%"], plot_vals[i, "90%"]), c(x_set[i] - bwd, x_set[i] + bwd), lwd = 1.5)
    
    lines(c(plot_vals[i, "50%"], plot_vals[i, "50%"]), c(x_set[i] - bwd, x_set[i] + bwd), lwd = 1.8)
    
  }
  
  lines(c(0, 0), c(0, max(x_set) + 1), lty = 2)
  
  if (is.null(labels))
    labels <- rownames(plot_vals)
  
  axis(2, cex.axis= 1.2, at = x_set, labels = labels, las = 1)
  
  mtext("Estimate", side = 1, adj = 0.5, line = 2.5, cex= 1)
  #mtext("Predictor", side = 2, adj = 0.5, line = 5.9)
  
}

#example usage
#old_mar <- par()$mar
#par(mar = c(5.1, 7.1, 1.1, 1.1))
#line_plot(m, pars = c("Fgen_std", "SexM", "Fgen_std:SexM"), digits = 3)
#par(mar = old_mar)

#old_mar <- par()$mar
#par(mar = c(5.1, 7.1, 1.1, 1.1))
#line_plot(m, pars = c("Fgen_std", "Fgen_std:SexM"), digits = 3)
#par(mar = old_mar)

#old_mar <- par()$mar
#par(mar = c(5.1, 7.1, 1.1, 1.1))
#box_plot(m, pars = c("Fgen_std", "Fgen_std:SexM"), digits = 3, bwd = 0.12)
#par(mar = old_mar)
