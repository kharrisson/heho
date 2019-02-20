# plotting functions for heho analysis
ind_plot <- function(mod, data, response, plot_settings = list()) {
  
  # set some general plot settings
  plot_set <- list(col_pal = viridis::inferno(256, alpha = 0.8)[c(40, 165, 120)],
                   col_pal_line = viridis::inferno(256, alpha = 1.0)[c(40, 165, 120)],
                   darkness = 0.45,
                   interval_width = 0.8,
                   num_pred = 200,
                   cex.axis = 1.2,
                   cex.lab = 1.2,
                   ylab = response,
                   xlab = "Heterozygosity",
                   variable = "MLH",
                   gs_value = 0,
                   hom_value = 0,
                   ylim_set = NULL)
  plot_set[names(plot_settings)] <- plot_settings

  # create a sequence of MLH values  
  pred_data <- data[[paste0(plot_set$variable, "_std")]]
  pred_seq <- data.frame(female = seq(min(pred_data[data$Sex == "F"], na.rm = TRUE),
                                      max(pred_data[data$Sex == "F"], na.rm = TRUE),
                                      length = plot_set$num_pred),
                         male = seq(min(pred_data[data$Sex == "M"], na.rm = TRUE),
                                    max(pred_data[data$Sex == "M"], na.rm = TRUE),
                                    length = plot_set$num_pred))
  
  # extract MCMC samples
  female_newdata <- data.frame(HOM_std = rep(plot_set$HOM_value, length(pred_seq$female)),
                               mean_gs_std = rep(plot_set$gs_value, length(pred_seq$female)),
                               Sex = rep("F", length(pred_seq$female)))
  female_newdata[[paste0(plot_set$variable, "_std")]] <- pred_seq$female
  female_fitted <- posterior_predict(mod,
                                     newdata = female_newdata,
                                     re.form = NA)
  female_fitted <- apply(female_fitted, 2, function(x) c("q10" = quantile(x, p = 0.1),
                                                         "mean" = mean(x),  
                                                         "q90" = quantile(x, p = 0.9)))
  male_newdata <- data.frame(HOM_std = rep(plot_set$HOM_value, length(pred_seq$male)),
                             mean_gs_std = rep(plot_set$gs_value, length(pred_seq$male)),
                             Sex = rep("M", length(pred_seq$male)))
  male_newdata[[paste0(plot_set$variable, "_std")]] <- pred_seq$male
  male_fitted <- posterior_predict(mod,
                                   newdata = male_newdata,
                                   re.form = NA)
  male_fitted <- apply(male_fitted, 2, function(x) c("q10" = quantile(x, p = 0.1),
                                                     "mean" = mean(x),
                                                     "q90" = quantile(x, p = 0.9)))
  ylim_set <- range(c(data[[response]], female_fitted, male_fitted), na.rm = TRUE)
  if (!is.null(plot_set$ylim_set))
    ylim_set <- plot_set$ylim_set
  plot(data[[response]] ~ pred_data,
       bty = "l",
       xlab = plot_set$xlab, ylab = plot_set$ylab,
       las = 1,
       cex.lab = plot_set$cex.lab,
       cex.axis = plot_set$cex.axis,
       type = "n",
       xaxt = "n",
       ylim = ylim_set,
       col = plot_set$col_pal[as.integer(data$Sex)])
  polygon(c(pred_seq$female, rev(pred_seq$female)),
          c(female_fitted[1, ], rev(female_fitted[3, ])),
          border = NA, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  polygon(c(pred_seq$male, rev(pred_seq$male)),
          c(male_fitted[1, ], rev(male_fitted[3, ])),
          border = NA, col = ggplot2::alpha(plot_set$col_pal[2], plot_set$darkness))
  lines(female_fitted[2, ] ~ pred_seq$female, lty = 1, lwd = 2,
        col = plot_set$col_pal_line[1])
  lines(male_fitted[2, ] ~ pred_seq$male, lty = 1, lwd = 2,
        col = plot_set$col_pal_line[2])
  points(data[[response]] ~ pred_data,
         pch = 16,
         col = plot_set$col_pal[as.integer(data$Sex)],
         cex = 1.2)
  axis(1, cex.axis = plot_set$cex.axis, at = seq(min(pred_seq), max(pred_seq), length = 5),
       labels = round(seq(min(data[[plot_set$variable]], na.rm = TRUE),
                          max(data[[plot_set$variable]], na.rm = TRUE), length = 5), 2))
  
}

nest_plot <- function(mod, data, response, plot_settings = list()) {
  
  # set some general plot settings
  plot_set <- list(col_pal = viridis::inferno(256, alpha = 0.8)[c(40, 200, 120)],
                   darkness = 0.35,
                   interval_width = 0.8,
                   num_pred = 150,
                   mtext_line = 0.5,
                   mtext_adj = 0,
                   mtext_cex = 1.5,
                   cex.axis = 1.3,
                   cex.lab = 1.3,
                   ylab = response,
                   xlab = "heterozygosity",
                   variable = "MLH",
                   num_trials = "n",
                   plot_prop = FALSE)
  plot_set[names(plot_settings)] <- plot_settings
  
  # extract samples from fitted model
  sims <- as.matrix(mod)
  
  # create sequences of MLH and GS
  female_pred_data <- data[[paste0("female_", plot_set$variable, "_std")]]
  male_pred_data <- data[[paste0("male_", plot_set$variable, "_std")]]
  pred_seq <- data.frame(female = seq(min(female_pred_data, na.rm = TRUE),
                                      max(female_pred_data, na.rm = TRUE),
                                      length = plot_set$num_pred),
                         male = seq(min(male_pred_data, na.rm = TRUE),
                                    max(male_pred_data, na.rm = TRUE),
                                    length = plot_set$num_pred),
  GS_seq <- seq(min(data$pair_GS_std, na.rm = TRUE),
                max(data$pair_GS_std, na.rm = TRUE),
                length = plot_set$num_pred),
  male_age_seq <- seq(min(data$male_age_std, na.rm = TRUE),
                max(data$male_age_std, na.rm = TRUE),
                length = plot_set$num_pred),
  female_age_seq <- seq(min(data$female_age_std, na.rm = TRUE),
                      max(data$female_age_std, na.rm = TRUE),
                      length = plot_set$num_pred))
  
  
  # plot female MLH effect
  if (!plot_set$plot_prop) {
    fitted <- exp(sweep(sims[, paste0("female_", plot_set$variable, "_std")] %o% pred_seq$female, 1, sims[, "(Intercept)"], "+"))
    plot_data <- data[[response]]
  } else {
    fitted <- plogis(sweep(sims[, paste0("female_", plot_set$variable, "_std")] %o% pred_seq$female, 1,
                           sims[, "(Intercept)"], "+")) 
    plot_data <- data[[response]] / data[[plot_set$num_trials]]
  }
  fitted <- apply(fitted, 2, quantile, p = c(0.5 - (plot_set$interval_width / 2),
                                             0.5,
                                             0.5 + (plot_set$interval_width / 2)))
  plot(plot_data ~ female_pred_data,
       bty = "l",
       xlab = expression(paste("Female ", "Homozygosity")), ylab = plot_set$ylab,
       las = 1,
       type = "n",
       xaxt = "n",
       cex.lab = plot_set$cex.lab,
       cex.axis = plot_set$cex.axis,
       ylim = range(c(plot_data, fitted), na.rm = TRUE),
       col = plot_set$col_pal[1])
  polygon(c(pred_seq$female, rev(pred_seq$female)),
          c(fitted[1, ], rev(fitted[3, ])),
          border = NA, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  lines(fitted[2, ] ~ pred_seq$female, lwd = 2, lty = 1, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  points(plot_data ~ female_pred_data,
         pch = 16,
         cex = 1.2)
  axis(1, cex.axis = plot_set$cex.axis, at = seq(min(female_pred_data, na.rm = TRUE),
                   max(female_pred_data, na.rm = TRUE), length = 5),
       labels = round(seq(min(data[[paste0("female_", plot_set$variable)]], na.rm = TRUE),
                          max(data[[paste0("female_", plot_set$variable)]], na.rm = TRUE), length = 5), 2))
  mtext("A", side = 3, line = plot_set$mtext_line, adj = plot_set$mtext_adj, cex = plot_set$mtext_cex)

  # plot male MLH effect
  if (!plot_set$plot_prop) {
    fitted <- exp(sweep(sims[, paste0("male_", plot_set$variable, "_std")] %o% pred_seq$male, 1, sims[, "(Intercept)"], "+"))
  } else {
    fitted <- plogis(sweep(sims[, paste0("male_", plot_set$variable, "_std")] %o% pred_seq$male, 1,
                           sims[, "(Intercept)"], "+")) 
  }
  fitted <- apply(fitted, 2, quantile, p = c(0.5 - (plot_set$interval_width / 2),
                                             0.5, 
                                             0.5 + (plot_set$interval_width / 2)))
  plot(plot_data ~ male_pred_data,
       bty = "l",
       xlab = expression(paste("Male ", "Homozygosity")), ylab = plot_set$ylab,
       las = 1,
       type = "n",
       xaxt = "n",
       cex.lab = plot_set$cex.lab,
       cex.axis = plot_set$cex.axis,
       ylim = range(c(plot_data, fitted), na.rm = TRUE),
       col = plot_set$col_pal[1])
  polygon(c(pred_seq$male, rev(pred_seq$male)),
          c(fitted[1, ], rev(fitted[3, ])),
          border = NA, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  lines(fitted[2, ] ~ pred_seq$male, lwd = 2, lty = 1, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  points(plot_data ~ male_pred_data,
         pch = 16,
         cex = 1.2)
  axis(1, cex.axis = plot_set$cex.axis, at = seq(min(male_pred_data, na.rm = TRUE),
                   max(male_pred_data, na.rm = TRUE), length = 5),
       labels = round(seq(min(data[[paste0("male_", plot_set$variable)]], na.rm = TRUE),
                          max(data[[paste0("male_", plot_set$variable)]], na.rm = TRUE), length = 5), 2))
  mtext("B", side = 3, line = plot_set$mtext_line, adj = plot_set$mtext_adj, cex = plot_set$mtext_cex)

  # plot GS effects
  if (!plot_set$plot_prop) {
    fitted <- exp(sweep(sims[, "pair_GS_std"] %o% GS_seq, 1, sims[, "(Intercept)"], "+"))
  } else {
    fitted <- plogis(sweep(sims[, "pair_GS_std"] %o% GS_seq, 1, sims[, "(Intercept)"], "+")) 
  }
  fitted <- apply(fitted, 2, quantile, p = c(0.5 - (plot_set$interval_width / 2),
                                             0.5,
                                             0.5 + (plot_set$interval_width / 2)))
  
  plot(plot_data ~ data$pair_GS_std,
       bty = "l",
       xlab = "Pairwise GS", ylab = plot_set$ylab,
       las = 1,
       type = "n",
       cex.lab = plot_set$cex.lab,
       cex.axis = plot_set$cex.axis,
       xaxt = "n",
       ylim = range(c(plot_data, fitted), na.rm = TRUE),
       col = plot_set$col_pal[1])
  polygon(c(GS_seq, rev(GS_seq)),
          c(fitted[1, ], rev(fitted[3, ])),
          border = NA, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  lines(fitted[2, ] ~ GS_seq, lwd = 2, lty = 1, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  points(plot_data ~ data$pair_GS_std,
         pch = 16,
         cex = 1.2)
  axis(1, cex.axis = plot_set$cex.axis, at = seq(min(data$pair_GS_std, na.rm = TRUE),
                   max(data$pair_GS_std, na.rm = TRUE), length = 5),
       labels = round(seq(min(data$pair_GS, na.rm = TRUE),
                          max(data$pair_GS, na.rm = TRUE), length = 5), 2))
  mtext("C", side = 3, line = plot_set$mtext_line, adj = plot_set$mtext_adj, cex = plot_set$mtext_cex)

}

# line plots of stan model parameters
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

}

# box plots of stan model parameters
box_plot <- function(model, pars, labels = NULL, bwd = 0.05, xlim = NULL, ...) {
  
  mod_summary <- summary(model, probs = c(0.025, 0.1, 0.5, 0.9, 0.975), ...)
  
  plot_vals <- mod_summary[, c(1, 3:8)]
  
  xlim_set <- xlim
  if (is.null(xlim_set))
    xlim_set <- range(plot_vals)
  
  if (!missing(pars))
    plot_vals <- plot_vals[match(pars, rownames(plot_vals)), ]
  
  x_set <- rev(seq_len(nrow(plot_vals)))
  
  plot(plot_vals[, "50%"], x_set, type = "n", bty = "l",
       xlab = "Estimate", ylab = "",
       yaxt = "n", xlim = xlim_set, ylim = c(0.5, max(x_set) + 0.5),
       cex.axis = 1.4, cex.lab = 1.4)
  
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

}
