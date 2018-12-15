# plotting functions for heho analysis

ind_plot <- function(mod, data, response, plot_settings = list()) {
  
  # set some general plot settings
  plot_set <- list(col_pal = viridis::inferno(256, alpha = 0.8)[c(40, 165, 120)],
                   col_pal_line = viridis::inferno(256, alpha = 1.0)[c(40, 165, 120)],
                   darkness = 0.45,
                   interval_width = 0.8,
                   num_pred = 150,
                   cex.axis = 1.2,
                   cex.lab = 1.2,
                   ylab = response,
                   xlab = "Heterozygosity",
                   variable = "MLH")
  plot_set[names(plot_settings)] <- plot_settings
  # plot_set$col_pal_line[2] <- "#883000"
  # plot_set$col_pal[2] <- "#883000"
  
  # create a sequence of MLH values  
  pred_data <- data[[paste0(plot_set$variable, "_std")]]
  pred_seq <- data.frame(female = seq(min(pred_data[data$Sex == "F"]),
                                      max(pred_data[data$Sex == "F"]),
                                      length = plot_set$num_pred),
                         male = seq(min(pred_data[data$Sex == "M"]),
                                    max(pred_data[data$Sex == "M"]),
                                    length = plot_set$num_pred))
  
  # extract MCMC samples
  sims <- as.matrix(mod)
  female_fitted <- exp(sweep(sims[, paste0(plot_set$variable, "_std")] %o% pred_seq$female, 1, sims[, "(Intercept)"], "+"))
  female_fitted <- apply(female_fitted, 2, quantile, p = c(0.5 - (plot_set$interval_width / 2),
                                                           0.5,
                                                           0.5 + (plot_set$interval_width / 2)))
  male_fitted <- exp(sweep((sims[, paste0(plot_set$variable, "_std")] + sims[, paste0(plot_set$variable, "_std:SexM")]) %o% pred_seq$male,
                           1,
                           (sims[, "(Intercept)"] + sims[, "SexM"]),
                           "+"))
  male_fitted <- apply(male_fitted, 2, quantile, p = c(0.5 - (plot_set$interval_width / 2),
                                                       0.5,
                                                       0.5 + (plot_set$interval_width / 2)))
  plot(data[[response]] ~ pred_data,
       bty = "l",
       xlab = plot_set$xlab, ylab = plot_set$ylab,
       las = 1,
       cex.lab = plot_set$cex.lab,
       cex.axis = plot_set$cex.axis,
       type = "n",
       xaxt = "n",
       ylim = range(c(data[[response]], female_fitted, male_fitted), na.rm = TRUE),
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
       labels = round(seq(min(data[[plot_set$variable]]), max(data[[plot_set$variable]]), length = 5), 2))
  
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
  
  # create sequences of MLH and GD
  female_pred_data <- data[[paste0("female_", plot_set$variable, "_std")]]
  male_pred_data <- data[[paste0("male_", plot_set$variable, "_std")]]
  pred_seq <- data.frame(female = seq(min(female_pred_data, na.rm = TRUE),
                                      max(female_pred_data, na.rm = TRUE),
                                      length = plot_set$num_pred),
                         male = seq(min(male_pred_data, na.rm = TRUE),
                                    max(male_pred_data, na.rm = TRUE),
                                    length = plot_set$num_pred),
  GD_seq <- seq(min(data$pair_GD_std, na.rm = TRUE),
                max(data$pair_GD_std, na.rm = TRUE),
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
       xlab = expression(paste("Female ", italic("F")["gen"])), ylab = plot_set$ylab,
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
  mtext("D", side = 3, line = plot_set$mtext_line, adj = plot_set$mtext_adj, cex = plot_set$mtext_cex)

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
       xlab = expression(paste("Male ", italic("F")["gen"])), ylab = plot_set$ylab,
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
  mtext("E", side = 3, line = plot_set$mtext_line, adj = plot_set$mtext_adj, cex = plot_set$mtext_cex)

  # plot GD effects
  if (!plot_set$plot_prop) {
    fitted <- exp(sweep(sims[, "pair_GD_std"] %o% GD_seq, 1, sims[, "(Intercept)"], "+"))
  } else {
    fitted <- plogis(sweep(sims[, "pair_GD_std"] %o% GD_seq, 1, sims[, "(Intercept)"], "+")) 
  }
  fitted <- apply(fitted, 2, quantile, p = c(0.5 - (plot_set$interval_width / 2),
                                             0.5,
                                             0.5 + (plot_set$interval_width / 2)))
  
  plot(plot_data ~ data$pair_GD_std,
       bty = "l",
       xlab = "Pairwise GD", ylab = plot_set$ylab,
       las = 1,
       type = "n",
       cex.lab = plot_set$cex.lab,
       cex.axis = plot_set$cex.axis,
       xaxt = "n",
       ylim = range(c(plot_data, fitted), na.rm = TRUE),
       col = plot_set$col_pal[1])
  polygon(c(GD_seq, rev(GD_seq)),
          c(fitted[1, ], rev(fitted[3, ])),
          border = NA, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  lines(fitted[2, ] ~ GD_seq, lwd = 2, lty = 1, col = ggplot2::alpha(plot_set$col_pal[1], plot_set$darkness))
  points(plot_data ~ data$pair_GD_std,
         pch = 16,
         cex = 1.2)
  axis(1, cex.axis = plot_set$cex.axis, at = seq(min(data$pair_GD_std, na.rm = TRUE),
                   max(data$pair_GD_std, na.rm = TRUE), length = 5),
       labels = round(seq(min(data$pair_GD, na.rm = TRUE),
                          max(data$pair_GD, na.rm = TRUE), length = 5), 2))
  mtext("F", side = 3, line = plot_set$mtext_line, adj = plot_set$mtext_adj, cex = plot_set$mtext_cex)

}