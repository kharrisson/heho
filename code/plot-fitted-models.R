# plot fitted helmeted honeyeater models

# optional: set working directory
# setwd("PATH/TO/DIR")

# load fitted models
mod_bayes_lifespan_HOM <- "outputs/mod_bayes_lifespan_HOM.rds"
mod_bayes_LRS_HOM <- "outputs/mod_bayes_LRS_HOM.rds"
mod_bayes_LRS_HOM_mean_gs <- "outputs/mod_bayes_LRS_HOM_mean_gs.rds"
mod_bayes_eggs_nest_HOM <- "outputs/mod_bayes_eggs_nest_HOM.rds"
mod_bayes_hatch_nest_HOM <- "outputs/mod_bayes_hatch_nest_HOM.rds"
mod_bayes_fledge_nest_HOM <- "outputs/mod_bayes_fledge_nest_HOM.rds"

# plot figure 1
pdf("figures/Figure1.pdf", width = 9, height = 4)
par(mfrow = c(1, 1),
    mar = c(5, 5, 2, 5))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.2
m_line <- 0
xlimits <- c(-0.6, 0.6)
old_mar <- par()$mar
par(mar = c(5.1, 12.1, 1.1, 1.1))
box_plot(mod_bayes_LRS_HOM_mean_gs, pars = c("HOM_std", "mean_gs_std", "HOM_std:SexM", "SexM:mean_gs_std", "HOM_std:mean_gs_std"), digits = 3, bwd = 0.12,
         xlim = xlimits)
par(mar = old_mar)
dev.off()

# plot supp figure model effects LRS full dataset
pdf("figures/supp_figure_LRS_full_dataset.pdf", width = 5, height = 5)
par(mfrow = c(1, 1),
    mar = c(5, 5, 2, 5))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.2
m_line <- 0
xlimits <- c(-0.6, 0.4)
old_mar <- par()$mar
par(mar = c(5.1, 9.1, 1.1, 1.1))
box_plot(mod_bayes_LRS_HOM, pars = c("HOM_std", "HOM_std:SexM"), digits = 3, bwd = 0.12,
         xlim = xlimits)
dev.off()

# plot fig 2
pdf("figures/Figure2.pdf", width = 10, height = 5)
par(mfrow = c(1, 2),
    mar = c(5, 5, 5, 5),
    mai = c(1, 1, 1, 0.5))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.4
m_line <- 0
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.4, cex.lab = 1.4,
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Mean similarity",
                              variable = "mean_gs", cex.axis = 1.4, cex.lab = 1.4,
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
dev.off()

# plot figure 3
pdf("figures/Figure3.pdf", width = 10, height = 5)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 5))
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "HOM",
                              variable = "HOM", cex.axis = 1.7, cex.lab = 1.7, ylim_set = c(0, 80),
                              gs_value = -1.3,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("A", side = 3, line = 0, adj = 0, cex = 1.4)
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "HOM",
                              variable = "HOM", cex.axis = 1.7, cex.lab = 1.7, ylim_set = c(0, 80),
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("B", side = 3, line = 0, adj = 0, cex = 1.4)

ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "HOM",
                              variable = "HOM", cex.axis = 1.7, cex.lab = 1.7, ylim_set = c(0, 80),
                              gs_value = 4,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("C", side = 3, line = 0, adj = 0, cex = 1.4)
dev.off()

# plot supp figure full LRS model
pdf("figures/supp_mat_LRS_full_model_plot.pdf", width = 5, height = 5)
par(mfrow = c(1, 1),
    mar = c(5, 5, 5, 5),
    mai = c(1, 1, 1, 0.5))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.4
m_line <- 0
ind_plot(mod = mod_bayes_LRS_HOM,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.4, cex.lab = 1.4,
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
dev.off()

# plot Fig 4 
pdf("figures/Figure4.pdf", width = 10, height = 5)
par(mfrow = c(1, 2),
    mar = c(5, 5, 2, 5))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.4
m_cex2 <- 1.4
m_line <- 0
xlimits <- c(-0.3, 0.3)
old_mar <- par()$mar
par(mar = c(5.1, 9.1, 2.1, 1.1))
box_plot(mod_bayes_lifespan_HOM, pars = c("HOM_std", "HOM_std:SexM"), digits = 3, bwd = 0.12,
         xlim = xlimits)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
old_mar <- par()$mar
par(mar = c(5.1, 6.1, 2.1, 1.1))
ind_plot(mod = mod_bayes_lifespan_HOM,
         data = inddata,
         response = "y_lifespan",
         plot_settings = list(ylab = "",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.4, cex.lab = 1.4,
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
title(ylab="Lifetime in days", line = 4, cex.lab = 1.4)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
dev.off()

# plot figure 5
pdf("figures/Figure5.pdf", width = 4, height = 10)
par(mfrow = c(3, 1), mar = c(5.1, 5.1, 2.1, 1.1))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1
m_line <- 0
xlimits <- c(-0.5, 0.3)
old_mar <- par()$mar
par(mar = c(5.1, 9.2, 1.1, 1.1))
box_plot(mod_bayes_eggs_nest_HOM, pars = c("female_HOM_std", "male_HOM_std", 
                                           "pair_GS_std", 
                                           "female_age_std", 
                                           "male_age_std"), digits = 3, bwd = 0.12, xlim = xlimits)
par(mar = old_mar)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Eggs", side = 3, line = m_line, adj = m_adj2, cex = m_cex2)
old_mar <- par()$mar
par(mar = c(5.1, 9.2, 1.1, 1.1))
box_plot(mod_bayes_hatch_nest_HOM, pars = c("female_HOM_std", 
                                            "male_HOM_std", "pair_GS_std",
                                            "female_age_std", "male_age_std"),
         digits = 3, bwd = 0.12, xlim = xlimits)
par(mar = old_mar)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Hatchlings", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
old_mar <- par()$mar
par(mar = c(5.1, 9.2, 1.1, 1.1))
box_plot(mod_bayes_fledge_nest_HOM, pars = c("female_HOM_std", 
                                             "male_HOM_std", "pair_GS_std",
                                             "female_age_std", "male_age_std"),
         digits = 3, bwd = 0.12, xlim = xlimits)
par(mar = old_mar)
mtext("C", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Fledglings", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
dev.off()

# plot supp figure
pdf("./figures/FigureS5.pdf", width = 4, height = 10)
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 2.1, 1.1))
m_adj <- -0.35
m_cex <- 1.5
m_line <- 0
par(mfrow = c(3, 1), mar = c(5.1, 5.1, 2.1, 1.1))
nest_plot(mod = mod_bayes_hatch_nest_HOM,
          data = nestdata,
          response = "y_hatch_nest",
          plot_settings = list(ylab = "Hatchlings per clutch",
                               xlab="HOM", col_pal = gray.colors(1, start=0.3, end = 0.3, gamma=2.2),
                               variable = c("HOM")))
dev.off()

# plot supp figure
pdf("figures/FigureS6.pdf", width = 4, height = 10)
par(mfrow = c(3, 2), mar = c(5.1, 5.1, 2.1, 1.1))
m_adj <- -0.35
m_cex <- 1.5
m_line <- 0
par(mfrow = c(3, 1), mar = c(5.1, 5.1, 2.1, 1.1))
nest_plot(mod = mod_bayes_fledge_nest_HOM,
          data = nestdata,
          response = "y_fledge_nest", 
          plot_settings = list(ylab = "Fledglings per clutch",
                               xlab="HOM", col_pal = gray.colors(1, start=0.3, end = 0.3, gamma=2.2),
                               variable = c("HOM")))
dev.off()
