# plot fitted helmeted honeyeater models

# optional: set working directory
# setwd("PATH/TO/DIR")

# load packages
library(rstanarm)
library(dplyr)

# source helper functions
source("code/plot_functions.R")

# read in data
inddata <- read.csv("data/S1_clean_data_used_in_ind_models120419.csv",
                    header = TRUE,
                    stringsAsFactors = FALSE)
nestdata <- read.csv("data/S2_clean_data_used_in_nest_models120419.csv",
                     header = TRUE, 
                     stringsAsFactors = FALSE)


# load fitted models
mod_bayes_lifespan_HOM <- readRDS("outputs/BBBmod_bayes_lifespan_HOM.rds")
mod_bayes_LRS_HOM_mean_gs <- readRDS("outputs/BBBmod_bayes_LRS_HOM_mean_gs.rds")
mod_bayes_eggs_nest_HOM <- readRDS("outputs/BBBmod_bayes_eggs_nest_HOM.rds")
mod_bayes_hatch_nest_HOM <- readRDS("outputs/BBBmod_bayes_hatch_nest_HOM.rds")
mod_bayes_fledge_nest_HOM <- readRDS("outputs/BBBmod_bayes_fledge_nest_HOM.rds")

# plot figure 1
pdf("figures/Figure1_revision.pdf", width = 11, height = 4)
par(mfrow = c(1, 3),
    mar = c(5, 5, 2, 5))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.2
m_line <- 0
xlimits <- c(-0.6, 0.6)
old_mar <- par()$mar
par(mar = c(5.1, 12.1, 1.7, 1.1))
box_plot(mod_bayes_LRS_HOM_mean_gs, pars = c("HOM_std", "mean_gs_std", "HOM_std:SexM", "SexM:mean_gs_std", "HOM_std:mean_gs_std"), digits = 3, bwd = 0.12,
         xlim = xlimits, labels=c("HOM", "Mean GS", "HOM:SexM", "Mean GS:SexM", "HOM:Mean GS"))
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
old_mar <- par()$mar
par(mar = c(5.1, 5, 1.7, 1))
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.4, cex.lab = 1.4,
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
old_mar <- par()$mar
par(mar = c(5.1, 5, 1.7, 5))
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Mean mate similarity",
                              variable = "mean_gs", cex.axis = 1.4, cex.lab = 1.4,
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("C", side = 3, line = m_line, adj = m_adj, cex = m_cex)

dev.off()

# plot figure 2
pdf("figures/Figure2.pdf", width = 10, height = 5)
par(mfrow = c(1, 3),
    mar = c(5, 5, 5, 2))
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.7, cex.lab = 1.7, ylim_set = c(0, 80),
                              gs_value = -1.3,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("A - Low mate GS", side = 3, line = 0, adj = 0, cex = 1.4)
ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.7, cex.lab = 1.7, ylim_set = c(0, 80),
                              gs_value = 0,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("B - Mean mate GS", side = 3, line = 0, adj = 0, cex = 1.4)

ind_plot(mod = mod_bayes_LRS_HOM_mean_gs,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Total fledglings produced",
                              xlab = "Homozygosity",
                              variable = "HOM", cex.axis = 1.7, cex.lab = 1.7, ylim_set = c(0, 80),
                              gs_value = 4,   # only matters if plotting HOM (high = what does HOM do when GS is high?)
                              HOM_value = 0)) # onlmy matters if plotting GS (high = what does GS do when HOM is high?)
mtext("C - High mate GS", side = 3, line = 0, adj = 0, cex = 1.4)
dev.off()

# plot figure 3
pdf("figures/Figure3_revised.pdf", width = 4, height = 10)
par(mfrow = c(3, 1), mar = c(5.1, 5.1, 5.1, 1.1))
m_adj <- 0
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1
m_line <- 2
xlimits <- c(-0.5, 0.3)
old_mar <- par()$mar
par(mar = c(5.1, 9.2, 3.1, 1.1))
box_plot(mod_bayes_eggs_nest_HOM, pars = c("female_HOM_std", "male_HOM_std", 
                                           "pair_GS_std", 
                                           "female_age_std", 
                                           "male_age_std"), digits = 3, bwd = 0.12, 
                                            xlim = xlimits, labels = c("Female HOM",
                                            "Male HOM", "Pair GS", "Female age", 
                                            "Male age"))
par(mar = old_mar)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Eggs", side = 3, line = m_line, adj = m_adj2, cex = m_cex2)
old_mar <- par()$mar
par(mar = c(5.1, 9.2, 3.1, 1.1))
box_plot(mod_bayes_hatch_nest_HOM, pars = c("female_HOM_std", 
                                            "male_HOM_std", "pair_GS_std",
                                            "female_age_std", "male_age_std"),
                                             digits = 3, bwd = 0.12, xlim = xlimits,
                                             labels = c("Female HOM",
                                             "Male HOM", "Pair GS", "Female age", 
                                             "Male age"))
par(mar = old_mar)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Hatchlings", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
old_mar <- par()$mar
par(mar = c(5.1, 9.2, 3.1, 1.1))
box_plot(mod_bayes_fledge_nest_HOM, pars = c("female_HOM_std", 
                                             "male_HOM_std", "pair_GS_std",
                                             "female_age_std", "male_age_std"),
                                             digits = 3, bwd = 0.12, xlim = xlimits,
                                             labels = c("Female HOM",
                                             "Male HOM", "Pair GS", "Female age", 
                                             "Male age"))
par(mar = old_mar)
mtext("C", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Fledglings", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
dev.off()

# plot Fig S2 
pdf("figures/FigureS2_revised.pdf", width = 10, height = 5)
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
         xlim = xlimits, labels=c("HOM", "HOM:SexM"))
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

# plot supp figure S3
pdf("./figures/FigureS3a-c.pdf", width = 4, height = 10)
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
pdf("figures/FigureS3d-e.pdf", width = 4, height = 10)
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
