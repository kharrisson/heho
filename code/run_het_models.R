#set working directory
setwd("~/Dropbox/heho/heho_clean/heho-osm/")

# load packages
library(rstanarm)

# source helper functions
source("./code/plot_functions_updated.R")
source("./code/plot_estimates.R")

# read in data
inddata <- read.csv("./data/clean_data_used_in_ind_models.csv",
                    header = TRUE,
                    stringsAsFactors = FALSE)
nestdata <- read.csv("./data/clean_data_used_in_nest_models.csv",
                     header = TRUE,
                     stringsAsFactors = FALSE)

#remove birds that are still alive from individual lifespan/LRS dataset
inddata <- subset(inddata, !(y_lifespan=="alive"))
inddata$y_lifespan <- as.numeric(as.character(inddata$y_lifespan))

# define model settings
#I ran in with 10000
num_iter <- 20000
num_chains <- 4

# fit models
# individual-level model - lifespan 
mod_bayes_lifespan_mlh <- stan_glmer(y_lifespan ~ MLH_std * Sex +
                                   (1 | Site) +
                                   (1 | Origin) +
                                   (1 | hatch),
                                 data = inddata,
                                 family = poisson,
                                 iter = num_iter,
                                 chains = num_chains,
                                 na.action = "na.omit",
                                 control = list(adapt_delta = 0.98, max_treedepth = 30))

# individual-level model - lifetime reproductive success (LRS) 
mod_bayes_LRS_mlh <- stan_glmer(y_LRS ~ MLH_std * Sex +
                              (1 | Site) +
                              (1 | Origin) +
                              (1 | hatch),
                            data = inddata,
                            family = poisson,
                            iter = num_iter,
                            chains = num_chains,
                            control = list(adapt_delta = 0.98, max_treedepth = 30))

# nest-level model - clutch_size - number of eggs per clutch
mod_bayes_eggs_nest_mlh <- stan_glmer(y_eggs_nest ~ male_MLH_std + female_MLH_std +
                                    + pair_GD_std + (1 | nest_site) +
                                    + male_age_std + female_age_std +
                                    (1 | season) + (1 | male_origin) +
                                    (1 | female_origin) + (1 | male_id) +
                                    (1 | female_id),
                                  data = nestdata,
                                  family = poisson,
                                  iter = num_iter,
                                  chains = num_chains,
                                  na.action = "na.omit",
                                  control = list(adapt_delta = 0.98, max_treedepth = 30))

# nest-level model - number of hatchlings per clutch
mod_bayes_hatch_nest_mlh_mlh_mlh <- stan_glmer(y_hatch_nest ~ male_MLH_std + female_MLH_std + 
                                     pair_GD_std + male_age_std + female_age_std +
                                     (1 | nest_site) +
                                     (1 | season) + (1 | male_origin) +
                                     (1 | female_origin) + (1 | male_id) +
                                     (1 | female_id),
                                   data = nestdata,
                                   family = poisson,
                                   iter = num_iter,
                                   chains = num_chains,
                                   na.action = "na.omit",
                                   control = list(adapt_delta = 0.98, max_treedepth = 30))

# nest-level model - number of fledglings per clutch
mod_bayes_fledge_nest_mlh <- stan_glmer(y_fledge_nest ~ male_MLH_std + female_MLH_std + 
                                      pair_GD_std + male_age_std + female_age_std +
                                      (1 | nest_site) +
                                      (1 | season) + (1 | male_origin) +
                                      (1 | female_origin) + (1 | male_id) +
                                      (1 | female_id),
                                    data = nestdata,
                                    family = poisson,
                                    iter = num_iter,
                                    chains = num_chains,
                                    na.action = "na.omit",
                                    control = list(adapt_delta = 0.98, max_treedepth = 30))

#trace plots 
plot(mod_bayes_lifespan_mlh, plotfun = "mcmc_trace")
plot(mod_bayes_LRS_mlh, plotfun = "mcmc_trace")
plot(mod_bayes_eggs_nest_mlh, plotfun = "mcmc_trace")
plot(mod_bayes_hatch_nest_mlh_mlh, plotfun = "mcmc_trace")
plot(mod_bayes_fledge_nest_mlh, plotfun = "mcmc_trace")

# calculate R2 for each model
r2_lifetime <- round(cor(fitted(mod_bayes_lifespan_mlh), inddata$y_lifespan) ** 2, 3)
r2_LRS <- round(cor(fitted(mod_bayes_LRS_mlh), inddata$y_LRS) ** 2, 3)
r2_eggs_nest <- round(cor(fitted(mod_bayes_eggs_nest_mlh), nestdata$y_eggs_nest) ** 2, 3)
r2_hatch_nest <- round(cor(fitted(mod_bayes_hatch_nest_mlh_mlh), nestdata$y_hatch_nest) ** 2, 3)
r2_fledge_nest <- round(cor(fitted(mod_bayes_fledge_nest_mlh), nestdata$y_fledge_nest) ** 2, 3)

#plot the model effects
#plot figure S2
pdf("./figures/ind_model_effects_mlh.pdf", width=9, height=5)
par(mfrow=c(1,2),
    mar = c(5, 5, 2, 5))
m_adj <- -0.4
m_cex <- 1.5
m_line <- 0
xlimits <- c(-0.5, 0.8)
old_mar <- par()$mar
par(mar = c(5.1, 7.1, 1.1, 1.1))
box_plot(mod_bayes_LRS_mlh, pars = c("MLH_std", "MLH_std:SexM"), digits = 3, bwd = 0.12,
         xlim = xlimits)
par(mar = old_mar)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
old_mar <- par()$mar
par(mar = c(5.1, 7.1, 1.1, 1.1))
box_plot(mod_bayes_lifespan_mlh, pars = c("MLH_std", "MLH_std:SexM"), digits = 3, bwd = 0.12,
         xlim = xlimits)
par(mar = old_mar)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)

dev.off()

#plot figure S3
pdf("./figures/nest_model_effects_mlh.pdf", width=4, height=10)
par(mfrow = c(3, 1), mar = c(5.1, 5.1, 2.1, 1.1))
m_adj <- -0.3
m_cex <- 1.5
m_line <- 0
xlimits <- c(-0.2, 0.6)
old_mar <- par()$mar
par(mar = c(5.1, 7.1, 1.1, 1.1))
box_plot(mod_bayes_eggs_nest_mlh, pars = c("female_MLH_std", "male_MLH_std", 
                                "pair_GD_std", 
                                "female_age_std", 
                                "male_age_std"), digits = 3, bwd = 0.12, xlim = xlimits)
par(mar = old_mar)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
old_mar <- par()$mar
par(mar = c(5.1, 7.1, 1.1, 1.1))
box_plot(mod_bayes_hatch_nest_mlh_mlh_mlh, pars = c("female_MLH_std", 
                                    "male_MLH_std", 
                                    "female_age_std", "male_age_std",
                                    "pair_GD_std"), digits = 3, bwd = 0.12, xlim = xlimits)
par(mar = old_mar)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
old_mar <- par()$mar
par(mar = c(5.1, 7.1, 1.1, 1.1))
box_plot(mod_bayes_fledge_nest_mlh, pars = c("female_MLH_std", 
                                  "male_MLH_std", 
                                  "female_age_std", "male_age_std",
                                  "pair_GD_std"), digits = 3, bwd = 0.12, xlim = xlimits)
par(mar = old_mar)
mtext("C", side = 3, line = m_line, adj = m_adj, cex = m_cex)
dev.off()

#obtain ”point estimate” (posterior median)
print(mod_bayes_lifespan_mlh, digits = 3)
print(mod_bayes_LRS_mlh, digits = 3)
print(mod_bayes_eggs_nest_mlh, digits = 3)
print(mod_bayes_hatch_nest_mlh_mlh, digits = 3)
print(mod_bayes_fledge_nest_mlh, digits = 3)
