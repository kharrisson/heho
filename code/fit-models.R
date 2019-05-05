# optional: set working directory
# setwd("/PATH/TO/DIR")

# load packages
library(rstanarm)
library(dplyr)

# source helper functions
source("code/plot_functions.R")
source("code/extract_pr_pos.R")

# read in data
inddata <- read.csv("data/S1_clean_data_used_in_ind_models120419.csv",
                    header = TRUE,
                    stringsAsFactors = FALSE)
nestdata <- read.csv("data/S2_clean_data_used_in_nest_models120419.csv",
                     header = TRUE, 
                     stringsAsFactors = FALSE)

# define model settings
num_iter <- 20000
num_chains <- 4

# fit models
# individual-level model - lifespan 
mod_bayes_lifespan_HOM <- stan_glmer(y_lifespan ~ HOM_std * Sex +
                                       (1 | Site) +
                                       (1 | Origin) +
                                       (1 | hatch),
                                     data = inddata,
                                     family = poisson,
                                     iter = num_iter,
                                     chains = num_chains,
                                     na.action = "na.omit",
                                     control = list(adapt_delta = 0.98, max_treedepth = 30),
                                     cores = 4)

# individual-level model - lifetime reproductive success (LRS) - all individuals N=94
mod_bayes_LRS_HOM <- stan_glmer(y_LRS ~ HOM_std * Sex +
                                  (1 | Site) +
                                  (1 | Origin) +
                                  (1 | hatch),
                                data = inddata,
                                family = poisson,
                                iter = num_iter,
                                chains = num_chains,
                                control = list(adapt_delta = 0.98, max_treedepth = 30),
                                cores = 4)

# individual-level model - lifetime reproductive success (LRS) - subset individuals 
#which have LRS and pair genetic information, N=64
mod_bayes_LRS_HOM_mean_gs <- stan_glmer(y_LRS ~ HOM_std * Sex + mean_gs_std * Sex + HOM_std * mean_gs_std +
                                           (1 | Site) +
                                           (1 | Origin) +
                                           (1 | hatch),
                                         data = inddata,
                                         family = poisson,
                                         iter = num_iter,
                                         chains = num_chains,
                                         control = list(adapt_delta = 0.98, max_treedepth = 30),
                                         cores = 4)

# nest-level model - clutch_size
mod_bayes_eggs_nest_HOM <- stan_glmer(y_eggs_nest ~ male_HOM_std + female_HOM_std +
                                        + pair_GS_std + (1 | nest_site) +
                                        + male_age_std + female_age_std +
                                        (1 | season) + (1 | male_origin) +
                                        (1 | female_origin) + (1 | male_id) +
                                        (1 | female_id),
                                      data = nestdata,
                                      family = poisson,
                                      iter = num_iter,
                                      chains = num_chains,
                                      na.action = "na.omit",
                                      control = list(adapt_delta = 0.98, max_treedepth = 30),
                                      cores = 4)

# nest-level model - number of hatchlings per clutch
mod_bayes_hatch_nest_HOM <- stan_glmer(y_hatch_nest ~ male_HOM_std + female_HOM_std + 
                                         pair_GS_std + male_age_std + female_age_std +
                                         (1 | nest_site) +
                                         (1 | season) + (1 | male_origin) +
                                         (1 | female_origin) + (1 | male_id) +
                                         (1 | female_id),
                                       data = nestdata,
                                       family = poisson,
                                       iter = num_iter,
                                       chains = num_chains,
                                       na.action = "na.omit",
                                       control = list(adapt_delta = 0.98, max_treedepth = 30),
                                       cores = 4)

# nest-level model - number of fledglings per clutch
mod_bayes_fledge_nest_HOM <- stan_glmer(y_fledge_nest ~ male_HOM_std + female_HOM_std + 
                                          pair_GS_std + male_age_std + female_age_std +
                                          (1 | nest_site) +
                                          (1 | season) + (1 | male_origin) +
                                          (1 | female_origin) + (1 | male_id) +
                                          (1 | female_id),
                                        data = nestdata,
                                        family = poisson,
                                        iter = num_iter,
                                        chains = num_chains,
                                        na.action = "na.omit",
                                        control = list(adapt_delta = 0.98, max_treedepth = 30),
                                        cores = 4)

# save outputs to use in plotting scripts
saveRDS(mod_bayes_lifespan_HOM, file = "outputs/BBBmod_bayes_lifespan_HOM.rds")
saveRDS(mod_bayes_LRS_HOM, file = "outputs/BBBmod_bayes_LRS_HOM.rds")
saveRDS(mod_bayes_LRS_HOM_mean_gs, file = "outputs/BBBmod_bayes_LRS_HOM_mean_gs.rds")
saveRDS(mod_bayes_eggs_nest_HOM, file = "outputs/BBBmod_bayes_eggs_nest_HOM.rds")
saveRDS(mod_bayes_hatch_nest_HOM, file = "outputs/BBBmod_bayes_hatch_nest_HOM.rds")
saveRDS(mod_bayes_fledge_nest_HOM, file = "outputs/BBBmod_bayes_fledge_nest_HOM.rds")
