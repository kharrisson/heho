# calculate summary statistics from fitted helmeted honeyeater models

# optional: set working directory
# setwd("PATH/TO/DIR")

# load fitted models
mod_bayes_lifespan_HOM <- "outputs/mod_bayes_lifespan_HOM.rds"
mod_bayes_LRS_HOM <- "outputs/mod_bayes_LRS_HOM.rds"
mod_bayes_LRS_HOM_mean_gs <- "outputs/mod_bayes_LRS_HOM_mean_gs.rds"
mod_bayes_eggs_nest_HOM <- "outputs/mod_bayes_eggs_nest_HOM.rds"
mod_bayes_hatch_nest_HOM <- "outputs/mod_bayes_hatch_nest_HOM.rds"
mod_bayes_fledge_nest_HOM <- "outputs/mod_bayes_fledge_nest_HOM.rds"

# calculate R2 for each model
r2_lifetime <- round(cor(fitted(mod_bayes_lifespan_HOM), inddata$y_lifespan) ** 2, 3)
r2_LRS <- round(cor(fitted(mod_bayes_LRS_HOM), inddata$y_LRS) ** 2, 3)
r2_LRS_mean_gs <- round(cor(fitted(mod_bayes_LRS_HOM_mean_gs), c(inddata$y_LRS)[-mod_bayes_LRS_HOM_mean_gs$na.action]) ** 2, 3)
r2_eggs_nest <- round(cor(fitted(mod_bayes_eggs_nest_HOM), c(nestdata$y_eggs_nest)[-mod_bayes_eggs_nest_HOM$na.action]) ** 2, 3)
r2_hatch_nest <- round(cor(fitted(mod_bayes_hatch_nest_HOM), c(nestdata$y_hatch_nest)[-mod_bayes_hatch_nest_HOM$na.action]) ** 2, 3)
r2_fledge_nest <- round(cor(fitted(mod_bayes_fledge_nest_HOM), c(nestdata$y_fledge_nest)[-mod_bayes_fledge_nest_HOM$na.action]) ** 2, 3)

# obtain ”point estimate” (posterior median)
print(mod_bayes_lifespan_HOM, digits = 3)
print(mod_bayes_LRS_HOM, digits = 3)
print(mod_bayes_LRS_HOM_mean_gs, digits = 3)
print(mod_bayes_eggs_nest_HOM, digits = 3)
print(mod_bayes_hatch_nest_HOM, digits = 3)
print(mod_bayes_fledge_nest_HOM, digits = 3)

# print the credible intervals
posterior_interval(mod_bayes_lifespan_HOM,
                   pars = c("(Intercept)", "HOM_std", "SexM", "HOM_std:SexM"),  
                   prob = 0.8) 
posterior_interval(mod_bayes_LRS_HOM,
                   pars = c("(Intercept)", "HOM_std", "SexM", "HOM_std:SexM"),  
                   prob = 0.8) 
posterior_interval(mod_bayes_eggs_nest_HOM,
                   pars = c("(Intercept)", "male_F_std", "female_F_std", 
                            "pair_GS_std"),  
                   prob = 0.8) 
posterior_interval(mod_bayes_hatch_nest_HOM,
                   pars = c("(Intercept)", "male_F_std", "female_F_std", 
                            "pair_GS_std"),  
                   prob = 0.8) 
posterior_interval(mod_bayes_fledge_nest_HOM,
                   pars = c("(Intercept)", "male_F_std", "female_F_std", 
                            "pair_GS_std"),  
                   prob = 0.8) 

# calculate posterior probabilites that effects are positive
extract_pr_pos(mod_bayes_LRS_HOM,
               pars = c("(Intercept)", "HOM_std", "HOM_std:SexM")) 

extract_pr_pos(mod_bayes_lifespan_HOM,
               pars = c("(Intercept)", "HOM_std", "HOM_std:SexM")) 

extract_pr_pos(mod_bayes_eggs_nest_HOM,
               pars = c("(Intercept)", "male_HOM_std", "female_HOM_std", 
                        "pair_GS_std")) 
extract_pr_pos(mod_bayes_hatch_nest_HOM,
               pars = c("(Intercept)", "male_HOM_std", "female_HOM_std", 
                        "pair_GS_std")) 
extract_pr_pos(mod_bayes_fledge_nest_HOM,
               pars = c("(Intercept)", "male_HOM_std", "female_HOM_std", 
                        "pair_GS_std")) 