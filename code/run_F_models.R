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
mod_bayes_lifespan <- stan_glmer(y_lifespan ~ Fgen_std * Sex +
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
mod_bayes_LRS <- stan_glmer(y_LRS ~ Fgen_std * Sex +
                              (1 | Site) +
                              (1 | Origin) + 
                              (1 | hatch),
                            data = inddata,
                            family = poisson,
                            iter = num_iter,
                            chains = num_chains,
                            control = list(adapt_delta = 0.98, max_treedepth = 30))

# nest-level model - clutch_size - number of eggs per clutch
mod_bayes_eggs_nest <- stan_glmer(y_eggs_nest ~ male_F_std + female_F_std + 
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

# nest-level model - number of hatchlings per clutch
mod_bayes_hatch_nest <- stan_glmer(y_hatch_nest ~ male_F_std + female_F_std + 
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
mod_bayes_fledge_nest <- stan_glmer(y_fledge_nest ~ male_F_std + female_F_std + 
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
plot(mod_bayes_lifespan, plotfun = "mcmc_trace")
plot(mod_bayes_LRS, plotfun = "mcmc_trace")
plot(mod_bayes_eggs_nest, plotfun = "mcmc_trace")
plot(mod_bayes_hatch_nest, plotfun = "mcmc_trace")
plot(mod_bayes_fledge_nest, plotfun = "mcmc_trace")

# calculate R2 for each model
r2_lifetime <- round(cor(fitted(mod_bayes_lifespan), inddata$y_lifespan) ** 2, 3)
r2_LRS <- round(cor(fitted(mod_bayes_LRS), inddata$y_LRS) ** 2, 3)
r2_eggs_nest <- round(cor(fitted(mod_bayes_eggs_nest), nestdata$y_eggs_nest[-mod_bayes_eggs_nest$na.action]) ** 2, 3)
r2_hatch_nest <- round(cor(fitted(mod_bayes_hatch_nest), nestdata$y_hatch_nest[-mod_bayes_hatch_nest$na.action]) ** 2, 3)
r2_fledge_nest <- round(cor(fitted(mod_bayes_fledge_nest), nestdata$y_fledge_nest[-mod_bayes_fledge_nest$na.action]) ** 2, 3)

#plot figure 1
pdf("./figures/Figure1.pdf", width=9, height=5)
par(mfrow = c(1, 2),
    mar = c(5, 5, 2, 5))
m_adj <- -0.4
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.2
m_line <- 0
xlimits_ind <- c(-0.7, 0.3)
old_mar <- par()$mar
par(mar = c(5.1, 9.1, 1.1, 1.1))
box_plot(mod_bayes_LRS, pars = c("Fgen_std", "Fgen_std:SexM"), digits = 3, bwd = 0.12,
         xlim = xlimits_ind)
par(mar = old_mar)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("LRS", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
old_mar <- par()$mar
par(mar = c(5.1, 9.1, 1.1, 1.1))
box_plot(mod_bayes_lifespan, pars = c("Fgen_std", "Fgen_std:SexM"), digits = 3, bwd = 0.12,
         xlim = xlimits_ind)
par(mar = old_mar)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Lifespan", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)

dev.off()

#plot figure 3
pdf("./figures/Figure3.pdf", width=4, height=10)
par(mfrow = c(3, 1), mar = c(5.1, 5.1, 2.1, 1.1))
m_adj <- -0.3
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.2
m_line <- 0
xlimits <- c(-0.4, 0.55)
old_mar <- par()$mar
par(mar = c(5.1, 9.1, 1.1, 1.1))
box_plot(mod_bayes_eggs_nest, pars = c("female_F_std", "male_F_std", "pair_GD_std",
                                       "female_age_std", "male_age_std"), digits = 3, bwd = 0.12,
         xlim = xlimits)
par(mar = old_mar)
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Eggs", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)

old_mar <- par()$mar
par(mar = c(5.1, 9.1, 1.1, 1.1))
box_plot(mod_bayes_hatch_nest, pars = c("female_F_std", "male_F_std", "pair_GD_std",
                                        "female_age_std", "male_age_std"), digits = 3, bwd = 0.12,
         xlim = xlimits)
par(mar = old_mar)
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Hatchlings", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)

old_mar <- par()$mar
par(mar = c(5.1, 9.1, 1.1, 1.1))
box_plot(mod_bayes_fledge_nest, pars = c("female_F_std", "male_F_std", "pair_GD_std",
                                         "female_age_std", "male_age_std"), digits = 3, bwd = 0.12,
         xlim = xlimits)
par(mar = old_mar)
mtext("C", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Fledglings", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)

dev.off()

#obtain ”point estimate” (posterior median)
print(mod_bayes_lifespan, digits = 3)
print(mod_bayes_LRS, digits = 3)
print(mod_bayes_eggs_nest, digits = 3)
print(mod_bayes_hatch_nest, digits = 3)
print(mod_bayes_fledge_nest, digits = 3)

#plot figure 2 
pdf("./figures/Figure2.pdf", width=10, height=5)
par(mfrow=c(1,2),
    mar = c(5, 5, 5, 5))
m_adj <- -0.33
m_adj2 <- 1
m_cex <- 1.5
m_cex2 <- 1.2
m_line <- 0
ind_plot(mod = mod_bayes_LRS,
         data = inddata,
         response = "y_LRS",
         plot_settings = list(ylab = "Lifetime reproductive success",
                              xlab = expression(italic("F")["gen"]),
                              variable = "Fgen"))
mtext("A", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("LRS", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
ind_plot(mod = mod_bayes_lifespan,
         data = inddata,
         response = "y_lifespan",
         plot_settings = list(ylab = "",
                              xlab = expression(italic("F")["gen"]),
                              variable = "Fgen"))
mtext("B", side = 3, line = m_line, adj = m_adj, cex = m_cex)
mtext("Lifespan", side = 3, line = m_line, adj = m_adj2, cex= m_cex2)
mtext("Lifetime in days", side =2, line=3.5, adj= 0.5, cex=1.2)

dev.off()

#predicting LRS for females
Fgen_values <- c(0.25, 0.125, 0.05, 0)
Fgen_std_values <- (Fgen_values - mean(inddata$Fgen)) / sd(inddata$Fgen)
predicted_LRS_F <- posterior_predict(mod_bayes_LRS,
                                     newdata = data.frame(Fgen_std = Fgen_std_values,
                                                          Sex = rep("F", length(Fgen_std_values))),
                                     re.form = NA)
predicted_LRS_F_values <- t(apply(predicted_LRS_F, 2, function(x) c("mean" = mean(x),  
                                                                    "median" = median(x),
                                                                    "q10" = quantile(x, p = 0.1),
                                                                    "q90" = quantile(x, p = 0.9))))

#predicting LRS for males
Fgen_values <- c(0.25, 0.125, 0.05, 0)
Fgen_std_values <- (Fgen_values - mean(inddata$Fgen)) / sd(inddata$Fgen)
predicted_LRS_M <- posterior_predict(mod_bayes_LRS,
                                     newdata = data.frame(Fgen_std = Fgen_std_values,
                                                          Sex = rep("M", length(Fgen_std_values))),
                                     re.form = NA)
predicted_LRS_M_values <- t(apply(predicted_LRS_M, 2, function(x) c("mean" = mean(x), 
                                                                    "median" = median(x),
                                                                    "q10" = quantile(x, p = 0.1),
                                                                    "q90" = quantile(x, p = 0.9))))

#predicting lifespan for females
Fgen_values <- c(0.25, 0.125, 0.05, 0)
Fgen_std_values <- (Fgen_values - mean(inddata$Fgen)) / sd(inddata$Fgen)
predicted_lifetime_f <- posterior_predict(mod_bayes_lifespan,
                                          newdata = data.frame(Fgen_std = Fgen_std_values,
                                                               Sex = rep("F", length(Fgen_std_values))),
                                          re.form = NA)
predicted_lifetime_f_values <- t(apply(predicted_lifetime_f, 2, function(x) c("mean" = mean(x), "median" = median(x),
                                                                              "q10" = quantile(x, p = 0.1),
                                                                              "q90" = quantile(x, p = 0.9))))

#predicting lifespan for males
Fgen_values <- c(0.25, 0.125, 0.05, 0)
Fgen_std_values <- (Fgen_values - mean(inddata$Fgen)) / sd(inddata$Fgen)
predicted_lifetime_m <- posterior_predict(mod_bayes_lifespan,
                                          newdata = data.frame(Fgen_std = Fgen_std_values,
                                                               Sex = rep("M", length(Fgen_std_values))),
                                          re.form = NA)
predicted_lifetime_m_values <- t(apply(predicted_lifetime_m, 2, function(x) c("mean" = mean(x), 
                                                                              "median" = median(x),
                                                                              "q10" = quantile(x, p = 0.1),
                                                                              "q90" = quantile(x, p = 0.9))))

#predicting fledglings per clutch for male F
Fgen_values <- c(0.25, 0.125, 0.05, 0)
male_F_std_values <- (Fgen_values - mean(nestdata$male_F)) / sd(nestdata$male_F)
predicted_fledgenest_male_F <- posterior_predict(mod_bayes_fledge_nest,
                                                 newdata = data.frame(male_F_std = male_F_std_values,
                                                                      female_F_std = rep(0, length(male_F_std_values)),
                                                                      pair_GD_std = rep(0, length(male_F_std_values)), 
                                                                      male_age_std =rep(0, length(male_F_std_values)),
                                                                      female_age_std = rep(0, length(male_F_std_values))),
                                                 re.form = NA)
predicted_fledgenest_male_F_values <- t(apply(predicted_fledgenest_male_F, 2, function(x) c("mean" = mean(x), 
                                                                                            "median" = median(x),
                                                                                            "q10" = quantile(x, p = 0.1),
                                                                                            "q90" = quantile(x, p = 0.9))))

#predicting fledglings per clutch for female F
Fgen_values <- c(0.25, 0.125, 0.05, 0)
female_F_std_values <- (Fgen_values - mean(nestdata$female_F)) / sd(nestdata$female_F)
predicted_fledgenest_female_F <- posterior_predict(mod_bayes_fledge_nest,
                                                   newdata = data.frame(female_F_std = female_F_std_values,
                                                                        male_F_std = rep(0, length(female_F_std_values)),
                                                                        pair_GD_std = rep(0, length(male_F_std_values)), 
                                                                        male_age_std =rep(0, length(male_F_std_values)),
                                                                        female_age_std = rep(0, length(male_F_std_values))),
                                                   re.form = NA)
predicted_fledgenest_female_F_values <- t(apply(predicted_fledgenest_female_F, 2, function(x) c("mean" = mean(x),
                                                                                                "median" = median(x),
                                                                                                "q10" = quantile(x, p = 0.1),
                                                                                                "q90" = quantile(x, p = 0.9))))

#predicting hatchlings per clutch for male F
Fgen_values <- c(0.25, 0.125, 0.05, 0)
male_F_std_values <- (Fgen_values - mean(nestdata$male_F)) / sd(nestdata$male_F)
predicted_hatchnest_male_F <- posterior_predict(mod_bayes_hatch_nest,
                                                newdata = data.frame(male_F_std = male_F_std_values,
                                                                     female_F_std = rep(0, length(male_F_std_values)),
                                                                     pair_GD_std = rep(0, length(male_F_std_values)), 
                                                                     male_age_std =rep(0, length(male_F_std_values)),
                                                                     female_age_std = rep(0, length(male_F_std_values))),
                                                re.form = NA)
predicted_hatchnest_male_F_values <- t(apply(predicted_hatchnest_male_F, 2, function(x) c("mean" = mean(x), 
                                                                                          "median" = median(x),
                                                                                          "q10" = quantile(x, p = 0.1),
                                                                                          "q90" = quantile(x, p = 0.9))))

#predicting hatchlings per clutch for female F
Fgen_values <- c(0.25, 0.125, 0.05, 0)
female_F_std_values <- (Fgen_values - mean(nestdata$female_F)) / sd(nestdata$female_F)
predicted_hatchnest_female_F <- posterior_predict(mod_bayes_hatch_nest,
                                                  newdata = data.frame(female_F_std = female_F_std_values,
                                                                       male_F_std = rep(0, length(female_F_std_values)),
                                                                       pair_GD_std = rep(0, length(female_F_std_values)), 
                                                                       male_age_std =rep(0, length(female_F_std_values)),
                                                                       female_age_std = rep(0, length(female_F_std_values))),
                                                  re.form = NA)
predicted_hatchnest_female_F_values <- t(apply(predicted_hatchnest_female_F, 2, function(x) c("mean" = mean(x), 
                                                                                              "median" = median(x),
                                                                                              "q10" = quantile(x, p = 0.1),
                                                                                              "q90" = quantile(x, p = 0.9))))

#predicting fledglings per clutch for pair GD
GD_values <- c(0.1, 0.18, 0.25)
GD_std_values <- (GD_values - mean(nestdata$pair_GD, na.rm=TRUE)) / sd(nestdata$pair_GD, na.rm=TRUE)
predicted_fledgenest_gd <- posterior_predict(mod_bayes_fledge_nest,
                                             newdata = data.frame(pair_GD_std = GD_std_values,
                                                                  female_F_std = rep(0, length(GD_std_values)),
                                                                  male_F_std = rep(0, length(GD_std_values)),
                                                                  male_age_std =rep(0, length(GD_std_values)),
                                                                  female_age_std = rep(0, length(GD_std_values))),
                                             re.form = NA)
predicted_fledgenest_gd_values <- t(apply(predicted_fledgenest_gd, 2, function(x) c("mean" = mean(x), 
                                                                                    "q10" = quantile(x, p = 0.1),
                                                                                    "q90" = quantile(x, p = 0.9))))

#predicting hatchlings per clutch for pair GD
GD_values <- c(0.1, 0.18, 0.25)
GD_std_values <- (GD_values - mean(nestdata$pair_GD, na.rm=TRUE)) / sd(nestdata$pair_GD, na.rm=TRUE)
predicted_hatchnest_gd <- posterior_predict(mod_bayes_hatch_nest,
                                            newdata = data.frame(pair_GD_std = GD_std_values,
                                                                 female_F_std = rep(0, length(GD_std_values)),
                                                                 male_F_std = rep(0, length(GD_std_values)),
                                                                 male_age_std =rep(0, length(GD_std_values)),
                                                                 female_age_std = rep(0, length(GD_std_values))),  
                                            re.form = NA)
predicted_hatchnest_gd_values <- t(apply(predicted_hatchnest_gd, 2, function(x) c("mean" = mean(x), "q10" = quantile(x, p = 0.1),
                                                                                  "q90" = quantile(x, p = 0.9))))
#predicting fledgings per clutch for male age
age_values <- c(2, 4, 6, 8, 10)
age_std_values <- (age_values - mean(nestdata$male_age, na.rm=TRUE)) / sd(nestdata$male_age, na.rm=TRUE)
predicted_fledgenest_age_m <- posterior_predict(mod_bayes_fledge_nest,
                                              newdata = data.frame(male_age_std = age_std_values,
                                                                   pair_GD_std = rep(0, length(age_std_values)),
                                                                   female_F_std = rep(0, length(age_std_values)),
                                                                   male_F_std = rep(0, length(age_std_values)),
                                                                   female_age_std = rep(0, length(age_std_values))),  
                                              re.form = NA)
predicted_fledgenest_age_m_values <- t(apply(predicted_fledgenest_age_m, 2, function(x) c("mean" = mean(x), "q10" = quantile(x, p = 0.1),
                                                                                      "q90" = quantile(x, p = 0.9))))

#predicting hatchlings per clutch for male age
age_values <- c(2, 4, 6, 8, 10)
age_std_values <- (age_values - mean(nestdata$male_age, na.rm=TRUE)) / sd(nestdata$male_age, na.rm=TRUE)
predicted_hatchnest_age_m <- posterior_predict(mod_bayes_hatch_nest,
                                             newdata = data.frame(male_age_std = age_std_values,
                                                                  pair_GD_std = rep(0, length(age_std_values)),
                                                                  female_F_std = rep(0, length(age_std_values)),
                                                                  male_F_std = rep(0, length(age_std_values)),
                                                                  female_age_std = rep(0, length(age_std_values))),  
                                             re.form = NA)
predicted_hatchnest_age_m_values <- t(apply(predicted_hatchnest_age_m, 2, function(x) c("mean" = mean(x), "q10" = quantile(x, p = 0.1),
                                                                                    "q90" = quantile(x, p = 0.9)
)))


#predicting fledglings per clutch for female age
age_values <- c(2, 4, 6, 8, 10)
age_std_values <- (age_values - mean(nestdata$female_age, na.rm=TRUE)) / sd(nestdata$female_age, na.rm=TRUE)
predicted_fledgenest_age_f <- posterior_predict(mod_bayes_fledge_nest,
                                                newdata = data.frame(female_age_std = age_std_values,
                                                                     pair_GD_std = rep(0, length(age_std_values)),
                                                                     female_F_std = rep(0, length(age_std_values)),
                                                                     male_F_std = rep(0, length(age_std_values)),
                                                                     male_age_std = rep(0, length(age_std_values))),  
                                                re.form = NA)
predicted_fledgenest_age_values_f <- t(apply(predicted_fledgenest_age_f, 2, function(x) c("mean" = mean(x),
                                                                                          "median" = median(x),
                                                                                          "q10" = quantile(x, p = 0.1),
                                                                                          "q90" = quantile(x, p = 0.9))))                                              
##predicting hatchlings per clutch for female age
age_values <- c(2, 4, 6, 8, 10)
age_std_values <- (age_values - mean(nestdata$female_age, na.rm=TRUE)) / sd(nestdata$female_age, na.rm=TRUE)
predicted_hatchnest_age_f <- posterior_predict(mod_bayes_hatch_nest,
                                               newdata = data.frame(female_age_std = age_std_values,
                                                                    pair_GD_std = rep(0, length(age_std_values)),
                                                                    female_F_std = rep(0, length(age_std_values)),
                                                                    male_F_std = rep(0, length(age_std_values)),
                                                                    male_age_std = rep(0, length(age_std_values))),  
                                               re.form = NA)
predicted_hatchnest_age_values_f <- t(apply(predicted_hatchnest_age_f, 2, function(x) c("mean" = mean(x), 
                                                                                        "median" = median(x),
                                                                                        "q10" = quantile(x, p = 0.1),
                                                                                        "q90" = quantile(x, p = 0.9))))                                              
