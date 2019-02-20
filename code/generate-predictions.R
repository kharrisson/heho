# predicting LRS for females
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_HOM <- posterior_predict(mod_bayes_LRS_HOM,
                                       newdata = data.frame(HOM_std = HOM_std_values,
                                                            Sex = rep("F", length(HOM_std_values))),
                                       re.form = NA)
predicted_LRS_HOM_values <- t(apply(predicted_LRS_HOM, 2, function(x) c("mean" = mean(x),  
                                                                        "median" = median(x),
                                                                        "q10" = quantile(x, p = 0.1),
                                                                        "q90" = quantile(x, p = 0.9))))
# predicting LRS for females - gs model - mean gs
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_HOM_mean_gs <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                               newdata = data.frame(HOM_std = HOM_std_values,
                                                                    mean_gs_std = rep(0, length(HOM_std_values)),
                                                                    Sex = rep("F", length(HOM_std_values))),
                                               re.form = NA)
predicted_LRS_HOM_values_mean_gs <- t(apply(predicted_LRS_HOM_mean_gs, 2, function(x) c("mean" = mean(x),  
                                                                                        "median" = median(x),
                                                                                        "q10" = quantile(x, p = 0.1),
                                                                                        "q90" = quantile(x, p = 0.9))))

# predicting LRS for females - gs model - low gs
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_HOM_mean_gs_low <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                                   newdata = data.frame(HOM_std = HOM_std_values,
                                                                        mean_gs_std = rep(-1.2, length(HOM_std_values)),
                                                                        Sex = rep("F", length(HOM_std_values))),
                                                   re.form = NA)
predicted_LRS_HOM_values_mean_gs_low <- t(apply(predicted_LRS_HOM_mean_gs_low, 2, function(x) c("mean" = mean(x),  
                                                                                                "median" = median(x),
                                                                                                "q10" = quantile(x, p = 0.1),
                                                                                                "q90" = quantile(x, p = 0.9))))
# predicting LRS for females - gs model - high gs
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_HOM_mean_gs_high <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                                    newdata = data.frame(HOM_std = HOM_std_values,
                                                                         mean_gs_std = rep(3.5, length(HOM_std_values)),
                                                                         Sex = rep("F", length(HOM_std_values))),
                                                    re.form = NA)
predicted_LRS_HOM_values_mean_gs_high <- t(apply(predicted_LRS_HOM_mean_gs_high, 2, function(x) c("mean" = mean(x),  
                                                                                                  "median" = median(x),
                                                                                                  "q10" = quantile(x, p = 0.1),
                                                                                                  "q90" = quantile(x, p = 0.9))))

# predicting LRS for females - gs effect
mean_GS_values <- c(0.74, 0.78, 0.83, 0.86, 0.90)
mean_GS_std_values <- (mean_GS_values - mean(inddata$mean_gs, na.rm=TRUE)) / sd(inddata$mean_gs, na.rm=TRUE)
predicted_LRS_F_mean_gs_2 <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                               newdata = data.frame(mean_gs_std = mean_GS_std_values,
                                                                    HOM_std = rep(0, length(mean_GS_std_values)),
                                                                    Sex = rep("F", length(mean_GS_std_values))),
                                               re.form = NA)

predicted_LRS_F_values_mean_gs_2 <- t(apply(predicted_LRS_F_mean_gs_2, 2, function(x) c("mean" = mean(x), 
                                                                                        "median" = median(x),
                                                                                        "q10" = quantile(x, p = 0.1),
                                                                                        "q90" = quantile(x, p = 0.9))))

# predicting LRS for males
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_M <- posterior_predict(mod_bayes_LRS_HOM,
                                     newdata = data.frame(HOM_std = HOM_std_values,
                                                          Sex = rep("M", length(HOM_std_values))),
                                     re.form = NA)
predicted_LRS_M_values <- t(apply(predicted_LRS_M, 2, function(x) c("mean" = mean(x), 
                                                                    "median" = median(x),
                                                                    "q10" = quantile(x, p = 0.1),
                                                                    "q90" = quantile(x, p = 0.9))))

# predicting LRS for males - gs model - mean gs
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_M_mean_gs <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                             newdata = data.frame(HOM_std = HOM_std_values,
                                                                  mean_gs_std = rep(0, length(HOM_std_values)),
                                                                  Sex = rep("M", length(HOM_std_values))),
                                             re.form = NA)
predicted_LRS_M_values_mean_gs <- t(apply(predicted_LRS_M_mean_gs, 2, function(x) c("mean" = mean(x), 
                                                                                    "median" = median(x),
                                                                                    "q10" = quantile(x, p = 0.1),
                                                                                    "q90" = quantile(x, p = 0.9))))
# predicting LRS for males - gs model - low gs
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_M_mean_gs_low <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                                 newdata = data.frame(HOM_std = HOM_std_values,
                                                                      mean_gs_std = rep(-1.2, length(HOM_std_values)),
                                                                      Sex = rep("M", length(HOM_std_values))),
                                                 re.form = NA)
predicted_LRS_M_values_mean_gs_low <- t(apply(predicted_LRS_M_mean_gs_low, 2, function(x) c("mean" = mean(x), 
                                                                                            "median" = median(x),
                                                                                            "q10" = quantile(x, p = 0.1),
                                                                                            "q90" = quantile(x, p = 0.9))))

# predicting LRS for males - gs model - high gs
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_LRS_M_mean_gs_high <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                                  newdata = data.frame(HOM_std = HOM_std_values,
                                                                       mean_gs_std = rep(3.5, length(HOM_std_values)),
                                                                       Sex = rep("M", length(HOM_std_values))),
                                                  re.form = NA)
predicted_LRS_M_values_mean_gs_high <- t(apply(predicted_LRS_M_mean_gs_high, 2, function(x) c("mean" = mean(x), 
                                                                                              "median" = median(x),
                                                                                              "q10" = quantile(x, p = 0.1),
                                                                                              "q90" = quantile(x, p = 0.9))))

# predicting LRS for males - gs effect
mean_GS_values <- c(0.74, 0.78, 0.83, 0.86, 0.90)
mean_GS_std_values <- (mean_GS_values - mean(inddata$mean_gs, na.rm=TRUE)) / sd(inddata$mean_gs, na.rm=TRUE)
predicted_LRS_M_mean_gs_2 <- posterior_predict(mod_bayes_LRS_HOM_mean_gs,
                                               newdata = data.frame(mean_gs_std = mean_GS_std_values,
                                                                    HOM_std = rep(0, length(mean_GS_std_values)),
                                                                    Sex = rep("M", length(mean_GS_std_values))),
                                               re.form = NA)

predicted_LRS_M_values_mean_gs_2 <- t(apply(predicted_LRS_M_mean_gs_2, 2, function(x) c("mean" = mean(x), 
                                                                                        "median" = median(x),
                                                                                        "q10" = quantile(x, p = 0.1),
                                                                                        "q90" = quantile(x, p = 0.9))))

# predicting lifespan for females
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_lifetime_f <- posterior_predict(mod_bayes_lifespan_HOM,
                                          newdata = data.frame(HOM_std = HOM_std_values,
                                                               Sex = rep("F", length(HOM_std_values))),
                                          re.form = NA)
predicted_lifetime_f_values <- t(apply(predicted_lifetime_f, 2, function(x) c("mean" = mean(x), "median" = median(x),
                                                                              "q10" = quantile(x, p = 0.1),
                                                                              "q90" = quantile(x, p = 0.9))))

# predicting lifespan for males
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
HOM_std_values <- (HOM_values - mean(inddata$HOM)) / sd(inddata$HOM)
predicted_lifetime_m <- posterior_predict(mod_bayes_lifespan_HOM,
                                          newdata = data.frame(HOM_std = HOM_std_values,
                                                               Sex = rep("M", length(HOM_std_values))),
                                          re.form = NA)
predicted_lifetime_m_values <- t(apply(predicted_lifetime_m, 2, function(x) c("mean" = mean(x), 
                                                                              "median" = median(x),
                                                                              "q10" = quantile(x, p = 0.1),
                                                                              "q90" = quantile(x, p = 0.9))))

# predicting fledglings per clutch for male HOM
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
male_HOM_std_values <- (HOM_values - mean(nestdata$male_HOM)) / sd(nestdata$male_HOM)
predicted_fledgenest_male_HOM <- posterior_predict(mod_bayes_fledge_nest_HOM,
                                                   newdata = data.frame(male_HOM_std = male_HOM_std_values,
                                                                        female_HOM_std = rep(0, length(male_HOM_std_values)),
                                                                        pair_GS_std = rep(0, length(male_HOM_std_values)), 
                                                                        male_age_std =rep(0, length(male_HOM_std_values)),
                                                                        female_age_std = rep(0, length(male_HOM_std_values))),
                                                   re.form = NA)
predicted_fledgenest_male_HOM_values <- t(apply(predicted_fledgenest_male_HOM, 2, function(x) c("mean" = mean(x), 
                                                                                                "median" = median(x),
                                                                                                "q10" = quantile(x, p = 0.1),
                                                                                                "q90" = quantile(x, p = 0.9))))

# predicting fledglings per clutch for female HOM
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
female_HOM_std_values <- (HOM_values - mean(nestdata$female_HOM)) / sd(nestdata$female_HOM)
predicted_fledgenest_female_HOM <- posterior_predict(mod_bayes_fledge_nest_HOM,
                                                     newdata = data.frame(female_HOM_std = female_HOM_std_values,
                                                                          male_HOM_std = rep(0, length(female_HOM_std_values)),
                                                                          pair_GS_std = rep(0, length(male_HOM_std_values)), 
                                                                          male_age_std =rep(0, length(male_HOM_std_values)),
                                                                          female_age_std = rep(0, length(male_HOM_std_values))),
                                                     re.form = NA)
predicted_fledgenest_female_HOM_values <- t(apply(predicted_fledgenest_female_HOM, 2, function(x) c("mean" = mean(x),
                                                                                                    "median" = median(x),
                                                                                                    "q10" = quantile(x, p = 0.1),
                                                                                                    "q90" = quantile(x, p = 0.9))))

# predicting hatchlings per clutch for male HOM
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
male_HOM_std_values <- (HOM_values - mean(nestdata$male_HOM)) / sd(nestdata$male_HOM)
predicted_hatchnest_male_HOM <- posterior_predict(mod_bayes_hatch_nest_HOM,
                                                  newdata = data.frame(male_HOM_std = male_HOM_std_values,
                                                                       female_HOM_std = rep(0, length(male_HOM_std_values)),
                                                                       pair_GS_std = rep(0, length(male_HOM_std_values)), 
                                                                       male_age_std =rep(0, length(male_HOM_std_values)),
                                                                       female_age_std = rep(0, length(male_HOM_std_values))),
                                                  re.form = NA)
predicted_hatchnest_male_HOM_values <- t(apply(predicted_hatchnest_male_HOM, 2, function(x) c("mean" = mean(x), 
                                                                                              "median" = median(x),
                                                                                              "q10" = quantile(x, p = 0.1),
                                                                                              "q90" = quantile(x, p = 0.9))))

# predicting hatchlings per clutch for female HOM
HOM_values <- c(0.73, 0.75, 0.77, 0.79, 0.81)
female_HOM_std_values <- (HOM_values - mean(nestdata$female_HOM)) / sd(nestdata$female_HOM)
predicted_hatchnest_female_HOM <- posterior_predict(mod_bayes_hatch_nest_HOM,
                                                    newdata = data.frame(female_HOM_std = female_HOM_std_values,
                                                                         male_HOM_std = rep(0, length(female_HOM_std_values)),
                                                                         pair_GS_std = rep(0, length(female_HOM_std_values)), 
                                                                         male_age_std =rep(0, length(female_HOM_std_values)),
                                                                         female_age_std = rep(0, length(female_HOM_std_values))),
                                                    re.form = NA)
predicted_hatchnest_female_HOM_values <- t(apply(predicted_hatchnest_female_HOM, 2, function(x) c("mean" = mean(x), 
                                                                                                  "median" = median(x),
                                                                                                  "q10" = quantile(x, p = 0.1),
                                                                                                  "q90" = quantile(x, p = 0.9))))

# predicting fledglings per clutch for pair GS
GS_values <- c(0.74, 0.78, 0.83, 0.86, 0.90)
GS_std_values <- (GS_values - mean(nestdata$pair_GS, na.rm=TRUE)) / sd(nestdata$pair_GS, na.rm=TRUE)
predicted_fledgenest_GS <- posterior_predict(mod_bayes_fledge_nest_HOM,
                                             newdata = data.frame(pair_GS_std = GS_std_values,
                                                                  female_HOM_std = rep(0, length(GS_std_values)),
                                                                  male_HOM_std = rep(0, length(GS_std_values)),
                                                                  male_age_std =rep(0, length(GS_std_values)),
                                                                  female_age_std = rep(0, length(GS_std_values))),
                                             re.form = NA)
predicted_fledgenest_GS_values <- t(apply(predicted_fledgenest_GS, 2, function(x) c("mean" = mean(x), 
                                                                                    "q10" = quantile(x, p = 0.1),
                                                                                    "q90" = quantile(x, p = 0.9))))

# predicting hatchlings per clutch for pair GS
GS_values <- c(0.74, 0.78, 0.83, 0.86, 0.90)
GS_std_values <- (GS_values - mean(nestdata$pair_GS, na.rm=TRUE)) / sd(nestdata$pair_GS, na.rm=TRUE)
predicted_hatchnest_GS <- posterior_predict(mod_bayes_hatch_nest_HOM,
                                            newdata = data.frame(pair_GS_std = GS_std_values,
                                                                 female_HOM_std = rep(0, length(GS_std_values)),
                                                                 male_HOM_std = rep(0, length(GS_std_values)),
                                                                 male_age_std =rep(0, length(GS_std_values)),
                                                                 female_age_std = rep(0, length(GS_std_values))),  
                                            re.form = NA)
predicted_hatchnest_GS_values <- t(apply(predicted_hatchnest_GS, 2, function(x) c("mean" = mean(x), 
                                                                                  "q10" = quantile(x, p = 0.1),
                                                                                  "q90" = quantile(x, p = 0.9))))
