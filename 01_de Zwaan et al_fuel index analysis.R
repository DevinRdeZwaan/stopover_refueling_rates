################################################################################
### (1) Fuel index analysis R code for:

### "Mass gain and stopover dynamics among migrating songbirds are linked to 
### seasonal, environmental, and life-history effects"

###  DR de Zwaan, A Huang, Q McCallum, K Owen, M Lamont, and W Easton (2022)

###  Ornithology	

###  Run with R version 4.1.1

################################################################################

###################

### Overview:

### The following code conducts the fuel index analysis and the association
### between mass gain, departure probability, and stopover duration using birds
### that have been captured multiple times within a season.
### This code is specifically for spring migration data, but fall migration
### analysis can be run on the same code by substituting the spring for fall data.

###################

### Set working directory

setwd("")

### Load required packages

library(rstan)
library(plyr)
library(doBy)
library(dplyr)
library(lme4)
library(lmerTest)
library(piecewiseSEM)
library(ggplot2)
library(igraph)
library(ggeffects)
library(RColorBrewer)
library(brms)
library(performance)
library(tidybayes)
library(tidyverse)
library(shinystan)
library(modelr)
library(PNWColors)


### Read in recapture and daily capture data

recap_data <- read.csv("de Zwaan et al_2022_recap_data.csv")
daily_cap_data <- read.csv("de Zwaan et al_2022_daily_cap_data.csv")

### Subset each dataset into spring and fall seasons

recap_spring <- subset(recap_data, season=="spring")
recap_fall <- subset(recap_data, season=="fall")

daily_cap_spring <- subset(daily_cap_data, season=="spring")
daily_cap_fall <- subset(daily_cap_data, season=="fall")


### Subset insectivores and omnivores into their own dataset

insectivore_recap <- subset(recap_spring, guild == "insectivore")
omnivore_recap <- subset(recap_spring, guild == "omnivore")

insectivore_daily <- subset(daily_cap_spring, guild == "insectivore")
omnivore_daily <- subset(daily_cap_spring, guild == "omnivore")


#####################
### Clean data to include likely migrants only
#####################

### Histograms of observed stopover to identify outliers

# Insectivores

hist(insectivore_recap$period) # Very low sample size above 10.

sum(table(insectivore_recap$period)[12:37]) # 54

sum(table(insectivore_recap$period)[12:37])/length(insectivore_recap$period) 
# Only 4.8% of insectivores recaptured stayed beyond 10 days.


# Omnivores
hist(omnivore_recap$period)

sum(table(omnivore_recap$period)[12:39]) # 63

sum(table(omnivore_recap$period)[12:39])/length(omnivore_recap$period) 
# 24.0% of omnivores recaptured stayed beyond 10 days.


### Visualize trend between mass gain and observed stopover length 

(split_plot <- ggplot(aes(period, fuel_index_change), data = insectivore_recap) + 
   geom_point() + 
   facet_wrap(~ species) + # create a facet for each mountain range
   xlab("period") + 
   ylab("prop mass"))

(split_plot <- ggplot(aes(period, fuel_index_change), data = omnivore_recap) + 
    geom_point() + 
    facet_wrap(~ species) + # create a facet for each mountain range
    xlab("period") + 
    ylab("prop mass"))

### Note - Mass gain clearly levels off after 7-10 days.


### Discard individuals that stayed longer than 10 days to avoid including residents.

# Note - use 14 days for fall migration

insectivore_mig <- subset(insectivore_recap, period <=10)
omnivore_mig <- subset(omnivore_recap, period <=10)

insectivore_daily_mig <- subset(insectivore_daily, period <=10)
omnivore_daily_mig <- subset(omnivore_daily, period <=10)


### Re-level factors

insectivore_mig$species <- factor(insectivore_mig$species)
omnivore_mig$species <- factor(omnivore_mig$species)

insectivore_mig$bird_ID <- factor(insectivore_mig$bird_ID)
omnivore_mig$bird_ID <- factor(omnivore_mig$bird_ID)

insectivore_daily_mig$species <- factor(insectivore_daily_mig$species)
omnivore_daily_mig$species <- factor(omnivore_daily_mig$species)

insectivore_daily_mig$bird_ID <- factor(insectivore_daily_mig$bird_ID)
omnivore_daily_mig$bird_ID <- factor(omnivore_daily_mig$bird_ID)


###########################
### Scale predictors 
###########################

### Insectivores

insectivore_mig_sc <- insectivore_mig %>%
  mutate(
    fuel_index = scale(fuel_index, center = TRUE, scale = TRUE),
    fuel_index_change = scale(fuel_index_change, center = TRUE, scale = TRUE),
    period_sc = scale(period, center = TRUE, scale = TRUE),
    first_ord_date = scale(first_ord_date, center = TRUE, scale = TRUE),
    first_time = scale(first_time, center = TRUE, scale = TRUE),
    last_time = scale(last_time, center = TRUE, scale = TRUE),
    
    mean_temp = scale(mean_temp, center = TRUE, scale = TRUE),
    avg_precip = scale(avg_precip, center = TRUE, scale = TRUE),
    EVI_sum_rel = scale(EVI_sum_rel, center = TRUE, scale = TRUE),
    
    guild_density = scale(insectivore_density, center = TRUE, scale = TRUE)
  )


### Omnivores

omnivore_mig_sc <- omnivore_mig %>%
  mutate(
    fuel_index = scale(fuel_index, center = TRUE, scale = TRUE),
    fuel_index_change = scale(fuel_index_change, center = TRUE, scale = TRUE),
    period_sc = scale(period, center = TRUE, scale = TRUE),
    first_ord_date = scale(first_ord_date, center = TRUE, scale = TRUE),
    first_time = scale(first_time, center = TRUE, scale = TRUE),
    last_time = scale(last_time, center = TRUE, scale = TRUE),
    
    mean_temp = scale(mean_temp, center = TRUE, scale = TRUE),
    avg_precip = scale(avg_precip, center = TRUE, scale = TRUE),
    EVI_sum_rel = scale(EVI_sum_rel, center = TRUE, scale = TRUE),
    
    guild_density = scale(omnivore_density, center = TRUE, scale = TRUE)
  )


### Make year into a factor
insectivore_mig_sc$year <- factor(insectivore_mig_sc$year)
omnivore_mig_sc$year <- factor(omnivore_mig_sc$year)


### Remove NAs

insectivore_mig_sc <- insectivore_mig_sc[-which(is.na(insectivore_mig_sc$last_time)),]

omnivore_mig_sc <- omnivore_mig_sc[-which(is.na(omnivore_mig_sc$last_time)),]


##########################
### Check correlations
##########################

### Environmental variables

cor.test(insectivore_mig_sc$mean_temp, insectivore_mig_sc$avg_precip) # -0.44
cor.test(insectivore_mig_sc$mean_temp, insectivore_mig_sc$EVI_sum_rel) # 0.46
cor.test(insectivore_mig_sc$EVI_sum_rel, insectivore_mig_sc$avg_precip) # -0.08


### Environmental variables with date

cor.test(insectivore_mig_sc$mean_temp, insectivore_mig_sc$first_ord_date) # 0.66***
cor.test(insectivore_mig_sc$avg_precip, insectivore_mig_sc$first_ord_date) # -0.03
cor.test(insectivore_mig_sc$EVI_sum_rel, insectivore_mig_sc$first_ord_date) # 0.37


### Calculate residual of temperature and ordinal date

insectivore_mig_sc$temp_res <- summary(lm(mean_temp ~ first_ord_date, data=insectivore_mig_sc))$residuals
omnivore_mig_sc$temp_res <- summary(lm(mean_temp ~ first_ord_date, data=omnivore_mig_sc))$residuals


###################
### (1) Models for proportional change in mass between first and last capture
##################


### Set priors

# For models with random slope (requires lkj correlational prior)
prior1 <- c(
  prior(normal(0, 5), class = Intercept),
  prior(normal(0, 5), class = b), 
  prior(cauchy(0, 2), class = sd), 
  prior(cauchy(0, 2), class = sigma), 
  prior(lkj(1), class = cor) 
)

# For models without random slope
prior2 <- c(
  prior(normal(0, 5), class = Intercept),
  prior(normal(0, 5), class = b), 
  prior(cauchy(0, 2), class = sd), 
  prior(cauchy(0, 2), class = sigma)
)

############
### Fit and run models
############

### A) Random slope

# Insectivore

insectivore_bm_rs <- brm(fuel_index_change ~ fuel_index + period_sc 
                              + first_ord_date + last_time
                              + temp_res
                              + EVI_sum_rel
                              + avg_precip
                              + guild_density
                              + (1 + period_sc|species) 
                              + (1|year), 
                      warmup = 1000,
                      iter = 3000,
                      thin = 1,
                      chains = 4,
                      data=insectivore_mig_sc, family=skew_normal(), 
                      prior = prior1,
                      control = list(adapt_delta = 0.999, max_treedepth = 15),
                      seed = 100,
                      cores = 4)

summary(insectivore_bm_rs)

### Save

#saveRDS(insectivore_bm_rs, "insectivore_bm_rs.rds")


# Omnivore

omnivore_bm_rs <- brm(fuel_index_change ~ fuel_index + period_sc 
                          + first_ord_date + last_time
                          + temp_res
                          + EVI_sum_rel
                          + avg_precip
                          + guild_density
                          + (1 + period_sc|species) 
                          + (1|year), 
                      warmup = 1000,
                      iter = 3000,
                      thin = 1,
                      chains = 4,
                      data= omnivore_mig_sc, family=skew_normal(), 
                      prior = prior1,
                      control = list(adapt_delta = 0.9999, max_treedepth = 15),
                      seed = 100,
                      cores = 4)

summary(omnivore_bm_rs)

### Save

#saveRDS(omnivore_bm_rs, "omnivore_bm_rs.rds")


### B) Sex/age effect

# Insectivores

insectivore_bm_rs_sex <- brm(fuel_index_change ~ fuel_index + period_sc 
                      + first_ord_date + last_time
                      + temp_res
                      + EVI_sum_rel
                      + avg_precip
                      + guild_density
                      + sex
                      + (1 + period_sc|species) 
                      + (1|year), 
                      warmup = 1000,
                      iter = 3000,
                      chains = 4,
                      thin = 1,
                      data=insectivore_mig_sc, family=skew_normal(), 
                      prior = prior1,
                      control = list(adapt_delta = 0.999, max_treedepth = 15),
                      seed = 100,
                      cores = 4)

summary(insectivore_bm_rs_sex)

### Save

#saveRDS(insectivore_bm_rs_sex, "insectivore_bm_rs_sex.rds")


# Note - no sex model for omnivores because too inconsistently measured
# Note - add age for fall models.


### C) Random intercept

# Insectivore

insectivore_bm_ri <- brm(fuel_index_change ~ fuel_index + period_sc 
                      + first_ord_date + last_time
                      + temp_res
                      + EVI_sum_rel
                      + avg_precip
                      + guild_density
                      + (1|species)
                      + (1|year), 
                      warmup = 1000,
                      iter = 3000,
                      thin = 1,
                      chains = 4,
                      data=insectivore_mig_sc, family=skew_normal(), 
                      prior = prior2,
                      control = list(adapt_delta = 0.999, max_treedepth = 15),
                      seed = 100,
                      cores = 4)

summary(insectivore_bm_ri)


### Save

#saveRDS(insectivore_bm_ri, "insectivore_bm_ri.rds")


# Omnivore

omnivore_bm_ri <- brm(fuel_index_change ~ fuel_index + period_sc 
                          + first_ord_date + last_time
                          + temp_res
                          + EVI_sum_rel
                          + avg_precip
                          + guild_density
                          + (1|species)
                          + (1|year), 
                          warmup = 1000,
                          iter = 3000,
                          thin = 1,
                          chains = 4,
                          data=omnivore_mig_sc, family=skew_normal(), 
                          prior = prior2,
                          control = list(adapt_delta = 0.999, max_treedepth = 15),
                          seed = 100,
                          cores = 4)

summary(omnivore_bm_ri)


### Save

#saveRDS(omnivore_bm_ri, "omnivore_bm_ri.rds")


############################
### Model comparisons
############################

### LOO (leave-one-out) comparison

loo(insectivore_bm_rs, insectivore_bm_rs_sex, insectivore_bm_ri)

loo(omnivore_bm_rs, omnivore_bm_ri)


### Model weights

model_weights(insectivore_bm_rs, insectivore_bm_rs_sex, insectivore_bm_ri,
              weights = "loo")

model_weights(omnivore_bm_rs, omnivore_bm_ri,
              weights = "loo")



################################################################################
### Long and short-distance migrant models
################################################################################

### Subset insectivores and omnivores into short and long distance migrants

insectivore_mig_long_sc <- subset(insectivore_mig_sc, species == "WIWA" | species == "YEWA" | species == "OCWA")
omnivore_mig_long_sc <- subset(omnivore_mig_sc, species == "LISP" | species == "GCSP" | species == "WCSP")

insectivore_mig_short_sc <- subset(insectivore_mig_sc, species == "COYE" | species == "YRWA")
omnivore_mig_short_sc <- subset(omnivore_mig_sc, species == "SOSP" | species == "FOSP")

### Combine

mig_long_sc <- bind_rows(insectivore_mig_long_sc, omnivore_mig_long_sc)
mig_short_sc <- bind_rows(insectivore_mig_short_sc, omnivore_mig_short_sc)


### A) Random slope

# Long distance 

long_bm_rs <- brm(fuel_index_change ~ fuel_index + period_sc 
                      + first_ord_date + last_time
                      + temp_res
                      + EVI_sum_rel
                      + avg_precip
                      + guild_density
                      + (1 + period_sc|species) 
                      + (1|year), 
                      warmup = 1000,
                      iter = 3000,
                      thin = 1,
                      chains = 4,
                      data=mig_long_sc, family=skew_normal(), 
                      prior = prior1,
                      control = list(adapt_delta = 0.999, max_treedepth = 15),
                      seed = 100,
                      cores = 4)

summary(long_bm_rs)


# Short distance

short_bm_rs <- brm(fuel_index_change ~ fuel_index + period_sc 
                   + first_ord_date + last_time
                   + temp_res
                   + EVI_sum_rel
                   + avg_precip
                   + guild_density
                   + (1 + period_sc|species) 
                   + (1|year), 
                   warmup = 1000,
                   iter = 3000,
                   thin = 1,
                   chains = 4,
                   data=mig_short_sc, family=skew_normal(), 
                   prior = prior1,
                   control = list(adapt_delta = 0.999, max_treedepth = 15),
                   seed = 100,
                   cores = 4)

summary(short_bm_rs)


### B) Random intercept

# Long distance

long_bm_ri <- brm(fuel_index_change ~ fuel_index + period_sc 
                   + first_ord_date + last_time
                   + temp_res
                   + EVI_sum_rel
                   + avg_precip
                   + guild_density
                   + (1|species)
                   + (1|year), 
                   warmup = 1000,
                   iter = 3000,
                   thin = 1,
                   chains = 4,
                   data=mig_long_sc, family=skew_normal(), 
                   prior = prior2,
                   control = list(adapt_delta = 0.999, max_treedepth = 15),
                   seed = 100,
                   cores = 4)

summary(long_bm_ri)


# Short distance

short_bm_ri <- brm(fuel_index_change ~ fuel_index + period_sc 
                   + first_ord_date + last_time
                   + temp_res
                   + EVI_sum_rel
                   + avg_precip
                   + guild_density
                   + (1|species)
                   + (1|year), 
                   warmup = 1000,
                   iter = 3000,
                   thin = 1,
                   chains = 4,
                   data=mig_short_sc, family=skew_normal(), 
                   prior = prior2,
                   control = list(adapt_delta = 0.999, max_treedepth = 15),
                   seed = 100,
                   cores = 4)

summary(short_bm_ri)



### Save

#saveRDS(long_bm_rs, "long_bm_rs.rds")
#saveRDS(short_bm_rs, "short_bm_rs.rds")

#saveRDS(long_bm_ri, "long_bm_ri.rds")
#saveRDS(short_bm_ri, "short_bm_ri.rds")


### LOO

loo(long_bm_rs, long_bm_ri)
loo(short_bm_rs, short_bm_ri)

### Model weights

model_weights(long_bm_rs, long_bm_ri, weights = "loo")

model_weights(short_bm_rs, short_bm_ri, weights = "loo")




################################################################################
### Figures
################################################################################

########
### Figure 1 and 2 - individual species slopes
########

# Insectivores
data_interim_w <- insectivore_bm_rs_sex %>%
                  spread_draws(b_period_sc, r_species[species,])
data_interim_w$species <- factor(data_interim_w$species)

# Omnivores
data_interim_s <- omnivore_bm_rs %>%
                  spread_draws(b_period_sc, r_species[species,])
data_interim_s$species <- factor(data_interim_s$species)


### Generate plot of random draw distribution

# Insectivores

insectivore_slopes <- insectivore_bm_rs_sex %>%
                      spread_draws(b_period_sc, r_species[species,]) %>%
                      mutate(species_mean = b_period_sc + r_species) %>%
                      ggplot(aes(y = species, x = species_mean, fill = stat(x > 0))) +
                      geom_vline(xintercept = 0, linetype = "dashed") +
                      stat_halfeye(.width = c(.93, .83)) +
                      scale_fill_manual(values = c("#d2848d", "#7bbcd5")) +
                      scale_y_discrete(limits = rev(levels(data_interim_w$species))) +
                      scale_x_continuous(limits=c(-2,3), breaks=c(-2,-1,0,1,2,3)) +
                      labs(x = "Estimates", y = "Species") +
                      theme_classic() +
                      theme(legend.position = "none",
                            axis.text=element_text(size=18, colour="black"),
                            axis.title=element_text(size=22, colour="black"))

insectivore_slopes

### Save
png(file="file_name.png",width=3000,height=2500, res=600)
insectivore_slopes
dev.off()


# Omnivores

omnivore_slopes <- omnivore_bm_rs %>%
                   spread_draws(b_period_sc, r_species[species,]) %>%
                   mutate(species_mean = b_period_sc + r_species) %>%
                   ggplot(aes(y = species, x = species_mean, fill = stat(x > 0))) +
                   geom_vline(xintercept = 0, linetype = "dashed") +
                   stat_halfeye(.width = c(.93, .83)) +
                   scale_fill_manual(values = c("#d2848d", "#7bbcd5")) +
                   scale_y_discrete(limits = rev(levels(data_interim_s$species))) +
                   scale_x_continuous(limits=c(-2,3), breaks=c(-2,-1,0,1,2,3)) +
                   labs(x = "Estimates", y = "Species") +
                   theme_classic() +
                   theme(legend.position = "none",
                         axis.text=element_text(size=18, colour="black"),
                         axis.title=element_text(size=22, colour="black"))

omnivore_slopes

### Save
png(file="file_name.png",width=3000,height=2500, res=600)
omnivore_slopes
dev.off()


########
### Plot individual species slopes
########

### Predict and attach to original data

# Insectivores
insectivore_mig_plot <- subset(insectivore_mig_sc, is.na(fuel_index_change) == FALSE &
                                                   is.na(fuel_index) == FALSE)

insectivore_mig_plot$predict <- predict(insectivore_bm_rs_sex)

insectivore_mig_plot$predict <- insectivore_mig_plot$predict[,1]

# Omnivores

omnivore_mig_plot <- subset(omnivore_mig_sc, is.na(fuel_index) == FALSE &
                                             is.na(fuel_index_change) == FALSE &
                                             is.na(last_time) == FALSE)

omnivore_mig_plot$predict <- predict(omnivore_bm_rs)

omnivore_mig_plot$predict <- omnivore_mig_plot$predict[,1]


### Determine the SD to reconvert back to real values

insectivore_mig_subset <- subset(insectivore_mig, is.na(fuel_index_change) == FALSE)
sd(insectivore_mig_subset$fuel_index_change) # 0.07858354

omnivore_mig_subset <- subset(omnivore_mig, is.na(fuel_index_change) == FALSE)
sd(omnivore_mig_subset$fuel_index_change) # 0.07036332


### Convert prop mass change to percentage (and multiply by standard deviation)

### Reconvert predictions

insectivore_mig_plot$fuel_index_change_use <- insectivore_mig_plot$fuel_index_change * 0.07858354 * 100
insectivore_mig_plot$predict_prop <- insectivore_mig_plot$predict * 0.07858354 * 100

omnivore_mig_plot$fuel_index_change_use <- omnivore_mig_plot$fuel_index_change * 0.07036332 * 100
omnivore_mig_plot$predict_prop <- omnivore_mig_plot$predict * 0.07036332 * 100


### Set colour ramp
colour_ramp <- pnw_palette("Sailboat",5,type="discrete") 

# Insectivores

insectivore_stacked_plot <- ggplot(insectivore_mig_plot, aes(x = period, y = fuel_index_change_use, colour = species)) +
  geom_jitter(width=0.2,height=0.3,alpha = 0.5) +
  theme_classic() +
  geom_smooth(aes(x = period, y = predict_prop, colour=species), method = "lm", se=TRUE, 
              size =1, alpha=0.3) +
  scale_colour_manual(values= colour_ramp) +
  scale_y_continuous(limits=c(-20,30), breaks=c(-20,-10,0,10,20,30)) +
  scale_x_continuous(limits=c(-0.5,10.5), breaks=c(0,2,4,6,8,10)) +
  labs(x = "Stopover length (days)", y = "Delta Fuel index (%)") + 
  theme(legend.position = "none",
        axis.text=element_text(size=18, colour="black"),
        axis.title=element_text(size=22, colour="black"))

insectivore_stacked_plot

### Export
png(file="file_name.png",width=3000,height=2500, res=600)
insectivore_stacked_plot
dev.off()


### Set colour ramp
colour_ramp2 <- pnw_palette("Bay",5,type="discrete") 

# Omnivore

omnivore_stacked_plot <- ggplot(omnivore_mig_plot, aes(x = period, y = fuel_index_change_use, colour = species)) +
  geom_jitter(width=0.2,height=0.3,alpha = 0.5) +
  theme_classic() +
  geom_smooth(aes(x = period, y = predict_prop, colour=species), method = "lm", se=TRUE, 
              size =1, alpha=0.3) +
  scale_colour_manual(values= colour_ramp2) +
  scale_y_continuous(limits=c(-20,20), breaks=c(-20,-15,-10,-5,0,5,10,15,20)) +
  scale_x_continuous(limits=c(-0.5,10.5), breaks=c(0,2,4,6,8,10)) +
  labs(x = "Stopover length (days)", y = "Delta Fuel index (%)") + 
  theme(legend.position = "none",
        axis.text=element_text(size=18, colour="black"),
        axis.title=element_text(size=22, colour="black"))

omnivore_stacked_plot

### Export
png(file="file_name.png",width=3000,height=2500, res=600)
omnivore_stacked_plot
dev.off()



################################################################################
### Models for probability of departure
################################################################################

### Note - for this analysis, use the daily capture dataset

#######################
### Scale variables
#######################

# Insectivores

insectivore_daily_mig_sc <- insectivore_daily_mig %>%
  mutate(
    ord_date = scale(ord_date, center = TRUE, scale = TRUE),
    time = scale(time_h, center = TRUE, scale = TRUE),
    fuel_index = scale(fuel_index, center = TRUE, scale = TRUE),
    avg_temp = scale(avg_temp, center = TRUE, scale = TRUE),
    precip_sum = scale(precip_sum, center = TRUE, scale = TRUE),
    period = scale(period, center = TRUE, scale = TRUE),
)


# Omnivores

omnivore_daily_mig_sc <- omnivore_daily_mig %>%
  mutate(
    ord_date = scale(ord_date, center = TRUE, scale = TRUE),
    time = scale(time_h, center = TRUE, scale = TRUE),
    fuel_index = scale(fuel_index, center = TRUE, scale = TRUE),
    avg_temp = scale(avg_temp, center = TRUE, scale = TRUE),
    precip_sum = scale(precip_sum, center = TRUE, scale = TRUE),
    period = scale(period, center = TRUE, scale = TRUE)
  )


### Make year into a factor
insectivore_daily_mig_sc$year <- factor(insectivore_daily_mig_sc$year)
omnivore_daily_mig_sc$year <- factor(omnivore_daily_mig_sc$year)

### Calculate temperature residuals with date

insectivore_daily_mig_sc$temp_res <- summary(lm(avg_temp ~ ord_date, data=insectivore_daily_mig_sc))$residuals
omnivore_daily_mig_sc$temp_res <- summary(lm(avg_temp ~ ord_date, data= omnivore_daily_mig_sc))$residuals


#############
### Set priors for logistic model
#############

prior1 <- c(
  prior(normal(0, 5), class = Intercept),
  prior(normal(0, 5), class = b),
  prior(cauchy(0, 5), class = sd),
  prior(lkj(1), class = cor)
)

prior2 <- c(
  prior(normal(0, 5), class = Intercept),
  prior(normal(0, 5), class = b),
  prior(cauchy(0, 5), class = sd)
)


##############
### Fit and run models
##############

# Insectivores

insectivore_depart_rs <- brm(depart ~ time + ord_date
                           + fuel_index
                           + (1 + fuel_index|species) 
                           + (1|bird_ID) 
                           + (1|year),
                      data=insectivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                      prior = prior1,
                      control = list(adapt_delta = 0.999, max_treedepth=15),
                      seed = 100,
                      warmup = 1000,
                      iter = 3000,
                      thin = 1,
                      chains = 4,
                      cores = 4)

summary(insectivore_depart_rs)


### Save

#saveRDS(insectivore_depart_rs, "insectivore_depart_rs.rds")


insectivore_depart_rs_weather <- brm(depart ~ time + ord_date
                                     + fuel_index
                                     + temp_res
                                     + precip_sum
                                     + (1 + fuel_index|species) 
                                     + (1|bird_ID) 
                                     + (1|year),
                                     data=insectivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                                     prior = prior1,
                                     control = list(adapt_delta = 0.999, max_treedepth=15),
                                     seed = 100,
                                     warmup = 1000,
                                     iter = 3000,
                                     thin = 1,
                                     chains = 4,
                                     cores = 4)

summary(insectivore_depart_rs_weather)

### Save

#saveRDS(insectivore_depart_rs_weather, "insectivore_depart_rs_weather.rds")


insectivore_ri <- brm(depart ~ time + ord_date
                        + fuel_index
                        + (1|species) 
                        + (1|bird_ID) 
                        + (1|year),
                        data=insectivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                        prior = prior2,
                        control = list(adapt_delta = 0.999, max_treedepth=15),
                        seed = 100,
                        warmup = 1000,
                        iter = 3000,
                        thin = 1,
                        chains = 4,
                        cores = 4)

summary(insectivore_ri)


### Save

#saveRDS(insectivore_ri, "insectivore_ri.rds")


insectivore_depart_ri_weather <- brm(depart ~ time + ord_date
                                   + fuel_index
                                   + temp_res
                                   + precip_sum
                                   + (1|species) 
                                   + (1|bird_ID) 
                                   + (1|year),
                                   data=insectivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                                   prior = prior2,
                                   control = list(adapt_delta = 0.999, max_treedepth=15),
                                   seed = 100,
                                   warmup = 1000,
                                   iter = 3000,
                                   thin = 1,
                                   chains = 4,
                                   cores = 4)

summary(insectivore_depart_ri_weather)

### Save

#saveRDS(insectivore_depart_ri_weather, "insectivore_depart_ri_weather.rds")


# Omnivores

omnivore_depart_rs <- brm(depart ~ time + ord_date
                        + fuel_index
                        + (1 + fuel_index|species) 
                        + (1|bird_ID) 
                        + (1|year),
                        data=omnivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                        prior = prior1,
                        control = list(adapt_delta = 0.999, max_treedepth=15),
                        seed = 100,
                        warmup = 1000,
                        iter = 3000,
                        thin = 1,
                        chains = 4,
                        cores = 4)

summary(omnivore_depart_rs)

### Save

#saveRDS(omnivore_depart_rs, "omnivore_depart_rs.rds")


omnivore_depart_rs_weather <- brm(depart ~ time + ord_date
                                     + fuel_index
                                     + temp_res
                                     + precip_sum
                                     + (1 + fuel_index|species) 
                                     + (1|bird_ID) 
                                     + (1|year),
                                     data=omnivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                                     prior = prior1,
                                     control = list(adapt_delta = 0.999, max_treedepth=15),
                                     seed = 100,
                                     warmup = 1000,
                                     iter = 3000,
                                     thin = 1,
                                     chains = 4,
                                     cores = 4)

summary(omnivore_depart_rs_weather)

### Save

#saveRDS(omnivore_depart_rs_weather, "omnivore_depart_rs_weather.rds")


omnivore_depart_ri <- brm(depart ~ time + ord_date
                        + fuel_index
                        + (1|species) 
                        + (1|bird_ID) 
                        + (1|year),
                        data=omnivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                        prior = prior2,
                        control = list(adapt_delta = 0.999, max_treedepth=15),
                        seed = 100,
                        warmup = 1000,
                        iter = 3000,
                        thin = 1,
                        chains = 4,
                        cores = 4)

summary(omnivore_depart_ri)

### Save
#saveRDS(omnivore_depart_ri, "omnivore_depart_ri.rds")


omnivore_depart_ri_weather <- brm(depart ~ time + ord_date
                                   + fuel_index
                                   + temp_res
                                   + precip_sum
                                   + (1|species) 
                                   + (1|bird_ID) 
                                   + (1|year),
                                   data=omnivore_daily_mig_sc, family = bernoulli(link = "logit"), 
                                   prior = prior2,
                                   control = list(adapt_delta = 0.999, max_treedepth=15),
                                   seed = 100,
                                   warmup = 1000,
                                   iter = 3000,
                                   thin = 1,
                                   chains = 4,
                                   cores = 4)

summary(omnivore_depart_ri_weather)

### Save

#saveRDS(omnivore_depart_ri_weather, "omnivore_depart_ri_weather.rds")


####################
### Model comparison
####################

### LOO
loo(insectivore_depart_rs, insectivore_depart_ri, insectivore_depart_rs_weather,
    insectivore_depart_ri_weather)

loo(omnivore_depart_rs, omnivore_depart_ri, omnivore_depart_rs_weather,
    omnivore_depart_ri_weather) 

### Model weights

model_weights(insectivore_depart_rs, insectivore_depart_ri, insectivore_depart_rs_weather,
              insectivore_depart_ri_weather,
              weights = "loo")

model_weights(omnivore_depart_rs, omnivore_depart_ri, omnivore_depart_rs_weather,
              omnivore_depart_ri_weather,
              weights = "loo")


################################################################################
### Fit departure models for long and short-distance migrants
#################################################################################

### Subset datasets into short and long distance migrants

insectivore_daily_mig_long_sc <- subset(insectivore_daily_mig_sc, species == "WIWA" | species == "YEWA" | species == "OCWA")
omnivore_daily_mig_long_sc <- subset(omnivore_daily_mig_sc, species == "LISP" | species == "GCSP" | species == "WCSP")

insectivore_daily_mig_short_sc <- subset(insectivore_daily_mig_sc, species == "COYE" | species == "YRWA")
omnivore_daily_mig_short_sc <- subset(omnivore_daily_mig_sc, species == "SOSP" | species == "FOSP")

### Combine
mig_daily_long_sc <- bind_rows(insectivore_daily_mig_long_sc, omnivore_daily_mig_long_sc)
mig_daily_short_sc <- bind_rows(insectivore_daily_mig_short_sc, omnivore_daily_mig_short_sc)

######
### Fit models
######

# Long distance

long_depart_rs <- brm(depart ~ time + ord_date
                        + fuel_index
                        + (1 + fuel_index|species) 
                        + (1|bird_ID) 
                        + (1|year),
                        data = mig_daily_long_sc, family = bernoulli(link = "logit"), 
                        prior = prior1,
                        control = list(adapt_delta = 0.999, max_treedepth=15),
                        seed = 100,
                        warmup = 1000,
                        iter = 3000,
                        thin = 1,
                        chains = 4,
                        cores = 4)

long_depart_rs_weather <- brm(depart ~ time + ord_date
                        + fuel_index
                        + temp_res
                        + precip_sum
                        + (1 + fuel_index|species) 
                        + (1|bird_ID) 
                        + (1|year),
                        data=mig_daily_long_sc, family = bernoulli(link = "logit"), 
                        prior = prior1,
                        control = list(adapt_delta = 0.999, max_treedepth=15),
                        seed = 100,
                        warmup = 1000,
                        iter = 3000,
                        thin = 1,
                        chains = 4,
                        cores = 4)

### Model comparison

loo(long_depart_rs,long_depart_rs_weather)

model_weights(long_depart_rs,long_depart_rs_weather, weights = "loo")


# Short distance

short_depart_rs <- brm(depart ~ time + ord_date
                        + fuel_index
                        + (1 + fuel_index|species) 
                        + (1|bird_ID) 
                        + (1|year),
                        data=mig_daily_short_sc, family = bernoulli(link = "logit"), 
                        prior = prior1,
                        control = list(adapt_delta = 0.999, max_treedepth=15),
                        seed = 100,
                        warmup = 1000,
                        iter = 3000,
                        thin = 1,
                        chains = 4,
                        cores = 4)

short_depart_rs_weather <- brm(depart ~ time + ord_date
                     + fuel_index
                     + temp_res
                     + precip_sum
                     + (1 + fuel_index|species) 
                     + (1|bird_ID) 
                     + (1|year),
                     data=mig_daily_short_sc, family = bernoulli(link = "logit"), 
                     prior = prior1,
                     control = list(adapt_delta = 0.999, max_treedepth=15),
                     seed = 100,
                     warmup = 1000,
                     iter = 3000,
                     thin = 1,
                     chains = 4,
                     cores = 4)


### Model comparison

loo(short_depart_rs, short_depart_rs_weather)

model_weights(short_depart_rs, short_depart_rs_weather, weights = "loo")


################################################################################
### Figures
################################################################################

#########
### Plot individual species relationships for departure probability
#########

### Remove NAs
insectivore_daily_mig_plot <- subset(insectivore_daily_mig_sc, is.na(fuel_index) == FALSE &
                                   is.na(time) == FALSE)

omnivore_daily_mig_plot <- subset(omnivore_daily_mig_sc, is.na(fuel_index) == FALSE &
                                   is.na(time) == FALSE)

### Predict
predict_in <- predict(insectivore_depart_rs)
insectivore_daily_mig_plot$predict <- predict_in[,1]

predict_omn <- predict(omnivore_depart_ri)
omnivore_daily_mig_plot$predict <- predict_omn[,1]

### Convert predictions to 0 and 1s
insectivore_daily_mig_plot$predict_new <- ifelse(insectivore_daily_mig_plot$predict < 0.5, 0, 1)

omnivore_daily_mig_plot$predict_new <- ifelse(omnivore_daily_mig_plot$predict < 0.5, 0, 1)


### Determine the SD to reconvert
insectivore_daily_mig_subset <- subset(insectivore_daily_mig, is.na(fuel_index) == FALSE)
sd(insectivore_daily_mig_subset$fuel_index) # 0.0898801

omnivore_daily_mig_subset <- subset(omnivore_daily_mig, is.na(fuel_index) == FALSE)
sd(omnivore_daily_mig_subset$fuel_index) # 0.09954909


### Convert prop mass to percentage (and multiply by standard deviation)

insectivore_daily_mig_plot$fuel_index_use <- insectivore_daily_mig_plot$fuel_index * 0.0898801 * 100
omnivore_daily_mig_plot$fuel_index_use <- omnivore_daily_mig_plot$fuel_index * 0.09954909 * 100


### Plot 

### Set colour ramp
colour_ramp <- pnw_palette("Sailboat",5,type="discrete") 

# Insectivores

insectivore_plot <- ggplot(insectivore_daily_mig_plot, aes(x = fuel_index_use, y = depart, colour=species)) +
  geom_jitter(width=0.0,height=0.0,alpha = 0.5) +
  theme_classic() +
  geom_smooth(aes(x = fuel_index_use, y = predict_new, colour=species), method = "glm", method.args = list(family = "binomial"), 
              se=TRUE, size =1, alpha = 0.2) +
  scale_colour_manual(values= colour_ramp) +
  labs(x = "Relative fuel index (%)", y = "Probability of departure", title = "") + 
  theme(legend.position = "none",
        axis.text=element_text(size=18, colour="black"),
        axis.title=element_text(size=22, colour="black"))

insectivore_plot


### Export
png(file="file_name.png",width=3000,height=2500, res=600)
insectivore_plot
dev.off()


### Plot omnivores

### Set colour ramp
colour_ramp2 <- pnw_palette("Bay",5,type="discrete") 

omnivore_plot <- ggplot(omnivore_daily_mig_plot, aes(x = fuel_index_use, y = depart, colour=species)) +
  geom_jitter(width=0.0,height=0.0,alpha = 0.5) +
  theme_classic() +
  geom_smooth(aes(x = fuel_index_use, y = predict_new, colour=species), method = "glm", method.args = list(family = "binomial"), 
              se=TRUE, size =1, alpha = 0.2) +
  scale_colour_manual(values= colour_ramp2) +
  labs(x = "Relative fuel index (%)", y = "Probability of departure", title = "") + 
  theme(legend.position = "none",
        axis.text=element_text(size=18, colour="black"),
        axis.title=element_text(size=22, colour="black"))

omnivore_plot


### Export
png(file="file_name.png",width=3000,height=2500, res=600)
omnivore_plot
dev.off()



################### End of script ##############################################