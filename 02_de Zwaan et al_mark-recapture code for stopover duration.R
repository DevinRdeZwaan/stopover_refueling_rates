################################################################################
### (2) Stopover duration analysis R code for:

### "Mass gain and stopover dynamics among migrating songbirds are linked to 
### seasonal, environmental, and life-history effects"

###  DR de Zwaan, A Huang, Q McCallum, K Owen, M Lamont, and W Easton (2022)

###  Ornithology	

###  Run with R version 4.1.1

################################################################################

###################

### Overview:

### The following code estimates the 'true' stopover duration of all birds
### using capture-recapture data and assesses how fuel index at first capture
### (a proxy for arrival condition) is associated with length of stay.

### This code also demonstrates calculation of LBM and the fuel index.

### Note - this code demonstrates the work flow for the spring migration analysis
### only. To conduct the fall migration analysis, simply substitute the spring 
### data with the fall data

###################


### Load R2admb and set up file path packages

library(R2admb)

prepare_admb=function()
{
  Sys.setenv(PATH = paste("C:/ADMB/bin;C:/ADMB/utilities;C:/MinGW/bin;",
                          Sys.getenv("PATH"), sep = ";"))
  Sys.setenv(ADMB_HOME = "C:/ADMB")
  invisible()
}

prepare_admb()

### Load remaining required packages

library(marked)
library(ggplot2)
library(reshape2)
library(plyr)
library(dplyr)
library(tidyverse)


### Set working directory

setwd("")


################################################################################
### Data pre-processing
################################################################################

### Read in bird data

bird_data <- read.csv("de Zwaan et al_2022_mark_recap_data.csv")

### Subset seasons

spring_data <- subset(bird_data, season == "spring")
fall_data <- subset(bird_data, season == "fall")

### Subset by guild

insectivores <- subset(spring_data, guild == "insectivore") 

omnivores <- subset(spring_data, guild == "omnivore") 

### Convert species to factor

insectivores$species <- factor(insectivores$species)

omnivores$species <- factor(omnivores$species)

### Inspect sample size
table(insectivores$species)
table(omnivores$species)

################################################################################
### Create a capture event variable & retain birds captured within 10 days
### (see methods and previous code)
################################################################################

### Create a column of captures with 1

insectivores$capture <- rep(1, length(insectivores$bird_ID))
omnivores$capture <- rep(1, length(omnivores$bird_ID))

### Group by bird_ID and use cumsum to create consecutive capture events

insectivores <- insectivores %>%
                group_by(bird_ID,year) %>%
                mutate(capture_event = cumsum(capture))

omnivores <- omnivores %>%
             group_by(bird_ID,year) %>%
             mutate(capture_event = cumsum(capture))

### Calculate days since first capture

insectivores <- insectivores %>%
                group_by(bird_ID, year) %>%
                mutate(days_since_first_cap = max(ord_date) - min(ord_date))

omnivores <- omnivores %>%
             group_by(bird_ID, year) %>%
             mutate(days_since_first_cap = max(ord_date) - min(ord_date))

### Visualize

hist(insectivores$days_since_first_cap)
hist(omnivores$days_since_first_cap)

### Restrict to 10 days or less which is when the majority of recaps occur

insectivores_subset <- subset(insectivores, days_since_first_cap < 11)
omnivores_subset <- subset(omnivores, days_since_first_cap < 11)

# Note - Cut-off is an attempt to differentiate between migrants and residents.
# See methods for more detail.


################################################################################
### Create lbm and fuel index variables
################################################################################

### Select individuals captured with zero fat

insectivores_fat_0 <- subset(insectivores_subset, fat == 0)
omnivores_fat_0 <- subset(omnivores_subset, fat == 0)

### Subset individual species

WIWA <- subset(insectivores_fat_0, species == "WIWA")
OCWA <- subset(insectivores_fat_0, species == "OCWA")
YRWA <- subset(insectivores_fat_0, species == "YRWA")
YEWA <- subset(insectivores_fat_0, species == "YEWA")
COYE <- subset(insectivores_fat_0, species == "COYE")

GCSP <- subset(omnivores_fat_0, species == "GCSP")
WCSP <- subset(omnivores_fat_0, species == "WCSP")
LISP <- subset(omnivores_fat_0, species == "LISP")
SOSP <- subset(omnivores_fat_0, species == "SOSP")
FOSP <- subset(omnivores_fat_0, species == "FOSP")

### Fit regressions of mass and wing length

WIWA_lm <- lm(mass ~ wing, data=WIWA)
OCWA_lm <- lm(mass ~ wing, data=OCWA)
YRWA_lm <- lm(mass ~ wing, data=YRWA)
YEWA_lm <- lm(mass ~ wing, data=YEWA)
COYE_lm <- lm(mass ~ wing, data=COYE)

GCSP_lm <- lm(mass ~ wing, data=GCSP)
WCSP_lm <- lm(mass ~ wing, data=WCSP)
LISP_lm <- lm(mass ~ wing, data=LISP)
SOSP_lm <- lm(mass ~ wing, data=SOSP)
FOSP_lm <- lm(mass ~ wing, data=FOSP)


### Using relationship, calculate LBM for full dataset

insectivores_subset$lbm <- ifelse(insectivores_subset$species == "WIWA", 
                                     WIWA_lm$coefficients[1] + 
                                     WIWA_lm$coefficients[2]*insectivores_subset$wing,
                              ifelse(insectivores_subset$species == "OCWA", 
                                     OCWA_lm$coefficients[1] + 
                                     OCWA_lm$coefficients[2]*insectivores_subset$wing,
                              ifelse(insectivores_subset$species == "YRWA", 
                                     YRWA_lm$coefficients[1] + 
                                     YRWA_lm$coefficients[2]*insectivores_subset$wing,
                              ifelse(insectivores_subset$species == "YEWA", 
                                     YEWA_lm$coefficients[1] + 
                                     YEWA_lm$coefficients[2]*insectivores_subset$wing,
                              COYE_lm$coefficients[1] + 
                              COYE_lm$coefficients[2]*insectivores_subset$wing))))


omnivores_subset$lbm <- ifelse(omnivores_subset$species == "GCSP", 
                                     GCSP_lm$coefficients[1] + 
                                     GCSP_lm$coefficients[2]*omnivores_subset$wing,
                              ifelse(omnivores_subset$species == "WCSP", 
                                     WCSP_lm$coefficients[1] + 
                                     WCSP_lm$coefficients[2]*omnivores_subset$wing,
                              ifelse(omnivores_subset$species == "LISP", 
                                     LISP_lm$coefficients[1] + 
                                     LISP_lm$coefficients[2]*omnivores_subset$wing,
                              ifelse(omnivores_subset$species == "SOSP", 
                                     SOSP_lm$coefficients[1] + 
                                     SOSP_lm$coefficients[2]*omnivores_subset$wing,
                              FOSP_lm$coefficients[1] + 
                              FOSP_lm$coefficients[2]*omnivores_subset$wing))))


### Calculate fuel index (mass as a proportion of LBM - see Methods)

insectivores_subset$fuel_index <- (insectivores_subset$mass - insectivores_subset$lbm)/insectivores_subset$lbm
omnivores_subset$fuel_index <- (omnivores_subset$mass - omnivores_subset$lbm)/omnivores_subset$lbm

### Rearrange by date

insectivores_sorted <- arrange(insectivores_subset, year, bird_ID, ord_date)
omnivores_sorted <- arrange(omnivores_subset, year, bird_ID, ord_date)

### Set ordinal date to start at 1.

insectivores_sorted$ord_date <- insectivores_sorted$ord_date - 105 # 233 for fall
omnivores_sorted$ord_date <- omnivores_sorted$ord_date - 105 # 233 for fall

### Identify all birds that were captured before cut-off date

insectivore_index <- insectivores_sorted[which(insectivores_sorted$ord_date < 1),5]
omnivore_index <- omnivores_sorted[which(omnivores_sorted$ord_date < 1),5]

# Note - removed first few days that only occurred in original year

### Remove all rows that correspond to the identified birds

insectivores_use <- filter(insectivores_sorted, !(bird_ID %in% insectivore_index$bird_ID))

omnivores_use <- filter(omnivores_sorted, !(bird_ID %in% omnivore_index$bird_ID))

### Select only first capture for additional data to be added in later

insectivores_first_capture <- subset(insectivores_use, capture_event == 1)
omnivores_first_capture <- subset(omnivores_use, capture_event == 1)

### Split first capture data into separate years

# Note - 2011 not included because sample size is too small

insectivores_1st_cap_2010 <- subset(insectivores_first_capture, year == 2010)
insectivores_1st_cap_2012 <- subset(insectivores_first_capture, year == 2012)
insectivores_1st_cap_2013 <- subset(insectivores_first_capture, year == 2013)
insectivores_1st_cap_2014 <- subset(insectivores_first_capture, year == 2014)
insectivores_1st_cap_2015 <- subset(insectivores_first_capture, year == 2015)
insectivores_1st_cap_2016 <- subset(insectivores_first_capture, year == 2016)
insectivores_1st_cap_2017 <- subset(insectivores_first_capture, year == 2017)
insectivores_1st_cap_2018 <- subset(insectivores_first_capture, year == 2018)
insectivores_1st_cap_2019 <- subset(insectivores_first_capture, year == 2019)


omnivores_1st_cap_2010 <- subset(omnivores_first_capture, year == 2010)
omnivores_1st_cap_2012 <- subset(omnivores_first_capture, year == 2012)
omnivores_1st_cap_2013 <- subset(omnivores_first_capture, year == 2013)
omnivores_1st_cap_2014 <- subset(omnivores_first_capture, year == 2014)
omnivores_1st_cap_2015 <- subset(omnivores_first_capture, year == 2015)
omnivores_1st_cap_2016 <- subset(omnivores_first_capture, year == 2016)
omnivores_1st_cap_2017 <- subset(omnivores_first_capture, year == 2017)
omnivores_1st_cap_2018 <- subset(omnivores_first_capture, year == 2018)
omnivores_1st_cap_2019 <- subset(omnivores_first_capture, year == 2019)

### From full data, extract only band, year, j_date and capture for now

insectivores_reduced <- insectivores_use %>% select(bird_ID,year,ord_date,capture)
omnivores_reduced <- omnivores_use %>% select(bird_ID,year,ord_date,capture)

### Split up by year as well 

insectivores_2010 <- subset(insectivores_reduced, year == 2010)
insectivores_2012 <- subset(insectivores_reduced, year == 2012)
insectivores_2013 <- subset(insectivores_reduced, year == 2013)
insectivores_2014 <- subset(insectivores_reduced, year == 2014)
insectivores_2015 <- subset(insectivores_reduced, year == 2015)
insectivores_2016 <- subset(insectivores_reduced, year == 2016)
insectivores_2017 <- subset(insectivores_reduced, year == 2017)
insectivores_2018 <- subset(insectivores_reduced, year == 2018)
insectivores_2019 <- subset(insectivores_reduced, year == 2019)

omnivores_2010 <- subset(omnivores_reduced, year == 2010)
omnivores_2012 <- subset(omnivores_reduced, year == 2012)
omnivores_2013 <- subset(omnivores_reduced, year == 2013)
omnivores_2014 <- subset(omnivores_reduced, year == 2014)
omnivores_2015 <- subset(omnivores_reduced, year == 2015)
omnivores_2016 <- subset(omnivores_reduced, year == 2016)
omnivores_2017 <- subset(omnivores_reduced, year == 2017)
omnivores_2018 <- subset(omnivores_reduced, year == 2018)
omnivores_2019 <- subset(omnivores_reduced, year == 2019)


################################################################################
### Create matrix of capture events (0&1) and replace NAs with zero
################################################################################

### Pivot tables

# Insectivores

interim_2010 <- melt(insectivores_2010, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2010 <- acast(interim_2010, bird_ID ~ ord_date)
insectivores_matrix_2010[is.na(insectivores_matrix_2010)] <- 0

interim_2012 <- melt(insectivores_2012, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2012 <- acast(interim_2012, bird_ID ~ ord_date)
insectivores_matrix_2012[is.na(insectivores_matrix_2012)] <- 0

interim_2013 <- melt(insectivores_2013, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2013 <- acast(interim_2013, bird_ID ~ ord_date)
insectivores_matrix_2013[is.na(insectivores_matrix_2013)] <- 0

interim_2014 <- melt(insectivores_2014, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2014 <- acast(interim_2014, bird_ID ~ ord_date)
insectivores_matrix_2014[is.na(insectivores_matrix_2014)] <- 0

interim_2015 <- melt(insectivores_2015, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2015 <- acast(interim_2015, bird_ID ~ ord_date)
insectivores_matrix_2015[is.na(insectivores_matrix_2015)] <- 0

interim_2016 <- melt(insectivores_2016, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2016 <- acast(interim_2016, bird_ID ~ ord_date)
insectivores_matrix_2016[is.na(insectivores_matrix_2016)] <- 0

interim_2017 <- melt(insectivores_2017, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2017 <- acast(interim_2017, bird_ID ~ ord_date)
insectivores_matrix_2017[is.na(insectivores_matrix_2017)] <- 0

interim_2018 <- melt(insectivores_2018, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2018 <- acast(interim_2018, bird_ID ~ ord_date)
insectivores_matrix_2018[is.na(insectivores_matrix_2018)] <- 0

interim_2019 <- melt(insectivores_2019, id.var=c("bird_ID","ord_date"), measure.var= "capture")
insectivores_matrix_2019 <- acast(interim_2019, bird_ID ~ ord_date)
insectivores_matrix_2019[is.na(insectivores_matrix_2019)] <- 0


# Omnivores

omn_interim_2010 <- melt(omnivores_2010, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2010 <- acast(omn_interim_2010, bird_ID ~ ord_date)
omnivores_matrix_2010[is.na(omnivores_matrix_2010)] <- 0

omn_interim_2012 <- melt(omnivores_2012, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2012 <- acast(omn_interim_2012, bird_ID ~ ord_date)
omnivores_matrix_2012[is.na(omnivores_matrix_2012)] <- 0

omn_interim_2013 <- melt(omnivores_2013, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2013 <- acast(omn_interim_2013, bird_ID ~ ord_date)
omnivores_matrix_2013[is.na(omnivores_matrix_2013)] <- 0

omn_interim_2014 <- melt(omnivores_2014, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2014 <- acast(omn_interim_2014, bird_ID ~ ord_date)
omnivores_matrix_2014[is.na(omnivores_matrix_2014)] <- 0

omn_interim_2015 <- melt(omnivores_2015, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2015 <- acast(omn_interim_2015, bird_ID ~ ord_date)
omnivores_matrix_2015[is.na(omnivores_matrix_2015)] <- 0

omn_interim_2016 <- melt(omnivores_2016, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2016 <- acast(omn_interim_2016, bird_ID ~ ord_date)
omnivores_matrix_2016[is.na(omnivores_matrix_2016)] <- 0

omn_interim_2017 <- melt(omnivores_2017, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2017 <- acast(omn_interim_2017, bird_ID ~ ord_date)
omnivores_matrix_2017[is.na(omnivores_matrix_2017)] <- 0

omn_interim_2018 <- melt(omnivores_2018, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2018 <- acast(omn_interim_2018, bird_ID ~ ord_date)
omnivores_matrix_2018[is.na(omnivores_matrix_2018)] <- 0

omn_interim_2019 <- melt(omnivores_2019, id.var=c("bird_ID","ord_date"), measure.var= "capture")
omnivores_matrix_2019 <- acast(omn_interim_2019, bird_ID ~ ord_date)
omnivores_matrix_2019[is.na(omnivores_matrix_2019)] <- 0


### For insectivores, remove last column of 2019 because it is the only one 
### that goes past 46 days (i.e., end on 46)

insectivores_matrix_2019 <- insectivores_matrix_2019[,-45]


### Fill in non-capture days (NAs) with zeros from 1-46

# Insectivores

# Create a vector of capture days
capture_day_vec <- seq(1,46,1)

# 2010
vec_2010 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2010))
insectivores_2010_new <- matrix(0,nrow=dim(insectivores_matrix_2010)[1],ncol=46)  
insectivores_2010_new[,-vec_2010] <- insectivores_matrix_2010
insectivores_2010_new2 <- insectivores_2010_new[which(rowSums(insectivores_2010_new) != 0),]

# 2012
vec_2012 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2012))
insectivores_2012_new <- matrix(0,nrow=dim(insectivores_matrix_2012)[1],ncol=46)  
insectivores_2012_new[,-vec_2012] <- insectivores_matrix_2012 
insectivores_2012_new2 <- insectivores_2012_new[which(rowSums(insectivores_2012_new) != 0),]

# 2013
vec_2013 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2013))
insectivores_2013_new <- matrix(0,nrow=dim(insectivores_matrix_2013)[1],ncol=46)  
insectivores_2013_new[,-vec_2013] <- insectivores_matrix_2013 
insectivores_2013_new2 <- insectivores_2013_new[which(rowSums(insectivores_2013_new) != 0),]

# 2014
vec_2014 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2014))
insectivores_2014_new <- matrix(0,nrow=dim(insectivores_matrix_2014)[1],ncol=46)  
insectivores_2014_new[,-vec_2014] <- insectivores_matrix_2014 
insectivores_2014_new2 <- insectivores_2014_new[which(rowSums(insectivores_2014_new) != 0),]

# 2015
vec_2015 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2015))
insectivores_2015_new <- matrix(0,nrow=dim(insectivores_matrix_2015)[1],ncol=46)  
insectivores_2015_new[,-vec_2015] <- insectivores_matrix_2015 
insectivores_2015_new2 <- insectivores_2015_new[which(rowSums(insectivores_2015_new) != 0),]

# 2016
vec_2016 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2016))
insectivores_2016_new <- matrix(0,nrow=dim(insectivores_matrix_2016)[1],ncol=46)  
insectivores_2016_new[,-vec_2016] <- insectivores_matrix_2016 
insectivores_2016_new2 <- insectivores_2016_new[which(rowSums(insectivores_2016_new) != 0),]

# 2017
vec_2017 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2017))
insectivores_2017_new <- matrix(0,nrow=dim(insectivores_matrix_2017)[1],ncol=46)  
insectivores_2017_new[,-vec_2017] <- insectivores_matrix_2017 
insectivores_2017_new2 <- insectivores_2017_new[which(rowSums(insectivores_2017_new) != 0),]

# 2018
vec_2018 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2018))
insectivores_2018_new <- matrix(0,nrow=dim(insectivores_matrix_2018)[1],ncol=46)  
insectivores_2018_new[,-vec_2018] <- insectivores_matrix_2018 
insectivores_2018_new2 <- insectivores_2018_new[which(rowSums(insectivores_2018_new) != 0),]

# 2019
vec_2019 <- setdiff(capture_day_vec, colnames(insectivores_matrix_2019))
insectivores_2019_new <- matrix(0,nrow=dim(insectivores_matrix_2019)[1],ncol=46)  
insectivores_2019_new[,-vec_2019] <- insectivores_matrix_2019 
insectivores_2019_new2 <- insectivores_2019_new[which(rowSums(insectivores_2019_new) != 0),]


# Omnivores

# 2010
vec_2010 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2010))
omnivores_2010_new <- matrix(0,nrow=dim(omnivores_matrix_2010)[1],ncol=46)  
omnivores_2010_new[,-vec_2010] <- omnivores_matrix_2010
omnivores_2010_new2 <- omnivores_2010_new[which(rowSums(omnivores_2010_new) != 0),]

# 2012
vec_2012 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2012))
omnivores_2012_new <- matrix(0,nrow=dim(omnivores_matrix_2012)[1],ncol=46)  
omnivores_2012_new[,-vec_2012] <- omnivores_matrix_2012 
omnivores_2012_new2 <- omnivores_2012_new[which(rowSums(omnivores_2012_new) != 0),]

# 2013
vec_2013 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2013))
omnivores_2013_new <- matrix(0,nrow=dim(omnivores_matrix_2013)[1],ncol=46)  
omnivores_2013_new[,-vec_2013] <- omnivores_matrix_2013 
omnivores_2013_new2 <- omnivores_2013_new[which(rowSums(omnivores_2013_new) != 0),]

# 2014
vec_2014 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2014))
omnivores_2014_new <- matrix(0,nrow=dim(omnivores_matrix_2014)[1],ncol=46)  
omnivores_2014_new[,-vec_2014] <- omnivores_matrix_2014 
omnivores_2014_new2 <- omnivores_2014_new[which(rowSums(omnivores_2014_new) != 0),]

# 2015
vec_2015 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2015))
omnivores_2015_new <- matrix(0,nrow=dim(omnivores_matrix_2015)[1],ncol=46)  
omnivores_2015_new[,-vec_2015] <- omnivores_matrix_2015 
omnivores_2015_new2 <- omnivores_2015_new[which(rowSums(omnivores_2015_new) != 0),]

# 2016
vec_2016 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2016))
omnivores_2016_new <- matrix(0,nrow=dim(omnivores_matrix_2016)[1],ncol=46)  
omnivores_2016_new[,-vec_2016] <- omnivores_matrix_2016 
omnivores_2016_new2 <- omnivores_2016_new[which(rowSums(omnivores_2016_new) != 0),]

# 2017
vec_2017 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2017))
omnivores_2017_new <- matrix(0,nrow=dim(omnivores_matrix_2017)[1],ncol=46)  
omnivores_2017_new[,-vec_2017] <- omnivores_matrix_2017 
omnivores_2017_new2 <- omnivores_2017_new[which(rowSums(omnivores_2017_new) != 0),]

# 2018
vec_2018 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2018))
omnivores_2018_new <- matrix(0,nrow=dim(omnivores_matrix_2018)[1],ncol=46)  
omnivores_2018_new[,-vec_2018] <- omnivores_matrix_2018 
omnivores_2018_new2 <- omnivores_2018_new[which(rowSums(omnivores_2018_new) != 0),]

# 2019
vec_2019 <- setdiff(capture_day_vec, colnames(omnivores_matrix_2019))
omnivores_2019_new <- matrix(0,nrow=dim(omnivores_matrix_2019)[1],ncol=46)  
omnivores_2019_new[,-vec_2019] <- omnivores_matrix_2019 
omnivores_2019_new2 <- omnivores_2019_new[which(rowSums(omnivores_2019_new) != 0),]


### Collapse all capture histories into one column

insectivores_2010_collapsed <- apply(insectivores_2010_new2,1,paste,collapse="")
insectivores_2012_collapsed <- apply(insectivores_2012_new2,1,paste,collapse="")
insectivores_2013_collapsed <- apply(insectivores_2013_new2,1,paste,collapse="")
insectivores_2014_collapsed <- apply(insectivores_2014_new2,1,paste,collapse="")
insectivores_2015_collapsed <- apply(insectivores_2015_new2,1,paste,collapse="")
insectivores_2016_collapsed <- apply(insectivores_2016_new2,1,paste,collapse="")
insectivores_2017_collapsed <- apply(insectivores_2017_new2,1,paste,collapse="")
insectivores_2018_collapsed <- apply(insectivores_2018_new2,1,paste,collapse="")
insectivores_2019_collapsed <- apply(insectivores_2019_new2,1,paste,collapse="")

omnivores_2010_collapsed <- apply(omnivores_2010_new2,1,paste,collapse="")
omnivores_2012_collapsed <- apply(omnivores_2012_new2,1,paste,collapse="")
omnivores_2013_collapsed <- apply(omnivores_2013_new2,1,paste,collapse="")
omnivores_2014_collapsed <- apply(omnivores_2014_new2,1,paste,collapse="")
omnivores_2015_collapsed <- apply(omnivores_2015_new2,1,paste,collapse="")
omnivores_2016_collapsed <- apply(omnivores_2016_new2,1,paste,collapse="")
omnivores_2017_collapsed <- apply(omnivores_2017_new2,1,paste,collapse="")
omnivores_2018_collapsed <- apply(omnivores_2018_new2,1,paste,collapse="")
omnivores_2019_collapsed <- apply(omnivores_2019_new2,1,paste,collapse="")


### Convert to dataframe, assign year, bird_ID and all other information

# Insectivores

#2010
insectivores_2010_df <- as.data.frame(insectivores_2010_collapsed)
colnames(insectivores_2010_df) <- "ch"
insectivores_2010_df$year <- rep(2010, length(insectivores_2010_df$ch))
insectivores_1st_cap_2010_new <- insectivores_1st_cap_2010[which(rowSums(insectivores_2010_new) != 0),]
### This is necessary to make the dataframes line up (same length)
insectivores_2010_df$bird_ID <- insectivores_1st_cap_2010_new$bird_ID
insectivores_2010_df$age <- insectivores_1st_cap_2010_new$age
insectivores_2010_df$sex <- insectivores_1st_cap_2010_new$sex
insectivores_2010_df$fuel_index <- insectivores_1st_cap_2010_new$fuel_index
insectivores_2010_df$species <- insectivores_1st_cap_2010_new$species
insectivores_2010_df$ord_date <- insectivores_1st_cap_2010_new$ord_date

#2012
insectivores_2012_df <- as.data.frame(insectivores_2012_collapsed)
colnames(insectivores_2012_df) <- "ch"
insectivores_2012_df$year <- rep(2012, length(insectivores_2012_df$ch))
insectivores_1st_cap_2012_new <- insectivores_1st_cap_2012[which(rowSums(insectivores_2012_new) != 0),]
insectivores_2012_df$bird_ID <- insectivores_1st_cap_2012_new$bird_ID
insectivores_2012_df$age <- insectivores_1st_cap_2012_new$age
insectivores_2012_df$sex <- insectivores_1st_cap_2012_new$sex
insectivores_2012_df$fuel_index <- insectivores_1st_cap_2012_new$fuel_index
insectivores_2012_df$species <- insectivores_1st_cap_2012_new$species
insectivores_2012_df$ord_date <- insectivores_1st_cap_2012_new$ord_date

#2013
insectivores_2013_df <- as.data.frame(insectivores_2013_collapsed)
colnames(insectivores_2013_df) <- "ch"
insectivores_2013_df$year <- rep(2013, length(insectivores_2013_df$ch))
insectivores_1st_cap_2013_new <- insectivores_1st_cap_2013[which(rowSums(insectivores_2013_new) != 0),]
insectivores_2013_df$bird_ID <- insectivores_1st_cap_2013_new$bird_ID
insectivores_2013_df$age <- insectivores_1st_cap_2013_new$age
insectivores_2013_df$sex <- insectivores_1st_cap_2013_new$sex
insectivores_2013_df$fuel_index <- insectivores_1st_cap_2013_new$fuel_index
insectivores_2013_df$species <- insectivores_1st_cap_2013_new$species
insectivores_2013_df$ord_date <- insectivores_1st_cap_2013_new$ord_date

#2014
insectivores_2014_df <- as.data.frame(insectivores_2014_collapsed)
colnames(insectivores_2014_df) <- "ch"
insectivores_2014_df$year <- rep(2014, length(insectivores_2014_df$ch))
insectivores_1st_cap_2014_new <- insectivores_1st_cap_2014[which(rowSums(insectivores_2014_new) != 0),]
insectivores_2014_df$bird_ID <- insectivores_1st_cap_2014_new$bird_ID
insectivores_2014_df$age <- insectivores_1st_cap_2014_new$age
insectivores_2014_df$sex <- insectivores_1st_cap_2014_new$sex
insectivores_2014_df$fuel_index <- insectivores_1st_cap_2014_new$fuel_index
insectivores_2014_df$species <- insectivores_1st_cap_2014_new$species
insectivores_2014_df$ord_date <- insectivores_1st_cap_2014_new$ord_date

#2015
insectivores_2015_df <- as.data.frame(insectivores_2015_collapsed)
colnames(insectivores_2015_df) <- "ch"
insectivores_2015_df$year <- rep(2015, length(insectivores_2015_df$ch))
insectivores_1st_cap_2015_new <- insectivores_1st_cap_2015[which(rowSums(insectivores_2015_new) != 0),]
insectivores_2015_df$bird_ID <- insectivores_1st_cap_2015_new$bird_ID
insectivores_2015_df$age <- insectivores_1st_cap_2015_new$age
insectivores_2015_df$sex <- insectivores_1st_cap_2015_new$sex
insectivores_2015_df$fuel_index <- insectivores_1st_cap_2015_new$fuel_index
insectivores_2015_df$species <- insectivores_1st_cap_2015_new$species
insectivores_2015_df$ord_date <- insectivores_1st_cap_2015_new$ord_date

#2016
insectivores_2016_df <- as.data.frame(insectivores_2016_collapsed)
colnames(insectivores_2016_df) <- "ch"
insectivores_2016_df$year <- rep(2016, length(insectivores_2016_df$ch))
insectivores_1st_cap_2016_new <- insectivores_1st_cap_2016[which(rowSums(insectivores_2016_new) != 0),]
insectivores_2016_df$bird_ID <- insectivores_1st_cap_2016_new$bird_ID
insectivores_2016_df$age <- insectivores_1st_cap_2016_new$age
insectivores_2016_df$sex <- insectivores_1st_cap_2016_new$sex
insectivores_2016_df$fuel_index <- insectivores_1st_cap_2016_new$fuel_index
insectivores_2016_df$species <- insectivores_1st_cap_2016_new$species
insectivores_2016_df$ord_date <- insectivores_1st_cap_2016_new$ord_date

#2017
insectivores_2017_df <- as.data.frame(insectivores_2017_collapsed)
colnames(insectivores_2017_df) <- "ch"
insectivores_2017_df$year <- rep(2017, length(insectivores_2017_df$ch))
insectivores_1st_cap_2017_new <- insectivores_1st_cap_2017[which(rowSums(insectivores_2017_new) != 0),]
insectivores_2017_df$bird_ID <- insectivores_1st_cap_2017_new$bird_ID
insectivores_2017_df$age <- insectivores_1st_cap_2017_new$age
insectivores_2017_df$sex <- insectivores_1st_cap_2017_new$sex
insectivores_2017_df$fuel_index <- insectivores_1st_cap_2017_new$fuel_index
insectivores_2017_df$species <- insectivores_1st_cap_2017_new$species
insectivores_2017_df$ord_date <- insectivores_1st_cap_2017_new$ord_date

#2018
insectivores_2018_df <- as.data.frame(insectivores_2018_collapsed)
colnames(insectivores_2018_df) <- "ch"
insectivores_2018_df$year <- rep(2018, length(insectivores_2018_df$ch))
insectivores_1st_cap_2018_new <- insectivores_1st_cap_2018[which(rowSums(insectivores_2018_new) != 0),]
insectivores_2018_df$bird_ID <- insectivores_1st_cap_2018_new$bird_ID
insectivores_2018_df$age <- insectivores_1st_cap_2018_new$age
insectivores_2018_df$sex <- insectivores_1st_cap_2018_new$sex
insectivores_2018_df$fuel_index <- insectivores_1st_cap_2018_new$fuel_index
insectivores_2018_df$species <- insectivores_1st_cap_2018_new$species
insectivores_2018_df$ord_date <- insectivores_1st_cap_2018_new$ord_date

#2019
insectivores_2019_df <- as.data.frame(insectivores_2019_collapsed)
colnames(insectivores_2019_df) <- "ch"
insectivores_2019_df$year <- rep(2019, length(insectivores_2019_df$ch))
insectivores_1st_cap_2019_new <- insectivores_1st_cap_2019[which(rowSums(insectivores_2019_new) != 0),]
insectivores_2019_df$bird_ID <- insectivores_1st_cap_2019_new$bird_ID
insectivores_2019_df$age <- insectivores_1st_cap_2019_new$age
insectivores_2019_df$sex <- insectivores_1st_cap_2019_new$sex
insectivores_2019_df$fuel_index <- insectivores_1st_cap_2019_new$fuel_index
insectivores_2019_df$species <- insectivores_1st_cap_2019_new$species
insectivores_2019_df$ord_date <- insectivores_1st_cap_2019_new$ord_date


# Omnivores

#2010
omnivores_2010_df <- as.data.frame(omnivores_2010_collapsed)
colnames(omnivores_2010_df) <- "ch"
omnivores_2010_df$year <- rep(2010, length(omnivores_2010_df$ch))
omnivores_2010_df$bird_ID <- omnivores_1st_cap_2010$bird_ID
omnivores_2010_df$age <- omnivores_1st_cap_2010$age
omnivores_2010_df$sex <- omnivores_1st_cap_2010$sex
omnivores_2010_df$fuel_index <- omnivores_1st_cap_2010$fuel_index
omnivores_2010_df$species <- omnivores_1st_cap_2010$species
omnivores_2010_df$ord_date <- omnivores_1st_cap_2010$ord_date

#2012
omnivores_2012_df <- as.data.frame(omnivores_2012_collapsed)
colnames(omnivores_2012_df) <- "ch"
omnivores_2012_df$year <- rep(2012, length(omnivores_2012_df$ch))
omnivores_2012_df$bird_ID <- omnivores_1st_cap_2012$bird_ID
omnivores_2012_df$age <- omnivores_1st_cap_2012$age
omnivores_2012_df$sex <- omnivores_1st_cap_2012$sex
omnivores_2012_df$fuel_index <- omnivores_1st_cap_2012$fuel_index
omnivores_2012_df$species <- omnivores_1st_cap_2012$species
omnivores_2012_df$ord_date <- omnivores_1st_cap_2012$ord_date

#2013
omnivores_2013_df <- as.data.frame(omnivores_2013_collapsed)
colnames(omnivores_2013_df) <- "ch"
omnivores_2013_df$year <- rep(2013, length(omnivores_2013_df$ch))
omnivores_2013_df$bird_ID <- omnivores_1st_cap_2013$bird_ID
omnivores_2013_df$age <- omnivores_1st_cap_2013$age
omnivores_2013_df$sex <- omnivores_1st_cap_2013$sex
omnivores_2013_df$fuel_index <- omnivores_1st_cap_2013$fuel_index
omnivores_2013_df$species <- omnivores_1st_cap_2013$species
omnivores_2013_df$ord_date <- omnivores_1st_cap_2013$ord_date

#2014
omnivores_2014_df <- as.data.frame(omnivores_2014_collapsed)
colnames(omnivores_2014_df) <- "ch"
omnivores_2014_df$year <- rep(2014, length(omnivores_2014_df$ch))
omnivores_2014_df$bird_ID <- omnivores_1st_cap_2014$bird_ID
omnivores_2014_df$age <- omnivores_1st_cap_2014$age
omnivores_2014_df$sex <- omnivores_1st_cap_2014$sex
omnivores_2014_df$fuel_index <- omnivores_1st_cap_2014$fuel_index
omnivores_2014_df$species <- omnivores_1st_cap_2014$species
omnivores_2014_df$ord_date <- omnivores_1st_cap_2014$ord_date

#2015
omnivores_2015_df <- as.data.frame(omnivores_2015_collapsed)
colnames(omnivores_2015_df) <- "ch"
omnivores_2015_df$year <- rep(2015, length(omnivores_2015_df$ch))
omnivores_2015_df$bird_ID <- omnivores_1st_cap_2015$bird_ID
omnivores_2015_df$age <- omnivores_1st_cap_2015$age
omnivores_2015_df$sex <- omnivores_1st_cap_2015$sex
omnivores_2015_df$fuel_index <- omnivores_1st_cap_2015$fuel_index
omnivores_2015_df$species <- omnivores_1st_cap_2015$species
omnivores_2015_df$ord_date <- omnivores_1st_cap_2015$ord_date

#2016
omnivores_2016_df <- as.data.frame(omnivores_2016_collapsed)
colnames(omnivores_2016_df) <- "ch"
omnivores_2016_df$year <- rep(2016, length(omnivores_2016_df$ch))
omnivores_2016_df$bird_ID <- omnivores_1st_cap_2016$bird_ID
omnivores_2016_df$age <- omnivores_1st_cap_2016$age
omnivores_2016_df$sex <- omnivores_1st_cap_2016$sex
omnivores_2016_df$fuel_index <- omnivores_1st_cap_2016$fuel_index
omnivores_2016_df$species <- omnivores_1st_cap_2016$species
omnivores_2016_df$ord_date <- omnivores_1st_cap_2016$ord_date

#2017
omnivores_2017_df <- as.data.frame(omnivores_2017_collapsed)
colnames(omnivores_2017_df) <- "ch"
omnivores_2017_df$year <- rep(2017, length(omnivores_2017_df$ch))
omnivores_2017_df$bird_ID <- omnivores_1st_cap_2017$bird_ID
omnivores_2017_df$age <- omnivores_1st_cap_2017$age
omnivores_2017_df$sex <- omnivores_1st_cap_2017$sex
omnivores_2017_df$fuel_index <- omnivores_1st_cap_2017$fuel_index
omnivores_2017_df$species <- omnivores_1st_cap_2017$species
omnivores_2017_df$ord_date <- omnivores_1st_cap_2017$ord_date

#2018
omnivores_2018_df <- as.data.frame(omnivores_2018_collapsed)
colnames(omnivores_2018_df) <- "ch"
omnivores_2018_df$year <- rep(2018, length(omnivores_2018_df$ch))
omnivores_2018_df$bird_ID <- omnivores_1st_cap_2018$bird_ID
omnivores_2018_df$age <- omnivores_1st_cap_2018$age
omnivores_2018_df$sex <- omnivores_1st_cap_2018$sex
omnivores_2018_df$fuel_index <- omnivores_1st_cap_2018$fuel_index
omnivores_2018_df$species <- omnivores_1st_cap_2018$species
omnivores_2018_df$ord_date <- omnivores_1st_cap_2018$ord_date

#2019
omnivores_2019_df <- as.data.frame(omnivores_2019_collapsed)
colnames(omnivores_2019_df) <- "ch"
omnivores_2019_df$year <- rep(2019, length(omnivores_2019_df$ch))
omnivores_2019_df$bird_ID <- omnivores_1st_cap_2019$bird_ID
omnivores_2019_df$age <- omnivores_1st_cap_2019$age
omnivores_2019_df$sex <- omnivores_1st_cap_2019$sex
omnivores_2019_df$fuel_index <- omnivores_1st_cap_2019$fuel_index
omnivores_2019_df$species <- omnivores_1st_cap_2019$species
omnivores_2019_df$ord_date <- omnivores_1st_cap_2019$ord_date


#### Combine all years together

insectivores_marked_complete <- bind_rows(insectivores_2010_df,
                                          insectivores_2012_df,
                                          insectivores_2013_df,
                                          insectivores_2014_df,
                                          insectivores_2015_df,
                                          insectivores_2016_df,
                                          insectivores_2017_df,
                                          insectivores_2018_df,
                                          insectivores_2019_df)


omnivores_marked_complete <- bind_rows(omnivores_2010_df,
                                       omnivores_2012_df,
                                       omnivores_2013_df,
                                       omnivores_2014_df,
                                       omnivores_2015_df,
                                       omnivores_2016_df,
                                       omnivores_2017_df,
                                       omnivores_2018_df,
                                       omnivores_2019_df)

################################################################################
### CJS survival models
################################################################################

# Note - the example below is for insectivores during spring migration only.
# For fall migration or for omnivores, use the same code but substitute the 
# required datasets


### Subset data to known sex categories (insectivores only)

insectivores_marked_use <- subset(insectivores_marked_complete, sex == "M" |
                                      sex == "F")

### Make sex a factor

insectivores_marked_use$sex <- factor(insectivores_marked_use$sex)

### Remove NAs in fuel index to allow models to fit properly

insectivores_marked_use <- subset(insectivores_marked_use, is.na(fuel_index) == FALSE)

### Convert year into a factor for random effect

insectivores_marked_use$year_f <- factor(insectivores_marked_use$year)

### Remove age and original year variables or model will not fit.

insectivores_marked_use <- insectivores_marked_use %>% select(-c(age,year))


########################
### Fit models
########################

### Designate Phi (survival) variables
design_Phi <- list(static = c("fuel_index", "sex", "species","ord_date", "year_f"))

### Designate p (detection probability) variables
design_p <- list(static=c("species", "fuel_index", "ord_date", "year_f"))

### Set formulas for full model to estimate beta coefficients
Phi_sam_r <- list(formula = ~ fuel_index + ord_date + species + (1|year_f))
p_sas_r <- list(formula = ~ species + (1|year_f))

### Set formulas for reduced model to estimated stopover duration 
### (i.e., species-specific averages across years)

Phi_sam <- list(formula = ~ species)
p_sas <- list(formula = ~ species)


### Full model - influence of arrival fuel index

insectivores_cjs_r <- crm(insectivores_marked_use,
                          begin.time = 1,
                          model = "probitCJS",
                          model.parameters = list(Phi = Phi_sam_r, p = p_sas_r),
                          design_parameters = list(Phi = design_Phi, p = design_p),
                          groups = "year_f",
                          use.admb = TRUE,
                          burnin = 1000,iter=3000,
                          seed = 100)
                         
insectivores_cjs_r


### Reduced model - stopover duration estimates

insectivores_cjs <- crm(insectivores_marked_use,
                        begin.time = 1,
                        model = "probitCJS",
                        model.parameters = list(Phi = Phi_sam, p = p_sas),
                        design_parameters = list(Phi = design_Phi, p = design),
                        use.admb = TRUE,
                        burnin = 100,iter=3000,
                        seed = 100)


### Extract real values and convert into stopover duration

insectivores_cjs$results$reals

### Back-calculate for estimates of length

# Mean duration
-1 / log(0.7373115) ## 3.28 days COYE
-1 / log(0.6172053) ## 2.07 days YRWA
-1 / log(0.5226010) ## 1.54 days OCWA
-1 / log(0.4882654) ## 1.39 days WIWA
-1 / log(0.5271789) ## 1.56 days YEWA

# Standard deviation
-1 / log(0.02968552) ## 0.28
-1 / log(0.02545338) ## 0.27
-1 / log(0.02533678) ## 0.27
-1 / log(0.01331329) ## 0.23
-1 / log(0.06326715) ## 0.36

# Lower CI
-1 / log(0.6619451) ## 2.42
-1 / log(0.5708838) ## 1.78
-1 / log(0.4683624) ## 1.32
-1 / log(0.4520903) ## 1.26
-1 / log(0.4492777) ## 1.25

# Upper CI
-1 / log(0.7817117) ## 4.06
-1 / log(0.6659627) ## 2.46
-1 / log(0.5660832) ## 1.76
-1 / log(0.5030081) ## 1.46
-1 / log(0.6943090) ## 2.74


################################################################################
### Reverse-time models
################################################################################

#####################
### Data pre-processing
#####################

### Reverse order of columns in original matrix

insectivores_2010_rev <- insectivores_2010_new2[,ncol(insectivores_2010_new2):1]
insectivores_2012_rev <- insectivores_2012_new2[,ncol(insectivores_2012_new2):1]
insectivores_2013_rev <- insectivores_2013_new2[,ncol(insectivores_2013_new2):1]
insectivores_2014_rev <- insectivores_2014_new2[,ncol(insectivores_2014_new2):1]
insectivores_2015_rev <- insectivores_2015_new2[,ncol(insectivores_2015_new2):1]
insectivores_2016_rev <- insectivores_2016_new2[,ncol(insectivores_2016_new2):1]
insectivores_2017_rev <- insectivores_2017_new2[,ncol(insectivores_2017_new2):1]
insectivores_2018_rev <- insectivores_2018_new2[,ncol(insectivores_2018_new2):1]
insectivores_2019_rev <- insectivores_2019_new2[,ncol(insectivores_2019_new2):1]

omnivores_2010_rev <- omnivores_2010_new2[,ncol(omnivores_2010_new2):1]
omnivores_2012_rev <- omnivores_2012_new2[,ncol(omnivores_2012_new2):1]
omnivores_2013_rev <- omnivores_2013_new2[,ncol(omnivores_2013_new2):1]
omnivores_2014_rev <- omnivores_2014_new2[,ncol(omnivores_2014_new2):1]
omnivores_2015_rev <- omnivores_2015_new2[,ncol(omnivores_2015_new2):1]
omnivores_2016_rev <- omnivores_2016_new2[,ncol(omnivores_2016_new2):1]
omnivores_2017_rev <- omnivores_2017_new2[,ncol(omnivores_2017_new2):1]
omnivores_2018_rev <- omnivores_2018_new2[,ncol(omnivores_2018_new2):1]
omnivores_2019_rev <- omnivores_2019_new2[,ncol(omnivores_2019_new2):1]


### Collapse all matrices so reversed capture history in one column

insectivores_2010_collapsed_rev <- apply(insectivores_2010_rev,1,paste,collapse="")
insectivores_2012_collapsed_rev <- apply(insectivores_2012_rev,1,paste,collapse="")
insectivores_2013_collapsed_rev <- apply(insectivores_2013_rev,1,paste,collapse="")
insectivores_2014_collapsed_rev <- apply(insectivores_2014_rev,1,paste,collapse="")
insectivores_2015_collapsed_rev <- apply(insectivores_2015_rev,1,paste,collapse="")
insectivores_2016_collapsed_rev <- apply(insectivores_2016_rev,1,paste,collapse="")
insectivores_2017_collapsed_rev <- apply(insectivores_2017_rev,1,paste,collapse="")
insectivores_2018_collapsed_rev <- apply(insectivores_2018_rev,1,paste,collapse="")
insectivores_2019_collapsed_rev <- apply(insectivores_2019_rev,1,paste,collapse="")

omnivores_2010_collapsed_rev <- apply(omnivores_2010_rev,1,paste,collapse="")
omnivores_2012_collapsed_rev <- apply(omnivores_2012_rev,1,paste,collapse="")
omnivores_2013_collapsed_rev <- apply(omnivores_2013_rev,1,paste,collapse="")
omnivores_2014_collapsed_rev <- apply(omnivores_2014_rev,1,paste,collapse="")
omnivores_2015_collapsed_rev <- apply(omnivores_2015_rev,1,paste,collapse="")
omnivores_2016_collapsed_rev <- apply(omnivores_2016_rev,1,paste,collapse="")
omnivores_2017_collapsed_rev <- apply(omnivores_2017_rev,1,paste,collapse="")
omnivores_2018_collapsed_rev <- apply(omnivores_2018_rev,1,paste,collapse="")
omnivores_2019_collapsed_rev <- apply(omnivores_2019_rev,1,paste,collapse="")


### Convert to data frame, assign year, bird_ID, and all other information

# Insectivores

#2010
insectivores_2010_r_df <- as.data.frame(insectivores_2010_collapsed_rev)
colnames(insectivores_2010_r_df) <- "ch"
insectivores_2010_r_df <- bind_cols(insectivores_2010_r_df, insectivores_2010_df[,2:8])

#2012
insectivores_2012_r_df <- as.data.frame(insectivores_2012_collapsed_rev)
colnames(insectivores_2012_r_df) <- "ch"
insectivores_2012_r_df <- bind_cols(insectivores_2012_r_df, insectivores_2012_df[,2:8])

#2013
insectivores_2013_r_df <- as.data.frame(insectivores_2013_collapsed_rev)
colnames(insectivores_2013_r_df) <- "ch"
insectivores_2013_r_df <- bind_cols(insectivores_2013_r_df, insectivores_2013_df[,2:8])

#2014
insectivores_2014_r_df <- as.data.frame(insectivores_2014_collapsed_rev)
colnames(insectivores_2014_r_df) <- "ch"
insectivores_2014_r_df <- bind_cols(insectivores_2014_r_df, insectivores_2014_df[,2:8])

#2015
insectivores_2015_r_df <- as.data.frame(insectivores_2015_collapsed_rev)
colnames(insectivores_2015_r_df) <- "ch"
insectivores_2015_r_df <- bind_cols(insectivores_2015_r_df, insectivores_2015_df[,2:8])

#2016
insectivores_2016_r_df <- as.data.frame(insectivores_2016_collapsed_rev)
colnames(insectivores_2016_r_df) <- "ch"
insectivores_2016_r_df <- bind_cols(insectivores_2016_r_df, insectivores_2016_df[,2:8])

#2017
insectivores_2017_r_df <- as.data.frame(insectivores_2017_collapsed_rev)
colnames(insectivores_2017_r_df) <- "ch"
insectivores_2017_r_df <- bind_cols(insectivores_2017_r_df, insectivores_2017_df[,2:8])

#2018
insectivores_2018_r_df <- as.data.frame(insectivores_2018_collapsed_rev)
colnames(insectivores_2018_r_df) <- "ch"
insectivores_2018_r_df <- bind_cols(insectivores_2018_r_df, insectivores_2018_df[,2:8])

#2019
insectivores_2019_r_df <- as.data.frame(insectivores_2019_collapsed_rev)
colnames(insectivores_2019_r_df) <- "ch"
insectivores_2019_r_df <- bind_cols(insectivores_2019_r_df, insectivores_2019_df[,2:8])


# Omnivores

#2010
omnivores_2010_r_df <- as.data.frame(omnivores_2010_collapsed_rev)
colnames(omnivores_2010_r_df) <- "ch"
omnivores_2010_r_df <- bind_cols(omnivores_2010_r_df, omnivores_2010_df[,2:8])

#2012
omnivores_2012_r_df <- as.data.frame(omnivores_2012_collapsed_rev)
colnames(omnivores_2012_r_df) <- "ch"
omnivores_2012_r_df <- bind_cols(omnivores_2012_r_df, omnivores_2012_df[,2:8])

#2013
omnivores_2013_r_df <- as.data.frame(omnivores_2013_collapsed_rev)
colnames(omnivores_2013_r_df) <- "ch"
omnivores_2013_r_df <- bind_cols(omnivores_2013_r_df, omnivores_2013_df[,2:8])

#2014
omnivores_2014_r_df <- as.data.frame(omnivores_2014_collapsed_rev)
colnames(omnivores_2014_r_df) <- "ch"
omnivores_2014_r_df <- bind_cols(omnivores_2014_r_df, omnivores_2014_df[,2:8])

#2015
omnivores_2015_r_df <- as.data.frame(omnivores_2015_collapsed_rev)
colnames(omnivores_2015_r_df) <- "ch"
omnivores_2015_r_df <- bind_cols(omnivores_2015_r_df, omnivores_2015_df[,2:8])

#2016
omnivores_2016_r_df <- as.data.frame(omnivores_2016_collapsed_rev)
colnames(omnivores_2016_r_df) <- "ch"
omnivores_2016_r_df <- bind_cols(omnivores_2016_r_df, omnivores_2016_df[,2:8])

#2017
omnivores_2017_r_df <- as.data.frame(omnivores_2017_collapsed_rev)
colnames(omnivores_2017_r_df) <- "ch"
omnivores_2017_r_df <- bind_cols(omnivores_2017_r_df, omnivores_2017_df[,2:8])

#2018
omnivores_2018_r_df <- as.data.frame(omnivores_2018_collapsed_rev)
colnames(omnivores_2018_r_df) <- "ch"
omnivores_2018_r_df <- bind_cols(omnivores_2018_r_df, omnivores_2018_df[,2:8])

#2019
omnivores_2019_r_df <- as.data.frame(omnivores_2019_collapsed_rev)
colnames(omnivores_2019_r_df) <- "ch"
omnivores_2019_r_df <- bind_cols(omnivores_2019_r_df, omnivores_2019_df[,2:8])


#### Combine all years together

insectivores_marked_complete_rev <- bind_rows(insectivores_2010_r_df,
                                              insectivores_2012_r_df,
                                              insectivores_2013_r_df,
                                              insectivores_2014_r_df,
                                              insectivores_2015_r_df,
                                              insectivores_2016_r_df,
                                              insectivores_2017_r_df,
                                              insectivores_2018_r_df,
                                              insectivores_2019_r_df)


omnivores_marked_complete_rev <- bind_rows(omnivores_2010_r_df,
                                           omnivores_2012_r_df,
                                           omnivores_2013_r_df,
                                           omnivores_2014_r_df,
                                           omnivores_2015_r_df,
                                           omnivores_2016_r_df,
                                           omnivores_2017_r_df,
                                           omnivores_2018_r_df,
                                           omnivores_2019_r_df)


################################################################################
### Reverse time CJS model fitting
################################################################################

# Note - the example below is for insectivores during spring migration only.
# For fall migration or for omnivores, use the same code but substitute the 
# required datasets


### Subset data to known sex categories (insectivores only)
insectivores_marked_rev <- subset(insectivores_marked_complete_rev, sex == "M" |
                                      sex == "F")
### Make sex a factor
insectivores_marked_rev$sex <- factor(insectivores_marked_rev$sex)

### Remove NAs in fuel index or models will not fit
insectivores_marked_rev <- subset(insectivores_marked_rev, is.na(fuel_index) == FALSE)

### Convert year into a factor for random effect
insectivores_marked_rev$year_f <- factor(insectivores_marked_rev$year)

### Remove age, and original year variables or model will not fit.
insectivores_marked_rev <- insectivores_marked_rev %>%
                             select(-c(age,year))


########################
### Fit models
########################


### Full model

insectivores_cjs_rev_r <- crm(insectivores_marked_rev,
                             begin.time = 1,
                             model = "probitCJS",
                             model.parameters = list(Phi = Phi_sam_r, p = p_sas_r),
                             design_parameters = list(Phi = design_Phi, p = design_p),
                             groups = "year_f",
                             use.admb = TRUE,
                             burnin = 1000,iter=3000,
                             seed = 100)

insectivores_cjs_rev_r


### Reduced model

insectivores_cjs_rev <- crm(insectivores_marked_rev,
                           begin.time = 1,
                           model = "hmmCJS",
                           model.parameters = list(Phi = Phi_sam, p = p_sas),
                           design_parameters = list(Phi = design_Phi, p = design),
                           use.admb = TRUE,
                           burnin = 1000,iter=3000,
                           seed = 100)



### Extract real values and convert into stopover duration

insectivores_cjs_rev$results$reals


# Mean duration
-1 / log(0.7469972) ## 3.43 days COYE
-1 / log(0.6330513) ## 2.19 days YRWA
-1 / log(0.5284006) ## 1.57 days OCWA
-1 / log(0.4814481) ## 1.37 days WIWA
-1 / log(0.4982643) ## 1.44 days YEWA



########################## End of code #########################################