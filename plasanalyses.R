## Plasticity analyses ##

library(tidyverse)
library(lubridate)
library(cowplot)
library(export)
library(glmmTMB)
library(DHARMa)
library(arm)
library(lme4)
library(lmerTest)
library(rptR)
library(mosaic)
library(broom)
library(reshape2)
library(rstatix)
library(common)
library(ggpubr)
library(climwin)
library(MuMIn)

mochhigh <- read.csv("mochhigh.csv")
mochlow <- read.csv("mochlow.csv")
snotel_data <- read.csv("snotel_data_2023.csv")

snotel_data$Date <- as.Date(snotel_data$Date, format = "%m/%d/%Y")
snotel_data$date <- format(snotel_data$Date, "%d/%m/%Y")
mochlow$FIRST.EGG <- format(as.Date(mochlow$FIRST.EGG, format = "%m/%d/%Y"), "%d/%m/%Y")
mochlow <- mochlow%>% filter(!is.na(FIRST.EGG))
mochhigh$FIRST.EGG <- format(as.Date(mochhigh$FIRST.EGG, format = "%m/%d/%Y"), "%d/%m/%Y")
mochhigh <- mochhigh%>% filter(!is.na(FIRST.EGG))

#### Climwin analyses ####

# Create a baseline model with no explanatory variables for first egg date at low elevation
baselow<-lm(J.FIRST.EGG ~ 1, data = mochlow) 

# Check residuals of baseline model using the 'DHARMa' package
simulationOutput <- simulateResiduals(baselow, plot = T) 

### Low Elevation
## Precipitation Accumulation
# x = Precipitation Accum
# Recalculate precip accum so that the year starts Sept 1
snotel_data$Date = as.Date(snotel_data$Date,"%m/%d/%Y")
for (j in 2:nrow(snotel_data)){
  if (str_detect(snotel_data$Date[j], regex(".-10-01", dotall = TRUE))) {
    snotel_data$low_precip_raw[j] = snotel_data$Low_Precip_Accum.mm.[j] 
  } else {
    snotel_data$low_precip_raw[j] = (snotel_data$Low_Precip_Accum.mm.[j] - snotel_data$Low_Precip_Accum.mm.[j-1])
  }
}
snotel_data$low_precip_accum_fixed<-NA
for (j in 2:nrow(snotel_data)){
  if (str_detect(snotel_data$Date[j], regex(".-09-02", dotall = TRUE))){
    snotel_data$low_precip_accum_fixed[j] = snotel_data$low_precip_raw[j]
  } else {
    snotel_data$low_precip_accum_fixed[j] = (snotel_data$low_precip_raw[j] + snotel_data$low_precip_accum_fixed[j-1])
  }
}


resultslow_precip_eggdate_mean_lin<-slidingwin(baseline = baselow, 
                                               exclude = c(14,-1), 
                                               xvar = list(snotel_data$low_precip_accum_fixed),
                                               type = c("absolute"), 
                                               range = c(244,0), 
                                               stat = c("mean"), 
                                               func = c("lin"),
                                               refday = c(3,5), 
                                               cmissing = "method1",
                                               cinterval = "day", 
                                               cdate = snotel_data$date, bdate = mochlow$FIRST.EGG)

# Save model file
saveRDS(resultslow_precip_eggdate_mean_lin,here::here("resultslow_precip_eggdate_mean_lin_11years.rds"))

# Read in model file
resultslow_precip_eggdate_mean_lin = readRDS(here::here("resultslow_precip_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultslow_precip_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultslow_precip_eggdate_mean_lin.dataset = resultslow_precip_eggdate_mean_lin[[1]]$Dataset
# Check summary
summary(resultslow_precip_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultslow_precip_eggdate_mean_lin.dataset)
# Plot all
plotall(resultslow_precip_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultslow_precip_eggdate_mean_lin.dataset,
         bestmodel = resultslow_precip_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultslow_precip_eggdate_mean_lin[[1]]$BestModelData)

## Average Daily Temperature
# x = Average Daily Temperature
resultslow_avgtemp_eggdate_mean_lin<-slidingwin(baseline = baselow, 
                                                exclude = c(14,-1), 
                                                xvar = list(snotel_data$Low_Air_Temp_Avg.degC.),
                                                type = c("absolute"), 
                                                range = c(91,0), 
                                                stat = c("mean"), 
                                                func = c("lin"),
                                                refday = c(3,5), 
                                                cmissing = "method1",
                                                cinterval = "day", 
                                                cdate = snotel_data$date, bdate = mochlow$FIRST.EGG)

# Save model file
saveRDS(resultslow_avgtemp_eggdate_mean_lin,here::here("resultslow_avgtemp_eggdate_mean_lin_11years.rds"))

# Read in model file
resultslow_avgtemp_eggdate_mean_lin = readRDS(here::here("resultslow_avgtemp_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultslow_avgtemp_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultslow_avgtemp_eggdate_mean_lin.dataset = resultslow_avgtemp_eggdate_mean_lin[[1]]$Dataset
# Check summary 
summary(resultslow_avgtemp_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultslow_avgtemp_eggdate_mean_lin.dataset)
# Plot all
plotall(resultslow_avgtemp_eggdate_mean_lin.dataset)
# Plot best
plotbest(dataset = resultslow_avgtemp_eggdate_mean_lin.dataset,
         bestmodel = resultslow_avgtemp_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultslow_avgtemp_eggdate_mean_lin[[1]]$BestModelData)

## Max Daily Temperature 
resultslow_maxtemp_eggdate_mean_lin<-slidingwin(baseline = baselow, 
                                                exclude = c(14,-1), 
                                                xvar = list(snotel_data$Low_Air_Temp_Max.degC.),
                                                type = c("absolute"), 
                                                range = c(91,0), 
                                                stat = c("mean"), 
                                                func = c("lin"),
                                                refday = c(3,5), 
                                                cmissing = "method1",
                                                cinterval = "day", 
                                                cdate = snotel_data$date, bdate = mochlow$FIRST.EGG)

# Save model file
saveRDS(resultslow_maxtemp_eggdate_mean_lin,here::here("resultslow_maxtemp_eggdate_mean_lin_11years.rds"))

# Read in model file
resultslow_maxtemp_eggdate_mean_lin = readRDS(here::here("resultslow_maxtemp_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultslow_maxtemp_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultslow_maxtemp_eggdate_mean_lin.dataset = resultslow_maxtemp_eggdate_mean_lin[[1]]$Dataset
# Check summary
summary(resultslow_maxtemp_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultslow_maxtemp_eggdate_mean_lin.dataset)
# Plot all
plotall(resultslow_maxtemp_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultslow_maxtemp_eggdate_mean_lin.dataset,
         bestmodel = resultslow_maxtemp_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultslow_maxtemp_eggdate_mean_lin[[1]]$BestModelData)

## Minimum Daily Temperature 
resultslow_mintemp_eggdate_mean_lin<-slidingwin(baseline = baselow, 
                                                exclude = c(14,-1), 
                                                xvar = list(snotel_data$Low_Air_Temp_Min.degC.),
                                                type = c("absolute"), 
                                                range = c(91,0), 
                                                stat = c("mean"), 
                                                func = c("lin"),
                                                refday = c(3,5), 
                                                cmissing = "method1",
                                                cinterval = "day", 
                                                cdate = snotel_data$date, bdate = mochlow$FIRST.EGG)

# Save model file
saveRDS(resultslow_mintemp_eggdate_mean_lin,here::here("resultslow_mintemp_eggdate_mean_lin_11years.rds"))

# Read in model file
resultslow_mintemp_eggdate_mean_lin = readRDS(here::here("resultslow_mintemp_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultslow_mintemp_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultslow_mintemp_eggdate_mean_lin.dataset = resultslow_mintemp_eggdate_mean_lin[[1]]$Dataset
# Check summary  
summary(resultslow_mintemp_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultslow_mintemp_eggdate_mean_lin.dataset)
# Plot all
plotall(resultslow_mintemp_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultslow_mintemp_eggdate_mean_lin.dataset,
         bestmodel = resultslow_mintemp_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultslow_mintemp_eggdate_mean_lin[[1]]$BestModelData)

 
### Randomizations - Low Elevation
#Create a list
randlist = seq(1,100,by=1)

## Low elevation

# Get month and day of first egg date
# mochlow$FIRST.EGG<-as.Date(mochlow$FIRST.EGG, "%d/%m/%Y")
mochlow$month = month(mochlow$FIRST.EGG)
mochlow$day = day(mochlow$FIRST.EGG)

# Get a list of the years
colnames(mochlow)[2] = "Year.real"

# Create a dataframe from list of unique years
Years = data.frame(mochlow[!duplicated(mochlow$Year.real),]$Year.real)
colnames(Years) = "Year.real"

# Avg temperature randomization
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochlow dataframe
    mochlow.rand = merge(mochlow,Years,by="Year.real")

    #Unite random years with month and day
    mochlow.rand2 = unite(mochlow.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F)

    #Get as date
    mochlow.rand2$Date.rand = as.Date(mochlow.rand2$Date.rand)

    #Run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(AvgTemp = snotel_data$Low_Air_Temp_Avg.degC.),
                      cdate = snotel_data$date,
                      bdate = mochlow.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochlow),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(91,1), #244, 1
                      refday = c(3,5), #julian day 123 min for low
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

# Apply the function to the whole list of 100 - get 100 randomizations
out<-lapply(randlist, randfunc)

# Export the data to a dataframe with the appropriate columns
rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

# Unlist the data from lapply to add to the export dataframe
for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

# Export to a csv
write.csv(rand.data,"randomizations_low_firstegg_mean_avgdailytemp_abs.csv")

# Max temp
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochlow dataframe
    mochlow.rand = merge(mochlow,Years,by="Year.real")

    #Unite random years with month and day
    mochlow.rand2 = unite(mochlow.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F)

    #Get as date
    mochlow.rand2$Date.rand = as.Date(mochlow.rand2$Date.rand)

    #Run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(MaxTemp = snotel_data$Low_Air_Temp_Max.degC.), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochlow.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochlow),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(91,1), #244, 1
                      refday = c(3,5), #julian day 123 min for low
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)
rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_low_firstegg_mean_maxdailytemp_abs.csv")

# Min temp
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochlow dataframe
    mochlow.rand = merge(mochlow,Years,by="Year.real")

    #Unite random years with month and day
    mochlow.rand2 = unite(mochlow.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F)

    #Get as date
    mochlow.rand2$Date.rand = as.Date(mochlow.rand2$Date.rand)

    #Run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(MinTemp = snotel_data$Low_Air_Temp_Min.degC.), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochlow.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochlow),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(91,1), #244, 1
                      refday = c(3,5), #julian day 123 min for low
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)

rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_low_firstegg_mean_mindailytemp_abs.csv")

## Precipitation

randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochlow dataframe
    mochlow.rand = merge(mochlow,Years,by="Year.real")

    #Unite random years with month and day
    mochlow.rand2 = unite(mochlow.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F)

    #Get as date
    mochlow.rand2$Date.rand = as.Date(mochlow.rand2$Date.rand)

    #Run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(Precip = snotel_data$Low_Precip_Accum.mm.), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochlow.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochlow),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(244,1), #244, 1
                      refday = c(3,5), #julian day 123 min for low
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)

rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_low_firstegg_mean_precip_abs.csv")

### High elevation
## Step 1. Baseline Model. 
basehigh<-lm(J.FIRST.EGG ~ 1, data = mochhigh) 

# Check residuals of model using the 'DHARMa' package
simulationOutput <- simulateResiduals(basehigh, plot = T) 

## Precipitation Accumulation
#High Elevation
#x = Precipitation Accum
#need to recalculate precip accum so that the year starts Sept 1
# snotel_data$Date = as.Date(snotel_data$Date,"%m/%d/%Y")
for (j in 2:nrow(snotel_data)){
  if (str_detect(snotel_data$Date[j], regex(".-10-01", dotall = TRUE))) {
    snotel_data$high_precip_raw[j] = snotel_data$X541_Precip_Accum.mm.[j] 
  } else {
    snotel_data$high_precip_raw[j] = (snotel_data$X541_Precip_Accum.mm.[j] - snotel_data$X541_Precip_Accum.mm.[j-1])
  }
}
snotel_data$high_precip_accum_fixed<-NA
for (j in 2:nrow(snotel_data)){
  if (str_detect(snotel_data$Date[j], regex(".-09-02", dotall = TRUE))){
    snotel_data$high_precip_accum_fixed[j] = snotel_data$high_precip_raw[j]
  } else {
    snotel_data$high_precip_accum_fixed[j] = (snotel_data$high_precip_raw[j] + snotel_data$high_precip_accum_fixed[j-1])
  }
}

resultshigh_precip_eggdate_mean_lin<-slidingwin(baseline = basehigh, 
                                                exclude = c(14,-1), 
                                                xvar = list(snotel_data$high_precip_accum_fixed),
                                                type = c("absolute"), 
                                                range = c(260,0), 
                                                stat = c("mean"), 
                                                func = c("lin"),
                                                refday = c(19,5), 
                                                cmissing = "method1",
                                                cinterval = "day", 
                                                cdate = snotel_data$date, bdate = mochhigh$FIRST.EGG)

# Save model file using the 'here' package
saveRDS(resultshigh_precip_eggdate_mean_lin,here::here("resultshigh_precip_eggdate_mean_lin_11years.rds"))

# Read in model file
resultshigh_precip_eggdate_mean_lin = readRDS(here::here("resultshigh_precip_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultshigh_precip_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultshigh_precip_eggdate_mean_lin.dataset = resultshigh_precip_eggdate_mean_lin[[1]]$Dataset
# Check summary 
summary(resultshigh_precip_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultshigh_precip_eggdate_mean_lin.dataset)
# Plotall
plotall(resultshigh_precip_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultshigh_precip_eggdate_mean_lin.dataset,
         bestmodel = resultshigh_precip_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultshigh_precip_eggdate_mean_lin[[1]]$BestModelData)

## Average Daily Temperature
# High Elevation
# x = Average Daily Temperature 
# linear function
resultshigh_avgtemp_eggdate_mean_lin<-slidingwin(baseline = basehigh, 
                                                 exclude = c(14,-1), 
                                                 xvar = list(snotel_data$X541_Air_Temp_Avg.degC.),
                                                 type = c("absolute"), 
                                                 range = c(107,0), 
                                                 stat = c("mean"), 
                                                 func = c("lin"),
                                                 refday = c(19,5), 
                                                 cmissing = "method1",
                                                 cinterval = "day", 
                                                 cdate = snotel_data$date, bdate = mochhigh$FIRST.EGG)

# Save model file
saveRDS(resultshigh_avgtemp_eggdate_mean_lin,here::here("resultshigh_avgtemp_eggdate_mean_lin_11years.rds"))

# Read in model file
resultshigh_avgtemp_eggdate_mean_lin = readRDS(here::here("resultshigh_avgtemp_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultshigh_avgtemp_eggdate_mean_lin[[1]]$Dataset,n=20)

# Get dataset
resultshigh_avgtemp_eggdate_mean_lin.dataset = resultshigh_avgtemp_eggdate_mean_lin[[1]]$Dataset
# Check summary
summary(resultshigh_avgtemp_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultshigh_avgtemp_eggdate_mean_lin.dataset)
# Plot all
plotall(resultshigh_avgtemp_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultshigh_avgtemp_eggdate_mean_lin.dataset,
         bestmodel = resultshigh_avgtemp_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultshigh_avgtemp_eggdate_mean_lin[[1]]$BestModelData)

## Maximum Daily Temperature
#High Elevation
#x = Maximum Daily Temperature
#linear function
resultshigh_maxtemp_eggdate_mean_lin<-slidingwin(baseline = basehigh, 
                                                 exclude = c(10,-1), 
                                                 xvar = list(snotel_data$X541_Air_Temp_Max.degC.),
                                                 type = c("absolute"), 
                                                 range = c(107,0), 
                                                 stat = c("mean"), 
                                                 func = c("lin"),
                                                 refday = c(19,5), 
                                                 cmissing = "method1",
                                                 cinterval = "day", 
                                                 cdate = snotel_data$date, bdate = mochhigh$FIRST.EGG)

# Save model file
saveRDS(resultshigh_maxtemp_eggdate_mean_lin,here::here("resultshigh_maxtemp_eggdate_mean_lin_11years.rds"))

# Read in model file
resultshigh_maxtemp_eggdate_mean_lin = readRDS(here::here("resultshigh_maxtemp_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultshigh_maxtemp_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultshigh_maxtemp_eggdate_mean_lin.dataset = resultshigh_maxtemp_eggdate_mean_lin[[1]]$Dataset
# Check summary 
summary(resultshigh_maxtemp_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultshigh_maxtemp_eggdate_mean_lin.dataset)
# Plot all
plotall(resultshigh_maxtemp_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultshigh_maxtemp_eggdate_mean_lin.dataset,
         bestmodel = resultshigh_maxtemp_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultshigh_maxtemp_eggdate_mean_lin[[1]]$BestModelData)

## Minimum Daily Temperature
# High Elevation
# x = Minimum Daily Temperature - probably want to look into year-round trends.
resultshigh_mintemp_eggdate_mean_lin<-slidingwin(baseline = basehigh, 
                                                 exclude = c(10,-1), 
                                                 xvar = list(snotel_data$X541_Air_Temp_Min.degC.),
                                                 type = c("absolute"), 
                                                 range = c(107,0), 
                                                 stat = c("mean"), 
                                                 func = c("lin"),
                                                 refday = c(19,5), 
                                                 cmissing = "method1",
                                                 cinterval = "day", 
                                                 cdate = snotel_data$date, bdate = mochhigh$FIRST.EGG)

# Save model file
saveRDS(resultshigh_mintemp_eggdate_mean_lin,here::here("resultshigh_mintemp_eggdate_mean_lin_11years.rds"))

# Read in model file
resultshigh_mintemp_eggdate_mean_lin = readRDS(here::here("resultshigh_mintemp_eggdate_mean_lin_11years.rds"))
# Look to see what windows are of best model 
head(resultshigh_mintemp_eggdate_mean_lin[[1]]$Dataset,n=20)
# Get dataset
resultshigh_mintemp_eggdate_mean_lin.dataset = resultshigh_mintemp_eggdate_mean_lin[[1]]$Dataset
# Check summary 
summary(resultshigh_mintemp_eggdate_mean_lin[[1]]$BestModel)
# Plot delta plot
plotdelta(dataset=resultshigh_mintemp_eggdate_mean_lin.dataset)
# Plot all
plotall(resultshigh_mintemp_eggdate_mean_lin.dataset)
# Plot best model
plotbest(dataset = resultshigh_mintemp_eggdate_mean_lin.dataset,
         bestmodel = resultshigh_mintemp_eggdate_mean_lin[[1]]$BestModel, 
         bestmodeldata = resultshigh_mintemp_eggdate_mean_lin[[1]]$BestModelData)


## Randomizations - High elevation
#Create a list
randlist = seq(1,100,by=1)

## High elevation, avg daily temp

# Get month and day of first egg date
mochhigh$FIRST.EGG<-as.Date(mochhigh$FIRST.EGG, "%d/%m/%Y")
mochhigh$month = month(mochhigh$FIRST.EGG)
mochhigh$day = day(mochhigh$FIRST.EGG)
mochhigh$FIRST.EGG<-format(as.Date(mochhigh$FIRST.EGG, "%d/%m/%Y"), "%d/%m/%Y")

#Get a list of the years
colnames(mochhigh)[2] = "Year.real"

Years = data.frame(mochhigh[!duplicated(mochhigh$Year.real),]$Year.real) #create a dataframe from list of unique years
colnames(Years) = "Year.real"

# Avg temp
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochhigh dataframe
    mochhigh.rand = merge(mochhigh,Years,by="Year.real")

    #Unite random years with month and day
    mochhigh.rand2 = unite(mochhigh.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F) #what does this do?

    #Get as date
    mochhigh.rand2$Date.rand = as.Date(mochhigh.rand2$Date.rand)

    #run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(AvgTemp = snotel_data$X541_Air_Temp_Avg.degC.), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochhigh.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochhigh),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(107,1), #230, 1
                      refday = c(19,5), #julian day 139 min first egg date at high
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    #dataset = dataset[which(dataset$ModelBeta<0),] #Get only negative relationships - observed relationship was negative
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)

rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_high_firstegg_mean_avgdailytemp_abs.csv")

# Maximum temp
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochhigh dataframe
    mochhigh.rand = merge(mochhigh,Years,by="Year.real")

    #Unite random years with month and day
    mochhigh.rand2 = unite(mochhigh.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F) #what does this do?

    #Get as date
    mochhigh.rand2$Date.rand = as.Date(mochhigh.rand2$Date.rand)

    #run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(MaxTemp = snotel_data$X541_Air_Temp_Max.degC.), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochhigh.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochhigh),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(107,1), #230, 1
                      refday = c(19,5), #julian day 139 min first egg date at high
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    #dataset = dataset[which(dataset$ModelBeta<0),] #Get only negative relationships - observed relationship was negative
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)

rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_high_firstegg_mean_maxdailytemp_abs.csv")

# Minumum temp
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochhigh dataframe
    mochhigh.rand = merge(mochhigh,Years,by="Year.real")

    #Unite random years with month and day
    mochhigh.rand2 = unite(mochhigh.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F) #what does this do?

    #Get as date
    mochhigh.rand2$Date.rand = as.Date(mochhigh.rand2$Date.rand)

    #run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(MinTemp = snotel_data$X541_Air_Temp_Min.degC.), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochhigh.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochhigh),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(107,1), #230, 1
                      refday = c(19,5), #julian day 139 min first egg date at high
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    #dataset = dataset[which(dataset$ModelBeta<0),] #Get only negative relationships - observed relationship was negative
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)

rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_high_firstegg_mean_mindailytemp_abs.csv")

# Precipitation
# Recalculate precipitation accumulation to begin on Sept 1
snotel_data$Date = as.Date(snotel_data$Date,"%m/%d/%Y")
for (j in 2:nrow(snotel_data)){
  if (str_detect(snotel_data$Date[j], regex(".-10-01", dotall = TRUE))) {
    snotel_data$high_precip_raw[j] = snotel_data$X541_Precip_Accum.mm.[j]
  } else {
    snotel_data$high_precip_raw[j] = (snotel_data$X541_Precip_Accum.mm.[j] - snotel_data$X541_Precip_Accum.mm.[j-1])
  }
}
snotel_data$high_precip_accum_fixed<-NA
for (j in 2:nrow(snotel_data)){
  if (str_detect(snotel_data$Date[j], regex(".-09-02", dotall = TRUE))){
    snotel_data$high_precip_accum_fixed[j] = snotel_data$high_precip_raw[j]
  } else {
    snotel_data$high_precip_accum_fixed[j] = (snotel_data$high_precip_raw[j] + snotel_data$high_precip_accum_fixed[j-1])
  }
}
randfunc = function(randlist) {
  for (i in 1:length(randlist)) {
    #Get a random year for each real year and don't let any year get assigned to itself
    repeat {Years$Year.rand = sample(Years$Year.real,replace=F) #Repeat year randomization until no year assigned to itself
    Years.test = Years[which(Years$Year.real==Years$Year.rand),]
    if(nrow(Years.test)==0) {break}
    }

    #Merge random year with mochhigh dataframe
    mochhigh.rand = merge(mochhigh,Years,by="Year.real")

    #Unite random years with month and day
    mochhigh.rand2 = unite(mochhigh.rand, Date.rand, c(Year.rand,month,day),sep="/",remove=F) #what does this do?

    #Get as date
    mochhigh.rand2$Date.rand = as.Date(mochhigh.rand2$Date.rand)

    #run a sliding window analysis and pick the best delta-AIC value and put into dataframe
    win <- slidingwin(xvar=list(Precip = snotel_data$high_precip_accum_fixed), #put an environmental variable here
                      cdate = snotel_data$date,
                      bdate = mochhigh.rand2$Date.rand,
                      baseline = lm(J.FIRST.EGG ~ 1, data = mochhigh),
                      type = "absolute",
                      cinterval = "day",
                      exclude = c(14, -1),
                      range = c(260,1), #230, 1
                      refday = c(19,5), #julian day 139 min first egg date at high
                      stat = "mean",
                      cmissing = "method1",
                      fun = "lin")

    dataset = win[[1]]$Dataset
    #dataset = dataset[which(dataset$ModelBeta<0),] #Get only negative relationships - observed relationship was negative
    dataset = dataset[order(dataset$deltaAICc),] #Make sure smallest AIC value is still first
    return(list(dataset$deltaAIC[1],dataset$WindowOpen[1],dataset$WindowClose[1],dataset$ModelBeta[1])) #return best deltaAIC value
  }}

out<-lapply(randlist, randfunc)

rand.data <- matrix(NA, ncol = 4, nrow = length(out))
colnames(rand.data) = c("deltaAICc","Open","Close","ModelBeta")

for (i in 1:length(out)) {
  rand.data[i,] = unlist(out[[i]])
}

write.csv(rand.data,"randomizations_high_firstegg_mean_precip_abs.csv")



#### Random Regression Models ####
## Extract climate variables and scale them
# Calculate the yday for snotel data
snotel_data$yDay <- yday(snotel_data$Date)

## High elevation
# Precipitation accumulation
preciphighwindow<-snotel_data[snotel_data$yDay > 47 & snotel_data$yDay < 82,] # Pull out the window from Julian day 48-81
preciphighwindow<-preciphighwindow%>%group_by(year)%>%summarise(precip = mean(high_precip_accum_fixed, na.rm = TRUE)) # Calculate the mean precipitation over that window
colnames(preciphighwindow)[1] <- "YEAR" # Rename column to match other dataframes

# Average temperature
avgtemphighwindow<- snotel_data[snotel_data$yDay > 85 & snotel_data$yDay < 122,] # Pull out the window from Julian day 86-121
avgtemphighwindow<-avgtemphighwindow %>% group_by(year) %>% summarise(avgtemp = mean(X541_Air_Temp_Avg.degC., na.rm = TRUE))
colnames(avgtemphighwindow)[1] <- "YEAR"

# Minimum temperature
mintemphighwindow<- snotel_data[snotel_data$yDay > 85 & snotel_data$yDay < 106,] # Pull out the window from Julian day 86-105
mintemphighwindow<-mintemphighwindow %>% group_by(year) %>% summarise(mintemp = mean(X541_Air_Temp_Min.degC., na.rm = TRUE))
colnames(mintemphighwindow)[1] <- "YEAR"

# Maximum temperature
maxtemphighwindow<- snotel_data[snotel_data$yDay > 85 & snotel_data$yDay < 117,] # Pull out the window from Julian day 86-116
maxtemphighwindow<-maxtemphighwindow %>% group_by(year) %>% summarise(maxtemp = mean(X541_Air_Temp_Max.degC., na.rm = TRUE))
colnames(maxtemphighwindow)[1] <- "YEAR"

# Merge dataframes
mochhigh1 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(mochhigh, preciphighwindow, avgtemphighwindow, mintemphighwindow, maxtemphighwindow))
mochhigh1$precip_scale <- zscore(mochhigh1$precip)
mochhigh1$avgtemp_scale <- zscore(mochhigh1$avgtemp)
mochhigh1$mintemp_scale <- zscore(mochhigh1$mintemp)
mochhigh1$maxtemp_scale <- zscore(mochhigh1$maxtemp)

## Low elevation
# Precipitation accumulation
preciplowwindow<-snotel_data[snotel_data$yDay > 14 & snotel_data$yDay < 39,] # Pull out the window from Julian day 15-38
preciplowwindow<-preciplowwindow%>%group_by(year)%>%summarise(precip = mean(low_precip_accum_fixed, na.rm = TRUE)) # Calculate the mean precipitation over that window
colnames(preciplowwindow)[1] <- "YEAR" # Rename column to match other dataframes

# Average temperature
avgtemplowwindow<- snotel_data[snotel_data$yDay > 86 & snotel_data$yDay < 124,] # Pull out the window from Julian day 87-123
avgtemplowwindow<-avgtemplowwindow %>% group_by(year) %>% summarise(avgtemp = mean(Low_Air_Temp_Avg.degC., na.rm = TRUE))
colnames(avgtemplowwindow)[1] <- "YEAR"

# Minimum temperature
mintemplowwindow<- snotel_data[snotel_data$yDay > 103 & snotel_data$yDay < 124,] # Pull out the window from Julian day 104-123
mintemplowwindow<-mintemplowwindow %>% group_by(year) %>% summarise(mintemp = mean(Low_Air_Temp_Min.degC., na.rm = TRUE))
colnames(mintemplowwindow)[1] <- "YEAR"

# Maximum temperature
maxtemplowwindow<- snotel_data[snotel_data$yDay > 86 & snotel_data$yDay < 124,] # Pull out the window from Julian day 87-123
maxtemplowwindow<-maxtemplowwindow %>% group_by(year) %>% summarise(maxtemp = mean(Low_Air_Temp_Max.degC., na.rm = TRUE))
colnames(maxtemplowwindow)[1] <- "YEAR"

# Merge dataframes
mochlow1 <- Reduce(function(x, y) merge(x, y, all=TRUE), list(mochlow, preciplowwindow, avgtemplowwindow, mintemplowwindow, maxtemplowwindow))
mochlow1$precip_scale <- zscore(mochlow1$precip)
mochlow1$avgtemp_scale <- zscore(mochlow1$avgtemp)
mochlow1$mintemp_scale <- zscore(mochlow1$mintemp)
mochlow1$maxtemp_scale <- zscore(mochlow1$maxtemp)

# Combine datasets 
mochdata <- rbind(mochhigh1, mochlow1)
mochdata<-mochdata%>%filter(!is.na(F.ID))%>%filter(!F.ID == "")%>%filter(!F.ID == "unknown")%>%group_by(F.ID) %>% filter(n()>2)


## Check for best temperature correlate
mtemp1 <- lmer(J.FIRST.EGG ~ ELEVATION*avgtemp_scale + (1|YEAR) + (avgtemp_scale|F.ID), data = mochdata)
summary(mtemp1)
Anova(mtemp1, Type = "3")
mtemp2 <- lmer(J.FIRST.EGG ~ ELEVATION*mintemp_scale + (1|YEAR) + (mintemp_scale|F.ID), data = mochdata)
summary(mtemp2)
Anova(mtemp2, Type = "3")
mtemp3 <- lmer(J.FIRST.EGG ~ ELEVATION*maxtemp_scale + (1|YEAR) + (maxtemp_scale|F.ID), data = mochdata)
summary(mtemp3)
Anova(mtemp3, Type = "3")
AICc(mtemp1, mtemp2, mtemp3) # min and max temp seem to be best fit

## Test whether precip or temp are better models
mprecip1 <- lmer(J.FIRST.EGG ~ ELEVATION*precip_scale + (1|YEAR) + (precip_scale|F.ID), data = mochdata)
summary(mprecip1)
Anova(mprecip1, Type = "3")
AICc(mtemp3, mprecip1) # precip is a better model than temperature - do we still want to report both variables? Probably...

### Extract individual slopes and intercepts 
## Maximum temperature models
temp_rn<-ranef(mtemp2)[["F.ID"]]
colnames(temp_rn)<-c("Intercept", "Slope")
temp_rn$F.ID<-rownames(temp_rn)

# Add the population-level slope and intercept to the model
temp_rn$Slope <- temp_rn$Slope + summary(mtemp2)$coef[3]
temp_rn$Intercept <- temp_rn$Intercept + summary(mtemp2)$coef[1]

# Extract standard errors using the arm package
temp_se<-se.ranef(mtemp2)

# Build dataframe
# Add standard errors to the temp_rn dataframe
temp_rn$Intercept_se <- temp_se$F.ID[,1]
temp_rn$Slope_se <- temp_se$F.ID[,2]

# Add elevation to the dataframe
temp_rn <- left_join(temp_rn, mochdata[, c(5,12)], multiple = "first")

# Plot the reaction norms
temp_rn_plot<-ggplot(temp_rn) +
  geom_abline(aes(slope = Slope, intercept = Intercept, colour = Slope)) +
  scale_color_continuous(type = "viridis")+
  xlim(c(min(mochdata$maxtemp_scale), max(mochdata$maxtemp_scale))) + ylim(c(142,165)) + 
  labs(x = "Scaled Mean Minimum Temperature (\u00B0C)", y = "First Egg Date (day of year)") + 
  theme_cowplot() + facet_wrap(factor(precip_rn$ELEVATION, levels = c("H", "L"), labels = c("High Elevation", "Low Elevation"))) + geom_vline(aes(xintercept = 0), color="black", linetype="dashed")
temp_rn_plot
ggsave("tempplasplot.jpg", temp_rn_plot, dpi = 600, width = 5, height = 4, bg = "white")

## Precip models
precip_rn<-ranef(mprecip1)[["F.ID"]]
colnames(precip_rn)<-c("Intercept", "Slope")
precip_rn$F.ID<-rownames(precip_rn)

# Add the population-level slope and intercept to the model
precip_rn$Slope <- precip_rn$Slope + summary(mprecip1)$coef[3]
precip_rn$Intercept <- precip_rn$Intercept + summary(mprecip1)$coef[1]

# Extract standard errors using the arm package
precip_se<-se.ranef(mprecip1)

# Build dataframe
# Add standard errors to the temp_rn dataframe
precip_rn$Intercept_se <- precip_se$F.ID[,1]
precip_rn$Slope_se <- precip_se$F.ID[,2]

# Add elevation to the dataframe
precip_rn <- left_join(precip_rn, mochdata[, c(5,12)], multiple = "first")

# Plot the reaction norms 
precip_rn_plot<-ggplot(precip_rn) +
  geom_abline(aes(slope = Slope, intercept = Intercept, colour = Slope)) +
  scale_color_continuous(type = "viridis")+
  xlim(c(min(mochdata$precip_scale), max(mochdata$precip_scale))) + ylim(c(142, 170)) + 
  labs(x = "Scaled Precipitation Accumulation (mm)", y = "First Egg Date (day of year)") + 
  theme_cowplot() + facet_wrap(factor(precip_rn$ELEVATION, levels = c("H", "L"), labels = c("High Elevation", "Low Elevation"))) + geom_vline(aes(xintercept = 0), color="black", linetype="dashed")
precip_rn_plot
ggsave("precipplasplot.jpg", precip_rn_plot, dpi = 600, width = 5, height = 4, bg = "white")

### Fitness proxies/measurements
## Lifetime breeding success
rm_birds<-unique(mochdata[mochdata$YEAR == 2023,]$F.ID)
lbs <- mochdata %>% filter(!is.na(F.ID)) %>% filter(F.ID != "") %>% filter(!F.ID%in%(rm_birds)) %>% 
  mutate(BROOD = ifelse(is.na(BROOD), 0, BROOD)) %>% group_by(F.ID) %>% filter(n()>1) %>% summarise(LBS = sum(BROOD))

# Models investigating relationship between LBS and slope of reaction norms
# Temperature
temp_data_lbs<-merge(temp_rn, lbs, by = "F.ID")
temp_data_lbs<-temp_data_lbs%>%filter(!is.na(Slope))
temp_lbs_ml <- lm(LBS ~ Slope*ELEVATION, data = temp_data_lbs)
summary(temp_lbs_ml) 
Anova(temp_lbs_ml, type = "III") # not significant
temp_lbs_ml_noelev <- lm(LBS ~ Slope, data = temp_data_lbs) # remove interaction
summary(temp_lbs_ml_noelev) # not significant
temp_lbs_mq<-lm(LBS ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_lbs)
summary(temp_lbs_mq) 
Anova(temp_lbs_mq, type = "III") # not significant
temp_lbs_mq_noelev<-lm(LBS ~ poly(Slope, degree = 2, raw = FALSE), data = temp_data_lbs) # remove interaction
summary(temp_lbs_mq_noelev) # not significant
AICc(temp_lbs_ml_noelev, temp_lbs_mq_noelev)# linear better


# Precip
precip_data_lbs<-merge(precip_rn, lbs, by = "F.ID")
precip_data_lbs<-precip_data_lbs%>%filter(!is.na(Slope))
precip_lbs_ml<- lm(LBS ~ Slope*ELEVATION, data = precip_data_lbs)
summary(precip_lbs_ml) 
Anova(precip_lbs_ml, type = "III") # no significant variables
precip_lbs_ml_noelev<- lm(LBS ~ Slope, data = precip_data_lbs) # remove interaction
summary(precip_lbs_ml_noelev) # not significant
precip_lbs_mq<-lm(LBS ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_lbs)
summary(precip_lbs_mq) 
Anova(precip_lbs_mq, type = "III") # no significant variables
precip_lbs_mq_noelev<-lm(LBS ~ poly(Slope, degree = 2, raw = FALSE), data = precip_data_lbs)
summary(precip_lbs_mq_noelev) # not significant
AICc(precip_lbs_ml_noelev, precip_lbs_mq_noelev)


# LBS and intercept of reaction norms
# Temperature
temp_lbs_intl<-lm(LBS ~ Intercept*ELEVATION, data = temp_data_lbs)
summary(temp_lbs_intl) 
Anova(temp_lbs_intl, type = "III") # not significant
temp_lbs_intl_noelev<-lm(LBS ~ Intercept, data = temp_data_lbs) # remove interaction
summary(temp_lbs_intl_noelev) # not significant
temp_lbs_intq<-lm(LBS ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_lbs)
summary(temp_lbs_intq) 
Anova(temp_lbs_intq, type = "III") # no significant variables
temp_lbs_intq_noelev<-lm(LBS ~ poly(Intercept, degree = 2, raw = FALSE), data = temp_data_lbs)
summary(temp_lbs_intq_noelev) # not significant
AICc(temp_lbs_intl_noelev, temp_lbs_intq_noelev)

# Precip
precip_lbs_intl<-lm(LBS ~ Intercept*ELEVATION, data = precip_data_lbs)
summary(precip_lbs_intl) 
Anova(precip_lbs_intl, type = "III") # no significant variables
precip_lbs_intl_noelev<-lm(LBS ~ Intercept, data = precip_data_lbs) # remove interaction
summary(precip_lbs_intl_noelev) # not significant
precip_lbs_intq<-lm(LBS ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_lbs)
summary(precip_lbs_intq) 
Anova(precip_lbs_intq, type = "III") # no significant variables
precip_lbs_intq_noelev<-lm(LBS ~ poly(Intercept, degree = 2, raw = FALSE), data = precip_data_lbs) # remove interaction
summary(precip_lbs_intq_noelev) # not significant
AICc(precip_lbs_intl_noelev, precip_lbs_intq_noelev)

## Relative clutch size
avgrelcs<-mochdata%>%group_by(YEAR, ELEVATION)%>%mutate(yearavgcs = mean(CLUTCH, na.rm = TRUE))%>%
  mutate(relcs = CLUTCH - yearavgcs) %>% group_by(F.ID) %>% summarise(avgrelcs = mean(relcs, na.rm = TRUE))

# Models investigating relationship between relative clutch size and slope of reaction norms
# Temperature
temp_data_acs<-merge(temp_rn, avgrelcs, by = "F.ID")
temp_data_acs<-temp_data_acs%>%filter(!is.na(Slope))
temp_acs_ml <- lm(avgrelcs ~ Slope*ELEVATION, data = temp_data_acs)
summary(temp_acs_ml) 
Anova(temp_acs_ml, type = "III") # interaction not significant
temp_acs_ml_noelev <- lm(avgrelcs ~ Slope, data = temp_data_acs) # remove interaction
summary(temp_acs_ml_noelev) # significant!
temp_acs_mq<-lm(avgrelcs ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_acs)
summary(temp_acs_mq) 
Anova(temp_acs_mq, type = "III") 
temp_acs_mq_noelev<-lm(avgrelcs ~ poly(Slope, degree = 2, raw = FALSE), data = temp_data_acs) # remove interaction
summary(temp_acs_mq_noelev) 
AICc(temp_acs_ml_noelev, temp_acs_mq_noelev)# neither model better

# Linear term significance testing - Slope
datresample <- temp_data_acs

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Slope[k], sd = datresample$Slope_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(temp_acs_ml_noelev)$coefficients[2,1] > 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist_temp_m_acs<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist_temp_m_acs

# Create a prediction dataframe for plotting & getting the confidence intervals
# For x, use a sequence from lowest value to highest value of relative timing in the model dataset
predict.df.temp_acs_ml = data.frame(Slope = seq(min(temp_data_acs$Slope), max(temp_data_acs$Slope), 0.005))
# Predict and get confidence intervals
model.p.temp_acs_ml = predict(temp_acs_ml_noelev,newdata=predict.df.temp_acs_ml,type="response",re.form=NA,interval="confidence", se.fit = TRUE)
# Combine prediction and model datasets
predict.plot.temp_acs_ml = data.frame(predict.df.temp_acs_ml,model.p.temp_acs_ml)
# Calculate 95% confidence interval
predict.plot.temp_acs_ml =  predict.plot.temp_acs_ml %>% mutate(ci.low = fit.fit-(se.fit*1.96),ci.high=fit.fit+(se.fit*1.96))

# Plot relationship between relative timing and brood reduction for first time breeding females
temp_avgcs_slope_mainplot<-ggplot(data = temp_data_acs, aes(x = Slope, y = avgrelcs)) +
  geom_point() + labs(x = "Reaction Norm Slope", y = "Average Relative Clutch Size") +
  geom_ribbon(data=predict.plot.temp_acs_ml, aes(x=Slope,y=fit.fit,ymin=ci.low,ymax=ci.high),fill="grey65",alpha=0.45) +
  geom_line(data=predict.plot.temp_acs_ml,aes(x=Slope,y=fit.fit), linewidth = 1, color = 'blue')+ theme_cowplot()
temp_avgcs_slope_mainplot
ggsave("temp_avgcs_slope_mainplot.jpg", temp_avgcs_slope_mainplot, dpi = 600, width = 5, height = 4, bg = "white")


# Outlier removal
temp_acs_ml_noelev_orm <- lm(avgrelcs ~ Slope, data = temp_data_acs[temp_data_acs$F.ID!="B_1530",]) # remove interaction
summary(temp_acs_ml_noelev_orm) # significant

# Linear term significance testing - Slope - outlier removed
datresample <- temp_data_acs[temp_data_acs$F.ID!="B_1530",]

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Slope[k], sd = datresample$Slope_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(temp_acs_ml_noelev_orm)$coefficients[2,1] > 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist
ggsave("temp_avgcs_orm_slope_samplehist.jpg", samplehist, dpi = 600, width = 6, height = 5, bg = "white")

# Precip
precip_data_acs<-merge(precip_rn, avgrelcs, by = "F.ID")
precip_data_acs<-precip_data_acs%>%filter(!is.na(Slope))
precip_acs_ml<- lm(avgrelcs ~ Slope*ELEVATION, data = precip_data_acs)
summary(precip_acs_ml) 
Anova(precip_acs_ml, type = "III") # interaction not significant
precip_acs_ml_noelev<- lm(avgrelcs ~ Slope, data = precip_data_acs) # remove interaction
summary(precip_acs_ml_noelev) # not significant
precip_acs_mq<-lm(avgrelcs ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_acs)
summary(precip_acs_mq) 
Anova(precip_acs_mq, type = "III") # not significant
precip_acs_mq_noelev<-lm(avgrelcs ~ poly(Slope, degree = 2, raw = FALSE), data = precip_data_acs) # not significant
summary(precip_acs_mq_noelev) 
AICc(precip_acs_ml_noelev, precip_acs_mq_noelev) # neither model better

# Relative clutch size and intercept of reaction norms
# Temperature
temp_acs_intl<-lm(avgrelcs ~ Intercept*ELEVATION, data = temp_data_acs)
summary(temp_acs_intl) 
Anova(temp_acs_intl, type = "III") # interaction with elevation not significant
temp_acs_intl_noelev<-lm(avgrelcs ~ Intercept, data = temp_data_acs) # remove interaction
summary(temp_acs_intl_noelev) # significant!
temp_acs_intq<-lm(avgrelcs ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_acs)
summary(temp_acs_intq) 
Anova(temp_acs_intq, type = "III") # interaction not significant
temp_acs_intq_noelev<-lm(avgrelcs ~ poly(Intercept, degree = 2, raw = FALSE), data = temp_data_acs) # remove interaction
summary(temp_acs_intq_noelev) 
AICc(temp_acs_intl_noelev, temp_acs_intq_noelev) # neither is better

# Linear term significance testing 
datresample <- temp_data_acs

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Intercept[k], sd = datresample$Intercept_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(temp_acs_intl_noelev)$coefficients[2,1] < 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist_temp_int_acs<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist_temp_int_acs
# ggsave("temp_avgcs_int_samplehist.jpg", samplehist_temp_int_acs, dpi = 600, width = 6, height = 5, bg = "white")

# Create a prediction dataframe for plotting & getting the confidence intervals
# For x, use a sequence from lowest value to highest value of relative timing in the model dataset
predict.df.temp_acs_intl = data.frame(Intercept = seq(min(datresample$Intercept), max(datresample$Intercept), 0.01))
# Predict and get confidence intervals
model.p.temp_acs_intl = predict(temp_acs_intl_noelev,newdata=predict.df.temp_acs_intl,type="response",re.form=NA,interval="confidence", se.fit = TRUE)
# Combine prediction and model datasets
predict.plot.temp_acs_intl = data.frame(predict.df.temp_acs_intl,model.p.temp_acs_intl)
# Calculate 95% confidence interval
predict.plot.temp_acs_intl =  predict.plot.temp_acs_intl %>% mutate(ci.low = fit.fit-(se.fit*1.96),ci.high=fit.fit+(se.fit*1.96))

# Plot relationship between relative timing and brood reduction for first time breeding females
temp_avgcs_int_mainplot<-ggplot(data = temp_data_acs, aes(x = Intercept, y = avgrelcs)) +
  geom_point() + labs(x = "Reaction Norm Intercept", y = "Average Relative Clutch Size") +
  geom_ribbon(data=predict.plot.temp_acs_intl, aes(x=Intercept,y=fit.fit,ymin=ci.low,ymax=ci.high),fill="grey65",alpha=0.45) +
  geom_line(data=predict.plot.temp_acs_intl,aes(x=Intercept,y=fit.fit), linewidth = 1, color = 'blue')+ theme_cowplot()
temp_avgcs_int_mainplot
ggsave("temp_avgcs_int_mainplot.jpg", temp_avgcs_int_mainplot, dpi = 600, width = 5, height = 4, bg = "white")


# Precip
precip_acs_intl<-lm(avgrelcs ~ Intercept*ELEVATION, data = precip_data_acs)
summary(precip_acs_intl) 
Anova(precip_acs_intl, type = "III") # intercept significant
precip_acs_intl_noelev<-lm(avgrelcs ~ Intercept, data = precip_data_acs) # remove interaction
summary(precip_acs_intl_noelev) # significant!
precip_acs_intq<-lm(avgrelcs ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_acs)
summary(precip_acs_intq) 
Anova(precip_acs_intq, type = "III") # interaction not significant
precip_acs_intq_noelev<-lm(avgrelcs ~ poly(Intercept, degree = 2, raw = FALSE), data = precip_data_acs) # remove interaction
summary(precip_acs_intq_noelev) 
AICc(precip_acs_intl_noelev, precip_acs_intq_noelev) # neither better

# Linear term significance testing 
datresample <- precip_data_acs

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Intercept[k], sd = datresample$Intercept_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelcs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(precip_acs_intl_noelev)$coefficients[2,1] < 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist_precip_int_acs<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist_precip_int_acs 
# ggsave("precip_avgcs_int_samplehist.jpg", samplehist_precip_int_acs, dpi = 600, width = 6, height = 5, bg = "white")


# Create a prediction dataframe for plotting & getting the confidence intervals
# For x, use a sequence from lowest value to highest value of relative timing in the model dataset
predict.df.precip_acs_intl = data.frame(Intercept = seq(min(precip_data_acs$Intercept), max(precip_data_acs$Intercept), 0.01))
# Predict and get confidence intervals
model.p.precip_acs_intl = predict(precip_acs_intl_noelev,newdata=predict.df.precip_acs_intl,type="response",re.form=NA,interval="confidence", se.fit = TRUE)
# Combine prediction and model datasets
predict.plot.precip_acs_intl = data.frame(predict.df.precip_acs_intl,model.p.precip_acs_intl)
# Calculate 95% confidence interval
predict.plot.precip_acs_intl =  predict.plot.precip_acs_intl %>% mutate(ci.low = fit.fit-(se.fit*1.96),ci.high=fit.fit+(se.fit*1.96))

# Plot relationship between relative timing and brood reduction for first time breeding females
precip_avgcs_int_mainplot<-ggplot(data = precip_data_acs, aes(x = Intercept, y = avgrelcs)) +
  geom_point() + labs(x = "Reaction Norm Intercept", y = "Average Relative Clutch Size") +
  geom_ribbon(data=predict.plot.precip_acs_intl, aes(x=Intercept,y=fit.fit,ymin=ci.low,ymax=ci.high),fill="grey65",alpha=0.45) +
  geom_line(data=predict.plot.precip_acs_intl,aes(x=Intercept,y=fit.fit), linewidth = 1, color = 'blue')+ theme_cowplot()
precip_avgcs_int_mainplot
ggsave("precip_avgcs_int_mainplot.jpg", precip_avgcs_int_mainplot, dpi = 600, width = 5, height = 4, bg = "white")


## Relative brood size
avgrelbs<-mochdata%>%group_by(YEAR, ELEVATION)%>%mutate(yearavgbs = mean(BROOD, na.rm = TRUE))%>%
  mutate(relbs = BROOD - yearavgbs) %>% group_by(F.ID) %>% summarise(avgrelbs = mean(relbs, na.rm = TRUE))

# Models investigating relationship between relative brood size and slope of reaction norms
# Temperature
temp_data_abs<-merge(temp_rn, avgrelbs, by = "F.ID")
temp_data_abs<-temp_data_abs%>%filter(!is.na(Slope))
temp_abs_ml <- lm(avgrelbs ~ Slope*ELEVATION, data = temp_data_abs)
summary(temp_abs_ml) 
Anova(temp_abs_ml, type = "III") # interaction not significant
temp_abs_ml_noelev <- lm(avgrelbs ~ Slope, data = temp_data_abs) # remove interaction
summary(temp_abs_ml_noelev) # significant!
plot(temp_data_abs$Slope, temp_data_abs$avgrelbs)
temp_abs_mq<-lm(avgrelbs ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_abs)
summary(temp_abs_mq) 
Anova(temp_abs_mq, type = "III") # interaction not significant
temp_abs_mq_noelev<-lm(avgrelbs ~ poly(Slope, degree = 2, raw = FALSE), data = temp_data_abs) # remove interaction
summary(temp_abs_mq_noelev) 
AICc(temp_abs_ml_noelev, temp_abs_mq_noelev)# perform similarly, select simpler model

# Linear term significance testing - Slope
datresample <- temp_data_abs

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Slope[k], sd = datresample$Slope_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(temp_abs_ml_noelev)$coefficients[2,1] > 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist_temp_m_abs<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist_temp_m_abs 
# ggsave("temp_avgbs_slope_samplehist.jpg", samplehist_temp_m_abs, dpi = 600, width = 6, height = 5, bg = "white")


# Create a prediction dataframe for plotting & getting the confidence intervals
# For x, use a sequence from lowest value to highest value of relative timing in the model dataset
predict.df.temp_abs_ml = data.frame(Slope = seq(min(temp_data_abs$Slope), max(temp_data_abs$Slope), 0.005))
# Predict and get confidence intervals
model.p.temp_abs_ml = predict(temp_abs_ml_noelev,newdata=predict.df.temp_abs_ml,type="response",re.form=NA,interval="confidence", se.fit = TRUE)
# Combine prediction and model datasets
predict.plot.temp_abs_ml = data.frame(predict.df.temp_abs_ml,model.p.temp_abs_ml)
# Calculate 95% confidence interval
predict.plot.temp_abs_ml =  predict.plot.temp_abs_ml %>% mutate(ci.low = fit.fit-(se.fit*1.96),ci.high=fit.fit+(se.fit*1.96))

# Plot relationship between relative timing and brood reduction for first time breeding females
temp_avgbs_slope_mainplot<-ggplot(data = temp_data_abs, aes(x = Slope, y = avgrelbs)) +
  geom_point() + labs(x = "Reaction Norm Slope", y = "Average Relative Brood Size") +
  geom_ribbon(data=predict.plot.temp_abs_ml, aes(x=Slope,y=fit.fit,ymin=ci.low,ymax=ci.high),fill="grey65",alpha=0.45) +
  geom_line(data=predict.plot.temp_abs_ml,aes(x=Slope,y=fit.fit), linewidth = 1, color = 'blue')+ theme_cowplot()
temp_avgbs_slope_mainplot
ggsave("temp_avgbs_slope_mainplot.jpg", temp_avgbs_slope_mainplot, dpi = 600, width = 5, height = 4, bg = "white")

# Outlier removal
temp_abs_ml_noelev_orm <- lm(avgrelbs ~ Slope, data = temp_data_abs[temp_data_abs$F.ID!="B_1530",]) # remove interaction
summary(temp_abs_ml_noelev_orm) # still significant!


# Linear term significance testing - Slope - outlier removed
datresample <- temp_data_abs[temp_data_abs$F.ID!="B_1530",]

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Slope[k], sd = datresample$Slope_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(temp_abs_ml_noelev_orm)$coefficients[2,1] > 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist 
ggsave("temp_avgbs_orm_slope_samplehist.jpg", samplehist, dpi = 600, width = 6, height = 5, bg = "white")

samplefitdis <- sample(fitness_dist, 1000)


# Precip
precip_data_abs<-merge(precip_rn, avgrelbs, by = "F.ID")
precip_data_abs<-precip_data_abs%>%filter(!is.na(Slope))
precip_abs_ml<- lm(avgrelbs ~ Slope*ELEVATION, data = precip_data_abs)
summary(precip_abs_ml) 
Anova(precip_abs_ml, type = "III") # interaction not significant
precip_abs_ml_noelev<- lm(avgrelbs ~ Slope, data = precip_data_abs) # remove interaction with elevation
summary(precip_abs_ml_noelev) # not significant
precip_abs_mq<-lm(avgrelbs ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_abs)
summary(precip_abs_mq) 
Anova(precip_abs_mq, type = "III")# interaction not significant
precip_abs_mq_noelev<-lm(avgrelbs ~ poly(Slope, degree = 2, raw = FALSE), data = precip_data_abs) # remove interaction 
summary(precip_abs_mq_noelev) # not significant
AICc(precip_abs_ml_noelev, precip_abs_mq_noelev) # neither model better

# Relative brood size and intercept of reaction norms
# Temperature
temp_abs_intl<-lm(avgrelbs ~ Intercept*ELEVATION, data = temp_data_abs)
summary(temp_abs_intl) 
Anova(temp_abs_intl, type = "III") # intercept significant
temp_abs_intl_noelev<-lm(avgrelbs ~ Intercept, data = temp_data_abs) # remove interaction
summary(temp_abs_intl_noelev) # significant!
plot(temp_data_abs$Intercept, temp_data_abs$avgrelbs)
temp_abs_intq<-lm(avgrelbs ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_abs)
summary(temp_abs_intq) 
Anova(temp_abs_intq, type = "III") # interaction not significant
temp_abs_intq_noelev<-lm(avgrelbs ~ poly(Intercept, degree = 2, raw = FALSE), data = temp_data_abs) # remove interaction
summary(temp_abs_intq_noelev) 
AICc(temp_abs_intl_noelev, temp_abs_intq_noelev) # neither is better


# Linear term significance testing 
datresample <- temp_data_abs

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Intercept[k], sd = datresample$Intercept_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(temp_abs_intl_noelev)$coefficients[2,1] < 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist_temp_int_abs<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist_temp_int_abs
# ggsave("temp_avgbs_int_samplehist.jpg", samplehist_temp_int_abs, dpi = 600, width = 6, height = 5, bg = "white")


# Create a prediction dataframe for plotting & getting the confidence intervals
# For x, use a sequence from lowest value to highest value of relative timing in the model dataset
predict.df.temp_abs_intl = data.frame(Intercept = seq(min(temp_data_abs$Intercept), max(temp_data_abs$Intercept), 0.01))
# Predict and get confidence intervals
model.p.temp_abs_intl = predict(temp_abs_intl_noelev,newdata=predict.df.temp_abs_intl,type="response",re.form=NA,interval="confidence", se.fit = TRUE)
# Combine prediction and model datasets
predict.plot.temp_abs_intl = data.frame(predict.df.temp_abs_intl,model.p.temp_abs_intl)
# Calculate 95% confidence interval
predict.plot.temp_abs_intl =  predict.plot.temp_abs_intl %>% mutate(ci.low = fit.fit-(se.fit*1.96),ci.high=fit.fit+(se.fit*1.96))

# Plot relationship between relative timing and brood reduction for first time breeding females
temp_avgbs_int_mainplot<-ggplot(data = temp_data_abs, aes(x = Intercept, y = avgrelbs)) +
  geom_point() + labs(x = "Reaction Norm Intercept", y = "Average Relative Brood Size") +
  geom_ribbon(data=predict.plot.temp_abs_intl, aes(x=Intercept,y=fit.fit,ymin=ci.low,ymax=ci.high),fill="grey65",alpha=0.45) +
  geom_line(data=predict.plot.temp_abs_intl,aes(x=Intercept,y=fit.fit), linewidth = 1, color = 'blue')+ theme_cowplot()
temp_avgbs_int_mainplot
ggsave("temp_avgbs_int_mainplot.jpg", temp_avgbs_int_mainplot, dpi = 600, width = 5, height = 4, bg = "white")


# Precip
precip_abs_intl<-lm(avgrelbs ~ Intercept*ELEVATION, data = precip_data_abs)
summary(precip_abs_intl) 
Anova(precip_abs_intl, type = "III") # interaction not significant
precip_abs_intl_noelev<-lm(avgrelbs ~ Intercept, data = precip_data_abs) # remove interaction
summary(precip_abs_intl_noelev) # significant
precip_abs_intq<-lm(avgrelbs ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_abs)
summary(precip_abs_intq) 
Anova(precip_abs_intq, type = "III") # interaction not significant
precip_abs_intq_noelev<-lm(avgrelbs ~ poly(Intercept, degree = 2, raw = FALSE), data = precip_data_abs) # remove interaction
summary(precip_abs_intq_noelev) 
AICc(precip_abs_intl_noelev, precip_abs_intq_noelev) # perform similarly 

# Linear term significance testing 
datresample <- precip_data_abs

# For loop for resampling
fitness_dist <- data.frame(Slope = rep(NA, 10000), Intercept = rep(NA, 10000))
saved <- NA

for (j in 1:10000){
  saved <- data.frame(F.ID = datresample$F.ID, sample = NA)
  for (k in 1:nrow(datresample)){
    saved$sample[k]<-rnorm(1, mean = datresample$Intercept[k], sd = datresample$Intercept_se[k])
  }
  saved<-merge(saved, datresample, by = "F.ID")
  fitness_dist$Slope[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[2]
  fitness_dist$Intercept[j]<-summary(lm(avgrelbs ~ sample, data = saved))$coefficients[1]
}

p <- signif(ifelse(summary(precip_abs_intl_noelev)$coefficients[2,1] < 0, count(fitness_dist$Intercept < 0)/nrow(fitness_dist), count(fitness_dist$Intercept > 0)/nrow(fitness_dist)), 2)
# "p-value": number of slopes above 0 divided by total number

mu <- round(mean(fitness_dist$Slope), 2)
sd <- round(sd(fitness_dist$Slope), 2)

samplehist_precip_int_abs<-ggplot(data = fitness_dist, aes(x = Slope)) + geom_histogram(color="grey", fill="white") + theme_cowplot() +
  labs(title = paste("Sample mean", mu, "\u00B1", sd, ", p =", p, sep = " "), x = "Samples", y = "Count") + geom_vline(xintercept = 0, color = "red")
samplehist_precip_int_abs
# ggsave("precip_avgbs_int_samplehist.jpg", samplehist_precip_int_abs, dpi = 600, width = 6, height = 5, bg = "white")


# Create a prediction dataframe for plotting & getting the confidence intervals
# For x, use a sequence from lowest value to highest value of relative timing in the model dataset
predict.df.precip_abs_intl = data.frame(Intercept = seq(min(precip_data_abs$Intercept), max(precip_data_abs$Intercept), 0.01))
# Predict and get confidence intervals
model.p.precip_abs_intl = predict(precip_abs_intl_noelev,newdata=predict.df.precip_abs_intl,type="response",re.form=NA,interval="confidence", se.fit = TRUE)
# Combine prediction and model datasets
predict.plot.precip_abs_intl = data.frame(predict.df.precip_abs_intl,model.p.precip_abs_intl)
# Calculate 95% confidence interval
predict.plot.precip_abs_intl =  predict.plot.precip_abs_intl %>% mutate(ci.low = fit.fit-(se.fit*1.96),ci.high=fit.fit+(se.fit*1.96))

# Plot relationship between relative timing and brood reduction for first time breeding females
precip_avgbs_int_mainplot<-ggplot(data = precip_data_abs, aes(x = Intercept, y = avgrelbs)) +
  geom_point() + labs(x = "Reaction Norm Intercept", y = "Average Relative Brood Size") +
  geom_ribbon(data=predict.plot.precip_abs_intl, aes(x=Intercept,y=fit.fit,ymin=ci.low,ymax=ci.high),fill="grey65",alpha=0.45) +
  geom_line(data=predict.plot.precip_abs_intl,aes(x=Intercept,y=fit.fit), linewidth = 1, color = 'blue')+ theme_cowplot()
precip_avgbs_int_mainplot
ggsave("precip_avgbs_int_mainplot.jpg", precip_avgbs_int_mainplot, dpi = 600, width = 5, height = 4, bg = "white")


## Relative mean nestling mass
avgrelmm<-mochdata%>%group_by(YEAR, ELEVATION)%>%mutate(yearavgmm = mean(MEANMASS, na.rm = TRUE))%>%
  mutate(relmm = MEANMASS - yearavgmm) %>% group_by(F.ID) %>% summarise(avgrelmm = mean(relmm, na.rm = TRUE))

# Models investigating relationship between relative mean nestling mass and slope of reaction norms
# Temperature
temp_data_amm<-merge(temp_rn, avgrelmm, by = "F.ID")
temp_data_amm<-temp_data_amm%>%filter(!is.na(Slope))
temp_amm_ml <- lm(avgrelmm ~ Slope*ELEVATION, data = temp_data_amm)
summary(temp_amm_ml) 
Anova(temp_amm_ml, type = "III")# interaction not significant
temp_amm_ml_noelev <- lm(avgrelmm ~ Slope, data = temp_data_amm) # remove interaction
summary(temp_amm_ml_noelev) # not significant
temp_amm_mq<-lm(avgrelmm ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_amm)
summary(temp_amm_mq) 
Anova(temp_amm_mq, type = "III") # interaction not significant
temp_amm_mq_noelev<-lm(avgrelmm ~ poly(Slope, degree = 2, raw = FALSE), data = temp_data_amm) # remove interaction
summary(temp_amm_mq_noelev) # not significant
AICc(temp_amm_ml_noelev, temp_amm_mq_noelev)# perform similarly

# Precip
precip_data_amm<-merge(precip_rn, avgrelmm, by = "F.ID")
precip_data_amm<-precip_data_amm%>%filter(!is.na(Slope))
precip_amm_ml<- lm(avgrelmm ~ Slope*ELEVATION, data = precip_data_amm)
summary(precip_amm_ml) 
Anova(precip_amm_ml, type = "III")# interaction not significant
precip_amm_ml_noelev<- lm(avgrelmm ~ Slope, data = precip_data_amm) # remove interaction
summary(precip_amm_ml_noelev) # not significant
precip_amm_mq<-lm(avgrelmm ~ poly(Slope, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_amm)
summary(precip_amm_mq) 
Anova(precip_amm_mq, type = "III") # not significant
precip_amm_mq_noelev<-lm(avgrelmm ~ poly(Slope, degree = 2, raw = FALSE), data = precip_data_amm) # remove interaction
summary(precip_amm_mq_noelev) 
AICc(precip_amm_ml_noelev, precip_amm_mq_noelev) # perform relatively equally


# Relative mean mass and intercept of reaction norms
# Temperature
temp_amm_intl<-lm(avgrelmm ~ Intercept*ELEVATION, data = temp_data_amm)
summary(temp_amm_intl) 
Anova(temp_amm_intl, type = "III") # interaction not significant
temp_amm_intl_noelev<-lm(avgrelmm ~ Intercept, data = temp_data_amm) # remove interaction
summary(temp_amm_intl_noelev) # not significant
temp_amm_intq<-lm(avgrelmm ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = temp_data_amm)
summary(temp_amm_intq) 
Anova(temp_amm_intq, type = "III") # interaction not significant
temp_amm_intq_noelev<-lm(avgrelmm ~ poly(Intercept, degree = 2, raw = FALSE), data = temp_data_amm) # remove interaction
summary(temp_amm_intq_noelev) # not significant
AICc(temp_amm_intl_noelev, temp_amm_intq_noelev) # perform similarly

# Precip
precip_amm_intl<-lm(avgrelmm ~ Intercept*ELEVATION, data = precip_data_amm)
summary(precip_amm_intl) 
Anova(precip_amm_intl, type = "III") # interaction not significant
precip_amm_intl_noelev<-lm(avgrelmm ~ Intercept, data = precip_data_amm) # remove interaction
summary(precip_amm_intl_noelev) # not significant
precip_amm_intq<-lm(avgrelmm ~ poly(Intercept, degree = 2, raw = FALSE)*ELEVATION, data = precip_data_amm)
summary(precip_amm_intq) 
Anova(precip_amm_intq, type = "III") # interaction not significant
precip_amm_intq_noelev<-lm(avgrelmm ~ poly(Intercept, degree = 2, raw = FALSE), data = precip_data_amm) # remove interaction
summary(precip_amm_intq_noelev) # not significant
AICc(precip_amm_intl_noelev, precip_amm_intq_noelev) # linear is better


# Combine graphs for main text
tempmainplotsarr <- ggarrange(temp_avgcs_slope_mainplot, temp_avgcs_int_mainplot, temp_avgbs_slope_mainplot, temp_avgbs_int_mainplot, nrow = 2, ncol = 2, labels = c("      (A)", "      (B)", "      (C)","      (D)"))
tempmainplotsarr
ggsave("tempmainplotsarr.jpg", tempmainplotsarr, dpi = 600, width = 10, height = 8, bg = "white")

# Combine main graphs for main text
precipmainplotsarr <- ggarrange(precip_avgcs_int_mainplot, precip_avgbs_int_mainplot, nrow = 1, ncol = 2, labels = c("      (A)", "      (B)"))
precipmainplotsarr
ggsave("precipmainplotsarr.jpg", precipmainplotsarr, dpi = 600, width = 10, height = 4, bg = "white")

# Combine sample histograms for the supplementary material
samplehistpreciparr <- ggarrange(samplehist_precip_int_acs, samplehist_precip_int_abs, nrow = 1, ncol = 2, labels = c("   (A)", "   (B)"))
samplehistpreciparr
ggsave("samplehistpreciparr.jpg", samplehistpreciparr, dpi = 600, width = 10, height = 4, bg = "white")

# Combine sample histograms for the supplementary material
samplehisttemparr <- ggarrange(samplehist_temp_m_acs, samplehist_temp_int_acs, samplehist_temp_m_abs, samplehist_temp_int_abs, nrow = 2, ncol = 2, labels = c("   (A)", "   (B)", "   (C)", "   (D)"))
samplehisttemparr
ggsave("samplehisttemparr.jpg", samplehisttemparr, dpi = 600, width = 10, height = 8, bg = "white")


#### Repeatability ####

# First egg date
# Combined elevations - need to include elevation as random effect, otherwise repeatability estimates are inflated

repcomb <- rpt(J.FIRST.EGG ~ ELEVATION + (1|YEAR) + (1|F.ID), 
                grname = "F.ID", data = mochdata, datatype = "Gaussian", nboot = 1000, npermut = 1000)
print(repcomb)

repcombt <- rpt(J.FIRST.EGG ~ ELEVATION + (1|F.ID), 
               grname = "F.ID", data = mochdata, datatype = "Gaussian", nboot = 1000, npermut = 1000)
print(repcombt)

repcomb1_a <- rpt(J.FIRST.EGG ~ ELEVATION + precip_scale + mintemp_scale + (1|YEAR) + (precip_scale + mintemp_scale|F.ID), 
                grname = "F.ID", data = mochdata, datatype = "Gaussian", nboot = 1000, npermut = 1000)
print(repcomb1_a) 


# Clutch size
# Combined elevations
simulationOutput <- simulateResiduals(glmer(CLUTCH ~ ELEVATION*precip_scale + (1|YEAR) + (precip_scale|F.ID), data = mochdata, family = poisson), plot = T, refit = F)
testDispersion(simulationOutput)

repcomb3 <- rpt(CLUTCH ~ ELEVATION + (1|YEAR) + (1|F.ID), 
                grname = "F.ID", data = mochdata, datatype = "Gaussian", nboot = 1000, npermut = 1000)
print(repcomb3)

repcomb3t <- rpt(CLUTCH ~ ELEVATION + (1|F.ID), 
                grname = "F.ID", data = mochdata, datatype = "Gaussian", nboot = 1000, npermut = 1000)
print(repcomb3t)

repcomb4 <- rpt(CLUTCH ~ ELEVATION + precip_scale + mintemp_scale + (1|YEAR) + (precip_scale + mintemp_scale|F.ID), 
                grname = "F.ID", data = mochdata, datatype = "Gaussian", nboot = 1000, npermut = 1000)
print(repcomb4)


# Repeatability figures
# Calculate the bird-year (year 1, 2, 3, 4, etc.)
mochdata<-mochdata%>%group_by(F.ID)%>%mutate(minYear = min(YEAR))
mochdata$BirdYear <- ifelse(mochdata$YEAR == mochdata$minYear, 1, mochdata$YEAR-mochdata$minYear + 1)

# High elevation, first egg
repeatplothfe<-ggplot(mochdata[mochdata$ELEVATION == "H",]) +
  geom_line(aes(x = BirdYear, y = J.FIRST.EGG, colour = F.ID), alpha = 0.3) +
  # xlim(c(min(model_data$clim_scale), max(model_data$clim_scale))) + ylim(c(min(model_data$J.FIRST.EGG), max(model_data$J.FIRST.EGG))) +
  labs(x = "Year", y = "First Egg Date (day of year)") + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_cowplot() + theme(legend.position="none")
repeatplothfe
# ggsave("mochhighferepeatplot.jpg", repeatplothfe, width = 6, height = 5, bg = "white")


# Low elevation, first egg
repeatplotlfe<-ggplot(mochdata[mochdata$ELEVATION == "L",]) +
  geom_line(aes(x = BirdYear, y = J.FIRST.EGG, colour = F.ID), alpha = 0.3) +
  # xlim(c(min(model_data$clim_scale), max(model_data$clim_scale))) + ylim(c(min(model_data$J.FIRST.EGG), max(model_data$J.FIRST.EGG))) +
  labs(x = "Year", y = "First Egg Date (day of year)") + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_cowplot() + theme(legend.position="none")
repeatplotlfe
# ggsave("mochlowferepeatplot.jpg", repeatplotlfe, width = 6, height = 5, bg = "white")

# High elevation, clutch
repeatplothcl<-ggplot(mochdata[mochdata$ELEVATION == "H",]) +
  geom_line(aes(x = BirdYear, y = CLUTCH, colour = F.ID), alpha = 0.3) +
  # xlim(c(min(model_data$clim_scale), max(model_data$clim_scale))) + ylim(c(min(model_data$J.FIRST.EGG), max(model_data$J.FIRST.EGG))) +
  labs(x = "Year", y = "Clutch Size") + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_cowplot() + theme(legend.position="none")
repeatplothcl
# ggsave("mochhighclrepeatplot.jpg", repeatplothcl, width = 6, height = 5, bg = "white")

# Low elevation, clutch
repeatplotlcl<-ggplot(mochdata[mochdata$ELEVATION == "L",]) +
  geom_line(aes(x = BirdYear, y = CLUTCH, colour = F.ID), alpha = 0.3) +
  # xlim(c(min(model_data$clim_scale), max(model_data$clim_scale))) + ylim(c(min(model_data$J.FIRST.EGG), max(model_data$J.FIRST.EGG))) +
  labs(x = "Year", y = "Clutch Size") + scale_x_continuous(breaks = c(1,2,3,4,5,6,7,8)) +
  theme_cowplot() + theme(legend.position="none")
repeatplotlcl
# ggsave("mochlowclrepeatplot.jpg", repeatplotlcl, width = 6, height = 5, bg = "white")


# Arrange plots
repeatarr<-ggarrange(repeatplothfe, repeatplotlfe, repeatplothcl, repeatplotlcl, nrow = 2, ncol = 2, labels = c("(A) High Elevation, First Egg Date", "(B) Low Elevation, First Egg Date", "(C) High Elevation, Clutch Size", "(D) Low Elevation, Clutch Size"), align = "hv", label.x = c(-0.1, -0.1, -0.07, -0.07))
repeatarr
ggsave("repeatplotarr.jpg", repeatarr, width = 12, height = 10, bg = "white")
