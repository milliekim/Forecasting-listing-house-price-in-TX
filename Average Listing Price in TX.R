# This file loads Average Listing House Price in Texas in a given market during the 
# specified month.
# We'll see how accurately the SARIMA model predicts for the last four observations.

#Clear all variables in my workspace
rm(list=ls())

#load packages
library(fpp3)
library(urca)
library(tidyverse)

#load data
temp = readr::read_csv("AVELISPRITX.csv")
temp
view(temp)

#check for NAS
temp%>%summarise(count=sum(is.na(AVELISPRITX)))


# convert tibble into tsibble and data is year-month-day
AveragePrice=temp%>%rename(AvePrice=AVELISPRITX)%>%
  mutate(DATE=yearmonth(DATE))%>%as_tsibble(index=DATE)
view(AveragePrice)

#Create training data from 7/1/2016 - 12/1/2021
#Hold out period: last four observations 9/1/2021 - 12/1/2021
TrData=AveragePrice%>%filter_index("2016 Jul"~"2021 Aug")
view(TrData)

# Visualize it
TrData%>%autoplot(AvePrice)
#1. potential non-linear trends
#2. definite changes in volatility
#3. this data is not seasonally adjusted

# Average Listing House Price obtained from: 
#https://fred.stlouisfed.org/series/AVELISPRITX


# visualize it with ggplot
# For seasonality
TrData%>%gg_subseries(AvePrice)

# gg_season
TrData%>%gg_season(AvePrice)


######Apply Box-Cox transformation###########
#The box-Cox transformation  transforms the data so that it closely resembles a normal
#distribution.In many statistical techniques, we assume that the errors are normally
#distrubuted. This assumption allows us to construct confidence intervals and conduct 
#hypothesis tests.

lambda=AveragePrice%>%features(AvePrice,features=guerrero)%>%pull(lambda_guerrero)
lambda
#lambda -0.89992 is below from 0 and far from 1. This gives the benefit to use 
#logarithmic transformation

###Comparing the Box-Cox and Log transformation plots

TrData%>%autoplot(box_cox(AvePrice,lambda))
TrData%>%autoplot(log(AvePrice))

# data shows potential time trend, so we should include one in our testing.

##Step1.Box-Jenkins Methodology###

# First,check Covariance Stationary

#The condition of stationarity
#1.mean has to be constant
#2.standard deviation(fluctuation) of time series is constant over all time
#3.No seasonality

# This data does not appear to be stationary. need to apply for unit root test.

#####Unit root test: ADF Test#####
#Checking for Seasonal differencing of AvePrice and Inflation
TrData%>%features(log(AvePrice), unitroot_nsdiffs)
#D=1, the test returns a 1; hence there is seasonal root and need seasonal differencing
#is needed.

### Now let's test with ADF test#########
## To run the test, convert the data series into a format accepted by the test. Elements
#of a tsibble are easily converted to the time series structure needed by our test

# unitroot test: ADF test on AvePrice for non-seasonal differencing
# Since ADF test doesn't allow us to run the unit root test in tsibble environment, so
# we need to extract the desired series from our tsibble and convert it to 
# a time-series object.
testData=as.ts(select(TrData,AvePrice))
summary(ur.df(log(testData),type="trend", lags=20,selectlags="AIC"))
# The ADF test result: test statistic is -3.0937 which is greater than 
# critical value is -3.45 with 5% test size, hence we fail to reject 
#the null hypothesis that has a unit root which means it's not a stationary.

#we need to difference

summary(ur.df(diff(log(testData),type="drift",lags=20, selectlags="AIC")))
#The ADF test result: test statistic is -4.4078 which is less than the critical value
#-1.95 with 5% test size, so we reject the null hypothesis and the data is stationary
#IN conclusion, D=1 and d=1


###Another way of testing for non-seasonal differencing using KPSS test#####
TrData%>%features(AvePrice,features=unitroot_kpss)
#P-value 0.01 is less than 0.05, we reject the null hypothesis that is trend
#stationary,and differencing is needed.  d=1


## we need to difference
TrData%>%features(difference(log(AvePrice)),unitroot_kpss)
#The test result: p-value 0.1 is greater than 0.05, we fail to reject the 
#null hypothesis,hence the data series is trend stationary


#######Step2:Box-Jenkins Methodology############

#### Now for selecting SARIMA model.############

#Make initial guesses about the most appropriate model
#visualize of ACF,PACF plots:
TrData%>%gg_tsdisplay(difference(log(AvePrice)),lag_max=48,
                      plot_type="partial")
#1.p and P both equal 1 (one significant positive spikes in ACF and PACF plots)
#2.d and D both equal 1(first difference and seasonal difference)
#3.q and Q both equal 0 (no significant negative spikes in PACF and one in ACF)

######Step 3: Box Jenkins Methodology######################
#Estimate the candidate models, and testing of various ARIMA/SARIMA models:
report(TrData%>% model(ARIMA(log(AvePrice)~pdq(1,1,0)+PDQ(1,1,0))))
#AIC=-285.62   BIC=-279.95
report(TrData%>% model(ARIMA(log(AvePrice)~pdq(2,1,0)+PDQ(1,1,0))))
#AIC=-283.76    BIC=-276.19

report(TrData%>% model(ARIMA(log(AvePrice))))
# R recommends ARIMA(1,1,0)(0,1,1)and AIC=-287.33     BIC=-281.65

report(TrData%>% model(ARIMA(log(AvePrice)~pdq(1,1,0)+PDQ(0,1,0))))
#AIC=-274.3  BIC=-270.52

report(TrData%>% model(ARIMA(log(AvePrice)~pdq(1,1,1)+PDQ(0,1,0))))
#AIC=-272.35   BIC=-266.67

report(TrData%>% model(ARIMA(log(AvePrice)~pdq(2,1,0)+PDQ(0,1,0))))
#AIC=-272.36   BIC=-266.68

report(TrData%>% model(ARIMA(log(AvePrice)~pdq(1,1,1)+PDQ(0,1,0))))
#AIC=-272.35   BIC=-266.67

#Estimate our final candidate with lowest AIC and BIC
EstMODEL=TrData%>%model(ARIMA(log(AvePrice)~pdq(2,1,0)+PDQ(1,1,0)))
gg_tsresiduals(EstMODEL, lag=24)

report(EstMODEL)
FinModel=EstMODEL%>%forecast(h=4)
FinModel%>%as.data.frame(select(Date, .mean))
FinModel%>%autoplot(AveragePrice%>%filter_index("2019 Jul"~"2021 Dec"))

###Analysis:
# With the SARIM model, the prediction was not very accurately predicted. The SARIMA
# model works well under normal circumstances,but here it seems like the covid was 
# affecting on the prediction. Overall, the price decreased after Jan 2022.
