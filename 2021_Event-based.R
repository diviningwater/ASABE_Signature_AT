##
## HYdrologic Signatures - event based

library(smoother)
library(lubridate)
library(insurancerating)
library(classInt)
library(dplyr)
library(ggplot2)
library(dplyr)
library(zoo)

rm(list=ls(all=TRUE))

setwd("G:/My Drive/Tidewater NUE/Data Dump/Data2021")
Rain <- read.csv("G:/Shared drives/Tidewater-HydSignatures/Rainfall/hourly2021.csv")
Rain$Hour <- round(Rain$Time * 24)
Rain$date<- as.POSIXct(paste(Rain$Year, Rain$Month, Rain$Day, Rain$Hour, sep = "-"), format = "%Y-%m-%d-%H")
Rain <- Rain[,c(8,6)]
Irrigation <- read.csv("G:/Shared drives/Tidewater-HydSignatures/Rainfall/Irrigation2021.csv")
## 1. Read in soil moisture data and do basic cleaning/formatting
#######################################################################
StnIds <- data.frame(
  ID=c("002060E6","00206FAD","00206F86","00206F61","0120A7B7",
       "01209C43","002060AE"),
  Name=c("Corn Full Irrigation","Corn No Irrigation","Corn Precision","Cotton Full Irrigation",
         "Cotton No Irrigation",
         "Cotton Precision","MetSt"))
DepthIds <- c("4 in","8 in","12 in","16 in","20 in","24 in","28 in","32 in","36 in")

sensor.data <- list()
for (s in 1:7){
  sensor.data[[s]] <- read.csv(file=list.files(pattern="csv")[s+157])
}
soil_moisture <- list()

for (s in 1:6){
  
  soil_moisture[[s]] <- sensor.data[[s]][,c(2,6:14)]  
}
#Let's keep only rows with SM value >=1
for (s in 1: 6){
  data <- soil_moisture[[s]]
  data <- data[data[,2] >= 1,]
  data <- na.omit(data)
  soil_moisture[[s]] <- data
}
## Convert datecolumns to posix
for (s in 1:length(soil_moisture)){  
  soil_moisture[[s]][,1] <- as.POSIXct(soil_moisture[[s]][,1],
                                       format="%Y-%m-%d %H:%M:%OS")
}

## 2. Functions to convert hourly SM data to smoothed value and generate plots
#####################################################################################

sensor.plot.full <- function(s, d){ 
  # s = probe id (1 through 6)
  # d = sensor depth (1 through 9)
  plot.title <<- paste(StnIds$Name[s]," ",DepthIds[d])
  data <- soil_moisture[[s]]
  #data <<- data[(data$date> "2021-08-01" & data$date< "2021-09-1"),]
  plot(data[,1],data[,(d+1)],type="l",main = plot.title,xlab="Date",ylab="Soil moisture (%)")
  #dev.off()
}



data.smooth <- function(s, t){
  # s is probe id
  # t is time interval for smoothing (e.g. t=4 leads to 4 hour moving average)
  # returns smoothed data frame (same format as soil_moisture) with all depths
  data <- soil_moisture[[s]]
  data.sm <- data
  X <- as.numeric(data$date)/3600 # as.numeric provides number of seconds since 1970
  for (d in 1:9){
    Y <- data[,(d+1)]
    Y.smth <- as.vector(stats::filter(Y,rep(1/t,t),circular=TRUE))
    data.sm[,(d+1)] <- Y.smth
  }
  return(data.sm)
}

test <- data.smooth(1,2)
# four-hour moving average for all sensors/depths
sm_smooth <- list()
for (s in 1:length(soil_moisture)){
  data.sm <- data.smooth(s,4)
  sm_smooth[[s]] <- data.sm
}

smooth.plot<- function(s,d, incl.orig=TRUE){
  # incl.orig is TRUE/FALSE - should original data be included in plot for comparison?
  plot.title <- paste(StnIds$Name[s]," ",DepthIds[d])
  data <- soil_moisture[[s]]
  data.sm <- sm_smooth[[s]]
  if (incl.orig == TRUE){
    plot(data[,1],data[,(d+1)],type="l",main = plot.title,xlab="Date",ylab="Soil moisture (%)", 
         col="grey")
    lines(data.sm[,1],data.sm[,(d+1)])
    legend("bottomright",legend=c("Smoothed","Original"),col=c("black","gray"),lty=1)
  }else{
    plot.title <- paste(StnIds$Name[s]," ",DepthIds[d]," smoothed")
    plot(data.sm[,1],data.sm[,(d+1)],type="l",main = plot.title,xlab="Date",ylab="Soil moisture (%)")
  }
}
#smooth.plot(1,3)


## 3. 4-zone classification 
##################################################################################

## First, check how heavy-tailed (skewed) the distribution of dY values are. 
## Heavy tails imply the need for head/tails classification
## If fewer than 40% of observations are beyond mean value, assume heavy tailed

# Function to calculate percentage of dY values greater than mean (for positive) and
#   less than mean (for negative)
dist.tail <- function(s,d){
  ## Determine if positive and negative dY values are heavy-tailed
  data <- soil_moisture[[s]]
  ts <- data.frame(date=data$date, Y=data[,(d+1)])
  dY <- diff(ts$Y)/diff(as.numeric(ts$date)/3600)# the derivative 
  ts$dY = c(dY[1],dY)  #making X and dY of same length
  Negative=ts$dY[ts$dY<0]
  Positive=ts$dY[ts$dY>0]
  pos.tail <- sum(Positive > mean(Positive))/length(Positive)
  neg.tail <- sum(Negative < mean(Negative))/length(Negative)
  return (list(pos.tail, neg.tail))
}

# Create matrix with percentage of dY values beyond mean
# Skewed distributions have a small percentage of values beyond mean
# Lower values means more skewed/heavy tailed
pos.tails <- matrix(data=NA, nrow=9, ncol=6)
neg.tails <- matrix(data=NA, nrow=9, ncol=6)

for (s in 1:6){
  for (d in 1:9){
    tails <- dist.tail(s,d)
    pos.tails[d,s] <- tails[[1]]
    neg.tails[d,s] <- tails[[2]]
    rm(tails)
  }
}

# Averages for each depth interval across sensors
apply(pos.tails, 1, mean) # lower percentages in positive values - more heavy-tailed
apply(neg.tails, 1, mean) # high percentages - less heavy tailed
# smaller numbers for positive values - distribution of dY is more heavy-tailed
# negative values are still below 0.4 - use head/tails classification


## Test below on single sensor and depth (need to functionalize)
#s <- 1
#d <- 3
# t is time interval for smoothing (e.g. t=4 leads to 4 hour moving average)
# returns smoothed data frame (same format as soil_moisture) with all depths

Main <- function(s,d,t){
  if (t>0){
    data <- data.smooth(s,t)
  }else{
    data <- soil_moisture[[s]]
  }
  class.data <- data.frame(date=data$date, Y=data[,(d+1)], dY=NA)
  
  #Add irrigation events to full and precision dataframe
  values_added <- FALSE
  Rain$Type <- "Rainfall"
  if (s %in% c(3) && !values_added) {
    Rain[4284,2]<-Rain[4284,2]+12.7
    Rain[4692,2]<-Rain[4692,2]+12.7
    Rain[4932,2]<-Rain[4932,2]+12.7
    values_added <- TRUE
    Rain[4284,2]  <- "Irrigation"
    Rain[4692,2]  <- "Irrigation"
    Rain[4932,2]  <- "Irrigation"
  }
  if (s %in% c(1,4) && !values_added) {
    Rain[4284,2]<-Rain[4284,2]+17.78
    Rain[4692,2]<-Rain[4692,2]+19.0
    Rain[4932,2]<-Rain[4932,2]+19.0
    values_added <- TRUE
    Rain[4284,2]  <- "Irrigation"
    Rain[4692,2]  <- "Irrigation"
    Rain[4932,2]  <- "Irrigation"
  }
  if (s==6 && !values_added) {
    Rain[4284,2]<-Rain[4284,2]+10.16
    Rain[4692,2]<-Rain[4692,2]+10.16
    Rain[4932,2]<-Rain[4932,2]+10.16
    values_added <- TRUE
    Rain[4284,2]  <- "Irrigation"
    Rain[4692,2]  <- "Irrigation"
    Rain[4932,2]  <- "Irrigation"
  }
  
  ## Find first and second derivatives
  dY <- diff(class.data$Y)/diff(as.numeric(class.data$date)/3600)# the derivative 
  class.data$dY = c(dY[1],dY)  #making X and dY of same length
  d2y <- diff(class.data$dY)/diff(as.numeric(class.data$date)/3600)
  class.data$d2Y = c(d2y[1],d2y)  #making X and dY of same length
  
  #hist(class.data$dY)
  summary(class.data$dY) # max values are high (e.g., 14, 9), but look valid from VWC time series
  #sensor.plot.full(s,d)
  ## Identify negative inflection points (dY <0 and d2Y goes from neg to pos)
  class.data$infl <- 0
  for (r in 1:(nrow(class.data)-1)){
    if (class.data$dY[r] < 0){
      if (class.data$d2Y[r] < 0 & class.data$d2Y[r+1] > 0){
        class.data$infl[r] <-1
      }
    }
  }
  
  ## Extract different classes based on dY heads/tails classification
  Negative=class.data$dY[class.data$dY<0]
  Positive=class.data$dY[class.data$dY>0]
  
  ## First do positive - data generally more heavy-tailed
  Levels.ht.pos <- classIntervals(Positive, style="headtails", thr=0.30)
  ## Choose uppermost level of head as Inf events
  brks.pos <- c(0,tail(Levels.ht.pos$brks,2)[1], max(Positive) )
  perc.head <- sum(Positive > brks.pos[2])/length(Positive)
  
  ## Now negative - capture at least as many head records as positive
  Levels.ht.neg <- classIntervals((-1*Negative), style = "headtails", thr=0.35)
  neg.brks <- -1*Levels.ht.neg$brks   #converting back to negative values.
  neg.brks <- sort(neg.brks)
  neg.class <- cut(Negative,breaks=neg.brks, include.lowest = TRUE)
  cum.perc <- data.frame(id=seq(1:length(levels(neg.class))),
                         class=levels(neg.class),
                         count=as.numeric(summary(neg.class)) )
  cum.perc$cum.perc <-NA
  cum.perc$cum.perc[1] <- cum.perc$count[1]/sum(cum.perc$count)
  for (c in 2:nrow(cum.perc)){
    cum.perc$cum.perc[c] <- cum.perc$cum.perc[c-1] + 
      cum.perc$count[c]/sum(cum.perc$count)
  }
  ## Select the minimum break that results in a cum.perc at least as high as perc.head
  neg.brk <- neg.brks[1+min(cum.perc$id[cum.perc$cum.perc > perc.head])] #1 for 2nd value
  brks.neg <- c(min(Negative),neg.brk)
  final.breaks <- c(min(Negative), neg.brk,brks.pos)
  
  #Extracting individual breakpoints to modify
  class.data$class<-cut(class.data$dY,breaks=final.breaks, include.lowest = TRUE,
                        labels=c("GD","ET","RW","INF"))
  
  class.data$class <- as.character(class.data$class)
  class.data$cum_Levels <- NA
  class.data$frequency <-NA
  class.data$cum_Levels[1]<-1
  
  for(i in 1:(nrow(class.data)-1)){
    if(class.data$class[i+1]==class.data$class[i]){
      class.data$cum_Levels[i+1]<-class.data$cum_Levels[i]+1
    }else {
      class.data$cum_Levels[i+1]<-1
    }
  }
  
  #calculating frequency of events
  a<-1
  for(i in 1:(nrow(class.data)-1)){
    if(class.data$class[i+1]!=class.data$class[i]){#wont work for the last row
      freq<-class.data$cum_Levels[i]
      class.data$frequency[a:i]<-rep(freq,times=i)
      a=i+1
    }
  }
  
  #For the last levels with NA values, replace them with the end frequency values
  class.data$frequency[is.na(class.data$frequency)]<-tail(class.data$cum_Levels,1)
  
  
  ## What percentage of inflection points occur in ET zones?
  ## What percentage of inflection points occur in FC zones?
  sum(class.data$infl == 1 & class.data$class == "ET")/sum(class.data$infl==1) #94.5%, 380 points
  sum(class.data$infl == 1 & class.data$class == "GD")/
    sum(class.data$infl==1) # 5.5%, 22 points
  
  ## Identify points where gravity transitions to ET: 
  class.data$valid_trans <- 0
  for (r in 1:(nrow(class.data)-1)){
    if(class.data$class[r] == "GD"){
      if(class.data$cum_Levels[r] == class.data$frequency[r]){ # last point of sequence
        if (class.data$class[r+1] == "ET"){
          class.data$valid_trans[r] <- 1
        }
      }
    }
  }
  
  ## Try-let's Identify points where Inf transitions to ET:
  ##This is the mostly case in deeper sensors
  ##We name this new column ad valid_trans2
  class.data$valid_trans2 <- 0
  for (r in 1:(nrow(class.data)-1)){
    if(class.data$class[r] == "INF"){
      if(class.data$cum_Levels[r] == class.data$frequency[r]){ # last point of sequence
        if (class.data$class[r+1] == "ET"){
          class.data$valid_trans2[r] <- 1
        }
      }
    }
  }
  
  
  sum(class.data$valid_trans2 == 1) # 19 points
  
  valid.data <- class.data[class.data$valid_trans==1,]
  #Lets create new dataframe valid.data2 with Inf to ET transitions
  
  valid.data2 <- class.data[class.data$valid_trans2==1,]
  #FC <- mean(valid.data$Y)
  
  #Moved this section up before classifiend FC because one RW point after INF
  #is creating issue in FC determination
  
  #below code to replace all RW points after Inf to Inf
  for (r in 1:(nrow(class.data)-2)){  #2 because 1 one giving issue while doing r+i
    if(class.data$class[r] == "INF"){
      n <- class.data$frequency[r+1]
      for (i in 1:n){
        if (class.data$dY[r+i]>0){
          class.data$class[r+i] <- "INF"  #Let's updatethe whole "class" column
        }else (class.data$class[r+i]<-class.data$class[r+i]) #not just "Events" column for plot
      }
    }
  }
  #Similarly grvity drainage is possible only after INF, so consider replacing RW events 
  #before GD to INF?
  ###There are lots of points where there is ET-GD-ET transitions, we will focus on
  ##Inf-GD-ET transition.
  ##We name this new column as valid_trans3
  class.data$valid_trans3 <- 0
  for (r in 1:(nrow(class.data)-1)){
    #if(class.data$class[r] == "RW" || class.data$class[r] == "INF" ){
    if(class.data$class[r] == "INF" ){
      if(class.data$class[r+1] == "GD"){  
        if(isTRUE(class.data$class[class.data$frequency[r+1]+r+1] == "ET")){
          class.data$valid_trans3[class.data$frequency[r+1]+r] <- 1
        }
      }
    }
  }  #check (6,3,0) (6,4,0) (6,5,0) )(3,6,0) (3,1,0) (1,4,0) (5,3,0)
  #Next only filter those with GD frequency>5?
  
  valid.data3 <- class.data[class.data$valid_trans3==1,]
  
  #Added October23rd.
  # Create a new column in valid.data3 to store the results
  if (nrow(valid.data3) > 0) {
    valid.data3$NonZeroRainWithin3Days <- FALSE
    
    #Check if there was rainfall/irrigation at least 3 days before FC date
    # Loop through each row in valid.data3 and check for non-zero Rain_inc value within 3 days before the date
    for (i in 1:nrow(valid.data3)) {
      date_to_check <- valid.data3$date[i]
      previous_dates <- Rain$date[Rain$date<= date_to_check & Rain$date>= (date_to_check - 3*24*60*60)]
      for (j in 1:length(previous_dates)) {
        if (any(Rain$Rain_inc[Rain$date== previous_dates[j]] > 0)) {
          valid.data3$NonZeroRainWithin3Days[i] <- TRUE
          break
        }
      }
    }
    
    # Filter valid.data3 based on the NonZeroRainWithin3Days column
    valid.data3 <- valid.data3[valid.data3$NonZeroRainWithin3Days, ]
  }else{
    print("Valid.data3 has no rows")
  }
  
  
  FC <- mean(valid.data3$Y)
  ## Possible characterizations of FC:
  # valid transitions from gravity to ET
  # inflection points within GD values
  # compare VWC readings for all data, valid transition points, and inflection points
  
  summary(class.data$Y) 
  summary(class.data$Y[class.data$valid_trans ==1])
  summary(class.data$Y[class.data$infl == 1 & class.data$class == "GD"])
  
  ## Compare 95% values
  quantile(class.data$Y,c(0.05, 0.95)) # 12.3 to 31.8%
  quantile(class.data$Y[class.data$valid_trans ==1],c(0.05, 0.95)) # 23.2 to 30.9%
  quantile(class.data$Y[class.data$infl == 1 & 
                          class.data$class == "GD"],c(0.05, 0.95)) #24.5 to 31.6%
  
  #Inf only has 1-2 data point to draw a continuous line. 
  #We will consider all RW points before and after Inf as Inf
  #Save new class in column Events
  class.data$Events <- class.data$class
  for (r in 1:(nrow(class.data)-1)){
    if(class.data$class[r] == "INF"){
      n <- class.data$frequency[r-1]
      for (i in 1:n){
        if (class.data$dY[r-i]>0){
          class.data$Events[r-i] <- "INF"
        }#else (class.data$Events[r-i]<-class.data$Events[r-i])
      }
    }
  }
  
  
  #Lastly, if ET transitions to Inf,the ET line connects upward with Inf.
  #Lets, replace every ET point just before Inf to Inf.
  for (r in 1:(nrow(class.data)-1)){
    if(class.data$Events[r] == "INF"){
      if (isTRUE(class.data$Events[r-1]=="ET")){
        class.data$Events[r-1] <- "INF"
      }
    }
  }
  #Let's do samething so that ET wont seem rising in the plot during RW.
  for (r in 1:(nrow(class.data)-1)){
    if(class.data$Events[r] == "RW"){
      if (isTRUE(class.data$Events[r-1]=="ET")){
        class.data$Events[r-1] <- "RW"
      }
    }
  }
  #Added in June 2023
  #Gravity drainage seems transitioning to RW sometimes
  
  #Added later for PWP
  sum_derivatives <- function(x) {
    sum(x)
  }
  max_derivatives <- function(x) {
    max(x)
  }
  Inf_count <- function(x) {
    sum(x == "INF", na.rm=TRUE)
  }
  GD_count <- function(x) {
    sum(x == "GD", na.rm=TRUE)
  }
  
  RW_perc <- function(x) {
    RW_count <- sum(x == "RW", na.rm = TRUE)
    total_count <- length(x)
    perc <- RW_count / total_count * 100
    perc
  }
  
  
  window_size <- 26  #To test stable soil moisture/negligible reduction 
  step_size <- 26 #To evaluate non-overlapping window
  
  #Instead of high, we will use Inf event and make sure there's none during
  #that window of 26 hrs; 6am-6pm for 2 days
  #(no Inf, no GD, RW <25%?)
  
  class.data$Hour <- hour(class.data$date)
  PEL.data <- class.data[class.data$Hour >= 6 & class.data$Hour <= 18, ]
  PEL.data$rf_count <- rollapplyr(PEL.data$class, window_size, Inf_count, 
                                  align="left",fill=NA, by=step_size)
  PEL.data$gd_count <- rollapplyr(PEL.data$class, window_size, Inf_count, 
                                  align="left",fill=NA, by=step_size)  
  
  PEL.data$sums <- rollapplyr(PEL.data$dY, window_size, sum_derivatives, 
                              align="left",fill=NA, by=step_size)
  PEL.data$RW_perc <- rollapplyr(PEL.data$class, window_size, 
                                 RW_perc, align = "left", fill = NA, by=step_size)
  #class.data$high <- rollapplyr(class.data$dY, window_size+10, max_derivatives, align="left",fill=NA)
  
  
  # Replace NA values with the same non-NA value throughout the window 
  #(fromLast ensures filling by most recent NonNA)
  PEL.data$rf_count <- na.locf(PEL.data$rf_count, fromLast = FALSE, na.rm = FALSE)
  PEL.data$gd_count <- na.locf(PEL.data$gd_count, fromLast = FALSE, na.rm = FALSE)
  PEL.data$sums    <- na.locf(PEL.data$sums, fromLast = FALSE, na.rm = FALSE)
  PEL.data$RW_perc    <- na.locf(PEL.data$RW_perc, fromLast = FALSE, na.rm = FALSE)
  PEL.data$abs_diffs <- abs(PEL.data$sums)
  
  PEL.data <- PEL.data[PEL.data$abs_diffs<0.1 & PEL.data$rf_count==0 &
                         PEL.data$gd_count==0 &PEL.data$RW_perc<25,]
  #PEL.data <- PEL.data[PEL.data$abs_diffs<0.05 & class.data$high<0.03,]
  PEL.data <- PEL.data[PEL.data$Y < mean(class.data$Y),]  #less than average wc
  PEL.data <- PEL.data[PEL.data$Y < min(valid.data3$Y),]  #less than FC estimate
  PEL.data <- na.omit(PEL.data)
  #Added October23rd.
  #Check if there was rainfall/irrigation at least 5 days before PWP Date
  # Loop through each row in PWP.data and check for zero Rain_inc value 
  #within 5 days before the date
  # Create a new column in PWP.data to store the results
  #I will look for rainfall days only for top 5 sensors??& s%in% c(1,2,3,4,5)
  if (nrow(PEL.data) > 0) {
    PEL.data$ZeroRainWithin5Days <- FALSE
    
    for (i in 1:nrow(PEL.data)) {
      date_to_check <- PEL.data$date[i] #check until 1 day before because it may take 1 day to arrive to deep sensors
      previous_dates <- Rain$date[Rain$date<= date_to_check - 1*24*60*60 & 
                                    Rain$date>= (date_to_check - 5*24*60*60)]
      all_zeros <- TRUE
      for (j in 1:length(previous_dates)) {
        if (any(Rain$Rain_inc[Rain$date== previous_dates[j]] >1)) {
          all_zeros <- FALSE
          break
        }
      }
      PEL.data$ZeroRainWithin5Days[i] <- all_zeros
    }
    # Filter valid.data3 based on the NonZeroRainWithin3Days column
    PEL.data <- PEL.data[PEL.data$ZeroRainWithin5Days, ]
  }else{
    print("PEL.data has no rows")
  }
  
  PEL <- mean(PEL.data$Y)
  
  ## plot results
  #plot.title <- paste(StnIds$Name[s])
  #depth <- paste("Volumetric soil water content (%)",DepthIds[d],sep=", ")
  depth <- paste("VWC (%)")
  #min_Y <- min(PEL.data$Y)
  # Group the data by dateand select a representative point for each day
  PEL.data$date_only <- as.Date(PEL.data$date)
  representative_points <- PEL.data %>%
    group_by(date_only) %>%
    slice(1)  # Select the first row for each day
  
  # Convert the datecolumns to the same format
  
  rain.data <- left_join(class.data, Rain, by = "date")
  rain.data <- rain.data[,c(1,14,15)]
  rain.data$Events <- NA
  
  
  ###########################################added in OCt2023 for rain
  y_lower_limit <- min(class.data$Y) * 0.8
  y_upper_limit <- max(class.data$Y) * 1.1
  rain.data$Rain_inc <- as.numeric(rain.data$Rain_inc)
  rain.data$Rain_inc1 <-(rain.data$Rain_inc/2.5)+y_lower_limit
  #To set secndary y axis maximum limit
  rain_ymax <- (max(rain.data$Rain_inc1)-y_lower_limit)*2.5
  
  if (nrow(valid.data3)>=1 && nrow(PEL.data)>=1){
    p <- ggplot(class.data, aes(date, Y,col=Events, group=1))+
      xlab("Time (2021 Growing season)") + ylab(depth)+
      geom_line(size=1)+
      geom_point(data=valid.data3,aes(date,Y),color='blue',size=1.5)+ 
      geom_point(data=representative_points,aes(date,Y),color='red',size=1.5)+
      #one more line for rainfall
      geom_bar(data = rain.data, aes(date, Rain_inc1, fill = Type), 
               alpha = 0.5, stat = "identity", width = 20000) +
      scale_fill_manual(values = c("Rainfall" = "grey", "Irrigation" = "red")) +
      scale_y_continuous(expand=c(0,0), sec.axis = sec_axis(~(.-y_lower_limit)*2.5, 
                                                            name = "Rain or Irrigation (mm/d)",
                                                            breaks = seq(0, 60,by=15))
      ) +
      coord_cartesian(ylim = c(y_lower_limit, y_upper_limit)) +
      #valid.data2 above for viewing Inf-et transitions point
      #theme(plot.title = element_text(hjust = 0.5))+
      #ggtitle(plot.title) +
      theme_bw()+
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        legend.position = c(0.90, 0.33),
        legend.title=element_text(size=10),
        axis.title.x = element_text(size =12),  # Adjust the x-axis label size
        axis.title.y = element_text(size =12),
        axis.text.x = element_text(size =12),  # Adjust the x-axis tick text size
        axis.text.y = element_text(size =12),  
        legend.text=element_text(size=10))+
      guides(col=guide_legend(ncol=2))+  #bcoz col in aesthetic, not fill
      scale_x_datetime(date_labels = "%b", breaks="1 month")
    # Format the x-axis labels to display only the month
    
    
  } else if (nrow(valid.data3) < 1 && nrow(PEL.data)>=1){
    p <- ggplot(class.data, aes(date, Y,col=Events, group=1))+
      xlab("Time (2021 Growing season)") + ylab(depth)+
      geom_line(size=1)+
      geom_point(data=representative_points,aes(date,Y),color='red',size=1.5)+ 
      #theme(plot.title = element_text(hjust = 0.5))+
      #one more line for rainfall
      geom_bar(data = rain.data, aes(date, Rain_inc1, fill = Type), 
               alpha = 0.5, stat = "identity", width = 20000) +
      scale_fill_manual(values = c("Rainfall" = "grey", "Irrigation" = "red")) +
      scale_y_continuous(expand=c(0,0), sec.axis = sec_axis(~(.-y_lower_limit)*2.5, 
                                                            name = "Rain or Irrigation (mm/d)",
                                                            breaks = seq(0, 60,by=15))
      ) +
      coord_cartesian(ylim = c(y_lower_limit, y_upper_limit)) +
      #ggtitle(plot.title) +
      theme_bw()+
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        
        #legend.position = "none",
        legend.position = c(0.90, 0.43),
        legend.title=element_text(size=10),
        axis.title.x = element_text(size =12),  # Adjust the x-axis label size
        axis.title.y = element_text(size =12),
        axis.text.x = element_text(size =12),  # Adjust the x-axis tick text size
        axis.text.y = element_text(size =12),  
        legend.text=element_text(size=10))+
      guides(col=guide_legend(ncol=2))+
      scale_x_datetime(date_labels = "%b", breaks="1 month")
    # Format the x-axis labels to display only the month
    
  } else if (nrow(valid.data3) >= 1 && nrow(PEL.data) < 1){
    p <- ggplot(class.data, aes(date, Y,col=Events, group=1))+
      xlab("Time (2021 Growing season)") + ylab(depth)+
      geom_line(size=1)+
      geom_point(data=valid.data3,aes(date,Y),color='blue',size=1.5)+ 
      #theme(plot.title = element_text(hjust = 0.5))+
      #one more line for rainfall
      geom_bar(data = rain.data, aes(date, Rain_inc1, fill = Type), 
               alpha = 0.5, stat = "identity", width = 20000) +
      scale_fill_manual(values = c("Rainfall" = "grey", "Irrigation" = "red")) +
      scale_y_continuous(expand=c(0,0), sec.axis = sec_axis(~(.-y_lower_limit)*2.5, 
                                                            name = "Rain or Irrigation (mm/d)",
                                                            breaks = seq(0, 75,by=15))
      ) +
      coord_cartesian(ylim = c(y_lower_limit, y_upper_limit)) +
      #ggtitle(plot.title) +
      theme_bw()+
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        #legend.position = c(0.10, 0.63),
        legend.title=element_text(size=10),
        axis.title.x = element_text(size =12),  # Adjust the x-axis label size
        axis.title.y = element_text(size =12),
        axis.text.x = element_text(size =12),  # Adjust the x-axis tick text size
        axis.text.y = element_text(size =12),  
        legend.text=element_text(size=10))+
      guides(col=guide_legend(ncol=2))+
      scale_x_datetime(date_labels = "%b", breaks="1 month")  
    # Format the x-axis labels to display only the month
    
    # Format the x-axis labels to display only the month
    
  } else {
    p <- ggplot(class.data, aes(date, Y,col=Events, group=1))+
      xlab("Time (2021 Growing season)") + ylab(depth)+
      geom_line(size=1)+
      #geom_point(data=PEL.data,aes(date,Y),color='red',size=1.5)+ 
      #theme(plot.title = element_text(hjust = 0.5))+
      #ggtitle(plot.title) +
      #one more line for rainfall
      geom_bar(data = rain.data, aes(date, Rain_inc1, fill = Type), 
               alpha = 0.5, stat = "identity", width = 20000) +
      scale_fill_manual(values = c("Rainfall" = "grey", "Irrigation" = "red")) +
      scale_y_continuous(expand=c(0,0), sec.axis = sec_axis(~(.-y_lower_limit)*2.5, 
                                                            name = "Rain or Irrigation (mm/d)",
                                                            breaks = seq(0, 60,by=15))
      ) +
      coord_cartesian(ylim = c(y_lower_limit, y_upper_limit)) +
      #ggtitle(plot.title) +
      theme_bw()+
      theme(
        plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.90, 0.53),
        legend.title=element_text(size=10),
        axis.title.x = element_text(size =12),  # Adjust the x-axis label size
        axis.title.y = element_text(size =12),
        axis.text.x = element_text(size =12),  # Adjust the x-axis tick text size
        axis.text.y = element_text(size =12),  
        legend.text=element_text(size=10))+
      guides(col=guide_legend(ncol=2))+
      scale_x_datetime(date_labels = "%b", breaks="1 month")
    # Format the x-axis labels to display only the month
    
  }
  return(list(p,FC,PEL,class.data,valid.data,valid.data2,valid.data3,PEL.data,representative_points))
  #valid.data2 with rf-et transition
  #valid.data3 with rf/RW-gd-et transition
}

#Main(s,d,t)
# s = probe id (1 through 6)
# d = sensor depth (1 through 9)
# t is time interval for smoothing (e.g. t=4 leads to 4 hour moving average)
dataframe<-Main(2,4,0) #For non-smoothed results t=0. better?
#Main(1,1,0)
Main(2,2,0)  #For smoothed results t>0


#Storing values
#create an empty dataframe to store average FC values for all probes
FCaverage <- data.frame(
  Depth=c("4-inch","8-inch","12-inch",
          "16-inch","20-inch","24-inch",
          "28-inch","32-inch","36-inch"),
  CornFull=c(""),
  CornNon=c(""),
  CornPrec=c(""),
  CottonFull=c(""),
  CottonNon=c(""),
  CottonPrec=c(""))

for (S in 1:6){
  for (depth in 1:9){
    FCaverage[depth,S+1] <- Main(S,depth,0)[[2]]
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "2021_soilsensordata","October2023",
                        paste("FCaverage.csv",sep=""))
    write.csv(FCaverage,mypath)
  }
}

#Saving all images
# Set desired dimensions in inches
plot_width <- 5.5
plot_height <- 2.5

# Convert dimensions to pixels assuming 300 pixels per inch (adjust as needed)
width_pixels <- plot_width * 300
height_pixels <- plot_height * 300

for (S in 1:6){
  for (depth in 1:9){
    image.title <- paste(StnIds$Name[S],DepthIds[depth],sep=",")
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "2021_soilsensordata","October2023",
                        paste("SM",image.title,".png",sep=" "))
    png(mypath,width = width_pixels, height = height_pixels, res = 180) 
    #only run this while saving, not viewing plots
    print(Main(S,depth,0)[[1]])
    dev.off()
  }
}

#for (i in dev.list()[1]:dev.list()[length(dev.list())]) {dev.off()}

#storing percentile values
FCpercentile <- data.frame(
  Depth=c("4-inch","[25%, 75%]","8-inch","[25%, 75%]","12-inch","[25%, 75%]",
          "16-inch","[25%, 75%]","20-inch","[25%, 75%]","24-inch","[25%, 75%]",
          "28-inch","[25%, 75%]","32-inch","[25%, 75%]","36-inch","[25%, 75%]"),
  CornFull=c(""),
  CornNon=c(""),
  CornPrec=c(""),
  CottonFull=c(""),
  CottonNon=c(""),
  CottonPrec=c(""))

for (S in 1:6){
  for (depth in 1:9){
    d <- 2*depth-1
    FCpercentile[d,S+1] <- Main(S,depth,0)[[2]]
    FCpercentile[d+1,S+1] <- paste(quantile(Main(S,depth,0)[[6]]$Y,c(0.25)),
                                   quantile(Main(S,depth,0)[[6]]$Y,c(0.75)),
                                   sep=c(", "))
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "2021_soilsensordata","March2023","June_Results",
                        paste("FCpercentile.csv",sep=""))
    #write.csv(FCpercentile,mypath)
  }
}
#quantile(Main(2,6,0)[[6]]$Y,c(0.25, 0.75))

##5/3/2023
#Storing PWP values
#create an empty dataframe to store average PWP values for all probes
PWPaverage <- data.frame(
  Depth=c("4-inch","8-inch","12-inch",
          "16-inch","20-inch","24-inch",
          "28-inch","32-inch","36-inch"),
  CornFull=c(""),
  CornNon=c(""),
  CornPrec=c(""),
  CottonFull=c(""),
  CottonNon=c(""),
  CottonPrec=c(""))

for (S in 1:6){
  for (depth in 1:9){
    PWPaverage[depth,S+1] <- Main(S,depth,0)[[3]]
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "2021_soilsensordata","October2023",
                        paste("PWPaverage.csv",sep=""))
    write.csv(PWPaverage,mypath)
  }
}
