rm(list=ls(all=TRUE))
library(dplyr)
setwd("G:/My Drive/Tidewater NUE/Data Dump/Data2021")
## Station IDs
StnIds <- data.frame(
  ID=c("002060E6","00206FAD","00206F86","00206F61","0120A7B7",
            "01209C43","002060AE"),
  Name=c("Corn_Full_irr","Corn_No_irr","Corn_Precision","Cotton_Full_irr","Cotton_No_irr",
              "Cotton_Precision","MetSt"))

#files <- list.files(pattern="csv")
#files <- files[158:164]

#Loading the last 7 csv files
sensor.data <- list()
for (s in 1:7){
  sensor.data[[s]] <- read.csv(file=list.files(pattern="csv")[s+157])
}

#Let's keep only the soil moisture data from all the dataframes
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


ColorPalette <- colorRampPalette(c('red','orange','green','blue','purple'))
plot.colors<-ColorPalette(5)

#Using function to plot the density plots
sensor.plot <- function(s){
  plot.title <- StnIds$Name[s]
  mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                      "Density_Images2021",paste("SM",plot.title,".png",sep=" "))
  png(mypath,res=100)
  layout(mat = matrix(c(1,2,3),
                      nrow=3,
                      ncol=1),
         heights=c(0.05,0.045,0.06),      #heights of 3 rows #add title and label to know relative heights
         widths=c(1.5,1.5,1.5))           #widths of 1 column
  layout.show(3)
  
  myYlim1 <- max(c(density(soil_moisture[[s]][,2])$y,density(soil_moisture[[s]][,3])$y,density(soil_moisture[[s]][,4])$y))
  myXlim1 <- range(c(soil_moisture[[s]][,2],soil_moisture[[s]][,3],soil_moisture[[s]][,4]))
  par(mar = c(0,4,2,2))   #bottom,left,top,right margins
  plot(density(soil_moisture[[s]][,2]),main = plot.title, xlim=myXlim1,ylim=c(0,myYlim1), xlab="",ylab="Density",col=plot.colors[1])
  lines(density(soil_moisture[[s]][,3]), col=plot.colors[2])
  lines(density(soil_moisture[[s]][,4]), col=plot.colors[3])
  legend("topright",                                 
         legend = c("SM_1", "SM_2", "SM_3"),
         col = c("red", "orange", "green"),
         lty = 1,cex=0.8)
  #dev.off()   #this completes image creation. otherwise the image won't update and show
  
  #Plot2
  myYlim2 <- max(c(density(soil_moisture[[s]][,5])$y,density(soil_moisture[[s]][,6])$y,density(soil_moisture[[s]][,7])$y))
  myXlim2 <- range(c(soil_moisture[[s]][,5],soil_moisture[[s]][,6],soil_moisture[[s]][,7]))
  par(mar = c(0,4,1,2))
  plot(density(soil_moisture[[s]][,5]), main = "", xlim=myXlim2, ylim=c(0,myYlim2), xlab="",col=plot.colors[1])#main = "SM 4-5-6"
  lines(density(soil_moisture[[s]][,6]), col=plot.colors[2])
  lines(density(soil_moisture[[s]][,7]), col=plot.colors[3])
  legend("topright",                                 
         legend = c("SM_4", "SM_5", "SM_6"),
         col = c("red", "orange", "green"),
         lty = 1,cex=0.8)
  
  #dev.off()
  #Plot3
  myYlim3 <- max(c(density(soil_moisture[[s]][,8])$y,density(soil_moisture[[s]][,9])$y,density(soil_moisture[[s]][,10])$y))
  myXlim3 <- range(c(soil_moisture[[s]][,8],soil_moisture[[s]][,9],soil_moisture[[s]][,10]))
  par(mar = c(4,4,1,2))
  plot(density(soil_moisture[[s]][,8]), main = "", xlim=myXlim3,ylim=c(0,myYlim3),
       xlab="Soil moisture (%)",col=plot.colors[1]) #main = "SM 7-8-9"
  lines(density(soil_moisture[[s]][,9]), col=plot.colors[2])
  lines(density(soil_moisture[[s]][,10]), col=plot.colors[3])
  legend("topright",                                 
         legend = c("SM_7", "SM_8", "SM_9"),
         col = c("red", "orange", "green"),
         lty = 1,cex=0.8)
  
  dev.off()  #just run this dev. do not run other dev.off to save the overall image
}  
for(s in 1:6){
  #sensor.plot(s)
}


#Checking peaks and valleys without any thresholds
SM_peaks <- function(s,D,thresh){
  dat <- soil_moisture[[s]][,D]
  d <- density(dat)
  probe <- StnIds$Name[s]
  dep <- as.character(D-1)
  #plot.title <- paste(probe,"depth",dep,sep="_")
  plot.title <- paste("")
  plot(d, type = "l",main=plot.title,xlab="VWC (%)",cex.lab=1.35,cex.axis=1.35)
  pks <- which(diff(sign(diff(d$y, na.pad = FALSE)), na.pad = FALSE)<0) +2 #ca                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              lculates peaks
  pks1 <- which(diff(sign(diff(d$y, na.pad = FALSE)), na.pad = FALSE)>0) +2 #calculates valleys
  a <<- pks[d$y[pks]>thresh*max(d$y[pks])] #thresholding based on height of max peak
  #let's plot
  points(d$x[a], d$y[a], col="red", pch=2,lwd=2, legend = "Peaks")# to check peak
  points(d$x[pks1], d$y[pks1], col="blue", pch=8, lwd=1.5, legend = "Valleys") #to check 
  # Add legend
  legend("topleft", legend = c("Peaks", "Valleys"), col = c("red", "blue"),
         pch = c(2, 8), lwd = c(2, 1.5), cex = 0.9,bty = "n")
  
  # rm(d,x,y)
}


#calling function 
SM_peaks(2,6,0)

#Saving density plots with initial peaks into library
for (S in 1:6){
  for (depth in 2:10){
    probe <- StnIds$Name[S]
    D <- as.character(depth-1)
    plot.title <- paste(probe,"depth",D,sep="_")
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "Plot_initialpeaks",paste(plot.title,".png",sep="_"))
    png(mypath,res=180)
   
    #name <- paste(probe,"depth",D,sep="_")
    plot <- SM_peaks(S, depth, 0)
    dev.off()
  }
} 

#Creating empty dataframe to save peak and valley coordinates and threshold result
Pktable <- data.frame(
  X_peak=numeric(),
  Y_peak=numeric(),
  X_valley=numeric(),
  Y_valley=numeric(),
  PK_thres_Valley=numeric(),#saving values that are higher than nearest valley
  PK_thres_PKmax=numeric()  #Saving values that are higher than some % of max peak
)

##Below we will try setting threshold
#Main function
SM_peaks <- function(S,depth,thresh1,thresh2){
  dat <- soil_moisture[[S]][,depth]
  d <- density(dat)
  plot(d, type = "l")
  pks <- which(diff(sign(diff(d$y, na.pad = FALSE)), na.pad = FALSE)<0) +2 
  pks1 <- which(diff(sign(diff(d$y, na.pad = FALSE)), na.pad = FALSE)>0) +2
  n.pks <- length(pks)
  n.valley <- length(pks1)
  
  Pktable <- data.frame(  
    X_peak=rep(NA, n.pks),
    Y_peak= NA,
    X_valley=NA,
    Y_valley=NA,
    PK_thres_Valley=NA,
    PK_thres_PKmax=NA  
  )
  
  valley_x <- d$x[pks1]
  valley_y <- d$y[pks1]
  Pktable$X_peak <- d$x[pks]
  Pktable$Y_peak <- d$y[pks]
  Pktable$Y_valley[length(pks)] <- d$y[length(d$x)] 
  Pktable$X_valley[length(pks)] <- d$x[length(d$x)] 
  
  for (i in 1:(length(pks)-1)) {
    #Let's save all peaks and valleys X and Y coordinates
    Pktable$X_valley[i] <- valley_x[i]
    Pktable$Y_valley[i] <- valley_y[i]
    if (Pktable$Y_valley[i] < thresh1*Pktable$Y_peak[i]){
      Pktable$PK_thres_Valley[i] <- Pktable$Y_peak[i]
    } 
    
  }
  if ((Pktable$Y_valley[length(pks)] < thresh1*Pktable$Y_peak[length(pks)])||
      (Pktable$Y_valley[length(pks)-1] < thresh1*Pktable$Y_peak[length(pks)])){
    Pktable$PK_thres_Valley[length(pks)] <- Pktable$Y_peak[length(pks)]
  }
  ###
  for (i in 1:length(pks)){
    if (!is.na(Pktable$PK_thres_Valley[i])){
      if (Pktable$PK_thres_Valley[i]>thresh2*(max(Pktable$PK_thres_Valley, na.rm=TRUE))){
        Pktable$PK_thres_PKmax[i] <- Pktable$PK_thres_Valley[i]
      }
    }
    
  }
  
  #### 
  return(Pktable) 
  #points(Pktable$X_peak, Pktable$PK_thres_PKmax, col="red", pch=2,lwd=2)
  
} # close function

#Saving Peaks values dataframe into library
for (S in 1:6){
  for (depth in 2:10){
    probe <- StnIds$Name[S]
    D <- as.character(depth-1)
    #name <- paste(probe,"depth",D,sep="_")
    data <- SM_peaks(S, depth, thresh1=0.75, thresh2=0.30)
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "Peaks_and_valleys",paste(probe,"depth",D,".csv",sep="_"))
    write.csv(data,mypath)
  }
} 
test <- SM_peaks(S=3, depth=7, thresh1=0.75, thresh2=0.30) #depth 2 to 10 #Probe(S=1 to 6)


#create an empty dataframe to store number of peaks for all probes (9 sensors each)
PK_number <- data.frame(
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
  for (depth in 2:10){
    data <- SM_peaks(S, depth, thresh1=0.75, thresh2=0.30)
    PK_number[depth-1,S+1] <- nrow(na.omit(data))
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "Peak_values","Updated",paste("PK_number.csv",sep=""))
    write.csv(PK_number,mypath)
  }
}

#Let's automate code to save values higher than mean as FC (plot with single peak)
#If one peak and lower than 1, it will be PWP.
#create an empty dataframe to store high peak values for all probes
PK_highvalues <- data.frame(
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
  for (depth in 2:10){
    data <- SM_peaks(S, depth, thresh1=0.75, thresh2=0.30)
    data1 <- na.omit(data)
    mean <- mean(soil_moisture[[S]][,depth])
    if (nrow(data1)>1){
      PK_highvalues[depth-1,S+1] <- tail(data1$X_peak,1)#high peaks
      mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                          "Peak_values","Updated",paste("PK_highvalues.csv",sep=""))
      write.csv(PK_highvalues,mypath)
    }
    if (nrow(data1)==1){
      if (tail(data1$X_peak,1)>mean){
        PK_highvalues[depth-1,S+1] <- tail(data1$X_peak,1)
      mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                          "Peak_values","Updated",paste("PK_highvalues.csv",sep=""))
      write.csv(PK_highvalues,mypath)
      }
    }
  }
}

#create an empty dataframe to store low peak values for all probes
PK_lowvalues <- data.frame(
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
  for (depth in 2:10){
    data <- SM_peaks(S, depth, thresh1=0.75, thresh2=0.30)
    data1 <- na.omit(data)
    mean <- mean(soil_moisture[[S]][,depth])
    if (nrow(data1)>1){
    PK_lowvalues[depth-1,S+1] <- head(data1$X_peak,1)#low peaks
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "Peak_values","Updated",paste("PK_lowvalues.csv",sep=""))
    write.csv(PK_lowvalues,mypath)
    }
    if (nrow(data1)==1){
      if (tail(data1$X_peak,1)<mean){
        PK_lowvalues[depth-1,S+1] <- tail(data1$X_peak,1)
    mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                        "Peak_values","Updated",paste("PK_lowvalues.csv",sep=""))
    write.csv(PK_lowvalues,mypath)
      }
    }
  }
}



##DEnsity signature ends here


#USing function to plot time series plots
## Convert date columns to posix
for (s in 1:length(soil_moisture)){  
  soil_moisture[[s]][,1] <- as.POSIXct(soil_moisture[[s]][,1],
                                     format="%Y-%m-%d %H:%M:%OS")
}
  sensor.plot <- function(s){
  plot.title <- paste(StnIds$Name[s])
  
  image.title <- StnIds$Name[s]
  mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures","SM_Images2021",paste("SM",image.title,".png",sep=" "))
  #png(mypath,res=160) #only run this while saving, not viewing plots
  
  #png("Images/plot.png", res=100) #setting width and height gives error of large margins, changed margins for each of the 3 plots (might be an issue)
  layout(mat = matrix(c(1,2,3),
                      nrow=3,
                      ncol=1),
         heights=c(0.05,0.045,0.06),      #heights of 3 rows #add title and label to know relative heights
         widths=c(1.5,1.5,1.5))           #widths of 1 column
  layout.show(3)
  
  
  myYlim1 <- range(c(soil_moisture[[s]][,2],soil_moisture[[s]][,3],soil_moisture[[s]][,4]))
  myXlim <- range(c(soil_moisture[[s]][,1]))
  par(mar = c(0,4,2,2))   #bottom,left,top,right margins
  plot(soil_moisture[[s]][,1],soil_moisture[[s]][,2],type="l",xaxt="n",main = plot.title,ylim=myYlim1, xlim=myXlim,xlab="",ylab="VWC (%)",col=plot.colors[1])
  lines(soil_moisture[[s]][,1],soil_moisture[[s]][,3], col=plot.colors[4])
  lines(soil_moisture[[s]][,1],soil_moisture[[s]][,4], col=plot.colors[3])
  legend("bottomright",                                 
         legend = c("4 inches", "8 inches", "12 inches"),
         col = c("red", "blue", "green"),
         lty = 1,cex=0.5)
  #box()
  #dev.off()   #this completes image creation. otherwise the image won't update and show
  
  #Plot2
  #png("SM_456.png")
  
  myYlim2 <- range(c(soil_moisture[[s]][,5],soil_moisture[[s]][,6],soil_moisture[[s]][,7]))
  #myXlim2 <- as.character(range(c(soil_moisture[[s]][,1])))
  par(mar = c(0,4,1,2))
  plot(soil_moisture[[s]][,1],soil_moisture[[s]][,5], type="l",xaxt="n",main = "", ylim=myYlim2,xlim=myXlim, xlab="",ylab="VWC (%)",col=plot.colors[1])#main = "SM 4-5-6"
  lines(soil_moisture[[s]][,1],soil_moisture[[s]][,6], col=plot.colors[4])
  lines(soil_moisture[[s]][,1],soil_moisture[[s]][,7], col=plot.colors[3])
  legend("bottomright",                                 
         legend = c("16 inches", "20 inches", "24 inches"),
         col = c("red", "blue", "green"),
         lty = 1,cex=0.5)
  #box()
  #dev.off()
  
  #Plot3
  #png("SM_789.png")
  
  myYlim3 <- range(c(soil_moisture[[s]][,8],soil_moisture[[s]][,9],soil_moisture[[s]][,10]))
  #myXlim3 <- range(c(soil_moisture[[s]][,1]))
  par(mar = c(4,4,1,2))
  plot(soil_moisture[[s]][,1],soil_moisture[[s]][,8], type="l",main = "", ylim=myYlim3,xlim=myXlim, xlab=" ",ylab="VWC (%)",col=plot.colors[1]) #main = "SM 7-8-9"
  lines(soil_moisture[[s]][,1],soil_moisture[[s]][,9], col=plot.colors[4])
  lines(soil_moisture[[s]][,1],soil_moisture[[s]][,10], col=plot.colors[3])
  legend("bottomright",                                 
         legend = c("28 inches", "32 inches", "36 inches"),
         col = c("red", "blue", "green"),
         lty = 1,cex=0.5)
  
  #dev.off()  #just run this dev. do not run other dev.off to save the overall image
}  
sensor.plot(6)


for(s in 1:6){
  sensor.plot(s)
}


#Final plot for timeseries in manuscript

for (s in 1:length(soil_moisture)){  
  soil_moisture[[s]][,1] <- as.POSIXct(soil_moisture[[s]][,1],
                                       format="%Y-%m-%d %H:%M:%OS")
}
layout(mat = matrix(c(1,2,3),
                    nrow=3,
                    ncol=1),
       heights=c(0.05,0.045,0.06),      #heights of 3 rows #add title and label to know relative heights
       widths=c(1.5,1.5,1.5))           #widths of 1 column
layout.show(3)


myYlim1 <- range(c(soil_moisture[[2]][,2],soil_moisture[[2]][,3],soil_moisture[[2]][,4]))
myXlim <- range(c(soil_moisture[[2]][,1]))
par(mar = c(0,6,4,2))   #bottom,left,top,right margins
plot(soil_moisture[[2]][,1],soil_moisture[[2]][,2],type="l",xaxt="n",
     main = "",ylim=myYlim1, xlim=myXlim,xlab="",ylab="",col=plot.colors[1])
lines(soil_moisture[[2]][,1],soil_moisture[[2]][,3], col=plot.colors[4])
lines(soil_moisture[[2]][,1],soil_moisture[[2]][,4], col=plot.colors[3])
#abline(h=25, lty=5,col="green")
#abline(h=27, lty=5, col="blue")
#abline(h=23, lty=5, col="red")

#target_date <- soil_moisture[[2]][1800, 1]
legend("bottomright",  y=50,                               
       legend = c("10 cm", "20 cm", "30 cm"),
       col = c("red", "blue", "green"),inset=c(0,0.8), #(+x,+y)from bottomright
       lty = 1,cex=0.9,xpd=TRUE,horiz=TRUE,bty="n",y.intersp=0.8)

myYlim2 <- range(c(soil_moisture[[2]][,5],soil_moisture[[2]][,6],soil_moisture[[2]][,7]))
#myXlim2 <- as.character(range(c(soil_moisture[[s]][,1])))
par(mar = c(0,6,1,2))
plot(soil_moisture[[2]][,1],soil_moisture[[2]][,5], type="l",xaxt="n",
     main = "", ylim=myYlim2,xlim=myXlim, xlab="",
     ylab="",col=plot.colors[1])
lines(soil_moisture[[2]][,1],soil_moisture[[2]][,6], col=plot.colors[4])
lines(soil_moisture[[2]][,1],soil_moisture[[2]][,7], col=plot.colors[3])
#abline(h=, lty=2,col="green")
#abline(h=31, lty=2, col="blue")
#abline(h=28, lty=2, col="red")

#PWP
#abline(h=30, lty=4,col="green")
#abline(h=25, lty=4, col="blue")
#abline(h=, lty=4, col="red")

#target_date <- soil_moisture[[2]][1800, 1]
legend("bottomright",  y=35,                                 
       legend = c("40 cm", "50 cm", "60 cm"),
       col = c("red", "blue", "green"),inset=c(0,0.9),
       lty = 1,cex=0.9,xpd=TRUE,horiz=TRUE,bty="n",y.intersp=0.8)

myYlim3 <- range(c(soil_moisture[[2]][,8],soil_moisture[[2]][,9],soil_moisture[[2]][,10]))
#myXlim3 <- range(c(soil_moisture[[s]][,1]))
par(mar = c(4,6,1,2))
plot(soil_moisture[[2]][,1],soil_moisture[[2]][,8], type="l",main = "", 
     ylim=myYlim3,xlim=myXlim, xlab=" ",ylab="",col=plot.colors[1])
lines(soil_moisture[[2]][,1],soil_moisture[[2]][,9], col=plot.colors[4])
lines(soil_moisture[[2]][,1],soil_moisture[[2]][,10], col=plot.colors[3])
#abline(h=, lty=2,col="green")
#abline(h=, lty=2, col="blue")
#abline(h=, lty=2, col="red")
#PWP
#abline(h=26, lty=4,col="green")
#abline(h=29.5, lty=4, col="blue")
#abline(h=33, lty=4, col="red")

#target_date <- soil_moisture[[2]][1800, 1]
legend("bottomright",  y=30.1,                                  
       legend = c("70 cm", "80 cm", "90 cm"),
       col = c("red", "blue", "green"),inset=c(0,0.83),
       lty = 1,cex=0.9,xpd=TRUE,horiz=TRUE,bty="n",y.intersp=0.6)
y_title_text <- "VWC (%)"
mtext(y_title_text, side = 2, line = -2, outer = TRUE)

