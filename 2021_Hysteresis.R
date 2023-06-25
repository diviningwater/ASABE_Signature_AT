load("G:/Shared drives/Tidewater-HydSignatures/2021_soilsensordata/Latest.R")
library(dplyr)

#Plotting scatterplot
#Let's plot for full irrigated corn
X <- soil_moisture[[1]]$SM_1   #just create vector with $ unlike [] that create list
Y <- soil_moisture[[1]]$SM_3
plot(X,Y, main="Hysteresis example for shallow sensors, corn full", xlab="Soil moisture, 4-inch",
     ylab="Soil moisture, 12-inch", pch=16)

Drying <- list()
for (s in 1:6){
  #data <- soil_moisture[[s]]
  #Creating empty dataframe
  Drying[[s]] <- data.frame(
    Date=soil_moisture[[s]][,1],
    Depth1_SM=NA,
    Depth2_SM=NA,
    Is_dry=NA,
    Color=NA)
  Drying[[s]]$Depth1_SM <- soil_moisture[[s]][,5]#2 for depth1(4inch)
  Drying[[s]]$Depth2_SM <- soil_moisture[[s]][,9]#6 for depth5(20inch)
  for (t in 1:(length(Drying[[s]]$Date)-1)){
    if (Drying[[s]]$Depth1_SM[t+1]<Drying[[s]]$Depth1_SM[t] 
        &&
        Drying[[s]]$Depth2_SM[t+1]<Drying[[s]]$Depth2_SM[t]){
      Drying[[s]]$Is_dry[t] <- "Yes"
      Drying[[s]]$Color[t] <- "Orange"
    }else {
      if(Drying[[s]]$Depth1_SM[t+1]>Drying[[s]]$Depth1_SM[t] 
         &&
         Drying[[s]]$Depth2_SM[t+1]>Drying[[s]]$Depth2_SM[t]){
        Drying[[s]]$Is_dry[t] <- "No"
        Drying[[s]]$Color[t] <- "blue"
      }else {
        Drying[[s]]$Is_dry[t] <- "Complicated"
        Drying[[s]]$Color[t] <- "grey"
      }
    }
  }
}

#Using different colors for drying and wetting and plotting again
#Function to plot hysteresis plots
#Updated on 3/29/2023
StnIds$Name[1] <- "Corn Full"
StnIds$Name[2] <- "Corn Non"
StnIds$Name[3] <- "Corn Precision"
StnIds$Name[4] <- "Cotton Full"
StnIds$Name[5] <- "Cotton Non"
StnIds$Name[6] <- "Cotton Precision"

df_grids <- list()
Hyst <- function(S){
  probe <- StnIds$Name[S]
  plot.title <- paste(probe, sep="")
  mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
  "Hysteresis_2023","16-32inch",paste(plot.title,".png"))
  #png(mypath,res=100)
  data <- Drying[[S]]
  #data <- data[(data$Date> "2021-06-10" & data$Date < "2021-06-20"),]
  X <- data[,2]
  Y <- data[,3]
  #D1 <- paste("\u03B8 (40 cm)")
  D1 <- paste("VWC (%) (40 cm)")
  D2 <- paste("VWC (%) (80 cm)")
  #arrowLine(X,Y,N=length(X),xlab=D1,ylab=D2,main=plot.title)
  plot(X,Y,xlab=D1,ylab=D2,main="",pch=20,cex.lab=1.2,
       cex.axis=1.2,
       col=ifelse(data$Color == "Orange","orange",
                  ifelse(data$Color == "blue","blue","grey")))
  #grid(nx = 40, ny = 40, col = "lightgray", lty = "dotted",
  #lwd = par("lwd"), equilogs = TRUE)
  
  legend("topleft", col=c("Orange","Blue","Grey","red"),
         c("Drying","Wetting","Mixed", "High-density cells"), pch= c(20,20,20,8),
         cex=0.8,bty="n",y.intersp=0.5)
  a <- max(X)
  b <- max(Y)
  df_grids[[S]] <<- data%>%    ###double <<- helped in saving df1 to global environment
    mutate(
      cut_x = cut(X, breaks = seq(0,a, by = 0.25), include.lowest = T),
      cut_y = cut(Y, breaks = seq(0,b, by = 0.25), include.lowest = T)
    ) %>%
    group_by(cut_x,cut_y) %>%
    dplyr::mutate(num = n())%>%
    slice(1)%>%  #this shows total count for each grid only once, else shows count
    # corresponding to each value present in a particular grid.
    
    #Time to filter top 5% high density grids
    arrange(desc(num))%>%
    #Filter number of points within mesh
    filter(num>15)
  # filter(Depth1_SM==max(Depth1_SM))
  dist<-sqrt(df_grids[[S]]$Depth1_SM^2+df_grids[[S]]$Depth2_SM^2)
  
  points(df_grids[[S]]$Depth1_SM[which.max(dist)],
         df_grids[[S]]$Depth2_SM[which.max(dist)],
         type="p",pch=8,col="red",cex=2)
  
  points(df_grids[[S]]$Depth1_SM[which.min(dist)],
         df_grids[[S]]$Depth2_SM[which.min(dist)],
         type="p",pch=8,col="red",cex=2)
  #return(df_grids)
}

Hyst(1)
#Save frm export with 545 and 450 dimensions
#dev.off()

for (i in 1:6){
  
  #Hyst(i)
  dev.off()
}

#
#S=1 to 6 # depth1 & depth2 = 2 to 10
#Try above function for other plots too and see how corn precision differs

#Let's save all wet and dry attractor values in a dataframe
Attractor <- data.frame(
  Site=c("CornFull","CornNon","CornPrec",
          "CottonFull","CottonNon","CottonPrec")
)

for (S in 1:6){
  dist<-sqrt(df_grids[[S]]$Depth1_SM^2+df_grids[[S]]$Depth2_SM^2)
  Attractor$Dry_att_D1[S] <- df_grids[[S]]$Depth1_SM[which.min(dist)]
  Attractor$Dry_att_D2[S] <- df_grids[[S]]$Depth2_SM[which.min(dist)]
  Attractor$Wet_att_D1[S] <- df_grids[[S]]$Depth1_SM[which.max(dist)]
  Attractor$Wet_att_D2[S] <- df_grids[[S]]$Depth2_SM[which.max(dist)]
  
  mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                      "Hysteresis_2023","Attractors",paste("16inch-32inch.csv",sep=""))
  write.csv(Attractor,mypath)
}

#Lets do D2vsD6(4inch vs 20inch);D3vsD7;D4vsD8;D5vsD9
  
  
  
  
  