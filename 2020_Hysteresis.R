rm(list=ls(all=TRUE))
library(dplyr)
load("soil_moisture2.RData")
for (s in 1: 10){
  data <- soil_moisture[[s]]
  data <- data[data[,2] >= 1,]
  data <- na.omit(data)
  soil_moisture[[s]] <- data
}
StnIds <- data.frame(
  ID=c("01209A5A","01209A4A","01209C43",
       "01209A50","01209C51","01209A63",
       "0120A68D","0120A68F","0120A7A8","0120A67D"),
  Name=c("Corn1","Corn2","Corn3","Cotton1",
         "Cotton2","Cotton3","CornNon","CornFull","CottonNon","CottonFull") )


#summary table for tracking drying vs wetting points
#n <- length(soil_moisture[[4]][,1])
#Drying$Date <- soil_moisture[[4]][,1]
#Let's create lists of 6 dataframes
Drying <- list()
for (s in 1:10){
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

df_grids <- list()
Hyst <- function(S){
  probe <- StnIds$Name[S]
  plot.title <- paste(probe, sep="")
  mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                      "Hysteresis_2020","16-32inch",paste(plot.title,".png"))
  png(mypath,res=100)
  data <- Drying[[S]]
  #data <- data[(data$Date> "2021-06-10" & data$Date < "2021-06-20"),]
  X <- data[,2]
  Y <- data[,3]
  D1 <- paste("SWC % (16inch)")
  D2 <- paste("SWC % (32inch)")
  #arrowLine(X,Y,N=length(X),xlab=D1,ylab=D2,main=plot.title)
  plot(X,Y,xlab=D1,ylab=D2,main=plot.title,pch=20,
       col=ifelse(data$Color == "Orange","orange",
                  ifelse(data$Color == "blue","blue","grey")))
  #grid(nx = 40, ny = 40, col = "lightgray", lty = "dotted",
  #lwd = par("lwd"), equilogs = TRUE)
  
  legend("bottomright", col=c("Orange","Blue","Grey"),
         c("Drying","Wetting","Complicated"), pch= c(20))
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

#Hyst(1)
#dev.off()

for (i in 1:10){
  
  Hyst(i)
  dev.off()
}

#
#S=1 to 6 # depth1 & depth2 = 2 to 10
#Try above function for other plots too and see how corn precision differs

#Let's save all wet and dry attractor values in a dataframe
Attractor <- data.frame(
  Site=c("Corn1","Corn2","Corn3","Cotton1",
         "Cotton2","Cotton3","CornNon","CornFull","CottonNon","CottonFull")
)

for (S in 1:10){
  dist<-sqrt(df_grids[[S]]$Depth1_SM^2+df_grids[[S]]$Depth2_SM^2)
  Attractor$Dry_att_D1[S] <- df_grids[[S]]$Depth1_SM[which.min(dist)]
  Attractor$Dry_att_D2[S] <- df_grids[[S]]$Depth2_SM[which.min(dist)]
  Attractor$Wet_att_D1[S] <- df_grids[[S]]$Depth1_SM[which.max(dist)]
  Attractor$Wet_att_D2[S] <- df_grids[[S]]$Depth2_SM[which.max(dist)]
  
  mypath <- file.path("G:","Shared drives","Tidewater-HydSignatures",
                      "Hysteresis_2020","Attractors",paste("16inch-32inch.csv",sep=""))
  write.csv(Attractor,mypath)
}

#Likewise do D2vsD6(4inch vs 20inch);D3vsD7;D4vsD8;D5vsD9
