#### Step III -Plotting results ###

library(dplyr)
library(ggplot2)
library(viridis)
library(rworldmap)
library(sf)

#Load data coordinates
cord50<-data.table::fread("./input_data/coord.50m.tutto.ipcc.csv")
cord50$pixel<-paste("p",c(1:nrow(cord50)),sep="") #add a "p" before the number of the pixel to allow the match

#Load csv with DEB model outputs
finalcsv<-data.table::fread("../example/DEB_Output_Wild_test.csv",sep=";",dec = ".",
                            header = T,data.table = F)

#Merge data
dati<-left_join(cord50,finalcsv,by=c("pixel"="pixel"))

world<-map_data("world")


#Plotting Length
sumtot<-summary(dati$Length.cm)

ggplot(data=world)+
  geom_tile(data=dati,mapping = aes(x=Longitude,y=Latitude,fill=Length.cm))+
  geom_map(map = world,mapping = aes(map_id=region),col="grey40")+ 
  expand_limits(x = c(11,17), y = c(36, 39))+
  scale_y_continuous(breaks = seq(36,39,1))+
  scale_x_continuous(breaks = seq(11,17,1))+
  coord_sf(xlim = c(11,17), ylim = c(36, 39), expand = FALSE) +
  scale_fill_viridis(discrete = F ,option = "C", direction = -1,
                     limits=c(sumtot[1],sumtot[6]))+
  theme_classic()+theme(panel.background = element_rect(fill = "lightblue"))+
  theme(text=element_text(size=15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+ 
  labs(fill = "Length, cm")+
  xlab(NULL)+ylab(NULL)



#Plotting Eggs
sumtot<-summary(dati$N.eggs.CUM)

ggplot(data=world)+
  geom_tile(data=dati,mapping = aes(x=Longitude,y=Latitude,fill=N.eggs.CUM))+
  geom_map(map = world,mapping = aes(map_id=region),col="grey40")+ 
  expand_limits(x = c(11,17), y = c(36, 39))+
  scale_y_continuous(breaks = seq(36,39,1))+
  scale_x_continuous(breaks = seq(11,17,1))+
  coord_sf(xlim = c(11,17), ylim = c(36, 39), expand = FALSE) +
  scale_fill_viridis(discrete = F ,option = "C", direction = -1,
                     limits=c(sumtot[1],sumtot[6]))+
  theme_classic()+theme(panel.background = element_rect(fill = "lightblue"))+
  theme(text=element_text(size=15),
        axis.text = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))+ 
  labs(fill = "Eggs.cum, n")+
  xlab(NULL)+ylab(NULL)


