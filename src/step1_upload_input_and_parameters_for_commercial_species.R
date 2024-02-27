## DEB model application ##
##This code compute Dynamic Energy Budget model from an input of temperature and food and saves multiple output to a single csv file##

##############################################################################################################################################################
#Load necessary libraries
library(zoo)
library(data.table)

rm(list=ls()) #clean R environment 
options(scipen = 999) #adopt scientific notation
##############################################################################################################################################################

### Step I - Upload Input data

#Set working directory and load inputs of temperature and food
#setwd() #set working directory

#Load temperature and food csv matrices:
#rows represent the steps of the time series (e.g. hours or days)
#columns represent the respective spatial point (e.g. pixels or coordinates) of input data
temperaturecsv<-data.table::fread("./input_data/DEB temp test.csv",sep=";",dec = ".", header = T,data.table = F)
foodcsv<-data.table::fread("./input_data/DEB chl test.csv",sep=";",dec = ".", header = T,data.table = F)

#For test the model subset the matrices accounting only the first two spatial points (e.g. pixels or coordinates)
# temperaturecsv<- temperaturecsv[,1:2] #Uncomment to subset
# foodcsv<-foodcsv[,1:2] #Uncomment to subset


#### Parameters ####
#Add all required parameters for running the DEB model
#For full parameters description see PDF file of the DEB model tutorial

#Life expectation of model species
year.v<-4
#Desired commercial or target size
ComSize <- 5


## Size, shape and volumes ##
#Length at seeding (cm) in case of aquaculture species
Ls <- 2
# Length at puberty (cm)
Lp <- 3.44 
# Shape coefficient of adult
del_M <- 0.2133
# Volume at seeding (cm^3) in case of aquaculture species
V_s <- (Ls * del_M)^3
# Volume at puberty (cm^3)
V_p <- (Lp * del_M)^3 


## Feeding ##
# Area-specific maximum ingestion rate (J cm^-2 h^-1)
J_Xm <- 10.62
# Assimilation efficiency (Adimensional, goes from 0 to 1)
ae <- 0.73
# Half saturation coefficient (ug L^-1)
K <- 0.91
# Food energy conversion factor (J mg^-1)
mu_X <- 1


## Energetic ##
# Volume specific cost for structure (J cm^-3)
E_G <- 5367 
# Maximum storage density or reserve capacity (J cm^-3)
E_m <- 2713
# Volume specific maintenance costs (J cm^-3 h^-1)
p_M <- 0.92
# Fraction of energy allocated to maintenance plus growth
kap <- 0.8
# Egg cost (J)
egg_cost <- 0.1
# Fraction of energy used for reproduction
kap_R <- 0.7
# Maximum temperature for spawning (?C)
Tmax.Spawn <- 18
# Minimum temperature for spawning (?C)
Tmin.Spawn <- 10
# Cumulative daydegree value after threshold for spawning
DDtrigger <- 420 


## Temperature Correction Parameters ##
# Arrhenius temperature (?K)
TA <- 3243
# Lower boundary of thermal tolerance range (?K)
TL <- 278
# Upper boundary of thermal tolerance range (?K)
TH <- 308
# Rate of decrease at lower boundary (?K)
TAL <- 4139
# Rate of decrease at upper boundary (?K)
TAH <- 1739 



## Mass balance ##
# Chemical indices for organics, nO
nO <- data.frame(Food = c(1,1.8,0.5,0.2),
                 Structure = c(1,1.8,0.5,0.2),
                 Reserves = c(1,1.8,0.5,0.2),
                 Products = c(1,1.8,0.5,0.2))



## Mass conversions ##
# Specific densities for structure (g cm^-3)
d_V <- 1
# Specific densities for reserve (g cm^-3)
d_E <- 0.1
#Molecular weights for structure (g mol^-1)
w_V <- 12*nO[1,2] + 1*nO[2,2] +16*nO[3,2] +14*nO[4,2]
#Molecular weights for reserve (g mol^-1)
w_E <- 12*nO[1,3] + 1*nO[2,3] +16*nO[3,3] +14*nO[4,3]
# Chemical potentials (Energy content in ash-free dry mass) for structure (J C-mol^-1)
mu_V <- 550000
# Chemical potentials (Energy content in ash-free dry mass) for reserve (J C-mol^-1)
mu_E <- 550000
# Molar mass of Structure (c-mols cm^-3)
M_v <- d_V/w_V
# Percentage of water
pct_h2o <- 10 
# Number of eggs per gram of wet weight
eggs_g <- 1500
# Energetic content of an egg (J egg^-1)
egg_energy <- 5 


#Save Parameters R Data
save.image(file="./input_data/DEB_Parameters_Commercial.RData") 

