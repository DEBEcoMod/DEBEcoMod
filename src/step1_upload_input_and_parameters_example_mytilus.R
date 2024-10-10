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

#Specify the model type "std" or "abj"
#model_type<-"std"
model_type<-"abj"
#Specify the type of species "wild" or "commercial"
species<-"wild"
#species<-"commercial"

#Main parameters
#Initial maximum surface area-specific assimilation rate
p_Am0<-7.13375
# Assimilation efficiency (Adimensional) 
kap_X<-0.7
# Initial Energy conductance (cm2 d^-1)
v0<-0.023785
# Fraction of energy allocated to maintenance plus growth (Adimensional)
kap<-0.47679
# Fraction of energy used for reproduction (Adimensional)
kap_R<-0.95
# Volume specific maintenance costs (J cm^-3 d^-1)
p_M <- 9.0477
#Surface specific maintenance costs (J cm^-2 d^-1)
p_T<-0 
#Maturity maintenance rate coefficient (1/d)
k_J<-0.002
# Volume specific cost for structure (J cm^-3)
E_G <-2380.6249 
#Maturity at birth (J)
E_Hb<-5.6300e-05
#Maturity at metamorphosis (J)
E_Hj<-0.0403
#Maturity at puberty (J)
E_Hp<-1587
# Shape coefficient of adult
del_M<-0.23035
#Shape coefficient at birth
del_Mb<-0.23287


## Temperature Correction Parameters ##
# Arrhenius temperature (K)
TA <- 14821.191
# Lower boundary of thermal tolerance range (K)
TL <- 293.15
# Upper boundary of thermal tolerance range (K)
TH <- 303.2494
# Rate of decrease at lower boundary (K)
TAL <- -0.001021
# Rate of decrease at upper boundary (K)
TAH <- 31292.8769   


####Additional parameters####

## Life-history ##
#Lifespan (d)
am<-4379
#Life expectation of model species
year.v<-round(am/365)
#Desired commercial or target size (e.g. minimum landing size)
ComSize <- 5


## Size, shape and volumes ##

#Length at seeding (cm) in case of aquaculture species or length at birth in case of wild species
Ls<-2
Lb<-0.01195 

if (species == "wild") {
  L0<- Lb * del_Mb
  
} else {
  L0 <- Ls * del_Mb  
}

# Volume at seeding/birth (cm^3) 
V0<-L0^3


## Feeding ##
# Half saturation coefficient (ug L^-1)
K <- 2.9864e-06 #2.973952

# Food energy conversion factor (J mg^-1)
mu_X <- 525000


## Energetic ##
#Zoom factor (Adimensional)
z<-0.37593
#Acceleration factor at f=1 (Adimensional)
s_M0<-8.93448
# Maximum storage density or reserve capacity (J cm^-3)
E_m <- 299.9265
# Initial scaled reserve (a value between 0 and 1, 0 = starved, 1 = top condition)
e0<-1
# Maximum temperature for spawning (°C)
Tmax.Spawn <- 22
# Minimum temperature for spawning (°C)
Tmin.Spawn <- 10
# Cumulative daydegree value after threshold for spawning
DDtrigger <- 420 


## Mass balance ##
# Chemical indices for organics, nO
nO <- data.frame(Food = c(1,1.8,0.5,0.15),
                 Structure = c(1,1.8,0.5,0.15),
                 Reserves = c(1,1.8,0.5,0.15),
                 Products = c(1,1.8,0.5,0.15))


## Mass conversions ##
# Specific densities for structure (g cm^-3)
d_V <- 0.09
# Specific densities for reserve (g cm^-3)
d_E <- 0.09
#Molecular weights for structure (g mol^-1)
w_V <- 12*nO[1,2] + 1*nO[2,2] +16*nO[3,2] +14*nO[4,2]
#Molecular weights for reserve (g mol^-1)
w_E <- 12*nO[1,3] + 1*nO[2,3] +16*nO[3,3] +14*nO[4,3]
# Chemical potentials (Energy content in ash-free dry mass) for structure (J C-mol^-1)
mu_V <- 500000
# Chemical potentials (Energy content in ash-free dry mass) for reserve (J C-mol^-1)
mu_E <- 550000
# Molar mass of Structure (c-mols cm^-3)
M_v <- d_V/w_V
# Percentage of water
pct_h2o <- 9
# Number of eggs per gram of wet weight
eggs_g <- 1500
# Energetic content of an egg (J egg^-1)
egg_energy <- 5 


#Save Parameters R Data
save.image(file="./input_data/DEB_Parameters_Mytilus.RData") 

