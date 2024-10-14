##############################################################################################################################################################
#Load necessary libraries
library(zoo)
library(data.table)

rm(list=ls()) #clean R environment 
options(scipen = 999) #add scientific notation

#Load Parameters R Data
load("./input_data/DEB_Parameters_Mytilus.RData")
##############################################################################################################################################################


### Step II - Run DEB Model
#from here on the code compute the model for each pixel through a for cycle, no need to change anything

#Create Empty 3 dimensional array in which upload DEB outputs
#(x= time series steps of the model (year), y= number of model ouputs, z= spatial point)
results<-array(NA, c(year.v, 9, ncol(temperaturecsv)))
colnames(results)<-c("pixel","year","Length.cm","WetWeight.g",
                     "N.eggs","Faeces.g/h",
                     "N.eggs.CUM","Faeces.g/h.CUM",
                     "TimeToCommercial.days")


timestep<-1/24 #Specify the desired timestep (1/24 for hourly, 1 for daily)

system.time(
  for (pix in 1:ncol(temperaturecsv)) {
    
    #### 2. Specify Equations for DEB fluxes, state variables, and outputs ####
    
    
    #### Create DEB vectors where specify parameter values ####
    
    #Environmental input
    #Body temperature (K and C)
    Tbody <- coredata(temperaturecsv[,pix]) + 273
    TbodyC <- coredata(temperaturecsv[,pix])
    
    #Reference temperature
    reftemp <- 273+20 
    #Temperature Correction Function
    TC<- exp(TA/ reftemp - TA/ Tbody) * (1 + exp(TAL/ reftemp - TAL/ TL) + exp(TAH/ TH - TAH/ reftemp))/(1 + exp(TAL/ Tbody - TAL/ TL) + exp(TAH/ TH - TAH/ Tbody))
    #Food density for adults (ugChl-a l^-1)
    Food <- coredata(foodcsv[,pix])
    
    #set length for the loop
    total.rows<-length(Tbody)
    
    #State variables
    #Initial conditions
    #Initial maximum surface area-specific assimilation rate
    #p_Am0<-z * p_Am0 #Uncomment this if p_am is not declared among AmP parameters.
    # Area-specific maximum ingestion rate (J cm^-2 h^-1)
    J_Xm <- p_Am0 / kap_X
    # Initial Energy conductance (cm2 h^-1)
    #v0<-p_Am0 / E_m #Uncomment this if v is not declared among AmP parameters.
    #Initial reserve (E0)
    E0<-e0 * p_Am0/ v0 * V0
    #Initial maturity
    E_H0<-E_Hb
    #Initial reproduction buffer
    E_R0<-0
    
    #Set energy volume vector and start with the initial volume
    V <- rep(NA,total.rows)
    V[1] <- V0
    #Set energy reserve vector and start with the initial reserve
    E <- rep(NA,total.rows)
    E[1] <- E0
    #Set reproduction buffer vector and start with the initial reproduction buffer
    ER <- rep(NA,total.rows)
    ER[1] <- E_R0
    
    #Set maturity level vector and start with the initial maturity
    EH <- rep(NA,total.rows)
    EH[1] <- E_H0
    
    #Gametes
    gam <- rep(NA,total.rows)
    gam[1] <- 0
    
    #Day degrees
    DD <- rep(NA,total.rows)
    DD[1] <- 0
    
    
    #Set empty rows of the variables to calculate inside the loop
    #Feeding response
    f <- rep(NA,total.rows)
    Jx <- rep(NA,total.rows)
    Jx_cells <- rep(NA,total.rows)
    
    #Compound parameters
    s_M<-rep(NA,total.rows)
    v <- rep(NA,total.rows)
    v[1]<-v0
    kM <- rep(NA,total.rows)
    vB <- rep(NA,total.rows)
    
    #Fluxes
    #Assimilation flux
    p_Am <- rep(NA,total.rows)
    p_Am[1] <- p_Am0
    pA <- rep(NA,total.rows)
    #Somatic maintenance flux
    pS <- rep(NA,total.rows)
    #Utilization flux 
    pC <- rep(NA,total.rows)
    #Growth flux
    pG <- rep(NA,total.rows)
    #Maturity maintenance flux
    pJ<-rep(NA,total.rows)
    #Reproduction buffer flux
    pR<-rep(NA,total.rows)
    #Maturation flux
    pH <- rep(NA,total.rows)
    #Dissipation flux
    pD <- rep(NA,total.rows)
    
    #Reserve dynamics
    dE.dt <- rep(NA,total.rows)
    dEH.dt<- rep(NA,total.rows)
    
    #Size dynamics
    dV.dt <- rep(NA,total.rows)
    L <- rep(NA,total.rows)
    
    #Maturity and Reproduction
    dER.dt <- rep(NA,total.rows)
    spawn <- rep(NA,total.rows)
    spawn[1]<-0
    DDtrigger<-365*24 #Spawning time (modify the hourly cycle of reproduction)
    dDD.dt <- rep(NA,total.rows)
    
    #Mass
    resWw <- rep(NA,total.rows)
    bodyWw <- rep(NA,total.rows)
    bodyDw <- rep(NA,total.rows)
    somaDw <- rep(NA,total.rows)
    
    #For plotting
    rep_events <- rep(NA,total.rows)
    mature <- rep(NA,total.rows)
    
    #Faeces g/h
    Faeces.g.h <- rep(NA,total.rows)
    
    
    #### LOOP to run DEB models for each individual ####
    
    ## Calculate the variables using a loop that integrates them all ####
    
    for (t in 1:total.rows) {
      
      ## 1. Feeding response ####
      f[t] <- Food[t]/ (Food[t] + K) # Scaled Functional Response (-)
      #f[t] <- 1 #Uncomment this line for testing purposes or sensitivity analyses
      Jx[t] <- J_Xm * TC[t] * f[t] * V[t]^(2/3) # Ingestion rate (mg h^-1)
      Jx_cells[t] <- Jx[t]/(20.3*0.57) # Cells ingestion rate (mg h^-1)
      
      ## 2. Compound parameters ####
      #Acceleration factor at f=1 (Adimensional)
      if (species=="wild") {
        
        if (model_type =="std") {
          s_M[t]<-1
          
        } else if (model_type =="abj"&EH[t]<=E_Hb){
          s_M[t]<-1
        }
        else if (model_type =="abj"&EH[t]>E_Hb&EH[t]<E_Hj){
          s_M[t]<-(V[t]^(1/3))/L0
        }
        else if (model_type =="abj"&EH[t]>=E_Hj){
          s_M[t]<-(V[which(EH>=E_Hj)[1]]^(1/3))/L0
        }
        else {
          s_M[t]<-(V[which(EH>=E_Hj)[1]]^(1/3))/L0
        }
      }
      else {
        s_M[t]<-s_M0
      }
      
      v[t] <- TC[t]*v0*s_M[t] # Energy conductance (cm2 h^-1)
      kM[t] <- p_M * TC[t] / E_G # Maintenance rate coefficient (h^-1)
      vB[t] <- ((1/3) * (p_M * TC[t]/(kap * f[t] * E_m + E_G))) / 24 # von Bertalanffy growth coefficient (h^-1)
      
      ## 3. Fluxes ####
      ## Assimilation ###
      p_Am[t] <- TC[t]*p_Am0*s_M[t] # Max surface-specific assimilation rate (J h^-1 cm^-2)
      pA[t] <- f[t] * p_Am[t] * V[t]^(2/3) # Assimilation rate (J h^-1)
      
      ## Somatic Maintenance ###
      pS[t] <- p_M * TC[t] * V[t]+p_T*TC[t]*V[t]^(2/3) # Somatic maintenance flux (J h^-1)
      
      ## Energy flow to Utilization, Growth, Reproduction and Dissipation ###
      pC[t] <- E[t] * ((E_G * v[t] * V[t]^(2/3) + pS[t] ) / (kap * E[t] + E_G*V[t])) # Utilization (J h^-1)
      pG[t] <- max(0,kap*pC[t] - pS[t])     #J/h growth flux
      pJ[t]<-EH[t]*k_J*TC[t] #J/h maturity maintenance flux
      pH[t]<-ifelse(EH[t]<E_Hp,max(0,(1-kap)*pC[t]-pJ[t]),0) #J/h maturation flux
      pR[t]<-ifelse(EH[t]<E_Hp,0,max(0,(1-kap)*pC[t]-pJ[t])) #new code
      pD[t] = pS[t] + pJ[t] + pH[t] + (1-kap_R) * pR[t] #J/h dissipation flux
      
      ## 4. Differential equations ####
      ## Reserve dynamics ####
      dE.dt[t] <- (pA[t]-pC[t])*timestep # Reserve energy flow (J cm^-3 h^-1)
      dEH.dt[t]<-max(0,pH[t],na.rm=T)*timestep
      
      # Change in Structure (cm3 h-1)
      dV.dt[t] <- (pG[t]/E_G)*timestep
      dV.dt[t] <- ifelse(dV.dt[t] < 0, ifelse(E[t] > 0, 0, dV.dt[t]), dV.dt[t])
      
      # Change in reproduction buffer J/h
      dDD.dt[t] <- ifelse(TbodyC[t] > Tmin.Spawn&TbodyC[t] < Tmax.Spawn, 1, 0) # Daydegrees, DD
      dER.dt[t] <- max(0,pR[t] * kap_R * timestep,na.rm=T)
     
      ## 5. Size dynamics ####
      L[t] <- V[t]^(1/3) / del_M # Structural length (cm)
      
      ## 6. Mass ####
      resWw[t] <- E[t] * w_E/mu_E # Storage (reserve) wet weight (g)
      bodyWw[t] <- V[t] + resWw[t] # Soft tissue wet mass (g)
      bodyDw[t] <- bodyWw[t] * pct_h2o/100 # Total dry mass (g)
      somaDw[t] <- (V[t] + resWw[t]) * pct_h2o/100 # Soma dry mass (g) - structure + reserves
      
      ## 7 For plotting ####
      mature[t] <- ifelse(E[t] > E_Hp, 0, 1) # Maturity
      
      ## 8 Step-forward variables (State variables) ####
      E[t+1] <- max(0,E[t] + dE.dt[t],na.rm=T) # Energy density (J cm^-3)
      V[t+1] <- V[t] + dV.dt[t] # Structure (cm^-3)
      ER[t+1] <- max(0,ER[t] + dER.dt[t]- spawn[t],na.rm=T) # Reproduction buffer, ER (J) - [ER] (J cm^-3)
      EH[t+1]<-max(0,EH[t]+dEH.dt[t],na.rm=T)
      spawn[t+1] <- ifelse(EH[t]>=E_Hp&DD[t] > DDtrigger&TbodyC[t]<Tmax.Spawn,ER[t]*kap_R, 0)
      DD[t+1] <- ifelse(spawn[t] > 0, 0, DD[t] + dDD.dt[t]) # Daydegrees, DD
      
      ## 9 Faeces g/h
      if (species=="wild") {
      Faeces.g.h[t] <- (Jx_cells[t] * (1-kap_X))/1000 #This one will compute for the whole life cycle
      }
        else {
      Faeces.g.h[t] <- ifelse(L[t]<=ComSize,(Jx_cells[t] * (1-kap_X))/1000,0) # to compute up to commercial size
      
    }
    }
    
    ## END OF MODEL LOOP ##
    
    ## END OF INDIVIDUAL LOOP ##
    
    
    #create empty vectors for results
    single.year=rep(NA,year.v)
    pixel=rep(NA,year.v)
    year=rep(NA,year.v)
    Length.cm=rep(NA,year.v)
    WetWeight.g=rep(NA,year.v)
    N.eggs=rep(NA,year.v)
    Faeces=rep(NA,year.v)
    N.eggs.CUM=rep(NA,year.v)
    Faeces.CUM=rep(NA,year.v)
    TimeToCommercial.days=rep(NA,year.v)
    
    #estimate 
    for(y in 1:year.v)
    {
      single.year[y]<-total.rows/year.v*y
      pixel[y]=colnames(temperaturecsv)[pix]
      year[y]=y
      Length.cm[y]= L[single.year[y]]
      WetWeight.g[y]=bodyWw[single.year[y]]
      N.eggs.CUM[y]=(sum(spawn[1:single.year[y]]))/0.0019
      Faeces.CUM[y]=sum(Faeces.g.h[1:single.year[y]])
      TimeToCommercial.days[y]= sum(ifelse(L<ComSize,1,0))/24}
    
    #FOR INDIVIDUAL YEAR, NOT CUMULATIVE 
    
    for (yy in 0:(year.v-1)){
      N.eggs[yy+1]=sum(spawn[((yy/year.v)*total.rows):(((yy+1)/year.v)*total.rows)])/0.0019  
      Faeces[yy+1]=sum(Faeces.g.h[((yy/year.v)*total.rows):(((yy+1)/year.v)*total.rows)])
    }
    
    #table of results
    results[,,pix]<-cbind(
      pixel,    
      year,
      round(Length.cm,2),
      round(WetWeight.g,2),
      round(N.eggs,0),
      round(Faeces,2),
      round(N.eggs.CUM,0),
      round(Faeces.CUM,2),
      round(TimeToCommercial.days,0))
    
    print(data.frame(Completion=c(pix,pix/ncol(temperaturecsv)*100),
                     row.names=c("Pixel","Percentage")))
    
  })    

#see first rows of the Output matrix for a single pixel (3 dimensional array: x= time series steps (year) y= model outputs z=spatial point)
head(results[,,1])#change the number inside the square brackets to see the output of a different pixel

#For model validation and testing purposes (comment if not needed)
#Store state variables and zero-variate predictions at different life-stages (birth,juvenile,puberty,maximum)
Test.df<-data.frame(Eb=E[which(EH>=E_Hb)[1]],ERb=ER[which(EH>=E_Hb)[1]],EHb=EH[which(EH>=E_Hb)[1]],Time_b= round(which(EH>=E_Hb)[1]/24),
                    Vb=V[which(EH>=E_Hb)[1]],Lb=L[which(EH>=E_Hb)[1]],BodyWwb=bodyWw[which(EH>=E_Hb)[1]],BodyDwb=bodyDw[which(EH>=E_Hb)[1]],
                    Ej=E[which(EH>=E_Hj)[1]],ERj=ER[which(EH>=E_Hj)[1]],EHj=EH[which(EH>=E_Hj)[1]],Time_j= round(which(EH>=E_Hj)[1]/24),
                    Vj=V[which(EH>=E_Hj)[1]],Lj=L[which(EH>=E_Hj)[1]],BodyWwj=bodyWw[which(EH>=E_Hj)[1]],BodyDwj=bodyDw[which(EH>=E_Hj)[1]],
                    Ep=E[which(EH>=E_Hp)[1]],ERp=ER[which(EH>=E_Hp)[1]],EHp=EH[which(EH>=E_Hp)[1]],Time_p= round(which(EH>=E_Hp)[1]/24),
                    Vp=V[which(EH>=E_Hp)[1]],Lp=L[which(EH>=E_Hp)[1]],BodyWwp=bodyWw[which(EH>=E_Hp)[1]],BodyDwp=bodyDw[which(EH>=E_Hp)[1]],
                    Ei=E[am*24],ERi=ER[am*24],EHi=EH[am*24],Time_i= am,
                    Vi=V[am*24],Li=L[am*24],BodyWwi=bodyWw[am*24],BodyDwi=bodyDw[am*24])

data.table::fwrite(Test.df,"../example/DEB_validation_outputs.csv",sep = ";",dec = ".") #write the csv file inside the specified working directory

#Save model output of year of interest
#Create a two dimensional matrix (rows= spatial point, columns= model outputs )
finalcsv<-(matrix(NA,nrow=ncol(temperaturecsv),ncol = 9))  
colnames(finalcsv)<-c("pixel","year","Length.cm","WetWeight.g",
                      "N.eggs","Faeces.g/h",
                      "N.eggs.CUM","Faeces.g/h.CUM",
                      "TimeToCommercial.days")

#Assign model outputs of the year of interest to the matrix
for (riga in 1:ncol(temperaturecsv)){
  finalcsv[riga,]<-results[year.v,,riga]}

finalcsv<-as.data.frame(finalcsv)

for (num in 2:(ncol(finalcsv)-1)){
  finalcsv[,num]<-as.numeric(as.character(finalcsv[,num]))}

#visualize the output file before writing it
View(finalcsv) 

#Save output matrix in a csv file
data.table::fwrite(finalcsv,"../example/DEB_Output.csv",sep = ";",dec = ".") #write the csv file inside the specified working directory

