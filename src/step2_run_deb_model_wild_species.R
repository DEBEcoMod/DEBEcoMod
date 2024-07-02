##############################################################################################################################################################
#Load necessary libraries
library(zoo)
library(data.table)

rm(list=ls()) #clean R environment 
options(scipen = 999) #adopt scientific notation

#Load Parameters R Data
load("C:/Users/gabri/Dropbox/Pacchetto DEB R/DEB Tutorial/Step1_Upload_Input_Data/DEB_Parameters_Wild.RData")
##############################################################################################################################################################

### Step II - Run DEB Model
#from here on the code compute the model for each pixel through a for cycle, no need to change anything

#Create Empty 3 dimensional array in which upload DEB outputs (x= time series steps of the model (year), y= number of model ouputs, z= spatial point)
results<-array(NA, c(year.v, 14, ncol(temperaturecsv)))
colnames(results)<-c("pixel","year","Length.cm","WetWeight.g",
                     "N.eggs","N.ReprEvents","TRO","Faeces.g/h",
                     "N.eggs.CUM","N.ReprEvents.CUM","TRO.CUM","Faeces.g/h.CUM",
                     "TimeToMaturity.days", "LengthAtMaturity.cm")

system.time(
  for (pix in 1:ncol(temperaturecsv)) {
    
    # Body temperature (K and C)
    Tbody <- coredata(temperaturecsv[,pix]) + 273
    TbodyC <- coredata(temperaturecsv[,pix])
    total.rows<-length(Tbody)
    # Time (1 HOUR time steps)
    time <- seq.int(1, total.rows, length=total.rows)
    # Food density for adults (ugChl-a l^-1)
    X <- coredata(foodcsv[,pix])
    #Reference temperature
    reftemp <- 273+20 
    
    # Temperature Correction Function
    TC <- exp(TA * (1/reftemp-1/Tbody)) / (1 + exp(TAL * (1/Tbody-1/TL)) + exp(TAH *(1/TH-1/Tbody)))
    
    
    #### 2. Specify Equations for DEB fluxes, state variables, and outputs ####
    
    
    #### Create DEB vectors where specify parameter values ####
    f <- rep(NA,total.rows)
    J_Xm_TC <- rep(NA,total.rows)
    Jx <- rep(NA,total.rows)
    Jx_cells <- rep(NA,total.rows)
    p_Am <- rep(NA,total.rows)
    pA <- rep(NA,total.rows)
    p_M_TC <- rep(NA,total.rows)
    pM <- rep(NA,total.rows)
    pC <- rep(NA,total.rows)
    pD <- rep(NA,total.rows)
    pG <- rep(NA,total.rows)
    dE.dt <- rep(NA,total.rows)
    E <- rep(NA,total.rows)
    E[1] <- 1500
    dV.dt <- rep(NA,total.rows)
    V <- rep(NA,total.rows)
    V[1] <- V_b  #For wild organisms.
    V_energy <- rep(NA,total.rows)
    L <- rep(NA,total.rows)
    V_mols <- rep(NA,total.rows)
    v <- rep(NA,total.rows)
    kM <- rep(NA,total.rows)
    vB <- rep(NA,total.rows)
    pJ <- rep(NA,total.rows)
    pR <- rep(NA,total.rows)
    dER.dt <- rep(NA,total.rows)
    ENrep<- rep(NA,total.rows)
    ENrep[1]<-0
    ER <- rep(NA,total.rows)
    ER[1] <- 0
    ER2 <- rep(NA,total.rows)
    ER2[1] <- 0
    dgam.dt <- rep(NA,total.rows)
    gam <- rep(NA,total.rows)
    gam[1] <- 0
    gam_J <- rep(NA,total.rows)
    ngt <- rep(NA,total.rows)
    E_scaled <- rep(NA,total.rows)
    lys <- rep(NA,total.rows)
    dDD.dt <- rep(NA,total.rows)
    DD <- rep(NA,total.rows)
    DD[1] <- 0
    spawn <- rep(NA,total.rows)
    gon <- rep(NA,total.rows)
    pJ_pR <- rep(NA,total.rows)
    gonWw <- rep(NA,total.rows)
    resWw <- rep(NA,total.rows)
    bodyWw <- rep(NA,total.rows)
    bodyDw <- rep(NA,total.rows)
    somaDw <- rep(NA,total.rows)
    rep_events <- rep(NA,total.rows)
    V_rescaled <- rep(NA,total.rows)
    E_rescaled <- rep(NA,total.rows)
    mature <- rep(NA,total.rows)
    Length <- rep(NA,total.rows)
    k_pC <- rep(NA,total.rows)
    S_M <- rep(NA,total.rows)
    ukpC <- rep(NA,total.rows)
    Faeces.g.h <- rep(NA,total.rows)
    
    
    #### LOOP to run DEB models for each individual ####
    
    ## Calculate the variables using a loop that integrates them all ####
    
    for (t in 1:total.rows) {
      
      ## 1. Feeding response ####
      f[t] <- X[t]/ (X[t] + K) # Scaled Functional Response (-)
      J_Xm_TC[t] <- J_Xm * TC[t] # Max surface-specific ingestion rate Temp Corr (mg h-1 cm^-2)
      Jx[t] <- J_Xm_TC[t] * f[t] * V[t]^(2/3) # Ingestion rate (mg h^-1)
      Jx_cells[t] <- Jx[t]/(20.3*0.57) # Cells ingestion rate (mg h^-1)
      
      ## 2. Assimilation ####
      p_Am[t] <- J_Xm_TC[t] * ae * mu_X # Max surface-specific assimilation rate (J h^-1 cm^-2)
      pA[t] <- f[t] * p_Am[t] * V[t]^(2/3) # Assimilation rate (J h^-1)
      
      ## 3. Maintenance ####
      p_M_TC[t] <- p_M * TC[t] # Volume-specific maintenance cost (J cm^-3 h^-1)
      pM[t] <- p_M_TC[t] * V[t] # Somatic maintenance rate (J h^-1)
      
      ## 4.a Energy flow to Utilization, Dissipation, and Growth ####
      pC[t] <- (E[t]/ (E_G + kap*E[t])) * (((E_G*p_Am[t]^(2/3))/E_m) +  p_M_TC[t]*V[t])  # Utilization (J h^-1)
      
      ## 5. Reserve dynamics ####
      dE.dt[t] <- p_Am[t] * V[t]^(-1/3) * (f[t] - E[t]/E_m) # Volume-specific Reserve energy flow (J cm^-3 h^-1)
      
      ## 6. Size dynamics ####
      # Structural volume (cm-3)
      dV.dt[t] <- (kap * p_Am[t] * (E[t]/E_m) * V[t]^(2/3) - p_M_TC[t] * V[t]) / (kap * E[t] + E_G) # Change in Structure (cm3 h-1)
      dV.dt[t] <- ifelse(dV.dt[t] < 0, ifelse(E[t] > 0, 0, dV.dt[t]), dV.dt[t])
      V_energy[t] <- V[t] * E_G # Structural energy (J)
      L[t] <- V[t]^(1/3) / del_M # Shell length (cm)
      V_mols[t] <- V[t] * M_v # Structure c-mols (C-mols)
      
      ## 7. Compound parameters ####
      v[t] <- p_Am[t] / E_m # Energy conductance (cm2 h^-1)
      kM[t] <- p_M_TC[t] / E_G # Maintenance rate coefficient (h^-1)
      vB[t] <- ((1/3) * (p_M_TC[t]/(kap * f[t] * E_m + E_G))) / 24 # von Bertalanffy growth coefficient (h^-1)
      
      ## 8. Maturity and Reproduction ##
      pJ[t] <- min(V[t], V_p) * p_M_TC[t] * ((1-kap)/ (kap)) # Maturity maintenance (J h^-1)
      pR[t] <- (1 - kap) * pC[t] - pJ[t] # Reproduction (J h^-1)
      
      # Reproduction buffer, ER (J) - [ER] (J cm^-3)
      dER.dt[t] <- ifelse(V[t] > V_p, pR[t], 0) # Energy flow to reproduction buffer, dER/dt (J) (flow to maintenance and reproduction)
      ENrep[t+1]<-ifelse(gam[t]==0,(dER.dt[t]+ENrep[t]),0) # Energy to reproduction
      dgam.dt[t] <- dER.dt[t] * kap_R # Gamete production, dgam/dt (J) - [gametes] (J cm^-3)
      gam_J[t] <- gam[t] * V[t] # Gametes (J) - [gametes] (J)
      ngt[t] <- (V[t]^(1/3)) * (p_M_TC[t]) / (kap * p_Am[t]) # No growth threshold, ngt (-)
      E_scaled[t] <- E[t] / V[t] / E_m # Scaled energy density (-)
      lys[t] <- ifelse(E_scaled[t] <= ngt[t], (pM[t]-kap*pC[t])/kap_R, 0) # Lysis
      spawn[t] <- ifelse(DD[t] > DDtrigger , gam[t] * V[t], 0)  # Spawning, spawn (J)
      dDD.dt[t] <- ifelse(TbodyC[t] > Tmin.Spawn, (TbodyC[t+1] - Tmin.Spawn) / 24, 0) # Daydegrees, DD
      gon[t] <- (ENrep[t] + gam[t]) * (1 - egg_cost) # Gonads, gon (J)
      pJ_pR[t] <- dER.dt[t] + pJ[t] # pJ + pR (i.e. (1-kap)*pC)), (J h^-1)
      
      ## 9. Mass ####
      gonWw[t] <- gon[t] / (eggs_g * egg_energy)  # Gonads wet weight (g)
      resWw[t] <- (E[t] * V[t]) / E_G # Storage (reserve) wet weight (g)
      bodyWw[t] <- V[t] + gonWw[t] + resWw[t] # Soft tissue wet mass (g)
      bodyDw[t] <- bodyWw[t] * pct_h2o/100 # Total dry mass (g)
      somaDw[t] <- (V[t] + resWw[t]) * pct_h2o/100 # Soma dry mass (g) - structure + reserves
      
      
      ## 12. For plotting ####
      rep_events[t] <- ifelse(spawn[t] > 0, 3.5, -10) # Reproductive events
      V_rescaled[t] <- V[t] * 10000 # Volume re-scaled for plotting
      E_rescaled[t] <- E[t] / 1000 # Energy density re-scaled for plotting
      mature[t] <- ifelse(V[t] > V_p, 0, 1) # Maturity
      
      
      ## 13. Step-forward variables (State variables) ####
      E[t+1] <- E[t] + dE.dt[t] # Energy density (J cm^-3)
      V[t+1] <- V[t] + dV.dt[t] # Structure (cm^-3)
      ER[t+1] <- ER[t] + dER.dt[t] - gam[t] # Reproduction buffer, ER (J) - [ER] (J cm^-3)
      ER[t+1] <- ifelse(V[t] >= V_p, ER[t], 0)
      ER2[t+1] <- ER2[t] + dER.dt[t] 
      gam[t+1] <- ifelse(spawn[t] > 0, 0, ifelse(TbodyC[t+1]<Tmax.Spawn, ((dER.dt[t] + ENrep[t])*kap_R+gam[t]), gam[t])) # Gametes
      DD[t+1] <- ifelse(spawn[t] > 0, 0, DD[t] + dDD.dt[t]) # Daydegrees, DD
      
      #14. Faeces g/h
      Faeces.g.h[t] <- (Jx_cells[t] * (1-ae))/1000 #Use this one to compute for the whole life cycle
      
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
    N.ReprEvents=rep(NA,year.v)
    TRO=rep(NA,year.v)
    Faeces=rep(NA,year.v)
    N.eggs.CUM=rep(NA,year.v)
    N.ReprEvents.CUM=rep(NA,year.v)
    TRO.CUM=rep(NA,year.v)
    Faeces.CUM=rep(NA,year.v)
    TimeToMaturity.days=rep(NA,year.v)
    LengthAtMaturity.cm=rep(NA,year.v)
    
    
    #estimate 
    for(y in 1:year.v)
    {
      single.year[y]<-total.rows/year.v*y
      
      pixel[y]=colnames(temperaturecsv)[pix]
      year[y]=y
      Length.cm[y]= L[single.year[y]]
      WetWeight.g[y]=bodyWw[single.year[y]]
      N.eggs.CUM[y]=(sum(spawn[1:single.year[y]]))/0.0019
      N.ReprEvents.CUM[y]=sum(ifelse(rep_events[1:single.year[y]]==3.5,1,0))
      TRO.CUM[y]=N.eggs.CUM[y]/N.ReprEvents.CUM[y]
      Faeces.CUM[y]=sum(Faeces.g.h[1:single.year[y]])
      TimeToMaturity.days[y]=sum(mature)/24
      LengthAtMaturity.cm[y]=L[sum(mature)]
      }
    
    #FOR INDIVIDUAL YEAR, NOT CUMULATIVE 
    
    for (yy in 0:(year.v-1)){
      N.eggs[yy+1]=sum(spawn[((yy/year.v)*total.rows):(((yy+1)/year.v)*total.rows)])/0.0019  
      N.ReprEvents[yy+1]=sum(ifelse(rep_events[((yy/year.v)*total.rows):(((yy+1)/year.v)*total.rows)]==3.5,1,0))
      TRO[yy+1]=N.eggs[yy+1]/N.ReprEvents[yy+1]
      Faeces[yy+1]=sum(Faeces.g.h[((yy/year.v)*total.rows):(((yy+1)/year.v)*total.rows)])
    }
    
    #table of results
    results[,,pix]<-cbind(
      pixel,    
      year,
      round(Length.cm,2),
      round(WetWeight.g,2),
      round(N.eggs,0),
      round(N.ReprEvents,0),
      round(TRO,0),
      round(Faeces,2),
      round(N.eggs.CUM,0),
      round(N.ReprEvents.CUM,0),
      round(TRO.CUM,0),
      round(Faeces.CUM,2),
      round(TimeToMaturity.days,0),
      round(LengthAtMaturity.cm,2))
    
    print(data.frame(Completion=c(pix,pix/ncol(temperaturecsv)*100),
                     row.names=c("Pixel","Percentage")))
    
  })    

#see first rows of the Output matrix for a single pixel (3 dimensional array: x= time series steps (year) y= model outputs z=spatial point)
head(results[,,1])#change the number inside the square brackets to see the output of a different pixel

#Save model output of year of interest
#Create a two dimensional matrix (rows= spatial point, columns= model outputs )
finalcsv<-(matrix(NA,nrow=ncol(temperaturecsv),ncol = 14))  
colnames(finalcsv)<-c("pixel","year","Length.cm","WetWeight.g",
                      "N.eggs","N.ReprEvents","TRO","Faeces.g/h",
                      "N.eggs.CUM","N.ReprEvents.CUM","TRO.CUM","Faeces.g/h.CUM",
                      "TimeToMaturity.days", "LengthAtMaturity.cm")

#Assign model outputs of the year of interest to the matrix
for (riga in 1:ncol(temperaturecsv)){
  finalcsv[riga,]<-results[year.v,,riga]}

finalcsv<-as.data.frame(finalcsv)

for (num in 2:(ncol(finalcsv)-1)){
  finalcsv[,num]<-as.numeric(as.character(finalcsv[,num]))}

#visualize the output file before writing it
View(finalcsv) 

#Save output matrix in a csv file
setwd("C:/Users/gabri/Dropbox/Pacchetto DEB R/DEB Tutorial/Step2_Run_DEB_Model/Model_Output")

data.table::fwrite(finalcsv,"DEB_Output_Wild_test.csv",sep = ";",dec = ".") #write the csv file inside the specified working directory
