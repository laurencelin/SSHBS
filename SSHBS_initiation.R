arg=commandArgs(T)
library(dplyr)
source('SSHBS_library.R')
source('https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_dailytimeseries3.R')

# https://www.researchgate.net/publication/277670582_R_function_to_convert_a_list_of_dataframes_into_a_3-D_array
list2ary = function(input.list){  #input a list of lists
  rows.cols <- dim(input.list[[1]])
  sheets <- length(input.list)
  output.ary <- array(unlist(input.list), dim = c(rows.cols, sheets))
  colnames(output.ary) <- colnames(input.list[[1]])
  row.names(output.ary) <- row.names(input.list[[1]])
  return(output.ary)    # output as a 3-D array
}
vectorMin = function(aa,bb){ ifelse(aa>bb,bb,aa) }#function
vectorMax = function(aa,bb){ ifelse(aa>bb, aa,bb) }#function

zoneInfo = read.csv('inputs/GIS_channel_model_setup.csv') %>% 
  group_by(zoneID) %>%
  summarise(
    len = sum(dist, na.rm=T),
    laterFrac = sum(laterFrac,na.rm=T),
    PARinputFile = PARinputFile[1],
    LateralNFlexFile = LateralNFlexFile[1], # mgNpL
    LateralPFlexFile = LateralPFlexFile[1], # ugPpL
    LateralQFlexFile = LateralQFlexFile[1], # Lps
    waterTempFile = waterTempFile[1],
    ChannelVelocityFile = ChannelVelocityFile[1],
    detritusInputFile = detritusInputFile[1]
  )
numZone = dim(zoneInfo)[1]
numZoneIndex = seq_len(numZone)
#network = read.csv('inputs/stream_network.csv') %>% arrange(zoneID) %>% mutate(zindex = match(upstreamZoneID,zoneInfo$zoneID))


## -- building channel zone -- traveling time distribution	
Morphology_matrix=list()	
maxChannelBenthicArea = rep(NA, numZone)
for(i in seq_len(numZone)){
	
	ratingTable = read.csv(zoneInfo$ChannelVelocityFile[i])
	maxChannelBenthicArea[i] = max(ratingTable[,'wettedBperimeter'])*zoneInfo[i,'len']
	ratingTable$Corrected_volume = ratingTable$volume*0.1*as.numeric(zoneInfo[i,'len'])
		
	check = (ratingTable$discharge/ratingTable$Corrected_volume)[-1]
	if(sum(check>0.95)>0) print(paste(i,check))
	
	Morphology_matrix[[ i ]] <- new('CrossSectionQrelation',
		Q2Vol = approxfun(ratingTable[,'discharge'], ratingTable[,'volume']),
		Vol2Z = approxfun(ratingTable[,'volume'], ratingTable[,'depth']),
		Vol2V = approxfun(ratingTable[,'volume'], ratingTable[,'mvelocity']),
		Vol2WP = approxfun(ratingTable[,'volume'], ratingTable[,'wettedBperimeter']),
		Vol2Q = approxfun(ratingTable[,'volume'], ratingTable[,'discharge']),
		Vol2CA = approxfun(ratingTable[,'volume'], ratingTable[,'crossarea']),
		maxQ = max(ratingTable[,'discharge']),
		maxVol = max(ratingTable[ratingTable[,'bank']==0,'volume'])
		)#list
}#

# lateralResource_matrix = zoneEnv_matrix = zoneRes_matrix
zoneRes_matrix = list2ary(lapply(seq_len(dim(zoneInfo)[1]),function(ii){
  
  # For Inner Join
  return <- Reduce(
    function(x, y, ...) merge(x, y, ...), 
    list(
      # usually daily
      Lps = read.csv(zoneInfo$LateralQFlexFile[ii]) %>% slice(rep(1:n(), each = 24)) %>% mutate(hour=rep(1:24,n()/24)),
      mgNpL = read.csv(zoneInfo$LateralNFlexFile[ii]) %>% slice(rep(1:n(), each = 24)) %>% mutate(hour=rep(1:24,n()/24)),
      ugPpL = read.csv(zoneInfo$LateralPFlexFile[ii]) %>% slice(rep(1:n(), each = 24)) %>% mutate(hour=rep(1:24,n()/24)),
      detritus = read.csv(zoneInfo$detritusInputFile[ii]) %>% slice(rep(1:n(), each = 24)) %>% mutate(hour=rep(1:24,n()/24)),
      # usually hourly
      PARps = read.csv(zoneInfo$PARinputFile[ii]),
      temp = read.csv(zoneInfo$waterTempFile[ii])
    )
  ) %>% mutate(
    baseflow = fivedayblockbaseflow(flow*0.01)*zoneInfo$laterFrac[ii],
    baseflowNO3 = no3_mgNpl * baseflow*1000,
    baseflowNH4 = 0,
    baseflowPO4 = po4_ugPpl * baseflow,
    stormQ = flow*0.01 - baseflow,
    stormQ = ifelse(stormQ>0,stormQ,0)*zoneInfo$laterFrac[ii],
    stormQNO3 = no3_mgNpl * stormQ*1000,
    stormQNH4 = 0,
    stormQPO4 = po4_ugPpl * stormQ,
    PARf = PAR2PARF(PARps)
  ) %>% arrange(
    year, month, day, hour
  )# return 

}))#
# zoneRes_matrix[1:48,'hour',] # change the time order


# time, variable, zone
STSstepOption = c(1, 2, 4, 8, 24, 48)	
timemax_BTS = dim(zoneRes_matrix)[1]
timemax_STS = 3600

####################################### initial ###########################################			
iniDischarge = cumsum(zoneRes_matrix[1,'baseflow',] + zoneRes_matrix[1,'stormQ',] )
meanDischarge = cumsum(colMeans(zoneRes_matrix[,'baseflow',])+colMeans(zoneRes_matrix[,'stormQ',]))
iniWater = sapply(numZoneIndex,function(i){ min(Morphology_matrix[[i]]@Q2Vol(iniDischarge[i]),Morphology_matrix[[i]]@maxVol,na.rm=T) })	
meanWater = sapply(numZoneIndex,function(i){ min(Morphology_matrix[[i]]@Q2Vol(meanDischarge[i]), Morphology_matrix[[i]]@maxVol,na.rm=T) }) # assume baseflow is 2 L/s

hydraulic_matrix = matrix(0, numZone,length(hydraulic_variable)); colnames(hydraulic_matrix)= hydraulic_variable
hydraulicHOLD=sapply(numZoneIndex,function(i){ c(
	depth = Morphology_matrix[[i]]@Vol2Z(iniWater[i]), #1
	velocity = Morphology_matrix[[i]]@Vol2V(iniWater[i]), #2
	waterPerimeter = Morphology_matrix[[i]]@Vol2WP(iniWater[i]), #3
	discharge = Morphology_matrix[[i]]@Vol2Q(iniWater[i]), #4
	crossArea = Morphology_matrix[[i]]@Vol2CA(iniWater[i]), # 5
	maxDischarge = Morphology_matrix[[i]]@maxQ, #5 -> 6
	maxWater = Morphology_matrix[[i]]@maxVol*0.99, # 7
	alphaA_As =  (Morphology_matrix[[i]]@Vol2CA(meanWater[i]) / (Morphology_matrix[[i]]@Vol2WP(meanWater[i])+3) * (1-exp(-p_d*1.5))*p0/p_d) * alpha_exchange,  
	longtermWaterDepth = Morphology_matrix[[i]]@Vol2Z(meanWater[i])
	) })
hydraulic_matrix[,'depth'] = hydraulicHOLD[1,]
hydraulic_matrix[,'velocity'] = hydraulicHOLD[2,] #m/s
hydraulic_matrix[,'wetBenthicArea'] = zoneInfo$len*hydraulicHOLD[3,]
hydraulic_matrix[,'wetBenthicArea_1'] = 1.0 / hydraulic_matrix[,'wetBenthicArea']
hydraulic_matrix[,'crossarea'] = hydraulicHOLD[5,]
hydraulic_matrix[,'maxVol'] = hydraulicHOLD[7,]
hydraulic_matrix[,'discharge'] = hydraulicHOLD[4,]
hydraulic_matrix[,'As'] = hydraulicHOLD[8,] 


####################################### in-channel process parameters ###########################################
bioprocess_matrix = matrix(0, numZone,length(bioprocess_variable)); colnames(bioprocess_matrix)= bioprocess_variable		
bioprocess_matrix[,'algalc'] = 50*1000; #50 gC/m2 as inintial; model unit is mgC/m2
bioprocess_matrix[,'algaln'] = bioprocess_matrix[,'algalc']*0.5*(algal_parameter['maxNC']+algal_parameter['minNC'])
bioprocess_matrix[,'algalp'] = bioprocess_matrix[,'algalc']*0.5*(algal_parameter['maxPC']+algal_parameter['minPC'])
bioprocess_matrix[,'algalCN'] = bioprocess_matrix[,'algalc']/bioprocess_matrix[,'algaln']
bioprocess_matrix[,'algalCP'] = bioprocess_matrix[,'algalc']/bioprocess_matrix[,'algalp']
bioprocess_matrix[,'algalselfF'] = 1
bioprocess_matrix[,'algalnutrientF'] = 0
bioprocess_matrix[,'algallightF'] = 1
bioprocess_matrix[,'algalNuptakeF'] = 1
bioprocess_matrix[,'algalPuptakeF'] = 1
bioprocess_matrix[,'colonizedBenthicArea']=hydraulic_matrix[,'wetBenthicArea']*0.9
bioprocess_matrix[,'colonizedBenthicArea_1'] = 1.0/bioprocess_matrix[,'colonizedBenthicArea']

hyportheicResource_matrix = matrix(0, numZone,length(hyportheicResource_variable)); colnames(hyportheicResource_matrix)= hyportheicResource_variable
benthicResource_matrix = matrix(0, numZone,length(benthicResource_variable)); colnames(benthicResource_matrix)= benthicResource_variable	
	

####################################### model computational settings	
columnResource_matrix = matrix(0,numZone,length(columnResource_variable)); colnames(columnResource_matrix)= columnResource_variable
columnResource_matrix[,'discharge'] = iniWater
columnResource_matrix[,'waterNO3'] = (zoneRes_matrix[1,'baseflowNO3',]+zoneRes_matrix[1,'stormQNO3',])/(zoneRes_matrix[1,'baseflow',] + zoneRes_matrix[1,'stormQ',])*iniWater
columnResource_matrix[,'waterNH4'] = 0
columnResource_matrix[,'waterPO4'] = (zoneRes_matrix[1,'baseflowPO4',]+zoneRes_matrix[1,'stormQPO4',])/(zoneRes_matrix[1,'baseflow',] + zoneRes_matrix[1,'stormQ',])*iniWater

hyporheicResource_matrix = matrix(0,numZone,length(columnResource_variable)+5); 
colnames(hyporheicResource_matrix)=c(columnResource_variable,'baseNO3','baseNH4','baseQ','baseNO3balance','baseNH4balance')
hyporheicResource_matrix[,'discharge'] = hydraulic_matrix[,'As'] * zoneInfo$len
hyporheicResource_matrix[,'waterNO3'] = columnResource_matrix[,'waterNO3']/columnResource_matrix[,'discharge'] * hyporheicResource_matrix[,'discharge']
hyporheicResource_matrix[,'waterNH4'] = 0 * hyporheicResource_matrix[,'discharge']
hyporheicResource_matrix[,'waterPO4'] = 200 * hyporheicResource_matrix[,'discharge']
hyporheicResource_matrix[,'baseNO3'] = hyporheicResource_matrix[,'waterNO3']
hyporheicResource_matrix[,'baseNH4'] = hyporheicResource_matrix[,'waterNH4']
hyporheicResource_matrix[,'baseNO3balance'] = 0 * hyporheicResource_matrix[,'discharge']
hyporheicResource_matrix[,'baseNH4balance'] = 0 * hyporheicResource_matrix[,'discharge']
		
STSexport = matrix(0,numZone,length(columnResource_variable)); colnames(STSexport)= columnResource_variable
STSaggregateVar = matrix(0,numZone,length(columnResource_variable)); colnames(STSaggregateVar)= columnResource_variable
STSinput = matrix(0,numZone,length(columnResource_variable)); colnames(STSinput)= columnResource_variable
STSinputTrack = matrix(0,numZone,length(columnResource_variable)); colnames(STSinputTrack)= columnResource_variable
STSLateralinputTrack = matrix(0,numZone,length(columnResource_variable)); colnames(STSLateralinputTrack)= columnResource_variable
exportProp = rep(0,numZone)
hyporheicDenitrification = rep(0,numZone)

RiparianResource_matrix = matrix(0,numZone,length(columnResource_variable)); colnames(RiparianResource_matrix)= columnResource_variable
cellspace = rep(0,numZone)
exRiparian = rep(0,numZone)
concChange = rep(0,numZone)
concMassChange = rep(0,numZone)
dailyAverageWetBenthicArea = rep(0, numZone);

detritusMatrix = matrix(0, numZone, length(detritusName)); colnames(detritusMatrix)=detritusName;
detritusMatrix[,'microbes'] = (zoneRes_matrix[1,'detritus1c',]+zoneRes_matrix[1,'detritus2c',]+zoneRes_matrix[1,'detritus3c',])*0.09
detritusMatrix[,'miner'] = (zoneRes_matrix[1,'detritus1c',]+zoneRes_matrix[1,'detritus2c',]+zoneRes_matrix[1,'detritus3c',])*0.09 * 0.1
detritusMatrix[,'detritus1c'] = zoneRes_matrix[1,'detritus1c',]
detritusMatrix[,'detritus2c'] = zoneRes_matrix[1,'detritus2c',]
detritusMatrix[,'detritus3c'] = zoneRes_matrix[1,'detritus3c',]
detritusMatrix[,'detritus1n'] = zoneRes_matrix[1,'detritus1n',]
detritusMatrix[,'detritus2n'] = zoneRes_matrix[1,'detritus2n',]
detritusMatrix[,'detritus3n'] = zoneRes_matrix[1,'detritus3n',]
detritusMatrix[,'detritus1p'] = zoneRes_matrix[1,'detritus1c',]/detritusCP
detritusMatrix[,'detritus2p'] = zoneRes_matrix[1,'detritus2c',]/detritusCP
detritusMatrix[,'detritus3p'] = zoneRes_matrix[1,'detritus3c',]/detritusCP

## -------- model holders / containers (do not touch)
downstreamCond = rep(T, numZone); downstreamCond[1]=F
upstreamCond = rep(T, numZone); upstreamCond[length(upstreamCond)]=F
upstreamIndex = (numZoneIndex)-1
tmpCond = rep(F,numZone)
hold = rep(0, numZone)
holdII = rep(0, numZone)
holdIII = rep(0, numZone)
holdIV = rep(0, numZone)
holdV = rep(0, numZone)
sumPotentialNO3= rep(0, numZone)
sumPotentialNH4= rep(0, numZone)
sumPotentialPO4= rep(0, numZone)
concNO3= rep(0, numZone)
concNH4= rep(0, numZone)
concPO4= rep(0, numZone)
hypoconcNO3= rep(0, numZone)
hypoconcNH4= rep(0, numZone)
hypoconcPO4= rep(0, numZone)
denitrify_concNO3= rep(0, numZone)
algalQ10f= rep(0, numZone)
Pot_algalupNO3= rep(0, numZone)
Pot_algalupNH4= rep(0, numZone)
Pot_algalupP= rep(0, numZone)
bioQ10f= rep(0, numZone)
algalhabitatDmgIndex= rep(0, numZone)
algalhabitatDmg= rep(0, numZone)
Pot_denitrification= rep(0, numZone)
Pot_denitrificationFrac = rep(0, numZone)
Pot_nitrification= rep(0, numZone)
Pot_nitrificationFrac= rep(0, numZone)
sumPotentialNO3= rep(0, numZone)
sumPotentialNH4= rep(0, numZone)
sumPotentialPO4= rep(0, numZone)
Pot_immobNO3= rep(0, numZone)
Pot_immobNH4= rep(0, numZone)
Pot_immobPO4= rep(0, numZone)
NO3limit_scalar= rep(0, numZone)
NH4limit_scalar= rep(0, numZone)
Plimit_scalar= rep(0, numZone)
miner_class2decayRatio = rep(0, numZone)
miner_class3decayRatio = rep(0, numZone)
immo_class2decayRatio = rep(0, numZone)
immo_class3decayRatio = rep(0, numZone)




