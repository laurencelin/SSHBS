library(methods)
source("https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_misc.r")
source("https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_dailytimeseries3.R")
## all mass unit in the model is "mgC" or "mgN"; all distance unit is metric, "m"

####################################### cross-section morphology
setClass('CrossSectionQrelation',slots=list(
    Q2Vol='function',
    Vol2Z='function',
    Vol2V='function',
    Vol2WP='function',
    Vol2Q='function',
    Vol2CA='function',
    maxQ='numeric',
    maxVol='numeric') )

ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)} # moving average function

yy=seq(0,3000); zz = 1/ ( 1+exp(  -(yy-300)/50 ) ); PAR2PARF = splinefun(yy,zz) # algal light response (Webster, Newbold, and Lin, 2016)

####################################### define constants
algal_param = list()
biochem_param = list()
detritus_param = list()
other_param = list()

# Mulholland et al. 2008 (further defined from below)
algal_param$Mulholland_pow_denitr = -0.493
algal_param$Mulholland_max_denitr = 1.059254e-05
biochem_param$Mulholland_pow_denitr = -0.493
biochem_param$Mulholland_max_denitr = 1.059254e-05
detritus_param$Mulholland_pow_denitr = -0.493
detritus_param$Mulholland_max_denitr = 1.059254e-05

# Lin and Webster 2014; ???
biochem_param$max_nitrif = 0.00002339181/17
detritus_param$growthRate_immob = 4.456975e-06 # log(2)/(1.8*3600*24)
detritus_param$growthRate_miner = 3.342729e-06 # log(2)/(2.4*3600*24)
detritus_param$baseR = 1.16e-7
detritus_param$gcoef = 0.5

# Webster et al. 2009
algal_param$umaxP = 8.460648e-08 * 2
algal_param$Phalf = 2 # mgP/m3 --> 2µgP/L

# Cross et al. 2005 # from molar ratio to mass ratio
algal_param$maxCN = (1/0.0606) *12/14
algal_param$maxNC = 1.0 / algal_param$maxCN
algal_param$maxCP = (1/0.00377) *12/31
algal_param$maxPC = 1.0 / algal_param$maxCP
algal_param$minCN = (10.404) *12/14
algal_param$minNC = 1.0 / algal_param$minCN
algal_param$minCP = (123.817) *12/31
algal_param$minPC = 1.0 / algal_param$minCP
detritus_param$immob_cn = 7.0
detritus_param$immob_nc = 1.0 / detritus_param$immob_cn
detritus_param$immob_cp = 188
detritus_param$immob_pc = 1.0 / detritus_param$immob_cp
detritus_param$miner_cn = 5.0
detritus_param$miner_nc = 1.0 / detritus_param$miner_cn
detritus_param$miner_cp = 20
detritus_param$miner_pc = 1.0 / detritus_param$miner_cp

# Webster, Newbold, and Lin 2016
algal_param$gcoef = 0.8
algal_param$selflimitcoef = 0.0015 # m2/mgC
algal_param$maxgrowthRate = 2*2*1.0/86400 # 1/s; daily to hourly; half day is in dark; convert NPP to GPP]
algal_param$mineralRate = 0.01/86400 # 1/s
algal_param$Q10 = 2
algal_param$moralityRate = 0.004*algal_param$maxgrowthRate
biochem_param$Q10 = 2

# from RHESSys, SSURGO, ...
other_param$detritusCN = 40
other_param$detritusCP = 444.5
other_param$p0 = 0.485 # soil 8 silt_loam
other_param$p_d = 1/4000

# other
STSstepOption = c(1, 2, 4, 8, 24, 48)

####################################### read in header file
arg=commandArgs(T)
arg=c('proj=./example','header=SSHBS_inputHeader.csv', 'st=2012-1-1', 'ed=2017-12-31', 'nitri=1.3', 'denitri=1', 'upt=1', 'immb=1','output=test.csv')

ModelArg = read.tcsv(text=arg, nfield =2, sep='=') 
ModelLibDenit = read.csv('https://raw.githubusercontent.com/laurencelin/SSHBS/master/lib_denitrification_ceof.csv')
ModelLibUptake = read.csv('https://raw.githubusercontent.com/laurencelin/SSHBS/master/lib_uptake_ceof.csv')

SimulationTime = seq.Date(from=as.Date(ModelArg$st), to=as.Date(ModelArg$ed) ,by="day")
SimulationTimeHour = do.call(c,lapply(SimulationTime,function(j){paste(j,0:23,sep='-')}))
timemax_BTS=24 * length(SimulationTime);
timemax_STS=3600

biochem_param$max_nitrif = max(5.341868e-10, min(biochem_param$max_nitrif*ModelArg$nitri, 2.444973e-05))
biochem_param$Mulholland_max_denitr = 10^ModelLibDenit[ModelArg$denitri,1]*0.01 # cm/s -> m/s
biochem_param$Mulholland_pow_denitr = ModelLibDenit[ModelArg$denitri,2] # per second

algal_param$Mulholland_max_denitr = 10^ModelLibUptake[ModelArg$upt,1]*0.01 # cm/s -> m/s
algal_param$Mulholland_pow_denitr = ModelLibUptake[ModelArg$upt,2] # per second
algal_param$v1ms_left = 0.15

detritus_param$Mulholland_max_denitr = 10^ModelLibUptake[ModelArg$immb,1]*0.01 # cm/s -> m/s
detritus_param$Mulholland_pow_denitr = ModelLibUptake[ModelArg$immb,2] # per second
detritus_param$v1ms_left = 0.15


header = read.csv(paste(ModelArg$proj, ModelArg$header ,sep='/'))
header$reachID = sapply(header$reachID,toString)
reachID = unique(header$reachID)
reachZone = list()
reachHyProf = list()
reachLbQ = list()
reachLsQ = list()
reachLNO3 = list()
reachLPO4 = list()
reachDetri = list()
reachPAR = list()
reachDeg = list()

for(ii in reachID){
    cond = header$reachID==ii
    reachZone[[ii]] = data.frame(
        zoneID = header$zoneID[cond],
        cqIndex = match(header$Cqinput[cond], unique(header$Cqinput[cond])),
        bqFrac = header$lateralBaseFrac[cond],
        sqFrac = header$lateralStormFrac[cond],
        len = header$len[cond],
        lenAdj = 0.1*header$len[cond],
        lenAdj_1 = 10/header$len[cond],
        hydrauIndex = match(header$hydraulicTable[cond], unique(header$hydraulicTable[cond])),
        parIndex = match(header$hourly_PAR[cond], unique(header$hourly_PAR[cond])),
        tempIndex = match(header$hourly_degree[cond], unique(header$hourly_degree[cond])),
        detriIndex = match(header$daily_detritus[cond], unique(header$daily_detritus[cond]))
    )# end of data.frame
    
    reachHyProf[[ii]] = lapply(unique(header$hydraulicTable[cond]),function(xx){
        ratingTable = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        ## rateTable assume zone length is 10m
        return <- new('CrossSectionQrelation',
            Q2Vol = approxfun(ratingTable[,'discharge'], ratingTable[,'volume']),
            Vol2Z = approxfun(ratingTable[,'volume'], ratingTable[,'depth']),
            Vol2V = approxfun(ratingTable[,'volume'], ratingTable[,'mvelocity']),
            Vol2WP = approxfun(ratingTable[,'volume'], ratingTable[,'wettedBperimeter']),
            Vol2Q = approxfun(ratingTable[,'volume'], ratingTable[,'discharge']),
            Vol2CA = approxfun(ratingTable[,'volume'], ratingTable[,'crossarea']),
            maxQ = max(ratingTable[,'discharge']),
            maxVol = max(ratingTable[ratingTable[,'bank']==0,'volume'])
            )#list
    })# end of lapply
    
    
    ## reach[[ii]] contains { ts1[t] , ts2[t] }; difficult to access
    ## reach[[ii]] contains { 1{ts1[1],ts2[1]}, 2{ts[2],ts[2]} }; easy to access
    ## reach[[ii]] contains { data.frame(row=profiles; col=time) }; can lump all reaches
    # -- assume daily input for cq
    reachLbQ[[ii]] = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
        baseflow = fivedayblockbaseflow(tmp[match(SimulationTime,tmp$date),'Qlps']*0.001)
        #stormQ = sapply(totalflow_m3s-baseflow,function(x){return <- ifelse(x>0,x,0)})
        return <- rep(baseflow,each=24)
    })))# end of lapply
    
    reachLsQ[[ii]] = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
        totalflow_m3s = tmp[match(SimulationTime,tmp$date),'Qlps']*0.001
        baseflow = fivedayblockbaseflow(totalflow_m3s)
        stormQ = sapply(totalflow_m3s-baseflow,function(x){return <- ifelse(x>0,x,0)})
        return <- rep(stormQ,each=24)
    })))# end of lapply
    
    reachLNO3[[ii]] = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
        return <- rep(tmp[match(SimulationTime,tmp$date),'Qlps']*tmp[match(SimulationTime,tmp$date),'Nmgpl'],each=24)
    })))# end of lapply
    
    reachLPO4[[ii]] = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
        return <- rep(tmp[match(SimulationTime,tmp$date),'Qlps']*tmp[match(SimulationTime,tmp$date),'Pmgpl'],each=24)
    })))# end of lapply
    
    ## this one is difficult (not transformed yet!)
    reachDetri[[ii]] = lapply(unique(header$daily_detritus[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
        return <- data.frame(
            labileC = rep(tmp$labileC[match(SimulationTime,tmp$date)],each=24),
            celluloseC = rep(tmp$celluloseC[match(SimulationTime,tmp$date)],each=24),
            ligninC = rep(tmp$ligninC[match(SimulationTime,tmp$date)],each=24),
            labileN = rep(tmp$labileN[match(SimulationTime,tmp$date)],each=24),
            celluloseN = rep(tmp$celluloseN[match(SimulationTime,tmp$date)],each=24),
            ligninN = rep(tmp$ligninN[match(SimulationTime,tmp$date)],each=24),
            labileP = rep(tmp$labileP[match(SimulationTime,tmp$date)],each=24),
            celluloseP = rep(tmp$celluloseP[match(SimulationTime,tmp$date)],each=24),
            ligninP = rep(tmp$ligninP[match(SimulationTime,tmp$date)],each=24)
            )# end of dataframe
    })# end of lapply
    
    reachPAR[[ii]] = as.data.frame(do.call(rbind,lapply(unique(header$hourly_PAR[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = paste(as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-')),tmp$hour,sep='-')
        return <- PAR2PARF(tmp$PAR[match(SimulationTimeHour,tmp$date)])
    })))# end of lapply
    
    reachDeg[[ii]] = as.data.frame(do.call(rbind,lapply(unique(header$hourly_degree[cond]),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = paste(as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-')),tmp$hour,sep='-')
        return <- tmp$temperature[match(SimulationTimeHour,tmp$date)]
    })))# end of lapply
    
}# end of for loop ii



## here
rm(ModelLibDenit)
rm(ModelLibUptake)

	
		
####################################### initial ###########################################
ReachHydraulic = list()
ReachCol = list()
ReachBen = list()
ReachStor = list()
ReachRipa = list()
ReachHold = list()
# lumping all reaches into one table seems a good idea


for(ii in reachID){
    numZone = dim(reachZone[[ii]])[1]
    numZoneIndex = seq_len(numZone)
    
    # need to reduce the size of this
    ReachHold[[ii]] = data.frame( tmpCond = rep(F,numZone) )
    ReachHold[[ii]]$downstreamCond = rep(T, numZone); ReachHold[[ii]]$downstreamCond[1]=F
    ReachHold[[ii]]$upstreamCond = rep(T, numZone); ReachHold[[ii]]$upstreamCond[numZone]=F
    ReachHold[[ii]]$hold = rep(0, numZone)
    ReachHold[[ii]]$holdII = rep(0, numZone)
    ReachHold[[ii]]$holdIII = rep(0, numZone)
    ReachHold[[ii]]$holdIV = rep(0, numZone)
    ReachHold[[ii]]$holdV = rep(0, numZone)
    ReachHold[[ii]]$sumPotentialNO3= rep(0, numZone)
    ReachHold[[ii]]$sumPotentialNH4= rep(0, numZone)
    ReachHold[[ii]]$sumPotentialPO4= rep(0, numZone)
    ReachHold[[ii]]$concNO3= rep(0, numZone)
    ReachHold[[ii]]$concNH4= rep(0, numZone)
    ReachHold[[ii]]$concPO4= rep(0, numZone)
    ReachHold[[ii]]$hypoconcNO3= rep(0, numZone)
    ReachHold[[ii]]$hypoconcNH4= rep(0, numZone)
    ReachHold[[ii]]$hypoconcPO4= rep(0, numZone)
    ReachHold[[ii]]$denitrify_concNO3= rep(0, numZone)
    ReachHold[[ii]]$algalQ10f= rep(0, numZone)
    ReachHold[[ii]]$Pot_algalupNO3= rep(0, numZone)
    ReachHold[[ii]]$Pot_algalupNH4= rep(0, numZone)
    ReachHold[[ii]]$Pot_algalupP= rep(0, numZone)
    ReachHold[[ii]]$bioQ10f= rep(0, numZone)
    ReachHold[[ii]]$algalhabitatDmgIndex= rep(0, numZone)
    ReachHold[[ii]]$algalhabitatDmg= rep(0, numZone)
    ReachHold[[ii]]$Pot_denitrification= rep(0, numZone)
    ReachHold[[ii]]$Pot_denitrificationFrac = rep(0, numZone)
    ReachHold[[ii]]$Pot_nitrification= rep(0, numZone)
    ReachHold[[ii]]$Pot_nitrificationFrac= rep(0, numZone)
    ReachHold[[ii]]$sumPotentialNO3= rep(0, numZone)
    ReachHold[[ii]]$sumPotentialNH4= rep(0, numZone)
    ReachHold[[ii]]$sumPotentialPO4= rep(0, numZone)
    ReachHold[[ii]]$Pot_immobNO3= rep(0, numZone)
    ReachHold[[ii]]$Pot_immobNH4= rep(0, numZone)
    ReachHold[[ii]]$Pot_immobPO4= rep(0, numZone)
    ReachHold[[ii]]$NO3limit_scalar= rep(0, numZone)
    ReachHold[[ii]]$NH4limit_scalar= rep(0, numZone)
    ReachHold[[ii]]$Plimit_scalar= rep(0, numZone)
    ReachHold[[ii]]$miner_class2decayRatio = rep(0, numZone)
    ReachHold[[ii]]$miner_class3decayRatio = rep(0, numZone)
    ReachHold[[ii]]$immo_class2decayRatio = rep(0, numZone)
    ReachHold[[ii]]$immo_class3decayRatio = rep(0, numZone)
    
    
    
    ReachHold[[ii]]$hold = cumsum(reachLbQ[[ii]][reachZone[[ii]]$cqIndex,1]*reachZone[[ii]]$bqFrac + reachLsQ[[ii]][reachZone[[ii]]$cqIndex,1]*reachZone[[ii]]$sqFrac) # ini discharge
    ReachHold[[ii]]$holdII = sapply(numZoneIndex,function(i){
        min(reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Q2Vol(2*0.001),
            reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@maxVol)*reachZone[[ii]]$lenAdj[i] # channel vol, not discharge
    }) # channel vol at mean discharge 2 L/s
    
    ReachCol[[ii]] = data.frame(q = sapply(numZoneIndex,function(i){
        min(reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Q2Vol(ReachHold[[ii]]$hold[i]),
            reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@maxVol)*reachZone[[ii]]$lenAdj[i] # channel vol, not discharge
    }))
    ReachCol[[ii]]$no3 = reachLNO3[[ii]][reachZone[[ii]]$cqIndex,1]/(reachLbQ[[ii]][reachZone[[ii]]$cqIndex,1]+reachLsQ[[ii]][reachZone[[ii]]$cqIndex,1])*ReachCol[[ii]]$q # mgN/m3 * m3
    ReachCol[[ii]]$nh4 = 0
    ReachCol[[ii]]$po4 = reachLPO4[[ii]][reachZone[[ii]]$cqIndex,1]/(reachLbQ[[ii]][reachZone[[ii]]$cqIndex,1]+reachLsQ[[ii]][reachZone[[ii]]$cqIndex,1])*ReachCol[[ii]]$q # mgP/m3 * m3
    ReachCol[[ii]]$POC = 0
    ReachCol[[ii]]$PON = 0
    ReachCol[[ii]]$POP = 0
    ReachCol[[ii]]$DOC = 0
    ReachCol[[ii]]$DON = 0
    ReachCol[[ii]]$DOP = 0
    
    ReachHydraulic[[ii]] = as.data.frame(do.call(rbind,lapply(numZoneIndex,function(i){ c(
        z = reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2Z(ReachCol[[ii]]$q[i]*reachZone[[ii]]$lenAdj_1[i]), #1
        v = reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2V(ReachCol[[ii]]$q[i]*reachZone[[ii]]$lenAdj_1[i]), #2
        wba = reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2WP(ReachCol[[ii]]$q[i]*reachZone[[ii]]$lenAdj_1[i])*reachZone[[ii]]$len[i], #3
        q = reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2Q(ReachCol[[ii]]$q[i]*reachZone[[ii]]$lenAdj_1[i]), #4
        A = reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2CA(ReachCol[[ii]]$q[i]*reachZone[[ii]]$lenAdj_1[i]), # 5
        maxv = reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@maxVol*reachZone[[ii]]$lenAdj[i] # 7
        ) })))
    ReachHydraulic[[ii]]$wba_1 = 1/ReachHydraulic[[ii]]$wba
    #cbind(ReachHydraulic[[ii]]$q, ReachCol[[ii]]$q, ReachHydraulic[[ii]]$q/ReachCol[[ii]]$q)
   
    ReachBen[[ii]] = data.frame(algC = rep(1000*10,max(numZoneIndex)) ) #gC/m2 -> mgC/m2
    ReachBen[[ii]]$algN = ReachBen[[ii]]$algC*0.5*(algal_param$maxNC+algal_param$minNC)
    ReachBen[[ii]]$algP = ReachBen[[ii]]$algC*0.5*(algal_param$maxPC+algal_param$minPC)
    ReachBen[[ii]]$algCN = ReachBen[[ii]]$algC/ReachBen[[ii]]$algN
    ReachBen[[ii]]$algCP = ReachBen[[ii]]$algC/ReachBen[[ii]]$algP
    ReachBen[[ii]]$algSelfF = 1
    ReachBen[[ii]]$algNutrF = 0
    ReachBen[[ii]]$algParF = 1
    ReachBen[[ii]]$algUpNF = 1
    ReachBen[[ii]]$algUpPF = 1
    ReachBen[[ii]]$algArea = 0.9 * sapply(numZoneIndex,function(i){reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2WP(ReachHold[[ii]]$holdII[i]*reachZone[[ii]]$lenAdj_1[i])*reachZone[[ii]]$len[i]})
    #cbind(ReachHydraulic[[ii]]$wba, ReachBen[[ii]]$algArea)
    
    ReachStor[[ii]] = data.frame(As = sapply(numZoneIndex,function(i){
        (3 + reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2CA(ReachHold[[ii]]$holdII[i]*reachZone[[ii]]$lenAdj_1[i])/reachHyProf[[ii]][[reachZone[[ii]]$hydrauIndex[i]]]@Vol2Z(ReachHold[[ii]]$holdII[i]*reachZone[[ii]]$lenAdj_1[i])) * (1-exp(-other_param$p_d*1.5))*other_param$p0/other_param$p_d
    }))
    ReachStor[[ii]]$q = ReachStor[[ii]]$As * reachZone[[ii]]$len
    ReachStor[[ii]]$no3 = reachLNO3[[ii]][reachZone[[ii]]$cqIndex,1]/(reachLbQ[[ii]][reachZone[[ii]]$cqIndex,1]+reachLsQ[[ii]][reachZone[[ii]]$cqIndex,1])*ReachStor[[ii]]$q
    ReachStor[[ii]]$nh4 = 0
    ReachStor[[ii]]$po4 = reachLPO4[[ii]][reachZone[[ii]]$cqIndex,1]/(reachLbQ[[ii]][reachZone[[ii]]$cqIndex,1]+reachLsQ[[ii]][reachZone[[ii]]$cqIndex,1])*ReachStor[[ii]]$q
    ReachStor[[ii]]$DOC = 0
    ReachStor[[ii]]$DON = 0
    ReachStor[[ii]]$DOP = 0
    ReachStor[[ii]]$baseNO3 = ReachStor[[ii]]$no3
    ReachStor[[ii]]$baseNH4 = ReachStor[[ii]]$nh4
    ReachStor[[ii]]$baseNO3balance = 0
    ReachStor[[ii]]$baseNH4balance = 0
    #cbind(ReachStor[[ii]]$q, ReachHydraulic[[ii]]$maxv, ReachStor[[ii]]$q/ReachHydraulic[[ii]]$maxv, ReachCol[[ii]]$q)
    
    # remove these later
    STSexport = matrix(0,numZone,length(columnResource_variable)); colnames(STSexport)= columnResource_variable
    STSaggregateVar = matrix(0,numZone,length(columnResource_variable)); colnames(STSaggregateVar)= columnResource_variable
    STSinput = matrix(0,numZone,length(columnResource_variable)); colnames(STSinput)= columnResource_variable
    STSinputTrack = matrix(0,numZone,length(columnResource_variable)); colnames(STSinputTrack)= columnResource_variable
    STSLateralinputTrack = matrix(0,numZone,length(columnResource_variable)); colnames(STSLateralinputTrack)= columnResource_variable
    exportProp = rep(0,numZone)
    hyporheicDenitrification = rep(0,numZone)
    
    
    
    
}# end of for loop ii


		
	

		
		
	
	
	RiparianResource_matrix = matrix(0,numZone,length(columnResource_variable)); colnames(RiparianResource_matrix)= columnResource_variable
	cellspace = rep(0,numZone)
	exRiparian = rep(0,numZone)
	concChange = rep(0,numZone)
	concMassChange = rep(0,numZone)
	dailyAverageWetBenthicArea = rep(0, numZone);
	
	
	## -------- holders / containers
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
	
	STSstepOption = c(1, 2, 4, 8, 24, 48)	
	lightAdjust = rep(1, numZone); lightAdjust[1:8]=0.5 # no use anymore July 1, 2019
	lightAdjustII = rep(1, numZone); lightAdjustII[1:8]=1.2 # no use anymore July 1, 2019
	
	detritusMatrix = matrix(0, numZone, length(detritusName)); colnames(detritusMatrix)=detritusName;
	detritusMatrix[,'microbes'] = (lateralResource_matrix[,1,'detritus1c']+lateralResource_matrix[,1,'detritus2c']+lateralResource_matrix[,1,'detritus3c'])*0.09
	detritusMatrix[,'miner'] = (lateralResource_matrix[,1,'detritus1c']+lateralResource_matrix[,1,'detritus2c']+lateralResource_matrix[,1,'detritus3c'])*0.09 * 0.1
	detritusMatrix[,'detritus1c'] = lateralResource_matrix[,1,'detritus1c']
	detritusMatrix[,'detritus2c'] = lateralResource_matrix[,1,'detritus2c']
	detritusMatrix[,'detritus3c'] = lateralResource_matrix[,1,'detritus3c']
    detritusMatrix[,'detritus1n'] = lateralResource_matrix[,1,'detritus1n']#/detritusCN
    detritusMatrix[,'detritus2n'] = lateralResource_matrix[,1,'detritus2n']#/detritusCN
    detritusMatrix[,'detritus3n'] = lateralResource_matrix[,1,'detritus3n']#/detritusCN
    detritusMatrix[,'detritus1p'] = lateralResource_matrix[,1,'detritus1c']/detritusCP
    detritusMatrix[,'detritus2p'] = lateralResource_matrix[,1,'detritus2c']/detritusCP
    detritusMatrix[,'detritus3p'] = lateralResource_matrix[,1,'detritus3c']/detritusCP

	## ----- output aggregaters
	BTSaggregateVar = matrix(0,nrow=numZone, ncol=length(propfluxTitle));
colnames(BTSaggregateVar)=propfluxTitle
	STSaggregateVar = matrix(0,nrow=numZone, ncol=length(outputTitleDownstreamExport));
		colnames(STSaggregateVar)= outputTitleDownstreamExport

    crossSection = rep(0, numZone)

	vectorMin = function(aa,bb){
		sapply(seq_along(aa),function(ii){min(aa[ii],bb[ii])})
	}#function
	vectorMax = function(aa,bb){
		sapply(seq_along(aa),function(ii){max(aa[ii],bb[ii])})
	}#function
	
	outputFile = sindex$outputFile 
	outputFile_buff <- file(outputFile,'w')
	cat( c(outputTitle, paste('no3out',1:numZone,sep='_')), '\n', file=outputFile_buff,sep=',')
	
	# set these here!
    downZs = 11:34 # max(zoneInfo[,'zoneID'])
	upZs = 1:8
	BTSaggregateVarPrint = c('nitrification','denitrification','algalNO3uptake','algalNH4uptake','algalNPP','mineralN','algalCN','algalCP','algalC','algalCloss','netmineral','detric','detriCN','detriCP','microbC','minerC','algalR','detritusR','Sw','tmperature','PARf','hdmg','algalNlimit','algalPlimit')
	
system.time({
	for(BTS in 1:timemax_BTS ){	 #1:timemax_BTS, 8760
		
		## debugging
		# hyporheicResource_matrix
		# columnResource_matrix
		# bioprocess_matrix
		# detritusMatrix
		
		
		if(sum(detritusMatrix[,'detritus1c']<0)>0){print(BTS);break}
		
		exportProp = hydraulic_matrix[,'discharge']/columnResource_matrix[,'discharge'];
		STSstep = max(1,apply(exportProp %o% STSstepOption<0.25,2,prod)* STSstepOption)	
		
		
        #algalQ10f = (algal_parameter['Q10']^(0.1*zoneEnv_matrix[,BTS,'temperature']-1.2) ) ## Jack's equation for algae # before verson 8
        algalQ10f = (algal_parameter['Q10']^(0.1*zoneEnv_matrix[,BTS,'temperature']-0.2) ) ## Jack's equation for algae # version 8 testing
		bioQ10f = (biochem_parameter['Q10']^(0.1*zoneEnv_matrix[,BTS,'temperature']-0.2) ) ## Jack's equation for microbes
		
         ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------##
		## algal death (biological with some velocity related)
        # xx = seq(0,3,0.01)
        # plot(xx, 2/(1+exp(-(xx-0.5)/0.02)), type='b'); # daily
        # hold[] = 1/(1+exp(-(hydraulic_matrix[,'velocity']-0.1)/0.02)) * algalQ10f * algal_parameter['moralityRate']* timemax_STS # moralityRate = death rate inthe modle now
        hold[] = algalQ10f * algal_parameter['moralityRate'] * timemax_STS # moralityRate = death rate inthe modle now
		bioprocess_matrix[,'algaldeadC'] <- bioprocess_matrix[,'algalc'] * hold[]  
		bioprocess_matrix[,'algaldeadN'] <- bioprocess_matrix[,'algaln'] * hold[] 
		bioprocess_matrix[,'algaldeadP'] <- bioprocess_matrix[,'algalp'] * hold[]
       
        ## wash out; could this be buried instead?
        # see below for document
        hold[] =  1 - (1-(1-detritus_param$v1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(1/24) # removal rate
        hold[hydraulic_matrix[,'velocity']<0.1]=0
		bioprocess_matrix[,'algalentC'] <- bioprocess_matrix[,'algalc'] * hold[]  
		bioprocess_matrix[,'algalentN'] <- bioprocess_matrix[,'algaln'] * hold[] 
		bioprocess_matrix[,'algalentP'] <- bioprocess_matrix[,'algalp'] * hold[] 
		#xx=seq(0,1.0,0.01); plot(xx, 1/(1+exp(-(xx-0.4)/0.07)), type='b')	
			
			if(BTS%%24==0){
				
				# new dmg
                tmpCond = algalhabitatDmgIndex < hold ## from above
				algalhabitatDmg[tmpCond] = 1-(1-hold[tmpCond])^24 # dmg inital 0
				
				# new dmg apply to old dmg
				tmpCond = tmpCond & (algalhabitatDmg > algalhabitatDmgIndex)
				algalhabitatDmgIndex[tmpCond] = 1 #algalhabitatDmg[tmpCond]
				
				# recover
				tmpCond = algalhabitatDmgIndex>0
                algalhabitatDmgIndex[tmpCond] = algalhabitatDmgIndex[tmpCond] * algalhabitatDmg[tmpCond] # algalhabitatDmgIndex is changing value; algalhabitatDmg is the initial hit
				
				# recover reset; max is 0.09148242
				tmpCond = algalhabitatDmgIndex<1e-3
				algalhabitatDmgIndex[tmpCond]=0
				algalhabitatDmg[tmpCond]=0
			}# daily
		
		## mineralization (it's greater than uptake!!) reset
		hold[]  = algal_parameter['mineralRate'] * algalQ10f * timemax_STS
		bioprocess_matrix[,'algalMineralC'] <- bioprocess_matrix[,'algalc'] * hold[]  #<<-- this respiration if we treat NPP as GPP
		bioprocess_matrix[,'algalMineralN'] <- bioprocess_matrix[,'algaln'] * hold[] 
		bioprocess_matrix[,'algalMineralP'] <- bioprocess_matrix[,'algalp'] * hold[] 
		
		bioprocess_matrix[,'algalc'] = bioprocess_matrix[,'algalc'] - bioprocess_matrix[,'algaldeadC'] - bioprocess_matrix[,'algalentC'] - bioprocess_matrix[,'algalMineralC']
		bioprocess_matrix[,'algaln'] = bioprocess_matrix[,'algaln'] - bioprocess_matrix[,'algaldeadN'] - bioprocess_matrix[,'algalentN'] - bioprocess_matrix[,'algalMineralN'] 
		bioprocess_matrix[,'algalp'] = bioprocess_matrix[,'algalp'] - bioprocess_matrix[,'algaldeadP'] - bioprocess_matrix[,'algalentP'] - bioprocess_matrix[,'algalMineralP']
		
        ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------##
		## microbes death, potential decay and immob; bioprocess_matrix[,'colonizedBenthicArea'], hydraulic_matrix[,'wetBenthicArea']
		## logistic growth break down rP*(1-P/K) = rP - r(P/K)[P] --> here  remaining % = 1-r(P/K)
        holdII[] = detritus_parameter['growthRate']*10*(detritusMatrix[,'microbes']+ detritusMatrix[,'miner'])/(detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']+detritusMatrix[,'detritus3c'])*bioQ10f*timemax_STS
        holdII[holdII>1]=0.999
        hold[] = holdII * detritusMatrix[,'microbes'] 
        holdIV[] = holdII * detritusMatrix[,'microbes']*detritus_parameter['nc']
        holdV[] = holdII * detritusMatrix[,'microbes']*detritus_parameter['pc']
        detritusMatrix[,'microbes'] = detritusMatrix[,'microbes'] * (1-holdII)
        
        holdII[] = detritus_parameter['growthRate_miner']*10*(detritusMatrix[,'microbes']+ detritusMatrix[,'miner'])/(detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']+detritusMatrix[,'detritus3c'])*bioQ10f*timemax_STS
        holdII[holdII>1]=0.999
        hold[] = hold + holdII * detritusMatrix[,'miner']
      	holdIV[] = holdIV[] + holdII[] * detritusMatrix[,'miner']*detritus_parameter['miner_nc']
        holdV[] = holdV[] + holdII[] * detritusMatrix[,'miner']*detritus_parameter['miner_pc']
        detritusMatrix[,'miner'] = detritusMatrix[,'miner'] * (1-holdII)
        
        
        holdIII[] = detritusMatrix[,'miner']*detritus_parameter['baseR']*bioQ10f*timemax_STS
        detritusMatrix[,'resp'] = detritusMatrix[,'microbes'] * detritus_parameter['baseR']*bioQ10f*timemax_STS
        detritusMatrix[,'respN'] = detritusMatrix[,'resp']*detritus_parameter['nc'] + holdIII*detritus_parameter['miner_nc']
        detritusMatrix[,'respP'] = detritusMatrix[,'resp']*detritus_parameter['pc'] + holdIII*detritus_parameter['miner_pc']
        detritusMatrix[,'resp'] = detritusMatrix[,'resp'] + holdIII
        detritusMatrix[,'microbes'] = detritusMatrix[,'microbes'] * (1-detritus_parameter['baseR']*bioQ10f*timemax_STS) # basal respiration
        detritusMatrix[,'miner'] = detritusMatrix[,'miner'] * (1-detritus_parameter['baseR']*bioQ10f*timemax_STS) # basal respiration
       
        
        ## -- detritus inputs
            ## lateralResource_matrix[, BTS,'discharge']*3600 # hourly input water m3; 5 L/s  -> 18 m3/hour
			## xx=seq(0,15,0.1); ww=xx/1000*3600; yy=dnorm(ww,mean=18,sd=3)/dnorm(18,mean=18,sd=3); plot(xx,yy,type='b'); abline(v=3); abline(v=8)
			holdII = detritusCinput*dnorm(lateralResource_matrix[, BTS,'discharge']*3600,mean=18,sd=3)/dnorm(18,mean=18,sd=3)
            if(lateralResource_matrix[1,BTS,'detritus1c']+lateralResource_matrix[1,BTS,'detritus2c']+lateralResource_matrix[1,BTS,'detritus3c'] > 0){
	            class1Frac = lateralResource_matrix[1,BTS,'detritus1c']/(lateralResource_matrix[1,BTS,'detritus1c']+lateralResource_matrix[1,BTS,'detritus2c']+lateralResource_matrix[1,BTS,'detritus3c'])
	            class2Frac = lateralResource_matrix[1,BTS,'detritus2c']/(lateralResource_matrix[1,BTS,'detritus1c']+lateralResource_matrix[1,BTS,'detritus2c']+lateralResource_matrix[1,BTS,'detritus3c'])
	            class3Frac = 1 - class1Frac - class2Frac
            }# couldb e input zero
		detritusMatrix[,'detritus1c'] = detritusMatrix[,'detritus1c'] + class1Frac*holdII[] +
            lateralResource_matrix[,BTS,'detritus1c']*0.04166667 +
            hold*0.5 + bioprocess_matrix[,'algaldeadC']*0.5 # mgC/m2 unit area;
            # hold is dead microbes at this point
            # class1Frac*holdII[] is lateral detritus input, which is zero
		detritusMatrix[,'detritus1n'] = detritusMatrix[,'detritus1n'] + class1Frac*holdII[]/detritusCN +
            lateralResource_matrix[,BTS,'detritus1n']*0.04166667 +
			holdIV*0.5 + bioprocess_matrix[,'algaldeadN']*0.5 # mgN/m2 unit area
		detritusMatrix[,'detritus1p'] = detritusMatrix[,'detritus1p'] + class1Frac*holdII[]/detritusCP +
            lateralResource_matrix[,BTS,'detritus1c']*0.04166667/detritusCP +
			holdV*0.5 + bioprocess_matrix[,'algaldeadP']*0.5 # mgN/m2 unit area	
            
		detritusMatrix[,'detritus2c'] = detritusMatrix[,'detritus2c'] + class2Frac*holdII[] +
            lateralResource_matrix[,BTS,'detritus2c']*0.04166667 +
			hold*0.5 + bioprocess_matrix[,'algaldeadC']*0.5
		detritusMatrix[,'detritus2n'] = detritusMatrix[,'detritus2n'] + class2Frac*holdII[]/detritusCN +
            lateralResource_matrix[,BTS,'detritus2n']*0.04166667 +
			holdIV*0.5 + bioprocess_matrix[,'algaldeadN']*0.5
		detritusMatrix[,'detritus2p'] = detritusMatrix[,'detritus2p'] + class2Frac*holdII[]/detritusCP +
            lateralResource_matrix[,BTS,'detritus2c']*0.04166667/detritusCP +
			holdV*0.5 + bioprocess_matrix[,'algaldeadP']*0.5

		detritusMatrix[,'detritus3c'] = detritusMatrix[,'detritus3c'] + class3Frac*holdII[] +
            lateralResource_matrix[,BTS,'detritus3c']*0.04166667
		detritusMatrix[,'detritus3n'] = detritusMatrix[,'detritus3n'] + class3Frac*holdII[]/detritusCN +
            lateralResource_matrix[,BTS,'detritus3n']*0.04166667
		detritusMatrix[,'detritus3p'] = detritusMatrix[,'detritus3p'] + class3Frac*holdII[]/detritusCP +
            lateralResource_matrix[,BTS,'detritus3c']*0.04166667/detritusCP
        
		## micorbial potential growth
		hold[] = (detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']*miner_class2decayRatio +detritusMatrix[,'detritus3c']*miner_class3decayRatio)
		miner_class2decayRatio[] = (10000+detritusMatrix[,'detritus2c'])/(100000+detritusMatrix[,'detritus1c']); # relative rate for class-2 materials
		miner_class3decayRatio[] = (1000+detritusMatrix[,'detritus3c'])/(100000+detritusMatrix[,'detritus1c']) # relative rate for class-3 materials
			detritusMatrix[,'foodNC_miner'] = (detritusMatrix[,'detritus1n']+detritusMatrix[,'detritus2n']*miner_class2decayRatio +detritusMatrix[,'detritus3n']*miner_class3decayRatio) / hold								    
        	detritusMatrix[,'foodPC_miner'] = (detritusMatrix[,'detritus1p']+detritusMatrix[,'detritus2p']*miner_class2decayRatio +detritusMatrix[,'detritus3p']*miner_class3decayRatio) / hold
			
			## growth_miner & decay_miner are correct but how to back calculate 
			detritusMatrix[,'growth_miner'] = detritusMatrix[,'miner']*detritus_parameter['growthRate_miner'] * bioQ10f * timemax_STS; # potential miner growth
			detritusMatrix[,'decay_miner'] = vectorMin(
				hold, 
				vectorMax(
					detritusMatrix[,'growth_miner'],
					vectorMax(
						detritus_parameter['gcoef']*detritusMatrix[,'growth_miner']*detritus_parameter['miner_nc']/detritusMatrix[,'foodNC_miner'],
						detritus_parameter['gcoef']*detritusMatrix[,'growth_miner']*detritus_parameter['miner_pc']/detritusMatrix[,'foodPC_miner']	
					)));
				# cbind(hold,detritusMatrix[,'growth_miner'], detritusMatrix[,'decay_miner'])
			
				tmpCond[] = detritusMatrix[,'growth_miner'] > detritusMatrix[,'decay_miner']
				if(sum(tmpCond)>0){
					detritusMatrix[tmpCond,'growth_miner'] = vectorMin(
						detritusMatrix[tmpCond,'decay_miner']*detritusMatrix[tmpCond,'foodNC_miner']/detritus_parameter['miner_nc'],
						detritusMatrix[tmpCond,'decay_miner']*detritusMatrix[tmpCond,'foodPC_miner']/detritus_parameter['miner_pc']
						)
				}#if
				detritusMatrix[,'growth_miner'] = detritusMatrix[,'growth_miner']*detritus_parameter['gcoef']
			detritusMatrix[,'decay_minerN'] = detritusMatrix[,'decay_miner']*detritusMatrix[,'foodNC_miner'] - detritusMatrix[,'growth_miner']*detritus_parameter['miner_nc']
			detritusMatrix[,'decay_minerP'] = detritusMatrix[,'decay_miner']*detritusMatrix[,'foodPC_miner'] - detritusMatrix[,'growth_miner']*detritus_parameter['miner_pc']
			detritusMatrix[detritusMatrix[,'decay_minerN']<0,'decay_minerN']=0
			detritusMatrix[detritusMatrix[,'decay_minerP']<0,'decay_minerP']=0
		
		hold[] = (detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']*immo_class2decayRatio +detritusMatrix[,'detritus3c']*immo_class3decayRatio)
		immo_class2decayRatio[] = (10000+detritusMatrix[,'detritus2c'])/(1000000+detritusMatrix[,'detritus1c']); # relative rate for class-2 materials
		immo_class3decayRatio[] = (1000+detritusMatrix[,'detritus3c'])/(1000000+detritusMatrix[,'detritus1c']) # relative rate for class-3 materials
			detritusMatrix[,'foodNC'] = (detritusMatrix[,'detritus1n']+detritusMatrix[,'detritus2n']*immo_class2decayRatio +detritusMatrix[,'detritus3n']*immo_class3decayRatio) / hold
	        detritusMatrix[,'foodPC'] = (detritusMatrix[,'detritus1p']+detritusMatrix[,'detritus2p']*immo_class2decayRatio +detritusMatrix[,'detritus3p']*immo_class3decayRatio) / hold
			detritusMatrix[,'Ncoef'] = detritus_parameter['gcoef']/detritus_parameter['cn'] - detritusMatrix[,'foodNC']
	        detritusMatrix[,'Pcoef'] = detritus_parameter['gcoef']/detritus_parameter['cp'] - detritusMatrix[,'foodPC']
			detritusMatrix[,'decay'] = detritusMatrix[,'microbes']*detritus_parameter['growthRate'] * bioQ10f * timemax_STS;
		
			
        detritusMatrix[,'netmineral'] = detritusMatrix[,'decay'] * detritusMatrix[,'Ncoef'] #- detritusMatrix[,'decay_minerN'] # hourly 
		detritusMatrix[,'netmineralP'] = detritusMatrix[,'decay'] * detritusMatrix[,'Pcoef'] #- detritusMatrix[,'decay_minerP']
            
            
        
			# debug = data.frame(zone=1:numZone)
			# debug$C = (detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']+detritusMatrix[,'detritus3c'])
			# debug$CN = debug$C/(detritusMatrix[,'detritus1n']+detritusMatrix[,'detritus2n']+detritusMatrix[,'detritus3n'])
			# debug$mic = detritusMatrix[,'microbes']/debug$C*100
			# debug$netmineral = detritusMatrix[,'netmineral']
			# debug$decay = detritusMatrix[,'decay']
			
			## debugging
			# hyporheicResource_matrix
			# columnResource_matrix
			# bioprocess_matrix
			# detritusMatrix
        ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------##
		STS = 0
		while(STS < timemax_STS){ #1:timemax_STS
			
				tmpCond[] = zoneEnv_matrix[,BTS,'temperature']>0 ## freezing issue
				## export
					exportProp = hydraulic_matrix[,'discharge']/columnResource_matrix[,'discharge']; exportProp[!tmpCond]=0
					STSexport = columnResource_matrix * exportProp * STSstep
					STSaggregateVar = STSaggregateVar + STSexport
					columnResource_matrix[] = columnResource_matrix[] * (1-exportProp * STSstep)
					
									
				## input (per second)
				STSinput[tmpCond, interested_variable] = #STSinput[tmpCond, interested_variable] + 
					lateralResource_matrix[tmpCond, BTS,interested_variable] * STSstep + 
					lateralResource_matrix[tmpCond, BTS,c('stormQ','stormNO3','stormNH4','stormPO4')] * STSstep #lateral/spring
					# STSinput[,'waterNO3']/STSinput[,'discharge'] #mg/m3 = ug/L
					STSinput[tmpCond,'waterNO3'] = STSinput[tmpCond,'waterNO3']
				STSLateralinputTrack[] = STSLateralinputTrack[] + STSinput[]
				STSinput[downstreamCond&tmpCond,] = STSinput[downstreamCond&tmpCond,] + STSexport[upstreamIndex[downstreamCond&tmpCond],] #upstream
                STSinputTrack[] = STSinputTrack[] + STSinput[] # tracking lcoal lateral and upstream inputs
					
					## mineralization (hourly to seconds)
				STSinput[,'waterNH4'] = STSinput[,'waterNH4'] + (bioprocess_matrix[,'algalMineralN']+detritusMatrix[,'mineral']+detritusMatrix[,'respN'])*hydraulic_matrix[,'wetBenthicArea'] * 0.0002777778 * STSstep
				STSinput[,'waterPO4'] = STSinput[,'waterPO4'] + (bioprocess_matrix[,'algalMineralP']+detritusMatrix[,'decay_minerP'])*hydraulic_matrix[,'wetBenthicArea']* 0.0002777778 * STSstep
				
					## entraintment & deposition
                # sigmoid function - logistic form:  1 / (1+exp(-x))
                # detritus_param$v1ms_left is the x% of material remain within a day under velocity V -> removal rate = f(V)
                # f(V) = (1-x)/(1+exp(-(V-0.4)/0.06))
                # e.g., x = 15% material remains within a day,
                # V = seq(0,5,0.001);  daily_removal_rate = (1-x)/(1+exp(-(V-0.4)/0.06))
                #
                # then the hourly removal rate is ((1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24)
                # e.g., at v = 1.0 m/s, hourly removal rate = ((1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24) = 0.9932494,
                # which does not make sense to remove 0.9932494 at the first hour;
                #
                # instead, we need to model, remaining rate, z^24 = 1-(1-x)/(1+exp(-(V-0.4)/0.06))
                # z = (1-(1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24) => 0.9240067
                # hourly removal rate = 1-z = 1 - (1-(1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24) = 0.0759933 @ V = 1.0 m/s
                # hourly removal rate = 1-z = 1 - (1-(1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24) = 5.319941e-05 @ V = 0.01 m/s
                #
                # for remaining suspending in water column, at v = 1.0 m/s, hourly removal rate = ((1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24) = 0.9932494
                # at v = 0.01 m/s, hourly removal rate = ((1-x)/(1+exp(-(V-0.4)/0.06)))^(1/24) = 0.7575486
                #
                holdIII[] = ((1-detritus_param$v1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(STSstep/86400)
                hold[] =  1 - (1-(1-detritus_param$v1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(STSstep/86400) # removal rate
                hold[] = hold[] * 10*(detritusMatrix[,'microbes']+ detritusMatrix[,'miner'])/(detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']+detritusMatrix[,'detritus3c']) # count for micorbial processing to fragmentation (new June 17)
                holdII[] = STSinput[,'FOC']*(1-holdIII) * hydraulic_matrix[,'wetBenthicArea_1'] / (detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']+detritusMatrix[,'detritus3c']) #
				STSinput[,'FOC'] =
                    STSinput[,'FOC']* holdIII + ## still questioning about this part
                    (detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']+detritusMatrix[,'detritus3c']) * hold * hydraulic_matrix[,'wetBenthicArea'] +
                    bioprocess_matrix[,'algalentC']*0.0002777778 * STSstep *  hydraulic_matrix[,'wetBenthicArea'] ## pre-calculated outside STS loop
                detritusMatrix[,'detritus1c'] = detritusMatrix[,'detritus1c'] * (1-hold) * (1+holdII)
                detritusMatrix[,'detritus2c'] = detritusMatrix[,'detritus2c'] * (1-hold) * (1+holdII)
                detritusMatrix[,'detritus3c'] = detritusMatrix[,'detritus3c'] * (1-hold) * (1+holdII)
				
				holdII[] = STSinput[,'FON']*(1-holdIII) * hydraulic_matrix[,'wetBenthicArea_1'] / (detritusMatrix[,'detritus1n']+detritusMatrix[,'detritus2n']+detritusMatrix[,'detritus3n'])
				STSinput[,'FON'] =
                    STSinput[,'FON'] * holdIII +
                    (detritusMatrix[,'detritus1n']+detritusMatrix[,'detritus2n']+detritusMatrix[,'detritus3n'])*hold * hydraulic_matrix[,'wetBenthicArea'] +
					bioprocess_matrix[,'algalentN']*0.0002777778 * STSstep *  hydraulic_matrix[,'wetBenthicArea']
                detritusMatrix[,'detritus1n'] = detritusMatrix[,'detritus1n'] * (1-hold) * (1+holdII)
                detritusMatrix[,'detritus2n'] = detritusMatrix[,'detritus2n'] * (1-hold) * (1+holdII)
                detritusMatrix[,'detritus3n'] = detritusMatrix[,'detritus3n'] * (1-hold) * (1+holdII)
				
				holdII[] = STSinput[,'FOP']*(1-holdIII) * hydraulic_matrix[,'wetBenthicArea_1'] / (detritusMatrix[,'detritus1p']+detritusMatrix[,'detritus2p']+detritusMatrix[,'detritus3p'])
				STSinput[,'FOP'] =
					 STSinput[,'FOP'] * holdIII +
					 (detritusMatrix[,'detritus1p']+detritusMatrix[,'detritus2p']+detritusMatrix[,'detritus3p'])*hold * hydraulic_matrix[,'wetBenthicArea'] +
					 bioprocess_matrix[,'algalentP']*0.0002777778 * STSstep *  hydraulic_matrix[,'wetBenthicArea']
				detritusMatrix[,'detritus1p'] = detritusMatrix[,'detritus1p'] * (1-hold) * (1+holdII)
                detritusMatrix[,'detritus2p'] = detritusMatrix[,'detritus2p'] * (1-hold) * (1+holdII)
                detritusMatrix[,'detritus3p'] = detritusMatrix[,'detritus3p'] * (1-hold) * (1+holdII)
			
					
				columnResource_matrix[] = columnResource_matrix[] + STSinput
					# cbind(zoneInfo, STSexport)
					# cbind(zoneInfo, lateralResource_matrix[, BTS,interested_variable])
					# cbind(hydraulic_matrix, columnResource_matrix[,'discharge'])
					# cbind(zoneInfo, exportProp)
					
					
				## extract to riparian	
				cellspace[ ] = columnResource_matrix[,'discharge'] -  hydraulic_matrix[,'maxVol'] # positive = to riparian; negative = potential from riparian
				tmpCond[ ] = cellspace>0 # is bounded by numZone
				if(sum(tmpCond)>0){
					## go to riparian
					exRiparian[ ] = 0 # is bounded by numZone
					exRiparian[ ] = cellspace[ ]/columnResource_matrix[,'discharge'] 
					RiparianResource_matrix[tmpCond,'discharge'] = RiparianResource_matrix[tmpCond,'discharge'] + cellspace[tmpCond]
					RiparianResource_matrix[tmpCond,'waterNO3'] = RiparianResource_matrix[tmpCond,'waterNO3'] + exRiparian[tmpCond]* columnResource_matrix[tmpCond,'waterNO3']
					RiparianResource_matrix[tmpCond,'waterNH4'] = RiparianResource_matrix[tmpCond,'waterNH4'] + exRiparian[tmpCond]* columnResource_matrix[tmpCond,'waterNH4']
					RiparianResource_matrix[tmpCond,'waterPO4'] = RiparianResource_matrix[tmpCond,'waterPO4'] + exRiparian[tmpCond]* columnResource_matrix[tmpCond,'waterPO4']
					RiparianResource_matrix[tmpCond,'FOC'] = RiparianResource_matrix[tmpCond,'FOC'] + exRiparian[tmpCond]* columnResource_matrix[tmpCond,'FOC']
					RiparianResource_matrix[tmpCond,'FON'] = RiparianResource_matrix[tmpCond,'FON'] + exRiparian[tmpCond]* columnResource_matrix[tmpCond,'FON']
					RiparianResource_matrix[tmpCond,'FOP'] = RiparianResource_matrix[tmpCond,'FOP'] + exRiparian[tmpCond]* columnResource_matrix[tmpCond,'FOP']
					
					exRiparian = 1-exRiparian # is bounded by numZone				
					columnResource_matrix[tmpCond,'discharge'] = columnResource_matrix[tmpCond,'discharge'] * exRiparian[tmpCond]
					columnResource_matrix[tmpCond,'waterNO3'] = columnResource_matrix[tmpCond,'waterNO3'] * exRiparian[tmpCond]
					columnResource_matrix[tmpCond,'waterNH4'] = columnResource_matrix[tmpCond,'waterNH4'] * exRiparian[tmpCond]
					columnResource_matrix[tmpCond,'waterPO4'] = columnResource_matrix[tmpCond,'waterPO4'] * exRiparian[tmpCond]
					columnResource_matrix[tmpCond,'FOC'] = columnResource_matrix[tmpCond,'FOC'] * exRiparian[tmpCond]
					columnResource_matrix[tmpCond,'FON'] = columnResource_matrix[tmpCond,'FON'] * exRiparian[tmpCond]
					columnResource_matrix[tmpCond,'FOP'] = columnResource_matrix[tmpCond,'FOP'] * exRiparian[tmpCond]
					
				}# to riparian
					
				## input from riparian (per STS) 
				tmpCond[ ] = cellspace<0 & RiparianResource_matrix[,'discharge']>0 # is bounded by numZone
				if(sum(tmpCond)>0){
					# get back from riparian
					exRiparian[ ] = -cellspace/RiparianResource_matrix[,'discharge'] # is bounded by numZone
					exRiparian[is.na(exRiparian)|is.infinite(exRiparian)] = 0;
					exRiparian[exRiparian>1] = 1 
					exRiparian[exRiparian<0] = 0	
					
					columnResource_matrix[tmpCond,'discharge'] = columnResource_matrix[tmpCond,'discharge'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'discharge']
					columnResource_matrix[tmpCond,'waterNO3'] = columnResource_matrix[tmpCond,'waterNO3'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'waterNO3']
					columnResource_matrix[tmpCond,'waterNH4'] = columnResource_matrix[tmpCond,'waterNH4'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'waterNH4']
					columnResource_matrix[tmpCond,'waterPO4'] = columnResource_matrix[tmpCond,'waterPO4'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'waterPO4']
					columnResource_matrix[tmpCond,'FOC'] = columnResource_matrix[tmpCond,'FOC'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'FOC']
					columnResource_matrix[tmpCond,'FON'] = columnResource_matrix[tmpCond,'FON'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'FON']
					columnResource_matrix[tmpCond,'FOP'] = columnResource_matrix[tmpCond,'FOP'] + exRiparian[tmpCond]*RiparianResource_matrix[tmpCond,'FOP']
					exRiparian = 1-exRiparian
					RiparianResource_matrix[tmpCond,'discharge'] = RiparianResource_matrix[tmpCond,'discharge']*exRiparian[tmpCond]
					RiparianResource_matrix[tmpCond,'waterNO3'] = RiparianResource_matrix[tmpCond,'waterNO3']*exRiparian[tmpCond]
					RiparianResource_matrix[tmpCond,'waterNH4'] = RiparianResource_matrix[tmpCond,'waterNH4']*exRiparian[tmpCond]
					RiparianResource_matrix[tmpCond,'waterPO4'] = RiparianResource_matrix[tmpCond,'waterPO4']*exRiparian[tmpCond]
					RiparianResource_matrix[tmpCond,'FOC'] = RiparianResource_matrix[tmpCond,'FOC']*exRiparian[tmpCond]
					RiparianResource_matrix[tmpCond,'FON'] = RiparianResource_matrix[tmpCond,'FON']*exRiparian[tmpCond]
					RiparianResource_matrix[tmpCond,'FOP'] = RiparianResource_matrix[tmpCond,'FOP']*exRiparian[tmpCond]
				}# to riparian
	
	
				###################################################################################### potential uptake
                # FOC, FON, FOP are not dissolved
                
                ##-------------------------------------------
	                # columnResource_matrix
	                tmpCond[] = columnResource_matrix[,'discharge']<=0
	                concNO3[] = columnResource_matrix[,'waterNO3']/columnResource_matrix[,'discharge']; concNO3[is.na(concNO3)|is.infinite(concNO3)]=0; concNO3[tmpCond]=0
	                concNH4[] = columnResource_matrix[,'waterNH4']/columnResource_matrix[,'discharge']; concNH4[is.na(concNH4)|is.infinite(concNH4)]=0; concNH4[tmpCond]=0
	                concPO4[] = columnResource_matrix[,'waterPO4']/columnResource_matrix[,'discharge']; concPO4[is.na(concPO4)|is.infinite(concPO4)]=0; concPO4[tmpCond]=0
	                
	                # hyporheicResource_matrix (vol and mg per unit area)
	                
	                
	                
	                	## ------- constructing hyporheic 
		                #hyporheicResource_matrix[,'discharge'] = zoneEnv_matrix[,BTS,'baseQ'] ## what is it? "m3/m2"
		                hyporheicResource_matrix[,'waterNO3'] = hyporheicResource_matrix[,'baseNO3'] + hyporheicResource_matrix[,'baseNO3balance']
		                	#hyporheicResource_matrix[,'baseNH4'] = 0; #hyporheicResource_matrix[,'discharge']* concNH4
		                hyporheicResource_matrix[,'waterNH4'] = hyporheicResource_matrix[,'baseNH4'] + hyporheicResource_matrix[,'baseNH4balance']
		                	hyporheicResource_matrix[,'baseNO3'] = hyporheicResource_matrix[,'waterNO3']
		                	hyporheicResource_matrix[,'baseNH4'] = hyporheicResource_matrix[,'waterNH4'] 
		                #hyporheicResource_matrix[,'waterPO4'] = hyporheicResource_matrix[,'discharge']* concPO4 ##<<----- come back here
	                
	                
	                hypoconcNO3 = hyporheicResource_matrix[,'waterNO3']/hyporheicResource_matrix[,'discharge'] ##<<----- come back here
	                hypoconcNH4 = hyporheicResource_matrix[,'waterNH4']/hyporheicResource_matrix[,'discharge'] ##<<----- come back here
	                hypoconcPO4 = hyporheicResource_matrix[,'waterPO4']/hyporheicResource_matrix[,'discharge'] ##<<----- come back here
						tmpCond[] = concNO3 > hypoconcNO3;
					denitrify_concNO3[tmpCond] = concNO3[tmpCond]; denitrify_concNO3[!tmpCond]=hypoconcNO3[!tmpCond]
					

                ## potential uptake
				#Pot_algalupNO3[] = algal_parameter['umaxN']*concNO3/(algal_parameter['Nhalf']+concNO3) * bioprocess_matrix[,'algalNuptakeF'] * algalQ10f * STSstep
				Pot_algalupNO3[] =  algal_param$Mulholland_max_denitr*(concNO3 ^(1+ algal_param$Mulholland_pow_denitr)) * bioprocess_matrix[,'algalNuptakeF'] * algalQ10f * STSstep *4e-05 # correct by 25gC/m2
					tmpCond[] = Pot_algalupNO3>columnResource_matrix[,'waterNO3']
					Pot_algalupNO3[tmpCond] = columnResource_matrix[tmpCond,'waterNO3']
					# mgN/mgC/s * mgC/m2 * s  * m2 = mgN; algalNuptakeF = mgC/m2 * m2 * fraction
				       
				#Pot_algalupNH4[] = algal_parameter['umaxN']*(concNH4+concNO3)/(algal_parameter['Nhalf']+concNH4+concNO3)*bioprocess_matrix[,'algalNuptakeF']*algalQ10f * STSstep # mgN
				Pot_algalupNH4[] = algal_param$Mulholland_max_denitr*((concNH4+concNO3)^(1+ algal_param$Mulholland_pow_denitr)) *bioprocess_matrix[,'algalNuptakeF']*algalQ10f * STSstep *4e-05
				Pot_algalupNH4[] = Pot_algalupNH4 - Pot_algalupNO3;
				Pot_algalupNH4[Pot_algalupNH4<0] = 0
					tmpCond[] = Pot_algalupNH4>columnResource_matrix[,'waterNH4']
					Pot_algalupNH4[tmpCond] = columnResource_matrix[tmpCond,'waterNH4']
				
				Pot_algalupP[] = algal_parameter['umaxP']*concPO4/(algal_parameter['Phalf']+concPO4)*bioprocess_matrix[,'algalPuptakeF']*algalQ10f * STSstep # mgP
					tmpCond[] = Pot_algalupP>columnResource_matrix[,'waterPO4']
					Pot_algalupP[tmpCond] = columnResource_matrix[tmpCond,'waterPO4']
				
				
				# Mulholland et al. 2008 Fig 2b f(ugN/L) = f(mg/m3);  mgN/m2/s* m2 * s = mgN;  
				bioQ10f = (biochem_parameter['Q10']^(0.1*zoneEnv_matrix[,BTS,'temperature']-1.2) ) ## Jack's equation
				
				# ## mgN *1/s *s = mgN
				#Pot_nitrification[] = biochem_parameter['max_nitrif'] * bioQ10f * 0.5*(columnResource_matrix[,'waterNH4']+hydraulic_matrix[,'wetBenthicArea']*hyporheicResource_matrix[,'waterNH4']) * STSstep 
				Pot_nitrification[] = biochem_parameter_max_nitrif * bioQ10f * 0.5*(columnResource_matrix[,'waterNH4']+hyporheicResource_matrix[,'waterNH4']) * STSstep
				Pot_nitrificationFrac = (0.5*Pot_nitrification - hyporheicResource_matrix[,'waterNH4'])/Pot_nitrification
				Pot_nitrificationFrac[Pot_nitrificationFrac<0] = 0
				Pot_nitrificationFrac = Pot_nitrificationFrac+0.5 
				 
				 
				## immobilization --  mgN/m2/s* m2 * s = mgN
				#Pot_immobNO3 = detritus_parameter['umaxN']*concNO3/(detritus_parameter['Nhalf']+concNO3) * detritusMatrix[,'microbes']*hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep;
				#Pot_immobNH4 = detritus_parameter['umaxN']*(concNH4+concNO3)/(detritus_parameter['Nhalf']+concNH4+concNO3) * detritusMatrix[,'microbes']*hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep - Pot_immobNO3;
				Pot_immobNO3[] = detritus_param$Mulholland_max_denitr*(concNO3 ^(1+ detritus_param$Mulholland_pow_denitr)) * hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep *detritusMatrix[,'microbes'] *4e-05;
					tmpCond[] = Pot_immobNO3>columnResource_matrix[,'waterNO3']
					Pot_immobNO3[tmpCond] = columnResource_matrix[tmpCond,'waterNO3']
				Pot_immobNH4[] = detritus_param$Mulholland_max_denitr*((concNH4+concNO3)^(1+ detritus_param$Mulholland_pow_denitr)) * hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep * detritusMatrix[,'microbes'] *4e-05 - Pot_immobNO3;
				Pot_immobNH4[Pot_immobNH4 <0] = 0
					tmpCond[] = Pot_immobNH4>columnResource_matrix[,'waterNH4']
					Pot_immobNH4[tmpCond] = columnResource_matrix[tmpCond,'waterNH4']
				Pot_immobPO4[] = algal_parameter['umaxP']*concPO4/(algal_parameter['Phalf']+concPO4) *hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep *detritusMatrix[,'microbes'] ## no need for 4e-05!!
                	tmpCond[] = Pot_immobPO4>columnResource_matrix[,'waterPO4']
					Pot_immobPO4[tmpCond] = columnResource_matrix[tmpCond,'waterPO4']
                #  Pot_algalupP[] = algal_parameter['umaxP']*concPO4/(algal_parameter['Phalf']+concPO4)*bioprocess_matrix[,'colonizedBenthicArea']*algalQ10f * STSstep * bioprocess_matrix[,'algalc'] # mgP
                #  Pot_algalupP[] = algal_parameter['umaxP']*concPO4/(algal_parameter['Phalf']+concPO4)*bioprocess_matrix[,'algalPuptakeF']*algalQ10f * STSstep # mgP
                #  bioprocess_matrix[,'algalPuptakeF'] = bioprocess_matrix[,'algalc']*(1-bioprocess_matrix[,'algalPF'])*bioprocess_matrix[,'colonizedBenthicArea']
                
                ## make sure microbes are not taking more than potential
                tmpCond[] = detritusMatrix[,'netmineral']>0 # 0.5/cn - 1/food_CN -> uptake [inside STS loop]: wN is cumsum in STS loop
                if(sum(tmpCond)>0){
                    hold[tmpCond] = (detritusMatrix[tmpCond,'netmineral'] - detritusMatrix[tmpCond,'wN']) * hydraulic_matrix[tmpCond,'wetBenthicArea'] # unit neededN vs uptakenN
                    Pot_immobNO3[tmpCond] = sapply(vectorMin( Pot_immobNO3[tmpCond], hold[tmpCond]),function(x){ifelse(x>0,x,0.0)})
                    Pot_immobNH4[tmpCond] = sapply(vectorMin( hold[tmpCond]-Pot_immobNO3[tmpCond], Pot_immobNH4[tmpCond]),function(x){ifelse(x>0,x,0.0)})
                }#if
                tmpCond[] = detritusMatrix[,'netmineralP']>0 # 0.5/cn - 1/food_CN -> uptake [inside STS loop]: wN is cumsum in STS loop
                if(sum(tmpCond)>0){
                    hold[tmpCond] = (detritusMatrix[tmpCond,'netmineralP'] - detritusMatrix[tmpCond,'wP']) * hydraulic_matrix[tmpCond,'wetBenthicArea'] # unit neededN vs uptakenN
                    Pot_immobPO4[tmpCond] = sapply(vectorMin( Pot_immobPO4[tmpCond], hold[tmpCond]),function(x){ifelse(x>0,x,0.0)})
                }#if
				 
                ## denitrification (new June 17)
                #Pot_denitrification[] = hydraulic_matrix[,'wetBenthicArea']*biochem_parameter['Mulholland_max_denitr']*bioQ10f*(denitrify_concNO3 ^(1+biochem_parameter['Mulholland_pow_denitr'])) * STSstep
                hold[] = detritusMatrix[,'resp'] + detritusMatrix[,'decay']*(1-detritus_parameter['gcoef']) + detritusMatrix[,'decay_miner']-detritusMatrix[,'growth_miner'];
                Pot_denitrification[] = hydraulic_matrix[,'wetBenthicArea']*biochem_parameter_Mulholland_max_denitr*(denitrify_concNO3 ^(1+biochem_parameter_Mulholland_pow_denitr)) * bioQ10f * STSstep
                Pot_denitrification[hold<=0] = 0 # no decay activity
                Pot_denitrificationFrac[] = (Pot_denitrification - hyporheicResource_matrix[,'waterNO3'])/Pot_denitrification;
                Pot_denitrificationFrac[is.na(Pot_denitrificationFrac)]=0
                Pot_denitrificationFrac[Pot_denitrificationFrac<0] = 0 # for water column
				 
				## ------ splitting resources to processes
				tmpCond[] = detritusMatrix[,'netmineral']>0
				sumPotentialNH4[] = Pot_algalupNH4 + Pot_nitrification*Pot_nitrificationFrac#<----above	*******
                sumPotentialNH4[tmpCond] = sumPotentialNH4[tmpCond] + Pot_immobNH4[tmpCond] ##<<--------------------- immob
				NH4limit_scalar[] = 1; 
					tmpCond =(columnResource_matrix[,'waterNH4'] < sumPotentialNH4); 
					NH4limit_scalar[tmpCond] = columnResource_matrix[tmpCond,'waterNH4']/sumPotentialNH4[tmpCond];
					
				sumPotentialNO3[] = Pot_algalupNO3 + Pot_denitrification*Pot_denitrificationFrac #<----above	
				sumPotentialNO3[tmpCond] = sumPotentialNO3[tmpCond] + Pot_immobNO3[tmpCond] ##<<--------------------- immob
				NO3limit_scalar[] = 1; 
					tmpCond =(columnResource_matrix[,'waterNO3'] < sumPotentialNO3); 
					NO3limit_scalar[tmpCond] = columnResource_matrix[tmpCond,'waterNO3']/sumPotentialNO3[tmpCond];
				
				tmpCond[] = detritusMatrix[,'netmineralP']>0
				sumPotentialPO4[] = Pot_algalupP 
				sumPotentialPO4[tmpCond] = sumPotentialPO4[tmpCond] + Pot_immobPO4[tmpCond]
				Plimit_scalar[] = 1; 
					tmpCond =(columnResource_matrix[,'waterPO4'] < sumPotentialPO4); 
					Plimit_scalar[tmpCond] = columnResource_matrix[tmpCond,'waterPO4']/sumPotentialPO4[tmpCond];
				
                ## remove uptaken solute from water column
                sumPotentialNO3[] = sumPotentialNO3 * NO3limit_scalar; 
                	sumPotentialNO3[] = sumPotentialNO3/columnResource_matrix[,'waterNO3']; 
                	sumPotentialNO3[is.na(sumPotentialNO3)]=0; 
                	sumPotentialNO3[sumPotentialNO3>1]=1
                sumPotentialNH4[] = sumPotentialNH4 * NH4limit_scalar; ##<<------------------------
                	sumPotentialNH4[] = sumPotentialNH4/columnResource_matrix[,'waterNH4']; 
                	sumPotentialNH4[is.na(sumPotentialNH4)]=0; 
                	sumPotentialNH4[sumPotentialNH4>1]=1
                sumPotentialPO4[] = sumPotentialPO4 * Plimit_scalar; 
                	sumPotentialPO4[] = sumPotentialPO4/columnResource_matrix[,'waterPO4']; 
                	sumPotentialPO4[is.na(sumPotentialPO4)]=0; 
                	sumPotentialPO4[sumPotentialPO4>1]=1
                    
                    ##
                    columnResource_matrix[,'waterNO3'] = columnResource_matrix[,'waterNO3'] * (1-sumPotentialNO3) + Pot_nitrificationFrac* Pot_nitrification
                    columnResource_matrix[,'waterNH4'] = columnResource_matrix[,'waterNH4'] * (1-sumPotentialNH4)
                    columnResource_matrix[,'waterPO4'] = columnResource_matrix[,'waterPO4'] * (1-sumPotentialPO4)
                    hyporheicResource_matrix[,'waterNO3'] = hyporheicResource_matrix[,'waterNO3'] - 
                        Pot_denitrification*(1.0-Pot_denitrificationFrac) +  # *bioprocess_matrix[,'colonizedBenthicArea_1']
                        Pot_nitrification*(1.0-Pot_nitrificationFrac); #*bioprocess_matrix[,'colonizedBenthicArea_1'];
                        hyporheicResource_matrix[hyporheicResource_matrix[,'waterNO3']<0,'waterNO3']=0
                        hyporheicResource_matrix[,'waterNH4'] = hyporheicResource_matrix[,'waterNH4'] - Pot_nitrification*(1.0-Pot_nitrificationFrac); #*bioprocess_matrix[,'colonizedBenthicArea_1'];
                        hyporheicResource_matrix[hyporheicResource_matrix[,'waterNH4']<0,'waterNH4']=0
                    
                    Pot_algalupNO3[] = Pot_algalupNO3 * NO3limit_scalar *bioprocess_matrix[,'colonizedBenthicArea_1']# unit area
                    Pot_algalupNH4[] = Pot_algalupNH4 * NH4limit_scalar *bioprocess_matrix[,'colonizedBenthicArea_1']# unit area
                    Pot_algalupP[] = Pot_algalupP * Plimit_scalar *bioprocess_matrix[,'colonizedBenthicArea_1']# unit area
                
                	detritusMatrix[,'wN'] = detritusMatrix[,'wN'] + (detritusMatrix[,'netmineral']>0)*Pot_immobNO3*NO3limit_scalar*bioprocess_matrix[,'colonizedBenthicArea_1']
                	detritusMatrix[,'wN'] = detritusMatrix[,'wN'] + (detritusMatrix[,'netmineral']>0)*Pot_immobNH4*NH4limit_scalar*bioprocess_matrix[,'colonizedBenthicArea_1']
                	detritusMatrix[,'wP'] = detritusMatrix[,'wP'] + (detritusMatrix[,'netmineralP']>0)*Pot_immobPO4*Plimit_scalar*bioprocess_matrix[,'colonizedBenthicArea_1']
                	
				Pot_denitrification[] = Pot_denitrification * (1.0 - Pot_denitrificationFrac + NO3limit_scalar*Pot_denitrificationFrac) *hydraulic_matrix[,'wetBenthicArea_1']# unit area
				Pot_nitrification[] = Pot_nitrification * (1.0 - Pot_nitrificationFrac + NH4limit_scalar*Pot_nitrificationFrac) *hydraulic_matrix[,'wetBenthicArea_1']# unit area
				
				bioprocess_matrix[,'algalupNO3'] <- bioprocess_matrix[,'algalupNO3'] + Pot_algalupNO3 # mgN/m2 # unit area
				bioprocess_matrix[,'algalupNH4'] <- bioprocess_matrix[,'algalupNH4'] + Pot_algalupNH4 # mgN/m2 # unit area
				bioprocess_matrix[,'algalupP'] <- bioprocess_matrix[,'algalupP'] + Pot_algalupP # mgP/m2 # unit area
				bioprocess_matrix[,'nitrification'] <- bioprocess_matrix[,'nitrification'] + Pot_nitrification # unit area
				bioprocess_matrix[,'denitrification'] <- bioprocess_matrix[,'denitrification'] + Pot_denitrification # unit area
                #print(paste(STS, Pot_denitrification,sep=":"))
                
                
                
                
                
                
                ## hyporheic exchange and process *******
                	# zoneEnv_matrix[,BTS,'exchangeProp2SRIP'] m3/m2/s
                	# zoneEnv_matrix[,BTS,'exchangeProp2STR'] m3/m2/s
                #hold = 1.5*zoneEnv_matrix[,BTS,'exchangeProp2SRIP']*hydraulic_matrix[,'wetBenthicArea']/columnResource_matrix[,'discharge'] * STSstep; hold[hold>1]=1; ## from str 2 rip
				#zoneEnv_matrix[,BTS,'exchange'] = 1.5*zoneEnv_matrix[,BTS,'exchangeProp2STR']*hydraulic_matrix[,'wetBenthicArea']/hyporheicResource_matrix[,'discharge'] * STSstep ## from rip 2 str
				
				hold[] = zoneEnv_matrix_exchange * STSstep; hold[hold>1]=1; 
				holdII[] = hold[] * hydraulic_matrix[,'crossarea'] / hydraulic_matrix[,'As'] 
				
				
				concMassChange[] = 0
				concChange[] = (hypoconcNO3[] - concNO3[]) #*zoneEnv_matrix[tmpCond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
					tmpCond = concChange>0 # positive "concChange" mean increase channel water[NO3]
					if(sum(tmpCond)>0) concMassChange[tmpCond] = vectorMin(
						concChange[tmpCond] * hold[tmpCond] * columnResource_matrix[tmpCond,'discharge'],# NO3 mass in channel increase due to exchange
						hyporheicResource_matrix[tmpCond,'waterNO3'] ) ## <<---- should time barea? it's total mass, not areal in the new model
					tmpCond = concChange<0 # negative "concChange" mean increase riparian water[NO3]
					if(sum(tmpCond)>0) concMassChange[tmpCond] = -vectorMin(
						-concChange[tmpCond] * holdII[tmpCond] *hyporheicResource_matrix[tmpCond,'discharge'], # NO3 mass in ztorage increase due to exchange
						columnResource_matrix[tmpCond,'waterNO3']  ) ##<<-----	
                        hyporheicResource_matrix[,'waterNO3'] = hyporheicResource_matrix[,'waterNO3'] - concMassChange[] ##<<---- correct with barea? it's total mass, not areal in the new model
					columnResource_matrix[,'waterNO3'] = columnResource_matrix[,'waterNO3'] + concMassChange[]
					bioprocess_matrix[,'exchangeNO3'] = bioprocess_matrix[,'exchangeNO3'] + concMassChange[]*hydraulic_matrix[,'wetBenthicArea_1'] # areal
					
									
				concMassChange[] = 0	
				concChange[] = (hypoconcNH4[] - concNH4[]) #*zoneEnv_matrix[tmpCond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
				tmpCond = concChange>0 # positive "concChange" mean increase channel water[NO3]
					if(sum(tmpCond)>0) concMassChange[tmpCond] = vectorMin(
						concChange[tmpCond] * hold[tmpCond] * columnResource_matrix[tmpCond,'discharge'],
						hyporheicResource_matrix[tmpCond,'waterNH4'] )
					tmpCond = concChange<0 # negative "concChange" mean increase riparian water[NO3]
					if(sum(tmpCond)>0) concMassChange[tmpCond] = -vectorMin(
						-concChange[tmpCond] * holdII[tmpCond] *hyporheicResource_matrix[tmpCond,'discharge'], 
						columnResource_matrix[tmpCond,'waterNH4'] )
					hyporheicResource_matrix[,'waterNH4'] = hyporheicResource_matrix[,'waterNH4'] - concMassChange[]
					columnResource_matrix[,'waterNH4'] = columnResource_matrix[,'waterNH4'] + concMassChange[]
					bioprocess_matrix[,'exchangeNH4'] = bioprocess_matrix[,'exchangeNH4'] + concMassChange[]*hydraulic_matrix[,'wetBenthicArea_1']
					
				concMassChange[] = 0	
				concChange[] = (hypoconcPO4[] - concPO4[]) #*zoneEnv_matrix[tmpCond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
				tmpCond = concChange>0 # positive "concChange" mean increase channel water[NO3]
					if(sum(tmpCond)>0) concMassChange[tmpCond] = vectorMin(
						concChange[tmpCond] * hold[tmpCond] * columnResource_matrix[tmpCond,'discharge'],
						hyporheicResource_matrix[tmpCond,'waterPO4'] )
					tmpCond = concChange<0 # negative "concChange" mean increase riparian water[NO3]
					if(sum(tmpCond)>0) concMassChange[tmpCond] = -vectorMin(
						-concChange[tmpCond] * holdII[tmpCond] *hyporheicResource_matrix[tmpCond,'discharge'], 
						columnResource_matrix[tmpCond,'waterPO4'] )
					hyporheicResource_matrix[,'waterPO4'] = hyporheicResource_matrix[,'waterPO4'] - concMassChange[]
					columnResource_matrix[,'waterPO4'] = columnResource_matrix[,'waterPO4'] + concMassChange[]
                	bioprocess_matrix[,'exchangePO4'] = bioprocess_matrix[,'exchangePO4'] + concMassChange[]*hydraulic_matrix[,'wetBenthicArea_1']
                	
  					# columnResource_matrix
  					# hyporheicResource_matrix
  					# concNO3 = columnResource_matrix[,'waterNO3']/columnResource_matrix[,'discharge']; concNO3[is.na(concNO3)|is.infinite(concNO3)]=0
  					# hypoconcNO3 = hyporheicResource_matrix[,'waterNO3']/hyporheicResource_matrix[,'discharge']
  					# cbind(concNO3, hypoconcNO3)
  					# plot(concNO3, type='b'); abline(v=8,col='red')
  					# plot(bioprocess_matrix[,'algalupNO3'], type='b')
  					
  				hyporheicResource_matrix[,'baseNH4balance'] = hyporheicResource_matrix[,'waterNH4']-hyporheicResource_matrix[,'baseNH4']	
                hyporheicResource_matrix[,'baseNO3balance'] = hyporheicResource_matrix[,'waterNO3']-hyporheicResource_matrix[,'baseNO3']  
          		tmpCond = hyporheicResource_matrix[,'baseNO3balance']>0
          		hold[] = 0
                if(sum(tmpCond)>0){
	                # hyporheicDenitrification[tmpCond] = hyporheicDenitrification[tmpCond] + vectorMin(
	                	# zoneEnv_matrix[tmpCond,BTS,'Nuptake']*zoneInfo[tmpCond,'len'],
	                	# hyporheicResource_matrix[tmpCond,'baseNO3balance'])
	                # hyporheicResource_matrix[tmpCond,'baseNO3balance'] = hyporheicResource_matrix[tmpCond,'baseNO3balance'] - vectorMin(
	                	# zoneEnv_matrix[tmpCond,BTS,'Nuptake']*zoneInfo[tmpCond,'len'],
	                	# hyporheicResource_matrix[tmpCond,'baseNO3balance']) 
	                	
	                hold[tmpCond] = vectorMin(
		                # zoneEnv_matrix[,,'Nuptake']
	                	#zoneEnv_matrix$landuptake[BTS, tmpCond]*zoneInfo[tmpCond,'len'], 
                        zoneEnv_matrix[tmpCond, BTS,'Nuptake']*zoneInfo[tmpCond,'len'], # unit for this uptake is mgN/m2/m, right?
	                	hyporheicResource_matrix[tmpCond,'baseNO3balance']*hydraulic_matrix[tmpCond,'wetBenthicArea_1'] )
                	
	                hyporheicDenitrification[tmpCond] = hyporheicDenitrification[tmpCond] + hold[tmpCond]
                    hyporheicResource_matrix[tmpCond,'baseNO3balance'] = hyporheicResource_matrix[tmpCond,'baseNO3balance'] - hold[tmpCond]*hydraulic_matrix[tmpCond,'wetBenthicArea'] # consistent to the total mass assumption
                }# if
                	
                	hold = Pot_algalupNO3*bioprocess_matrix[,'colonizedBenthicArea']*hydraulic_matrix[,'wetBenthicArea_1'] + Pot_denitrification + Pot_immobNO3 - Pot_nitrification # <<----------this is all NO3
                    ## here Pot_algalupNO3 was corrected for colonizedBenthicArea, not all BenthicArea!
				bioprocess_matrix[,'totaluptake'] = bioprocess_matrix[,'totaluptake'] + hold #areal uptake
					hold = columnResource_matrix[,'waterNO3']/columnResource_matrix[,'discharge']*hydraulic_matrix[,'depth']*hydraulic_matrix[,'velocity']/hold
                bioprocess_matrix[,'Sw'] = bioprocess_matrix[,'Sw'] + hold*STSstep # uptake length
                
                ###################################################################################### update hydraulic parameters
                crossSection = crossSection + hydraulic_matrix[,'crossarea']
                
                hydraulicHOLD=sapply(numZoneIndex,function(i){ c(
					Morphology_matrix[[i]]@Vol2Z(columnResource_matrix[i,'discharge']), #1
		            Morphology_matrix[[i]]@Vol2V(columnResource_matrix[i,'discharge']), #2
		            Morphology_matrix[[i]]@Vol2WP(columnResource_matrix[i,'discharge']), #3
		            Morphology_matrix[[i]]@Vol2Q(columnResource_matrix[i,'discharge']), #4
                    Morphology_matrix[[i]]@Vol2CA(columnResource_matrix[i,'discharge']), # 5
                    Morphology_matrix[[i]]@maxQ, #5 -> 6
                    Morphology_matrix[[i]]@maxVol ) }) #6 -> 7
				hydraulic_matrix[,'depth'] = hydraulicHOLD[1,]
				hydraulic_matrix[,'velocity'] = hydraulicHOLD[2,] #m/s
				hydraulic_matrix[,'wetBenthicArea'] = zoneInfo[,'len']*hydraulicHOLD[3,]
				hydraulic_matrix[,'wetBenthicArea_1'] = 1.0 / hydraulic_matrix[,'wetBenthicArea']
                hydraulic_matrix[,'crossarea'] = hydraulicHOLD[5,]
				hydraulic_matrix[,'discharge'] = hydraulicHOLD[4,]

                #print(paste(c(BTS,STS, columnResource_matrix[,'discharge']),collapse=' '))
                if(is.na(sum(columnResource_matrix[,'discharge']))|sum(columnResource_matrix[,'discharge']<0)>0) 
                	print(paste(c(BTS,STS, STSstep,round(abs(hydraulic_matrix[,'discharge'])/columnResource_matrix[,'discharge']*STSstep,6)),collapse=' '))
                	
                if(sum(hyporheicResource_matrix[,'waterNH4']<0)>0)
                	print(paste(c(BTS,STS, STSstep, hyporheicResource_matrix[,'waterNH4']),collapse=' ')) 
                # debug
                # print(	paste(
                	# STS, '|',
                	# paste(round(columnResource_matrix[1:3,'discharge'],3),collapse=','),'|',
                	# paste(round(STSexport[1:3,'discharge'],6),collapse=','),'|',
                	# paste(round(STSinput[1:3,'discharge'],6),collapse=',')
                	# ))		
                STSinput[ ] = 0	
             
             	STS = STS + STSstep
             	exportProp = hydraulic_matrix[,'discharge']/columnResource_matrix[,'discharge'];
             	STSstep = min(timemax_STS-STS,max(1,apply(exportProp %o% STSstepOption<0.25,2,prod)* STSstepOption))
		}#STS
		
		bioprocess_matrix[,'Sw'] = bioprocess_matrix[,'Sw']* 0.0002777778 ## correcting aggregating
			
		## debugging
		# hyporheicResource_matrix
		# columnResource_matrix
		# bioprocess_matrix
		# detritusMatrix[,'requiredN']
		# detritusMatrix[,'dN']+detritusMatrix[,'wN']
		# mean(bioprocess_matrix[,'denitrification']*24) # 17.34134
		
		
		# end.time <- Sys.time()
		# time.taken <- end.time - start.time
		# time.taken	
		# RiparianResource_matrix
		# columnResource_matrix
		# columnResource_matrix[,'waterNO3']/columnResource_matrix[,'discharge']
		# hyporheicResource_matrix
		# STSaggregateVar/3600
		# hydraulic_matrix
		# cbind(concNO3, concNH4, concPO4, hyporheicResource_matrix)
		# cbind(concNO3,hyporheicResource_matrix[,'waterNO3'],hyporheicResource_matrix[,'waterNO3'] + hyporheic_nitrification)
		# 
			
		## algal process (hourly) # unit area
		
			## NPP
			bioprocess_matrix[,'algalNF'] = 1-(bioprocess_matrix[,'algalCN']-algal_parameter['minCN'])/(algal_parameter['maxCN']-algal_parameter['minCN']); 
				bioprocess_matrix[bioprocess_matrix[,'algalNF']<0,'algalNF']=0
				bioprocess_matrix[bioprocess_matrix[,'algalNF']>1,'algalNF']=1
				
			bioprocess_matrix[,'algalPF'] = 1-(bioprocess_matrix[,'algalCP']-algal_parameter['minCP'])/(algal_parameter['maxCP']-algal_parameter['minCP'])
				bioprocess_matrix[bioprocess_matrix[,'algalPF']<0,'algalPF']=0
				bioprocess_matrix[bioprocess_matrix[,'algalPF']>1,'algalPF']=1
			bioprocess_matrix[,'algalnutrientF'] <- apply(bioprocess_matrix[,c('algalNF','algalPF')],1,min)
			
			bioprocess_matrix[,'algalselfF'] <- 1 / (1 + algal_parameter['selflimitcoef']* bioprocess_matrix[,'algalc']) # algalc mgC/m2; selflimitcoef m2/mgC
			
			# exp(-(zoneEnv_matrix[,'PAR',BTS]-b)/a)
			# month = as.numeric(format(SimulationTime[as.integer(BTS/24)],'%m'))
            # if( as.numeric(format(SimulationTime[ceiling(BTS/24)],'%m')) %in% c(6,7,8,9,10)){
            #    bioprocess_matrix[,'algallightF'] <- zoneEnv_matrix[,BTS,'PARf'] #* lightAdjust # 6,7,8,9,10 month
            # }else{
            # 	bioprocess_matrix[,'algallightF'] <- zoneEnv_matrix[,BTS,'PARf'] #* lightAdjustII
            # }#
            bioprocess_matrix[,'algallightF'] <- zoneEnv_matrix[,BTS,'PARf'] #*(0.9+0.1*exp(-0.1666667*columnResource_matrix[,'FOC']/columnResource_matrix[,'discharge']));
            
            
			bioprocess_matrix[,'algalNPP'] <- bioprocess_matrix[,'algalc']* algal_parameter['maxgrowthRate'] * bioprocess_matrix[,'algalnutrientF'] * bioprocess_matrix[,'algalselfF'] * bioprocess_matrix[,'algallightF'] * algalQ10f * timemax_STS * (1-algalhabitatDmgIndex)
					
	
		bioprocess_matrix[,'resp'] <- bioprocess_matrix[,'algalNPP']*(1-algal_parameter['gcoef']) + bioprocess_matrix[,'algalMineralC']
		bioprocess_matrix[,'algalc'] <- bioprocess_matrix[,'algalc'] + algal_parameter['gcoef']*bioprocess_matrix[,'algalNPP'] 
		bioprocess_matrix[,'algaln'] <- bioprocess_matrix[,'algaln'] + bioprocess_matrix[,'algalupNO3'] + bioprocess_matrix[,'algalupNH4']
		bioprocess_matrix[,'algalp'] <- bioprocess_matrix[,'algalp'] + bioprocess_matrix[,'algalupP']
		
		## convert areal flux to total flux to exchange with water column resources
        #bioprocess_matrix[,'algalMineralC'] = bioprocess_matrix[,'algalMineralC'] * bioprocess_matrix[,'colonizedBenthicArea'] # from unit area to full area
        #bioprocess_matrix[,'algalMineralN'] = bioprocess_matrix[,'algalMineralN'] * bioprocess_matrix[,'colonizedBenthicArea']
        #bioprocess_matrix[,'algalMineralP'] = bioprocess_matrix[,'algalMineralP'] * bioprocess_matrix[,'colonizedBenthicArea']
        #bioprocess_matrix[,'algalentC'] = bioprocess_matrix[,'algalentC'] * bioprocess_matrix[,'colonizedBenthicArea'] # from unit area to full area
        #bioprocess_matrix[,'algalentN'] = bioprocess_matrix[,'algalentN'] * bioprocess_matrix[,'colonizedBenthicArea']
        #bioprocess_matrix[,'algalentP'] = bioprocess_matrix[,'algalentP'] * bioprocess_matrix[,'colonizedBenthicArea']
		
		## update C:N and C:P
		bioprocess_matrix[,'algalCN'] <- bioprocess_matrix[,'algalc']/bioprocess_matrix[,'algaln']
		bioprocess_matrix[,'algalCP'] <- bioprocess_matrix[,'algalc']/bioprocess_matrix[,'algalp']
        bioprocess_matrix[,'algalNuptakeF'] <- bioprocess_matrix[,'algalc']*(1-bioprocess_matrix[,'algalNF'])*bioprocess_matrix[,'colonizedBenthicArea']
        bioprocess_matrix[,'algalPuptakeF'] <- bioprocess_matrix[,'algalc']*(1-bioprocess_matrix[,'algalPF'])*bioprocess_matrix[,'colonizedBenthicArea']
        
		dailyAverageWetBenthicArea[] = dailyAverageWetBenthicArea  + hydraulic_matrix[,'wetBenthicArea']* 0.04166667
		
	
		# decomposition process (hourly ... add later)	
		tmpCond[] = detritusMatrix[,'netmineral']>0 # 0.5/cn - 1/CN (uptake)
        if(sum(tmpCond)>0) detritusMatrix[tmpCond,'decay'] = vectorMin(detritusMatrix[tmpCond,'decay'], detritusMatrix[tmpCond,'wN']/detritusMatrix[tmpCond,'Ncoef']) ## N can limit decay
        tmpCond[] = detritusMatrix[,'netmineralP']>0 # 0.5/cn - 1/CN (uptake)
        if(sum(tmpCond)>0) detritusMatrix[tmpCond,'decay'] = vectorMin(detritusMatrix[tmpCond,'decay'], detritusMatrix[tmpCond,'wP']/detritusMatrix[tmpCond,'Pcoef']) ## P can limit decay
        
		hold[] = detritusMatrix[,'decay']/(detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']* immo_class2decayRatio +detritusMatrix[,'detritus3c']* immo_class3decayRatio) # decay by immo and miner
		holdII[] = detritusMatrix[,'decay_miner']/(detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']* miner_class2decayRatio +detritusMatrix[,'detritus3c']* miner_class3decayRatio)
		holdII[holdII>1] = 1.0
		hold[hold>1]=1.0
			# cbind(
				# detritusMatrix[,'detritus1c']+detritusMatrix[,'detritus2c']* miner_class2decayRatio +detritusMatrix[,'detritus3c']* miner_class3decayRatio,
				# detritusMatrix[,'decay_miner']
			# )
		detritusMatrix[,'detritus1c'] = detritusMatrix[,'detritus1c'] * (1 - hold) * (1 - holdII) # mgC/m2 unit area
		detritusMatrix[,'detritus1n'] = detritusMatrix[,'detritus1n'] * (1 - hold) * (1 - holdII)# mgN/m2 unit area
		detritusMatrix[,'detritus1p'] = detritusMatrix[,'detritus1p'] * (1 - hold) * (1 - holdII)
		detritusMatrix[,'detritus2c'] = detritusMatrix[,'detritus2c'] * (1 - hold* immo_class2decayRatio) * (1 - holdII* miner_class2decayRatio)
		detritusMatrix[,'detritus2n'] = detritusMatrix[,'detritus2n'] * (1 - hold* immo_class2decayRatio) * (1 - holdII* miner_class2decayRatio)
		detritusMatrix[,'detritus2p'] = detritusMatrix[,'detritus2p'] * (1 - hold* immo_class2decayRatio) * (1 - holdII* miner_class2decayRatio)
		detritusMatrix[,'detritus3c'] = detritusMatrix[,'detritus3c'] * (1 - hold* immo_class3decayRatio) * (1 - holdII* miner_class3decayRatio)
		detritusMatrix[,'detritus3n'] = detritusMatrix[,'detritus3n'] * (1 - hold* immo_class3decayRatio) * (1 - holdII* miner_class3decayRatio)
		detritusMatrix[,'detritus3p'] = detritusMatrix[,'detritus3p'] * (1 - hold* immo_class3decayRatio) * (1 - holdII* miner_class3decayRatio)
		
		detritusMatrix[,'microbes'] = detritusMatrix[,'microbes'] + detritusMatrix[,'decay']*detritus_parameter['gcoef']
		detritusMatrix[,'miner'] = detritusMatrix[,'miner'] + detritusMatrix[,'decay']*detritus_parameter['gcoef']
        detritusMatrix[,'resp'] = detritusMatrix[,'resp'] + detritusMatrix[,'decay']*(1-detritus_parameter['gcoef']) + detritusMatrix[,'decay_miner']-detritusMatrix[,'growth_miner']
		detritusMatrix[,'mineral'] = 0
		tmpCond = detritusMatrix[,'netmineral']<0
		if(sum(tmpCond)>0) detritusMatrix[tmpCond,'mineral'] = -detritusMatrix[tmpCond,'netmineral']
		detritusMatrix[,'mineral'] = detritusMatrix[,'mineral'] + detritusMatrix[,'decay_minerN'] ## for miner to release
		
		# hourly aggregate
        BTSaggregateVar[,'velocity'] = BTSaggregateVar[,'velocity'] + hydraulic_matrix[,'velocity'] ## need to be /24
        BTSaggregateVar[,'depth'] = BTSaggregateVar[,'depth'] + hydraulic_matrix[,'depth'] ## need to be /24
        BTSaggregateVar[,'benthicArea'] = BTSaggregateVar[,'benthicArea'] + hydraulic_matrix[,'wetBenthicArea'] ## need to be /24
        BTSaggregateVar[,'denitrification'] = BTSaggregateVar[,'denitrification'] + bioprocess_matrix[,'denitrification'] # unit area
        BTSaggregateVar[,'nitrification'] = BTSaggregateVar[,'nitrification'] + bioprocess_matrix[,'nitrification'] # unit area
        BTSaggregateVar[,'algalNO3uptake'] = BTSaggregateVar[,'algalNO3uptake'] + bioprocess_matrix[,'algalupNO3'] # unit area
        BTSaggregateVar[,'algalNH4uptake'] = BTSaggregateVar[,'algalNH4uptake'] + bioprocess_matrix[,'algalupNH4'] # unit area
        BTSaggregateVar[,'algalNPP'] = BTSaggregateVar[,'algalNPP'] + bioprocess_matrix[,'algalNPP'] # unit area
        BTSaggregateVar[,'mineralN'] = BTSaggregateVar[,'mineralN'] + bioprocess_matrix[,'algalMineralN']
        BTSaggregateVar[,'algalCN'] = BTSaggregateVar[,'algalCN'] + bioprocess_matrix[,'algalCN'] ## need to be /24
        BTSaggregateVar[,'algalCP'] = BTSaggregateVar[,'algalCP'] + bioprocess_matrix[,'algalCP'] ## need to be /24
        BTSaggregateVar[,'algalC'] = BTSaggregateVar[,'algalC'] + bioprocess_matrix[,'algalc'] ## need to be /24
        BTSaggregateVar[,'algalCloss'] = BTSaggregateVar[,'algalCloss'] + bioprocess_matrix[,'algalentC'] + bioprocess_matrix[,'algaldeadC']
        BTSaggregateVar[,'algalR'] = BTSaggregateVar[,'algalR'] + bioprocess_matrix[,'resp']
        BTSaggregateVar[,'detritusR'] = BTSaggregateVar[,'detritusR'] + detritusMatrix[,'resp']
        BTSaggregateVar[,'Sw'] = BTSaggregateVar[,'Sw'] + bioprocess_matrix[,'Sw']## need to be /24
        
		BTSaggregateVar[,'netmineral'] = BTSaggregateVar[,'netmineral'] + detritusMatrix[,'mineral'] - detritusMatrix[,'wN'] # unit area 
		BTSaggregateVar[,'detric'] = BTSaggregateVar[,'detric'] + (detritusMatrix[,'detritus1c']+ detritusMatrix[,'detritus2c']+ detritusMatrix[,'detritus3c'])  ## need to be /24
		
        BTSaggregateVar[,'detriCN'] = BTSaggregateVar[,'detriCN'] + 1/detritusMatrix[,'foodNC'] ## need to be /24 <<----------------------- debugging
        BTSaggregateVar[,'detriCP'] = BTSaggregateVar[,'detriCP'] + 1/detritusMatrix[,'foodPC'] ## need to be /24 <<----------------------- debugging
		
		BTSaggregateVar[,'microbC'] = BTSaggregateVar[,'microbC'] + detritusMatrix[,'microbes']  ## need to be /24 <<----- ***
		BTSaggregateVar[,'minerC'] = BTSaggregateVar[,'minerC'] + detritusMatrix[,'miner']
        
		BTSaggregateVar[,'tmperature'] = BTSaggregateVar[,'tmperature'] + zoneEnv_matrix[,BTS,'temperature'] ## need to be /24
		BTSaggregateVar[,'PARf'] = BTSaggregateVar[,'PARf'] + zoneEnv_matrix[,BTS,'PARf'] ## need to be /24
		BTSaggregateVar[,'hdmg'] = BTSaggregateVar[,'hdmg'] + algalhabitatDmgIndex
        
        BTSaggregateVar[,'algalNlimit'] = BTSaggregateVar[,'algalNlimit'] + bioprocess_matrix[,'algalNF']
        BTSaggregateVar[,'algalPlimit'] = BTSaggregateVar[,'algalPlimit'] + bioprocess_matrix[,'algalPF']
        
		## reset
		bioprocess_matrix[,'algalupP']=0;
		bioprocess_matrix[,'algalupNO3']=0;
		bioprocess_matrix[,'algalupNH4']=0;
		# bioprocess_matrix[,'resp'] # calculated hourly
		# bioprocess_matrix[,'algalMineralC'] # calculated hourly
		# bioprocess_matrix[,'algalMineralN'] # calculated hourly
		# bioprocess_matrix[,'algalMineralP'] # calculated hourly
		# bioprocess_matrix[,'algalentC'] # calculated hourly
		# bioprocess_matrix[,'algalentN'] # calculated hourly
		# bioprocess_matrix[,'algalentP'] # calculated hourly
		bioprocess_matrix[,'nitrification']=0
		bioprocess_matrix[,'denitrification']=0
		detritusMatrix[,'wN'] = 0 ## N immobilization
        detritusMatrix[,'wP'] = 0
		# detritusMatrix[,'resp'] # calculated hourly
		# detritusMatrix[,'Ncoef'] # calculated hourly
		# detritusMatrix[,'netmineral'] # calculated hourly
		# detritusMatrix[,'decay'] # calculated hourly
		# detritusMatrix[,'mineral'] # calculated hourly
		
		if(BTS%%24==0){

			## adjust algal colonized area
				# cond = bioprocess_matrix[,'colonizedBenthicArea'] < dailyAverageWetBenthicArea
				# # crowded, need expansion; bioprocess_matrix[cond,'algalselfF']=low
				# bioprocess_matrix[cond,'colonizedBenthicArea'] = 
					# (0.9666667 - 0.03333333*(1-bioprocess_matrix[cond,'algalselfF']) )* bioprocess_matrix[cond,'colonizedBenthicArea'] + 
					# (0.03333333 + 0.03333333*(1-bioprocess_matrix[cond,'algalselfF']) )* dailyAverageWetBenthicArea[cond] #m2
					
				# cond = bioprocess_matrix[,'colonizedBenthicArea'] > 1.5*dailyAverageWetBenthicArea & dailyAverageWetBenthicArea>0
				# # no density effect
				# bioprocess_matrix[cond,'colonizedBenthicArea'] =
					# 0.99666667*bioprocess_matrix[cond,'colonizedBenthicArea'] + 
					# 0.003333333* dailyAverageWetBenthicArea[cond] #m2
					
				# cond = bioprocess_matrix[,'colonizedBenthicArea'] > maxChannelBenthicArea	
				# bioprocess_matrix[cond,'colonizedBenthicArea'] = maxChannelBenthicArea[cond]
				# bioprocess_matrix[,'colonizedBenthicArea_1'] = 1.0/bioprocess_matrix[,'colonizedBenthicArea']
			
			# bioprocess_matrix[,'algalc'] = bioprocess_matrix[,'algalc'] * bioprocess_matrix[,'colonizedBenthicArea_1']
			# bioprocess_matrix[,'algaln'] = bioprocess_matrix[,'algaln'] * bioprocess_matrix[,'colonizedBenthicArea_1']
			# bioprocess_matrix[,'algalp'] = bioprocess_matrix[,'algalp'] * bioprocess_matrix[,'colonizedBenthicArea_1']
			
			
            
            ### new output , after AGU
            BTSaggregateVar[,'velocity'] = BTSaggregateVar[,'velocity'] * 0.04166667 # over by 24
            BTSaggregateVar[,'depth'] = BTSaggregateVar[,'depth'] * 0.04166667
            BTSaggregateVar[,'benthicArea'] = BTSaggregateVar[,'benthicArea'] * 0.04166667
            BTSaggregateVar[,'algalCN'] = BTSaggregateVar[,'algalCN'] * 0.04166667
            BTSaggregateVar[,'algalCP'] = BTSaggregateVar[,'algalCP'] * 0.04166667
            BTSaggregateVar[,'algalC'] = BTSaggregateVar[,'algalC'] * 0.04166667
            BTSaggregateVar[,'detric'] = BTSaggregateVar[,'detric']* 0.04166667
            BTSaggregateVar[,'detriCN'] = BTSaggregateVar[,'detriCN']* 0.04166667
            BTSaggregateVar[,'detriCP'] = BTSaggregateVar[,'detriCP']* 0.04166667
            BTSaggregateVar[,'microbC'] = BTSaggregateVar[,'microbC']* 0.04166667
            BTSaggregateVar[,'minerC'] = BTSaggregateVar[,'minerC']* 0.04166667
            BTSaggregateVar[,'Sw'] = BTSaggregateVar[,'Sw'] * 0.04166667
            BTSaggregateVar[,'tmperature'] = BTSaggregateVar[,'tmperature'] * 0.04166667
			BTSaggregateVar[,'PARf'] = BTSaggregateVar[,'PARf'] * 0.04166667
            BTSaggregateVar[,'hdmg'] = BTSaggregateVar[,'hdmg'] * 0.04166667
            BTSaggregateVar[,'algalNlimit'] = BTSaggregateVar[,'algalNlimit'] * 0.04166667
            BTSaggregateVar[,'algalPlimit'] = BTSaggregateVar[,'algalPlimit'] * 0.04166667
            cat( c(
                as.numeric(format(SimulationTime[as.integer(BTS/24)],'%Y')),
                as.numeric(format(SimulationTime[as.integer(BTS/24)],'%m')),
                as.numeric(format(SimulationTime[as.integer(BTS/24)],'%d')), #3
                
                #downZs = 11:60
                #upZs = 1:8
                # upstream (daily) [1-8]
                STSinputTrack[upZs[1],] + colSums(STSLateralinputTrack[upZs[2]:upZs[length(upZs)],]), # inputs: discharge waterNO3 waterNH4 waterPO4 FOC FON FOP (flux/load) from upsteram & lateral IN
                STSaggregateVar[upZs[length(upZs)],], # outputs: discharge waterNO3 waterNH4 waterPO4 FOC FON FOP (flux/load)
                sum(BTSaggregateVar[upZs,'benthicArea']),
                (zoneInfo[upZs,'len']/sum(zoneInfo[upZs,'len'])) %*% BTSaggregateVar[upZs,c('velocity','depth')],
                (BTSaggregateVar[upZs,'benthicArea']/sum(BTSaggregateVar[upZs,'benthicArea'])) %*% BTSaggregateVar[upZs, BTSaggregateVarPrint],
               	100*(1-(STSaggregateVar[upZs[length(upZs)],'nitrate'])/(STSinputTrack[upZs[1],'waterNO3']+sum(STSLateralinputTrack[upZs[2]:upZs[length(upZs)],'waterNO3'])) )/sum(zoneInfo[upZs,'len']),
                
                #downstream (daily) full 11 - 34; odd = pool
                STSinputTrack[downZs[1],] + colSums(STSLateralinputTrack[downZs[2]:downZs[length(downZs)],]), # inputs: discharge waterNO3 waterNH4 waterPO4 FOC FON FOP (flux/load) from upsteram & lateral IN
                STSaggregateVar[downZs[length(downZs)],], # outputs: discharge waterNO3 waterNH4 waterPO4 FOC FON FOP (flux/load)
                sum(BTSaggregateVar[downZs,'benthicArea']),
                (zoneInfo[downZs,'len']/sum(zoneInfo[downZs,'len'])) %*% BTSaggregateVar[downZs,c('velocity','depth')],
                (BTSaggregateVar[downZs,'benthicArea']/sum(BTSaggregateVar[downZs,'benthicArea'])) %*% BTSaggregateVar[downZs, BTSaggregateVarPrint],
                100*(1-(STSaggregateVar[downZs[length(downZs)],'nitrate'])/(STSinputTrack[downZs[1],'waterNO3']+sum(STSLateralinputTrack[downZs[2]:downZs[length(downZs)],'waterNO3'])))/sum(zoneInfo[downZs,'len']),
                
                sum(bioprocess_matrix[upZs,'exchangeNO3']* (BTSaggregateVar[upZs,'benthicArea']/sum(BTSaggregateVar[upZs,'benthicArea'])) ), # it's areal
                sum(bioprocess_matrix[downZs,'exchangeNO3']* (BTSaggregateVar[downZs,'benthicArea']/sum(BTSaggregateVar[downZs,'benthicArea'])) ),
                
                ## input nitrate
                #columnResource_matrix[,'waterNO3'] /columnResource_matrix[,'discharge']
                #STSaggregateVar[,'nitrate'] /STSaggregateVar[,'discharge'] ## export
                #STSaggregateVar # tracking local output
                #STSinputTrack # tracking lcoal lateral and upstream inputs <-----
                #STSLateralinputTrack # tracking lcoal lateral only  <-----
                #downZs = 11:60 
				#upZs = 1:8
                sapply(upZs,function(i){
                    totalinput = STSinputTrack[upZs[1],'waterNO3']
                    if(i>upZs[1]) totalinput = totalinput + sum(STSLateralinputTrack[upZs[2]:i,'waterNO3'])
                    #return <- 100*(1-STSaggregateVar[i,'nitrate']/totalinput)
                    return <- totalinput
                }), # upstream
                sapply(9:10,function(i){
                    totalinput = STSinputTrack[9,'waterNO3']
                    if(i>9) totalinput = totalinput + sum(STSLateralinputTrack[10:i,'waterNO3'])
                    #return <- 100*(1-STSaggregateVar[i,'nitrate']/totalinput)
                    return <- totalinput
                }), # mid transition
                sapply(downZs[1]:max(zoneInfo[,'zoneID']),function(i){
                    totalinput = STSinputTrack[downZs[1],'waterNO3']
                    if(i>downZs[1]) totalinput = totalinput + sum(STSLateralinputTrack[downZs[2]:i,'waterNO3'])
                    #return <- 100*(1-STSaggregateVar[i,'nitrate']/totalinput)
                    return <- totalinput
                }), # downstream
                
                ## output nitrate
                STSaggregateVar[,'nitrate']
                ), '\n', file=outputFile_buff,sep=',')

                
                # downstream (daily) -- last 60 m [27,29,31,33]
                # STSinputTrack[27,] + colSums(STSLateralinputTrack[28:34,]), # inputs: discharge waterNO3 waterNH4 waterPO4 FOC FON FOP (flux/load) from upsteram & lateral IN
                # STSaggregateVar[34,], # outputs: discharge waterNO3 waterNH4 waterPO4 FOC FON FOP (flux/load)
                # sum(BTSaggregateVar[c(27,29,31,33),'benthicArea']),
                # (zoneInfo[c(27,29,31,33),'len']/sum(zoneInfo[c(27,29,31,33),'len'])) %*% BTSaggregateVar[c(27,29,31,33),c('velocity','depth')],
                # (BTSaggregateVar[c(27,29,31,33),'benthicArea']/sum(BTSaggregateVar[c(27,29,31,33),'benthicArea'])) %*% BTSaggregateVar[c(27,29,31,33),c(
                	# 'nitrification','denitrification','algalNO3uptake','algalNH4uptake','algalNPP','mineralN','algalCN','algalCP','algalC',
                	# 'algalCloss','netmineral','detric','detriCN','microbC','algalR','detritusR','Sw')],
                # 100*(1-(STSaggregateVar[34,'nitrate'])/(STSinputTrack[27,'waterNO3']+sum(STSLateralinputTrack[c(27,29,31,33),'waterNO3'])))/sum(zoneInfo[27:34,'len'])
                # ), '\n', file=outputFile_buff,sep=',')
				
			BTSaggregateVar[ ] = 0
			STSaggregateVar[ ] = 0
			dailyAverageWetBenthicArea[ ] = 0
			bioprocess_matrix[,'totaluptake'] = 0
			bioprocess_matrix[,'Sw'] = 0
			bioprocess_matrix[,'exchangeNO3'] = 0
			bioprocess_matrix[,'exchangeNH4'] = 0
			bioprocess_matrix[,'exchangePO4'] = 0
			hyporheicDenitrification[] =0
			STSinputTrack[] = 0
			STSLateralinputTrack[] = 0
		}## daily output -- BTS%%24==0	
		# RiparianResource_matrix
		# cbind(zoneInfo,columnResource_matrix)
		# STSaggregateVar/3600
		# hydraulic_matrix
		# check = data.frame(
			# detritus1 = detritusMatrix[,'detritus1c'],
			# detritus2 = detritusMatrix[,'detritus2c'],
			# detritus3 = detritusMatrix[,'detritus3c'],
			# detritus1cn = detritusMatrix[,'detritus1c']/detritusMatrix[,'detritus1n'],
			# detritus2cn = detritusMatrix[,'detritus2c']/detritusMatrix[,'detritus2n'],
			# detritus3cn = detritusMatrix[,'detritus3c']/detritusMatrix[,'detritus3n'],
			# microbes = detritusMatrix[,'microbes'],
			# miner = detritusMatrix[,'miner'],
			# growth_miner = detritusMatrix[,'growth_miner'],
			# decay_miner = detritusMatrix[,'decay_miner'],
			# decay_minerN = detritusMatrix[,'decay_minerN'],
			# decay_minerP = detritusMatrix[,'decay_minerP'],
			# netmineral = detritusMatrix[,'netmineral'],
			# netmineralP = detritusMatrix[,'netmineralP'],
			# decay = detritusMatrix[,'decay'],
			# foodNC = detritusMatrix[,'foodNC'],
			# foodPC = detritusMatrix[,'foodPC'],
			# foodNC_miner = detritusMatrix[,'foodNC_miner'],
			# foodPC_miner = detritusMatrix[,'foodPC_miner'],
			# resp = detritusMatrix[,'resp'],
			# algaldeadc = bioprocess_matrix[,'algaldeadC']
		# )
		
				
		if(BTS%%720==0) print(paste(BTS,"is done"))#monthly
	}#BTS
	close(outputFile_buff)
})  
##--------------------------------------------------------------- afterward

# # check = read.csv(outputFile)
# daily = read.csv('inputs/daily_totallateralinputs.csv')

# keyTable = cbind(1:dim(check)[2], colnames(check)); 
# distance = zoneInfo[,'len']
# sw = as.matrix(check[,paste('sw',zoneInfo[,1],sep='')])

# result = data.frame(month=check[,'month'])
# result$date = as.Date(paste(check[,'day'], check[,'month'], check[,'year'],sep="-"),format="%d-%m-%Y");
# result$unrestored_sw = (sw[,1:8] %*% distance[1:8])/sum(distance[1:8])
# result$restored_sw = (sw[,27:33] %*% distance[27:33])/sum(distance[27:33])
# result$inputQ = daily[,'discharge']
# result$inputNO3 = daily[,'no3']
# result$inputNO3c = daily[,'no3c']
# result$inputNH4c = daily[,'nh4c']
# result$outputNO3c = check[,'nitrate']/(check[,'discharge']*1000) # (mgN/s) / (L/s) = mgN/L
# result$outputNH4c = check[,'ammonium']/(check[,'discharge']*1000) # (mgN/s) / (L/s) = mgN/L

# plot(result$date,result$unrestored_sw, type='l', log='y', ylim=c(1,60000),ylab='Sw',xlab=''); 
# lines(result$date,result$restored_sw, col='red'); abline(h=25,lty=2,col='gray'); abline(h=100,lty=2,col='gray')

# plot(result$date, result$outputNO3c, type='l', xlab='', ylab='[NO3] mgN/L', log='y', ylim=c(1/1000/1000,2))
# lines(result$date,result$outputNO3c, col='red')

# plot(result$date, result$inputNH4c, type='l', xlab='', ylab='[NO3] mgN/L', log='y', ylim=c(1/1000/1000,2))
# lines(result$date,result$outputNH4c, col='red')









# layout(matrix(1:6,nrow=2)); par(mar=c(4,4,1,1))
# plot(result$inputQ,result$unrestored_sw)
# plot(result$inputQ,result$restored_sw)

# plot(result$inputNc,result$unrestored_sw, log='y')
# plot(result$inputNc,result$restored_sw, log='y')

# plot(result$temp,result$unrestored_sw, log='y')
# plot(result$temp,result$restored_sw, log='y')


