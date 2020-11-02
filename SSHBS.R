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

vectorMin = function(aa,bb){
    sapply(seq_along(aa),function(ii){min(aa[ii],bb[ii])})
}#function

vectorMax = function(aa,bb){
    sapply(seq_along(aa),function(ii){max(aa[ii],bb[ii])})
}#function

yy=seq(0,3000); zz = 1/ ( 1+exp(  -(yy-300)/50 ) ); PAR2PARF = splinefun(yy,zz) # algal light response (Webster, Newbold, and Lin, 2016)

####################################### define constants
algPar = list()
micPar = list()
other_param = list()

# Mulholland et al. 2008 (further defined from below)
algPar$MulhDenB = -0.493
algPar$MulhDenA = 1.059254e-05
micPar$MulhUpB = -0.493
micPar$MulhUpA = 1.059254e-05
micPar$MulhUpB = -0.493
micPar$MulhUpA = 1.059254e-05

# Lin and Webster 2014; ???
micPar$max_nitrif = 0.00002339181/17
micPar$growthRate_im = 4.456975e-06 # log(2)/(1.8*3600*24)
micPar$growthRate_mi = 3.342729e-06 # log(2)/(2.4*3600*24)
micPar$baseR = 1.16e-7
micPar$gcoef = 0.5

# Webster et al. 2009
algPar$umaxP = 8.460648e-08 * 2
algPar$Phalf = 2 # mgP/m3 --> 2µgP/L

# Cross et al. 2005 # from molar ratio to mass ratio
algPar$maxCN = (1/0.0606) *12/14
algPar$maxNC = 1.0 / algPar$maxCN
algPar$maxCP = (1/0.00377) *12/31
algPar$maxPC = 1.0 / algPar$maxCP
algPar$minCN = (10.404) *12/14
algPar$minNC = 1.0 / algPar$minCN
algPar$minCP = (123.817) *12/31
algPar$minPC = 1.0 / algPar$minCP
micPar$im_cn = 7.0
micPar$im_nc = 1.0 / micPar$im_cn
micPar$im_cp = 188
micPar$im_pc = 1.0 / micPar$im_cp
micPar$mi_cn = 5.0
micPar$mi_nc = 1.0 / micPar$mi_cn
micPar$mi_cp = 20
micPar$mi_pc = 1.0 / micPar$mi_cp

# Webster, Newbold, and Lin 2016
algPar$gcoef = 0.8
algPar$selflimitcoef = 0.0015 # m2/mgC
algPar$maxgrowthRate = 2*2*1.0/86400 # 1/s; daily to hourly; half day is in dark; convert NPP to GPP]
algPar$mialRate = 0.01/86400 # 1/s
algPar$Q10 = 2
algPar$moralityRate = 0.004*algPar$maxgrowthRate
micPar$Q10 = 2

# from RHESSys, SSURGO, ...
other_param$detritusCN = 40
other_param$detritusCP = 444.5
other_param$p0 = 0.485 # soil 8 silt_loam
other_param$p_d = 1/4000

# other
STSstepOption = c(1, 2, 4, 8, 24, 48)

####################################### read in header file and other input files
arg=commandArgs(T)
arg=c('proj=./example','header=SSHBS_inputHeader.csv', 'st=2012-1-1', 'ed=2017-12-31', 'nitri=1.3', 'denitri=1', 'upt=1', 'immb=1','output=test.csv')

ModelArg = read.tcsv(text=arg, nfield =2, sep='=') 
ModelLibDenit = read.csv('https://raw.githubusercontent.com/laurencelin/SSHBS/master/lib_denitrification_ceof.csv')
ModelLibUptake = read.csv('https://raw.githubusercontent.com/laurencelin/SSHBS/master/lib_uptake_ceof.csv')

SimulationTime = seq.Date(from=as.Date(ModelArg$st), to=as.Date(ModelArg$ed) ,by="day")
SimulationTimeHour = do.call(c,lapply(SimulationTime,function(j){paste(j,0:23,sep='-')}))
timemax_BTS=24 * length(SimulationTime);
timemax_STS=3600

micPar$max_nitrif = max(5.341868e-10, min(micPar$max_nitrif*ModelArg$nitri, 2.444973e-05))
micPar$MulhDenA = 10^ModelLibDenit[ModelArg$denitri,1]*0.01 # cm/s -> m/s
micPar$MulhDenB = ModelLibDenit[ModelArg$denitri,2] # per second

algPar$MulhUpA = 10^ModelLibUptake[ModelArg$upt,1]*0.01 # cm/s -> m/s
algPar$MulhUpB = ModelLibUptake[ModelArg$upt,2] # per second
algPar$v1ms_left = 0.15

micPar$MulhUpA = 10^ModelLibUptake[ModelArg$immb,1]*0.01 # cm/s -> m/s
micPar$MulhUpB = ModelLibUptake[ModelArg$immb,2] # per second
micPar$v1ms_left = 0.15
rm(ModelLibDenit)
rm(ModelLibUptake)

header = read.csv(paste(ModelArg$proj, ModelArg$header ,sep='/'))
header$reachID = sapply(header$reachID,toString)
reachID = unique(header$reachID)

ZoneConn = lapply(seq_len(dim(header)[1]),function(i){ return <- which(header$downstreamZoneID == header$zoneID[i]) })

Zone = data.frame(
    zoneID = header$zoneID,
    cqIndex = match(header$Cqinput, unique(header$Cqinput)),
    bqFrac = header$lateralBaseFrac,
    sqFrac = header$lateralStormFrac,
    len = header$len,
    lenAdj = 0.1*header$len,
    lenAdj_1 = 10/header$len,
    hydIndex = match(header$hydraulicTable, unique(header$hydraulicTable)),
    parIndex = match(header$hourly_PAR, unique(header$hourly_PAR)),
    degIndex = match(header$hourly_degree, unique(header$hourly_degree)),
    detriIndex = match(header$daily_detritus, unique(header$daily_detritus))
)# end of data.frame
    
HydPrf = lapply(unique(header$hydraulicTable),function(xx){
    ratingTable = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    ## rateTable assume zone length is 10m
    return <- new('CrossSectionQrelation',
        Q2Vol = approxfun(ratingTable$discharge, ratingTable[,'volume']),
        Vol2Z = approxfun(ratingTable[,'volume'], ratingTable[,'depth']),
        Vol2V = approxfun(ratingTable[,'volume'], ratingTable[,'mvelocity']),
        Vol2WP = approxfun(ratingTable[,'volume'], ratingTable[,'wettedBperimeter']),
        Vol2Q = approxfun(ratingTable[,'volume'], ratingTable$discharge),
        Vol2CA = approxfun(ratingTable[,'volume'], ratingTable[,'crossarea']),
        maxQ = max(ratingTable$discharge),
        maxVol = max(ratingTable[ratingTable[,'bank']==0,'volume'])
        )#list
})# end of lapply
    
    
## reach contains { ts1[t] , ts2[t] }; difficult to access
## reach contains { 1{ts1[1],ts2[1]}, 2{ts[2],ts[2]} }; easy to access
## reach contains { data.frame(row=profiles; col=time) }; can lump all reaches
# -- assume daily input for cq
TSbaseQ = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput),function(xx){
    tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
    baseflow = fivedayblockbaseflow(tmp[match(SimulationTime,tmp$date),'Qlps']*0.001)
    #stormQ = sapply(totalflow_m3s-baseflow,function(x){return <- ifelse(x>0,x,0)})
    return <- rep(baseflow,each=24)
})))# end of lapply
    
TSquickQ = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput),function(xx){
    tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
    totalflow_m3s = tmp[match(SimulationTime,tmp$date),'Qlps']*0.001
    baseflow = fivedayblockbaseflow(totalflow_m3s)
    stormQ = sapply(totalflow_m3s-baseflow,function(x){return <- ifelse(x>0,x,0)})
    return <- rep(stormQ,each=24)
})))# end of lapply
    
TSno3c = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput),function(xx){
    tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
    return <- rep(1000*tmp[match(SimulationTime,tmp$date),'Nmgpl'],each=24)
})))# end of lapply; concentrations mg/m3

TSpo4c = as.data.frame(do.call(rbind,lapply(unique(header$Cqinput),function(xx){
    tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
    return <- rep(1000*tmp[match(SimulationTime,tmp$date),'Pmgpl'],each=24)
})))# end of lapply; concentrations
    
    ## this one is difficult (not transformed yet!)
    tmp2 = lapply(unique(header$daily_detritus),function(xx){
        tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
        tmp$date = as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-'))
        return <- data.frame(
            LabC = tmp$labileC[match(SimulationTime,tmp$date)]*0.04166667,
            IntC = tmp$celluloseC[match(SimulationTime,tmp$date)]*0.04166667,
            RecC = tmp$ligninC[match(SimulationTime,tmp$date)]*0.04166667,
            LabN = tmp$labileN[match(SimulationTime,tmp$date)]*0.04166667,
            IntN = tmp$celluloseN[match(SimulationTime,tmp$date)]*0.04166667,
            RecN = tmp$ligninN[match(SimulationTime,tmp$date)]*0.04166667,
            LabP = tmp$labileP[match(SimulationTime,tmp$date)]*0.04166667,
            IntP = tmp$celluloseP[match(SimulationTime,tmp$date)]*0.04166667,
            RecP = tmp$ligninP[match(SimulationTime,tmp$date)]*0.04166667
            )# end of dataframe
    })# end of lapply

TSdetri = list()
TSdetri$LabC = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$LabC })))
TSdetri$LabN = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$LabN })))
TSdetri$LabP = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$LabP })))
TSdetri$IntC = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$IntC })))
TSdetri$IntN = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$IntN })))
TSdetri$IntP = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$IntP })))
TSdetri$RecC = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$RecC })))
TSdetri$RecN = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$RecN })))
TSdetri$RecP = as.data.frame(do.call(rbind,lapply(seq_along(tmp2),function(i){return <- tmp2[[i]]$RecP })))
rm(tmp2)
    
TSpar = as.data.frame(do.call(rbind,lapply(unique(header$hourly_PAR),function(xx){
    tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    tmp$date = paste(as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-')),tmp$hour,sep='-')
    return <- PAR2PARF(tmp$PAR[match(SimulationTimeHour,tmp$date)])
})))# end of lapply

TSdeg = as.data.frame(do.call(rbind,lapply(unique(header$hourly_degree),function(xx){
    tmp = read.csv(paste(ModelArg$proj, xx ,sep='/'))
    tmp$date = paste(as.Date(paste(tmp$year,tmp$month,tmp$day,sep='-')),tmp$hour,sep='-')
    return <- tmp$temperature[match(SimulationTimeHour,tmp$date)]
})))# end of lapply
    

	
################################## initial model variables ###################################

# lumping all reaches into one table seems a good idea
numZone = dim(Zone)[1]
numZoneIndex = seq_len(numZone)
TMP = data.frame( tmpCond = rep(F,numZone) )
TMP$STSstep = 1; ## yes
TMP$cond = rep(F, numZone)
TMP$h = rep(0, numZone) ## yes
TMP$h2 = rep(0, numZone) ## yes
TMP$h3 = rep(0, numZone) #<- hypoconcNO3
TMP$h4 = rep(0, numZone) #<- hypoconcNH4
TMP$h5 = rep(0, numZone) #<- hypoconcPO4
TMP$concNO3 = rep(0, numZone)
TMP$concNH4 = rep(0, numZone)
TMP$concPO4 = rep(0, numZone)



TMP$h = cumsum(TSbaseQ[Zone$cqIndex,1]*Zone$bqFrac + TSquickQ[Zone$cqIndex,1]*Zone$sqFrac) # ini discharge
TMP$h2 = sapply(numZoneIndex,function(i){
    min(HydPrf[[Zone$hydIndex[i]]]@Q2Vol(2*0.001),
        HydPrf[[Zone$hydIndex[i]]]@maxVol)*Zone$lenAdj[i] # channel vol, not discharge
}) # channel vol at mean discharge 2 L/s
    
ZoneCol = data.frame(q = sapply(numZoneIndex,function(i){
    min(HydPrf[[Zone$hydIndex[i]]]@Q2Vol(TMP$h[i]),
        HydPrf[[Zone$hydIndex[i]]]@maxVol)*Zone$lenAdj[i] # channel vol, not discharge
}))
ZoneCol$no3 = TSno3c[Zone$cqIndex,1]*ZoneCol$q # mgN/m3 * m3
ZoneCol$nh4 = 0
ZoneCol$po4 = TSpo4c[Zone$cqIndex,1]*ZoneCol$q # mgP/m3 * m3
ZoneCol$POC = 0
ZoneCol$PON = 0
ZoneCol$POP = 0
ZoneCol$DOC = 0
ZoneCol$DON = 0
ZoneCol$DOP = 0
      
ZoneHyd = as.data.frame(do.call(rbind,lapply(numZoneIndex,function(i){ c(
    z = HydPrf[[Zone$hydIndex[i]]]@Vol2Z(ZoneCol$q[i]*Zone$lenAdj_1[i]), #1
    v = HydPrf[[Zone$hydIndex[i]]]@Vol2V(ZoneCol$q[i]*Zone$lenAdj_1[i]), #2
    wba = HydPrf[[Zone$hydIndex[i]]]@Vol2WP(ZoneCol$q[i]*Zone$lenAdj_1[i])*Zone$len[i], #3
    q = HydPrf[[Zone$hydIndex[i]]]@Vol2Q(ZoneCol$q[i]*Zone$lenAdj_1[i]), #4 <<---------------- wrong
    A = HydPrf[[Zone$hydIndex[i]]]@Vol2CA(ZoneCol$q[i]*Zone$lenAdj_1[i]), # 5
    maxv = HydPrf[[Zone$hydIndex[i]]]@maxVol*Zone$lenAdj[i] # 7
    ) })))
ZoneHyd$wba_1 = 1/ZoneHyd$wba
#cbind(ZoneHyd$v, ZoneHyd$q, ZoneCol$q, ZoneHyd$q/ZoneCol$q)
  
   
   
   
   
ZoneBen = data.frame(algC = rep(1000*10,max(numZoneIndex)) ) #gC/m2 -> mgC/m2
ZoneBen$algN = ZoneBen$algC*0.5*(algPar$maxNC+algPar$minNC)
ZoneBen$algP = ZoneBen$algC*0.5*(algPar$maxPC+algPar$minPC)
ZoneBen$algDC = 0
ZoneBen$algDN = 0
ZoneBen$algDP = 0
ZoneBen$algEntC = 0
ZoneBen$algEntN = 0
ZoneBen$algEntP = 0
ZoneBen$alghaDmg = 0 ## check
ZoneBen$alghaDmgIndex = 0 ## check
ZoneBen$algMC = 0
ZoneBen$algMN = 0
ZoneBen$algMP = 0
ZoneBen$algCN = ZoneBen$algC/ZoneBen$algN
ZoneBen$algCP = ZoneBen$algC/ZoneBen$algP
ZoneBen$algSelfF = 1
ZoneBen$algNutrF = 0
ZoneBen$algParF = 1
ZoneBen$algUpNF = 1
ZoneBen$algUpPF = 1
ZoneBen$algPUpNO3 = 1
ZoneBen$algPUpNH4 = 1
ZoneBen$algPUpPO4 = 1
ZoneBen$algArea = 0.9 * sapply(numZoneIndex,function(i){HydPrf[[Zone$hydIndex[i]]]@Vol2WP(TMP$h2[i]*Zone$lenAdj_1[i])*Zone$len[i]})
#cbind(ZoneHyd$wba, ZoneBen$algArea)
    
ZoneBen$LabC = TSdetri$LabC[Zone$detriIndex,1]
ZoneBen$LabN = TSdetri$LabN[Zone$detriIndex,1]
ZoneBen$LabP = TSdetri$LabP[Zone$detriIndex,1]
ZoneBen$IntC = TSdetri$IntC[Zone$detriIndex,1]
ZoneBen$IntN = TSdetri$IntN[Zone$detriIndex,1]
ZoneBen$IntP = TSdetri$IntP[Zone$detriIndex,1]
ZoneBen$RecC = TSdetri$RecC[Zone$detriIndex,1]
ZoneBen$RecN = TSdetri$RecN[Zone$detriIndex,1]
ZoneBen$RecP = TSdetri$RecP[Zone$detriIndex,1]
ZoneBen$im = (ZoneBen$LabC+ZoneBen$IntC+ZoneBen$RecC)*0.09
ZoneBen$mi = ZoneBen$im * 0.1
ZoneBen$mi_class2decayRatio = 0
ZoneBen$mi_class3decayRatio = 0
ZoneBen$mi_resNC = 0
ZoneBen$mi_resPC = 0
ZoneBen$mi_grow = 0
ZoneBen$mi_decayC = 0
ZoneBen$mi_decayN = 0
ZoneBen$mi_decayP = 0
ZoneBen$immo_class2decayRatio = 0
ZoneBen$immo_class3decayRatio = 0
ZoneBen$im_resNC = 0
ZoneBen$im_resPC = 0
ZoneBen$Ncoef_im = 0
ZoneBen$Pcoef_im = 0
ZoneBen$im_decay = 0
ZoneBen$netmineralN = 0
ZoneBen$netmineralP = 0
ZoneBen$im_PUpNO3 = 0
ZoneBen$im_PUpNH4 = 0
ZoneBen$im_PUpPO4 = 0
ZoneBen$totMicRespC = 0
ZoneBen$totMicRespN = 0
ZoneBen$totMicRespP = 0
ZoneBen$Pnitri = 0
ZoneBen$PnitriFrac = 0
ZoneBen$Pdenitri = 0
ZoneBen$PdenitriFrac = 0


ZoneStor = data.frame(As = sapply(numZoneIndex,function(i){
    (3 + HydPrf[[Zone$hydIndex[i]]]@Vol2CA(TMP$h2[i]*Zone$lenAdj_1[i])/HydPrf[[Zone$hydIndex[i]]]@Vol2Z(TMP$h2[i]*Zone$lenAdj_1[i])) * (1-exp(-other_param$p_d*1.5))*other_param$p0/other_param$p_d
}))
ZoneStor$q = ZoneStor$As * Zone$len
ZoneStor$no3 = TSno3c[Zone$cqIndex,1]*ZoneStor$q
ZoneStor$nh4 = 0
ZoneStor$po4 = TSpo4c[Zone$cqIndex,1]*ZoneStor$q
ZoneStor$DOC = 0
ZoneStor$DON = 0
ZoneStor$DOP = 0
ZoneStor$baseNO3 = ZoneStor$no3
ZoneStor$baseNH4 = ZoneStor$nh4
ZoneStor$baseNO3balance = 0
ZoneStor$baseNH4balance = 0
#cbind(ZoneStor$q, ZoneHyd$maxv, ZoneStor$q/ZoneHyd$maxv, ZoneCol$q)
    
ZoneInput = ZoneCol; ZoneInput[]=0
ZoneExport = ZoneCol; ZoneExport[]=0
ZoneBank = ZoneCol; ZoneBank[]=0

system.time({
    outputFile = ifelse(grepl('/',ModelArg$output,fixed=T),ModelArg$output, paste(ModelArg$proj,ModelArg$output,sep='/'))
    outputFile_buff <- file(outputFile,'w')
    cat( c(outputTitle, paste('no3out',1:numZone,sep='_')), '\n', file=outputFile_buff,sep=',')
    
	for(BTS in 1:timemax_BTS ){	 #1:timemax_BTS, 8760
		
		# debug
		#if(sum(ZoneBen$LabC<0)>0){print(BTS);break}
		
        # set STS runs
        BTSDAY = floor((BTS-1)/24)+1
		TMP$h = ZoneHyd$q/ZoneCol$q;
		TMP$STSstep = max(1,apply(TMP$h %o% STSstepOption<0.25,2,prod)* STSstepOption)
		
		
        # inital hourly processes
        algalQ10f = (algPar$Q10^(0.1*TSdeg[Zone$degIndex,BTS]-1.2))
		bioQ10f = (micPar$Q10^(0.1*TSdeg[Zone$degIndex,BTS]-1.2))
		
        # algal death (biological with some velocity related)
        TMP$h = algalQ10f * algPar$moralityRate * timemax_STS
		ZoneBen$algDC <- ZoneBen$algC * TMP$h[]
		ZoneBen$algDN <- ZoneBen$algN * TMP$h[]
		ZoneBen$algDP <- ZoneBen$algP * TMP$h[]
       
        ## wash out;
        TMP$h[] =  1 - (1-(1-micPar$v1ms_left)/(1+exp(-(ZoneHyd$v-0.4)/0.06)))^(1/24) # removal rate
        TMP$h[ZoneHyd$v<0.1]=0
		ZoneBen$algEntC <- ZoneBen$algC * TMP$h[]
		ZoneBen$algEntN <- ZoneBen$algN * TMP$h[]
		ZoneBen$algEntP <- ZoneBen$algP * TMP$h[]
		
        # algal recolonization
        if(BTS%%24==0){
            # new dmg
            TMP$cond = ZoneBen$alghaDmgIndex < TMP$h
            ZoneBen$alghaDmg[TMP$cond] = 1-(1-TMP$h[TMP$cond])^24
            
            # new dmg apply to old dmg
            TMP$cond = TMP$cond & (ZoneBen$alghaDmg > ZoneBen$alghaDmgIndex)
            ZoneBen$alghaDmgIndex[TMP$cond] = 1
            
            # recover
            TMP$cond = ZoneBen$alghaDmgIndex>0
            ZoneBen$alghaDmgIndex[TMP$cond] = ZoneBen$alghaDmgIndex[TMP$cond] * ZoneBen$alghaDmg[TMP$cond] #
            
            # recover reset; max is 0.09148242
            TMP$cond = ZoneBen$alghaDmgIndex<1e-3
            ZoneBen$alghaDmgIndex[TMP$cond]=0
            ZoneBen$alghaDmg[TMP$cond]=0
        }# daily
		
		## mialization
		TMP$h[]  = algPar$mialRate * algalQ10f * timemax_STS
		ZoneBen[,'algMC'] <- ZoneBen[,'algC'] * TMP$h[]
		ZoneBen[,'algMN'] <- ZoneBen[,'algN'] * TMP$h[]
		ZoneBen[,'algMP'] <- ZoneBen[,'algP'] * TMP$h[]
		
		ZoneBen[,'algC'] = ZoneBen[,'algC'] - ZoneBen[,'algDC'] - ZoneBen[,'algEntC'] - ZoneBen[,'algMC']
		ZoneBen[,'algN'] = ZoneBen[,'algN'] - ZoneBen[,'algDN'] - ZoneBen[,'algEntN'] - ZoneBen[,'algMN']
		ZoneBen[,'algP'] = ZoneBen[,'algP'] - ZoneBen[,'algDP'] - ZoneBen[,'algEntP'] - ZoneBen[,'algMP']
		
        ##----------------------------------------------------------------------------------------------------------##
		## microbes death, potential decay and im; ZoneBen[,'colonizedBenthicArea'], ZoneHyd$wba
		## logistic growth break down rP*(1-P/K) = rP - r(P/K)[P] --> here  remaining % = 1-r(P/K)
        TMP$h2[] = micPar$growthRate_im*10*(ZoneBen$im+ ZoneBen$mi)/(ZoneBen$LabC+ZoneBen$IntC+ZoneBen$RecC)*bioQ10f*timemax_STS
        TMP$h2[TMP$h2>1]=0.999
        TMP$h[] = TMP$h2 * ZoneBen[,'im']
        TMP$h4[] = TMP$h2 * ZoneBen[,'im']*micPar$im_nc
        TMP$h5[] = TMP$h2 * ZoneBen[,'im']*micPar$im_pc
        ZoneBen[,'im'] = ZoneBen[,'im'] * (1-TMP$h2)
        
        TMP$h2[] = micPar$growthRate_mi*10*(ZoneBen$im+ ZoneBen$mi)/(ZoneBen$LabC+ZoneBen$IntC+ZoneBen$RecC)*bioQ10f*timemax_STS
        TMP$h2[TMP$h2>1]=0.999
        TMP$h[] = TMP$h + TMP$h2 * ZoneBen[,'mi']
      	TMP$h4[] = TMP$h4[] + TMP$h2[] * ZoneBen[,'mi']*micPar$mi_nc
        TMP$h5[] = TMP$h5[] + TMP$h2[] * ZoneBen[,'mi']*micPar$mi_pc
        ZoneBen[,'mi'] = ZoneBen[,'mi'] * (1-TMP$h2)
        
        TMP$h3[] = ZoneBen[,'mi']*micPar$baseR*bioQ10f*timemax_STS
        ZoneBen[,'totMicRespC'] = ZoneBen[,'im']*micPar$baseR*bioQ10f*timemax_STS
        ZoneBen[,'totMicRespN'] = ZoneBen[,'totMicRespC']*micPar$im_nc + TMP$h3*micPar$mi_nc
        ZoneBen[,'totMicRespP'] = ZoneBen[,'totMicRespC']*micPar$im_pc + TMP$h3*micPar$mi_pc
        ZoneBen[,'totMicRespC'] = ZoneBen[,'totMicRespC'] + TMP$h3[]
        ZoneBen[,'im'] = ZoneBen[,'im'] * (1-micPar$baseR*bioQ10f*timemax_STS) # basal respiration
        ZoneBen[,'mi'] = ZoneBen[,'mi'] * (1-micPar$baseR*bioQ10f*timemax_STS) #
       
        
        ## -- detritus inputs (1/24 = 0.04166667)
		ZoneBen[,'LabC'] = ZoneBen[,'LabC'] + TSdetri$LabC[Zone$detriIndex,BTSDAY] +
            TMP$h*0.5 + ZoneBen[,'algDC']*0.5 # mgC/m2 unit area;
            # hold is dead microbes at this point
		ZoneBen[,'LabN'] = ZoneBen[,'LabN'] + TSdetri$LabN[Zone$detriIndex,BTSDAY] +
			TMP$h4*0.5 + ZoneBen[,'algDN']*0.5 # mgN/m2 unit area
		ZoneBen[,'LabP'] = ZoneBen[,'LabP'] + TSdetri$LabN[Zone$detriIndex,BTSDAY] +
			TMP$h5*0.5 + ZoneBen[,'algDP']*0.5 # mgN/m2 unit area
            
        ZoneBen[,'IntC'] = ZoneBen[,'IntC'] + TSdetri$IntC[Zone$detriIndex,BTSDAY] +
            TMP$h*0.5 + ZoneBen[,'algDC']*0.5 # mgC/m2 unit area;
            # hold is dead microbes at this point
        ZoneBen[,'IntN'] = ZoneBen[,'IntN'] + TSdetri$IntN[Zone$detriIndex,BTSDAY] +
            TMP$h4*0.5 + ZoneBen[,'algDN']*0.5 # mgN/m2 unit area
        ZoneBen[,'IntP'] = ZoneBen[,'IntP'] + TSdetri$IntP[Zone$detriIndex,BTSDAY] +
            TMP$h5*0.5 + ZoneBen[,'algDP']*0.5 # mgN/m2 unit area
            
        ZoneBen[,'RecC'] = ZoneBen[,'RecC'] + TSdetri$IntC[Zone$detriIndex,BTSDAY]
        ZoneBen[,'RecN'] = ZoneBen[,'RecN'] + TSdetri$IntN[Zone$detriIndex,BTSDAY]
        ZoneBen[,'RecP'] = ZoneBen[,'RecP'] + TSdetri$IntP[Zone$detriIndex,BTSDAY]
	

		## micorbial potential growth -- mi
        ZoneBen$mi_class2decayRatio = (1000+ZoneBen$LabC)/(10000+ZoneBen$IntC);
        ZoneBen$mi_class3decayRatio = (1000+ZoneBen$LabC)/(100000+ZoneBen$RecC)
		TMP$h[] = (ZoneBen$LabC+ZoneBen$IntC*ZoneBen$mi_class2decayRatio +ZoneBen$RecC*ZoneBen$mi_class3decayRatio)
        ZoneBen[,'mi_resNC'] = (ZoneBen[,'LabN']+ZoneBen[,'IntN']*ZoneBen$mi_class2decayRatio +ZoneBen[,'RecN']*ZoneBen$mi_class3decayRatio) / TMP$h[]
        ZoneBen[,'mi_resPC'] = (ZoneBen[,'LabP']+ZoneBen[,'IntP']*ZoneBen$mi_class2decayRatio +ZoneBen[,'RecP']*ZoneBen$mi_class3decayRatio) / TMP$h[]
        ZoneBen[,'mi_grow'] = ZoneBen[,'mi']*micPar$growthRate_mi * bioQ10f * timemax_STS;
        ZoneBen[,'mi_decayC'] = vectorMin(
            TMP$h[], # availC
            vectorMax(
                ZoneBen[,'mi_grow'], # bio growth
                vectorMax(
                    micPar$gcoef*ZoneBen[,'mi_grow']*micPar$mi_nc/ZoneBen[,'mi_resNC'],
                    micPar$gcoef*ZoneBen[,'mi_grow']*micPar$mi_pc/ZoneBen[,'mi_resPC']
                    # nutrient resources
                ))); # actual decay
            TMP$cond[] = ZoneBen[,'mi_grow'] > ZoneBen[,'mi_decayC']
            if(sum(TMP$cond)>0){
                ZoneBen[TMP$cond,'mi_grow'] = vectorMin(
                    ZoneBen[TMP$cond,'mi_decayC']*ZoneBen[TMP$cond,'mi_resNC']/micPar$mi_nc,
                    ZoneBen[TMP$cond,'mi_decayC']*ZoneBen[TMP$cond,'mi_resPC']/micPar$mi_pc
            )}#if # use actual decay to adjust bio growth
        ZoneBen[,'mi_grow'] = ZoneBen[,'mi_grow']*micPar$gcoef
        ZoneBen[,'mi_decayN'] = ZoneBen[,'mi_decayC']*ZoneBen[,'mi_resNC'] - ZoneBen[,'mi_grow']*micPar$mi_nc
        ZoneBen[,'mi_decayP'] = ZoneBen[,'mi_decayC']*ZoneBen[,'mi_resPC'] - ZoneBen[,'mi_grow']*micPar$mi_pc
        ZoneBen[ZoneBen[,'mi_decayN']<0,'mi_decayN']=0
        ZoneBen[ZoneBen[,'mi_decayP']<0,'mi_decayP']=0
        
        
        ZoneBen$immo_class2decayRatio = (1000+ZoneBen$LabC)/(10000+ZoneBen$IntC);
        ZoneBen$immo_class3decayRatio = (1000+ZoneBen$LabC)/(1000000+ZoneBen$RecC)
        TMP$h[] = (ZoneBen$LabC+ZoneBen$IntC*ZoneBen$immo_class2decayRatio +ZoneBen$RecC*ZoneBen$immo_class3decayRatio)
        ZoneBen[,'im_resNC'] = (ZoneBen[,'LabN']+ZoneBen[,'IntN']*ZoneBen$immo_class2decayRatio +ZoneBen[,'RecN']*ZoneBen$immo_class3decayRatio) / TMP$h[]
        ZoneBen[,'im_resPC'] = (ZoneBen[,'LabP']+ZoneBen[,'IntP']*ZoneBen$immo_class2decayRatio +ZoneBen[,'RecP']*ZoneBen$immo_class3decayRatio) / TMP$h[]
        ZoneBen[,'im_decay'] = ZoneBen[,'im']*micPar$growthRate_im * bioQ10f * timemax_STS;
		ZoneBen[,'netmineralN'] = ZoneBen[,'im_decay'] * (micPar$gcoef/micPar$im_cn - ZoneBen[,'im_resNC'])
		ZoneBen[,'netmineralP'] = ZoneBen[,'im_decay'] * (micPar$gcoef/micPar$im_cp - ZoneBen[,'im_resPC'])
            

        ##-------------------------------------------------------------------------------------------##
		STS = 0; while(STS < timemax_STS){ #1:timemax_STS
                
                ## export
				TMP$cond = TSdeg[Zone$degIndex,BTS]>0 ## freezing issue
                TMP$h[] = ZoneHyd$q/ZoneCol$q;
                TMP$h[!TMP$cond] = 0
                ZoneExport = ZoneCol[] * (TMP$h * TMP$STSstep)
                #STSaggregateVar = STSaggregateVar + ZoneExport
                ZoneCol[] = ZoneCol[] * (1-TMP$h * TMP$STSstep)
					
									
				## lateral input (per second)
                TMP$h[] = TSquickQ[Zone$cqIndex,BTS]*Zone$sqFrac
                ZoneInput$q = TSbaseQ[Zone$cqIndex,BTS]*Zone$bqFrac*TMP$STSstep
                ZoneInput$no3 = TSno3c[Zone$cqIndex,BTS]*(TMP$h+ZoneInput$q)
                ZoneInput$po4 = TSpo4c[Zone$cqIndex,BTS]*(TMP$h+ZoneInput$q)
                ZoneInput$q = ZoneInput$q + TMP$h
                
                
                ## upstream input
                for(ii in numZoneIndex){
                    ZoneInput$q[ii] = ZoneInput$q[ii] + sum(ZoneExport$q[ZoneConn[[ii]]])
                    ZoneInput$no3[ii] = ZoneInput$no3[ii] + sum(ZoneExport$no3[ZoneConn[[ii]]])
                    ZoneInput$nh4[ii] = ZoneInput$nh4[ii] + sum(ZoneExport$nh4[ZoneConn[[ii]]])
                    ZoneInput$po4[ii] = ZoneInput$po4[ii] + sum(ZoneExport$po4[ZoneConn[[ii]]])
                    ZoneInput$POC[ii] = ZoneInput$POC[ii] + sum(ZoneExport$POC[ZoneConn[[ii]]])
                    ZoneInput$PON[ii] = ZoneInput$PON[ii] + sum(ZoneExport$PON[ZoneConn[[ii]]])
                    ZoneInput$POP[ii] = ZoneInput$POP[ii] + sum(ZoneExport$POP[ZoneConn[[ii]]])
                    ZoneInput$DOC[ii] = ZoneInput$DOC[ii] + sum(ZoneExport$DOC[ZoneConn[[ii]]])
                    ZoneInput$DON[ii] = ZoneInput$DON[ii] + sum(ZoneExport$DON[ZoneConn[[ii]]])
                    ZoneInput$DOP[ii] = ZoneInput$DOP[ii] + sum(ZoneExport$DOP[ZoneConn[[ii]]])
                }# end of for ii
                
			  	## (basal resp.) mial input (from hourly to seconds) # stop here
				ZoneInput$nh4 = ZoneInput$nh4 + (ZoneBen$algMN + ZoneBen$totMicRespN) * ZoneHyd$wba * 0.0002777778 * TMP$STSstep #  ZoneBen$mial not sure
				ZoneInput$po4 = ZoneInput$po4 + (ZoneBen$algMP + ZoneBen$totMicRespP) * ZoneHyd$wba * 0.0002777778 * TMP$STSstep
				
                
                
                ## entraintment & deposition
                TMP$h3[] = ((1-micPar$v1ms_left)/(1+exp(-(ZoneHyd$v-0.4)/0.06)))^(TMP$STSstep*1.157407e-05)
                TMP$h[] =  1 - (1-(1-micPar$v1ms_left)/(1+exp(-(ZoneHyd$v-0.4)/0.06)))^(TMP$STSstep*1.157407e-05) # removal rate
                TMP$h[] = TMP$h[] * 10*(ZoneBen$im+ ZoneBen$mi)/(ZoneBen$LabC+ZoneBen$IntC+ZoneBen$RecC) # count for micorbial processing to fragmentation (new June 17)
                TMP$h2[] = ZoneInput[,'POC']*(1-TMP$h3) * ZoneHyd$wba_1 / (ZoneBen$LabC+ZoneBen$IntC+ZoneBen$RecC) #
				ZoneInput[,'POC'] =
                    ZoneInput[,'POC']* TMP$h3 +
                    (ZoneBen$LabC+ZoneBen$IntC+ZoneBen$RecC) * TMP$h * ZoneHyd$wba +
                    ZoneBen$algEntC*0.0002777778 * TMP$STSstep *  ZoneHyd$wba ## pre-calculated outside STS loop
                ZoneBen$LabC = ZoneBen$LabC * (1-TMP$h) * (1+TMP$h2)
                ZoneBen$IntC = ZoneBen$IntC * (1-TMP$h) * (1+TMP$h2)
                ZoneBen$RecC = ZoneBen$RecC * (1-TMP$h) * (1+TMP$h2)
				
				TMP$h2[] = ZoneInput[,'PON']*(1-TMP$h3) * ZoneHyd$wba_1 / (ZoneBen$LabN+ZoneBen$IntN+ZoneBen$RecN)
				ZoneInput[,'PON'] =
                    ZoneInput[,'PON'] * TMP$h3 +
                    (ZoneBen$LabN+ZoneBen$IntN+ZoneBen$RecN)*TMP$h * ZoneHyd$wba +
					ZoneBen[,'algEntN']*0.0002777778 * TMP$STSstep *  ZoneHyd$wba
                ZoneBen$LabN = ZoneBen$LabN * (1-TMP$h) * (1+TMP$h2)
                ZoneBen$IntN = ZoneBen$IntN * (1-TMP$h) * (1+TMP$h2)
                ZoneBen$RecN = ZoneBen$RecN * (1-TMP$h) * (1+TMP$h2)
				
				TMP$h2[] = ZoneInput[,'POP']*(1-TMP$h3) * ZoneHyd$wba_1 / (ZoneBen$LabP+ZoneBen$IntP+ZoneBen$RecP)
				ZoneInput[,'POP'] =
					 ZoneInput[,'POP'] * TMP$h3 +
					 (ZoneBen$LabP+ZoneBen$IntP+ZoneBen$RecP)*TMP$h * ZoneHyd$wba +
					 ZoneBen[,'algEntP']*0.0002777778 * TMP$STSstep *  ZoneHyd$wba
				ZoneBen$LabP = ZoneBen$LabP * (1-TMP$h) * (1+TMP$h2)
                ZoneBen$IntP = ZoneBen$IntP * (1-TMP$h) * (1+TMP$h2)
                ZoneBen$RecP = ZoneBen$RecP * (1-TMP$h) * (1+TMP$h2)
			
				ZoneCol[] = ZoneCol[] + ZoneInput
					
                    
					
				## extract to riparian
                ## positive = to riparian; negative = potential from riparian
				TMP$h[ ] = ZoneCol$q -  ZoneHyd$maxv
				TMP$cond[ ] = TMP$h>0 # is bounded by numZone
				if(sum(TMP$cond)>0){
					## go to riparian
					TMP$h2[ ] = 0 # is bounded by numZone
					TMP$h2[ ] = TMP$h[ ]/ZoneCol$q
					ZoneBank[TMP$cond,'q'] = ZoneBank[TMP$cond,'q'] + TMP$h[TMP$cond]
					ZoneBank[TMP$cond,'no3'] = ZoneBank[TMP$cond,'no3'] + TMP$h2[TMP$cond]* ZoneCol[TMP$cond,'no3']
					ZoneBank[TMP$cond,'nh4'] = ZoneBank[TMP$cond,'nh4'] + TMP$h2[TMP$cond]* ZoneCol[TMP$cond,'nh4']
					ZoneBank[TMP$cond,'po4'] = ZoneBank[TMP$cond,'po4'] + TMP$h2[TMP$cond]* ZoneCol[TMP$cond,'po4']
					ZoneBank[TMP$cond,'POC'] = ZoneBank[TMP$cond,'POC'] + TMP$h2[TMP$cond]* ZoneCol[TMP$cond,'POC']
					ZoneBank[TMP$cond,'PON'] = ZoneBank[TMP$cond,'PON'] + TMP$h2[TMP$cond]* ZoneCol[TMP$cond,'PON']
					ZoneBank[TMP$cond,'POP'] = ZoneBank[TMP$cond,'POP'] + TMP$h2[TMP$cond]* ZoneCol[TMP$cond,'POP']
					
					TMP$h2 = 1-TMP$h2 # is bounded by numZone
					ZoneCol[TMP$cond,'q'] = ZoneCol[TMP$cond,'q'] * TMP$h2[TMP$cond]
					ZoneCol[TMP$cond,'no3'] = ZoneCol[TMP$cond,'no3'] * TMP$h2[TMP$cond]
					ZoneCol[TMP$cond,'nh4'] = ZoneCol[TMP$cond,'nh4'] * TMP$h2[TMP$cond]
					ZoneCol[TMP$cond,'po4'] = ZoneCol[TMP$cond,'po4'] * TMP$h2[TMP$cond]
					ZoneCol[TMP$cond,'POC'] = ZoneCol[TMP$cond,'POC'] * TMP$h2[TMP$cond]
					ZoneCol[TMP$cond,'PON'] = ZoneCol[TMP$cond,'PON'] * TMP$h2[TMP$cond]
					ZoneCol[TMP$cond,'POP'] = ZoneCol[TMP$cond,'POP'] * TMP$h2[TMP$cond]
					
				}# to riparian
					
				## input from riparian (per STS) 
				TMP$cond[ ] = TMP$h<0 & ZoneBank$q>0 # is bounded by numZone
				if(sum(TMP$cond)>0){
					# get back from riparian
					TMP$h2[ ] = -TMP$h/ZoneBank$q # is bounded by numZone
					TMP$h2[is.na(TMP$h2)|is.infinite(TMP$h2)] = 0;
					TMP$h2[TMP$h2>1] = 1
					TMP$h2[TMP$h2<0] = 0
					
					ZoneCol[TMP$cond,'q'] = ZoneCol[TMP$cond,'q'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'q']
					ZoneCol[TMP$cond,'no3'] = ZoneCol[TMP$cond,'no3'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'no3']
					ZoneCol[TMP$cond,'nh4'] = ZoneCol[TMP$cond,'nh4'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'nh4']
					ZoneCol[TMP$cond,'po4'] = ZoneCol[TMP$cond,'po4'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'po4']
					ZoneCol[TMP$cond,'POC'] = ZoneCol[TMP$cond,'POC'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'POC']
					ZoneCol[TMP$cond,'PON'] = ZoneCol[TMP$cond,'PON'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'PON']
					ZoneCol[TMP$cond,'POP'] = ZoneCol[TMP$cond,'POP'] + TMP$h2[TMP$cond]*ZoneBank[TMP$cond,'POP']
					TMP$h2 = 1-TMP$h2
					ZoneBank[TMP$cond,'q'] = ZoneBank[TMP$cond,'q']*TMP$h2[TMP$cond]
					ZoneBank[TMP$cond,'no3'] = ZoneBank[TMP$cond,'no3']*TMP$h2[TMP$cond]
					ZoneBank[TMP$cond,'nh4'] = ZoneBank[TMP$cond,'nh4']*TMP$h2[TMP$cond]
					ZoneBank[TMP$cond,'po4'] = ZoneBank[TMP$cond,'po4']*TMP$h2[TMP$cond]
					ZoneBank[TMP$cond,'POC'] = ZoneBank[TMP$cond,'POC']*TMP$h2[TMP$cond]
					ZoneBank[TMP$cond,'PON'] = ZoneBank[TMP$cond,'PON']*TMP$h2[TMP$cond]
					ZoneBank[TMP$cond,'POP'] = ZoneBank[TMP$cond,'POP']*TMP$h2[TMP$cond]
				}# to riparian
	
            
				###################################################################################### potential uptake
                # POC, PON, POP are not dissolved
                
                ##-------------------------------------------
                # ZoneCol
                TMP$cond[] = ZoneCol$q<=0
                TMP$concNO3[] = ZoneCol[,'no3']/ZoneCol$q; TMP$concNO3[is.na(TMP$concNO3)|is.infinite(TMP$concNO3)]=0; TMP$concNO3[TMP$cond]=0
                TMP$concNH4[] = ZoneCol[,'nh4']/ZoneCol$q; TMP$concNH4[is.na(TMP$concNH4)|is.infinite(TMP$concNH4)]=0; TMP$concNH4[TMP$cond]=0
                TMP$concPO4[] = ZoneCol[,'po4']/ZoneCol$q; TMP$concPO4[is.na(TMP$concPO4)|is.infinite(TMP$concPO4)]=0; TMP$concPO4[TMP$cond]=0
                
                
                    ## ------- constructing hyporheic
                    # ZoneStor (vol and mg per unit area)
                    # ZoneStor$q = zoneEnv_matrix[,BTS,'baseQ'] ## what is it? "m3/m2"
                    ZoneStor[,'no3'] = ZoneStor[,'baseNO3'] + ZoneStor[,'baseNO3balance']
                    ZoneStor[,'nh4'] = ZoneStor[,'baseNH4'] + ZoneStor[,'baseNH4balance']
                        ZoneStor[,'baseNO3'] = ZoneStor[,'no3']
                        ZoneStor[,'baseNH4'] = ZoneStor[,'nh4']
                    #ZoneStor[,'po4'] = ZoneStor$q* TMP$concPO4 ##<<----- come back here
                
                
                TMP$h3 = ZoneStor[,'no3']/ZoneStor$q
                TMP$h4 = ZoneStor[,'nh4']/ZoneStor$q
                TMP$h5 = ZoneStor[,'po4']/ZoneStor$q
                TMP$cond[] = TMP$concNO3 > TMP$h3;
				TMP$h2[TMP$cond] = TMP$concNO3[TMP$cond];
                TMP$h2[!TMP$cond]=TMP$h3[!TMP$cond]
					

                ## potential uptake
                # what is *4e-05?
				ZoneBen$algPUpNO3[] =  algPar$MulhUpA*(TMP$concNO3 ^(1+algPar$MulhUpB)) * ZoneBen[,'algUpNF'] * algalQ10f * TMP$STSstep *4e-05 # correct by 25gC/m2
                TMP$cond[] = ZoneBen$algPUpNO3>ZoneCol[,'no3']
                ZoneBen$algPUpNO3[TMP$cond] = ZoneCol[TMP$cond,'no3']
                # mgN/mgC/s * mgC/m2 * s  * m2 = mgN; algUpNF = mgC/m2 * m2 * fraction
				       
				ZoneBen$algPUpNH4[] = algPar$MulhUpA*((TMP$concNH4+TMP$concNO3)^(1+ algPar$MulhUpB)) *ZoneBen[,'algUpNF']*algalQ10f * TMP$STSstep *4e-05
				ZoneBen$algPUpNH4[] = ZoneBen$algPUpNH4 - ZoneBen$algPUpNO3;
				ZoneBen$algPUpNH4[ZoneBen$algPUpNH4<0] = 0
					TMP$cond[] = ZoneBen$algPUpNH4>ZoneCol[,'nh4']
					ZoneBen$algPUpNH4[TMP$cond] = ZoneCol[TMP$cond,'nh4']
				
				ZoneBen$algPUpPO4[] = algPar$umaxP*TMP$concPO4/(algPar$Phalf+TMP$concPO4)*ZoneBen[,'algUpPF']*algalQ10f * TMP$STSstep # mgP
					TMP$cond[] = ZoneBen$algPUpPO4>ZoneCol[,'po4']
					ZoneBen$algPUpPO4[TMP$cond] = ZoneCol[TMP$cond,'po4']
				
				## nitrification
				ZoneBen$Pnitri[] = micPar$max_nitrif * bioQ10f * 0.5*(ZoneCol[,'nh4']+ZoneStor[,'nh4']) * TMP$STSstep
				ZoneBen$PnitriFrac = (0.5*ZoneBen$Pnitri - ZoneStor[,'nh4'])/ZoneBen$Pnitri
				ZoneBen$PnitriFrac[ZoneBen$PnitriFrac<0] = 0
				ZoneBen$PnitriFrac = ZoneBen$PnitriFrac + 0.5
				 
				 
				## immobilization --  mgN/m2/s* m2 * s = mgN
				ZoneBen$im_PUpNO3[] = micPar$MulhUpA*(TMP$concNO3 ^(1+micPar$MulhUpB)) * ZoneHyd$wba*bioQ10f*TMP$STSstep *ZoneBen[,'im'] *4e-05;
					TMP$cond[] = ZoneBen$im_PUpNO3>ZoneCol[,'no3']
					ZoneBen$im_PUpNO3[TMP$cond] = ZoneCol[TMP$cond,'no3']
				ZoneBen$im_PUpNH4[] = micPar$MulhUpA*((TMP$concNH4+TMP$concNO3)^(1+ micPar$MulhUpB)) * ZoneHyd$wba*bioQ10f*TMP$STSstep * ZoneBen[,'im'] *4e-05 - ZoneBen$im_PUpNO3;
				ZoneBen$im_PUpNH4[ZoneBen$im_PUpNH4 <0] = 0
					TMP$cond[] = ZoneBen$im_PUpNH4>ZoneCol[,'nh4']
					ZoneBen$im_PUpNH4[TMP$cond] = ZoneCol[TMP$cond,'nh4']
				ZoneBen$im_PUpPO4[] = algPar$umaxP*TMP$concPO4/(algPar$Phalf+TMP$concPO4) *ZoneHyd$wba*bioQ10f*TMP$STSstep *ZoneBen[,'im'] ## no need for 4e-05!!
                	TMP$cond[] = ZoneBen$im_PUpPO4>ZoneCol[,'po4']
					ZoneBen$im_PUpPO4[TMP$cond] = ZoneCol[TMP$cond,'po4']
                
                
                ## make sure microbes are not taking more than potential
                ## 0.5/cn - 1/food_CN -> uptake [inside STS loop]: wN is cumsum in STS loop
                TMP$cond[] = ZoneBen[,'netmineralN']>0
                if(sum(TMP$cond)>0){
                    TMP$h[TMP$cond] = (ZoneBen[TMP$cond,'netmineral'] - ZoneBen[TMP$cond,'wN']) * ZoneHyd[TMP$cond,'wba']
                    # unit neededN vs uptakenN
                    ZoneBen$im_PUpNO3[TMP$cond] = sapply(vectorMin( ZoneBen$im_PUpNO3[TMP$cond], TMP$h[TMP$cond]),function(x){ifelse(x>0,x,0.0)})
                    ZoneBen$im_PUpNH4[TMP$cond] = sapply(vectorMin( TMP$h[TMP$cond]-ZoneBen$im_PUpNO3[TMP$cond], ZoneBen$im_PUpNH4[TMP$cond]),function(x){ifelse(x>0,x,0.0)})
                }#if
                TMP$cond[] = ZoneBen[,'netmineralP']>0 # 0.5/cn - 1/food_CN -> uptake [inside STS loop]:
                if(sum(TMP$cond)>0){
                    TMP$h[TMP$cond] = (ZoneBen[TMP$cond,'netmineralP'] - ZoneBen[TMP$cond,'wP']) * ZoneHyd[TMP$cond,'wba']
                    ZoneBen$im_PUpPO4[TMP$cond] = sapply(vectorMin( ZoneBen$im_PUpPO4[TMP$cond], TMP$h[TMP$cond]),function(x){ifelse(x>0,x,0.0)})
                }#if
				 
                 
                ## denitrification (new June 17)
                TMP$h[] = ZoneBen[,'totMicRespC'] + ZoneBen[,'im_decay']*(1-micPar$gcoef) + ZoneBen[,'mi_decayC']-ZoneBen[,'mi_grow'];
                ZoneBen$Pdenitri[] = ZoneHyd$wba * micPar$MulhDenA *(TMP$h2^(1+micPar$MulhDenB)) * bioQ10f * TMP$STSstep
                ZoneBen$Pdenitri[TMP$h<=0] = 0 # no decay activity
                ZoneBen$PdenitriFrac[] = (ZoneBen$Pdenitri - ZoneStor[,'no3'])/ZoneBen$Pdenitri;
                ZoneBen$PdenitriFrac[is.na(ZoneBen$PdenitriFrac)]=0
                ZoneBen$PdenitriFrac[ZoneBen$PdenitriFrac<0] = 0 # for water column
				 
                ## stop here
				## ------ splitting resources to processes
				TMP$cond[] = ZoneBen[,'netmineral']>0
				sumPotentialNH4[] = ZoneBen$algPUpNH4 + ZoneBen$Pnitri*ZoneBen$PnitriFrac
                sumPotentialNH4[TMP$cond] = sumPotentialNH4[TMP$cond] + ZoneBen$im_PUpNH4[TMP$cond]
				NH4limit_scalar[] = 1; 
					TMP$cond =(ZoneCol[,'nh4'] < sumPotentialNH4);
					NH4limit_scalar[TMP$cond] = ZoneCol[TMP$cond,'nh4']/sumPotentialNH4[TMP$cond];
					
				sumPotentialNO3[] = ZoneBen$algPUpNO3 + ZoneBen$Pdenitri*ZoneBen$PdenitriFrac
				sumPotentialNO3[TMP$cond] = sumPotentialNO3[TMP$cond] + ZoneBen$im_PUpNO3[TMP$cond]
				NO3limit_scalar[] = 1; 
					TMP$cond =(ZoneCol[,'no3'] < sumPotentialNO3);
					NO3limit_scalar[TMP$cond] = ZoneCol[TMP$cond,'no3']/sumPotentialNO3[TMP$cond];
				
				TMP$cond[] = ZoneBen[,'netmineralP']>0
				sumPotentialPO4[] = ZoneBen$algPUpPO4
				sumPotentialPO4[TMP$cond] = sumPotentialPO4[TMP$cond] + ZoneBen$im_PUpPO4[TMP$cond]
				Plimit_scalar[] = 1; 
					TMP$cond =(ZoneCol[,'po4'] < sumPotentialPO4);
					Plimit_scalar[TMP$cond] = ZoneCol[TMP$cond,'po4']/sumPotentialPO4[TMP$cond];
				
                ## remove uptaken solute from water column
                sumPotentialNO3[] = sumPotentialNO3 * NO3limit_scalar; 
                	sumPotentialNO3[] = sumPotentialNO3/ZoneCol[,'no3'];
                	sumPotentialNO3[is.na(sumPotentialNO3)]=0; 
                	sumPotentialNO3[sumPotentialNO3>1]=1
                sumPotentialNH4[] = sumPotentialNH4 * NH4limit_scalar; ##<<------------------------
                	sumPotentialNH4[] = sumPotentialNH4/ZoneCol[,'nh4'];
                	sumPotentialNH4[is.na(sumPotentialNH4)]=0; 
                	sumPotentialNH4[sumPotentialNH4>1]=1
                sumPotentialPO4[] = sumPotentialPO4 * Plimit_scalar; 
                	sumPotentialPO4[] = sumPotentialPO4/ZoneCol[,'po4'];
                	sumPotentialPO4[is.na(sumPotentialPO4)]=0; 
                	sumPotentialPO4[sumPotentialPO4>1]=1
                    
                    ##
                    ZoneCol[,'no3'] = ZoneCol[,'no3'] * (1-sumPotentialNO3) + ZoneBen$PnitriFrac* ZoneBen$Pnitri
                    ZoneCol[,'nh4'] = ZoneCol[,'nh4'] * (1-sumPotentialNH4)
                    ZoneCol[,'po4'] = ZoneCol[,'po4'] * (1-sumPotentialPO4)
                    ZoneStor[,'no3'] = ZoneStor[,'no3'] -
                        ZoneBen$Pdenitri*(1.0-ZoneBen$PdenitriFrac) +  # *ZoneBen[,'colonizedBenthicArea_1']
                        ZoneBen$Pnitri*(1.0-ZoneBen$PnitriFrac); #*ZoneBen[,'colonizedBenthicArea_1'];
                        ZoneStor[ZoneStor[,'no3']<0,'no3']=0
                        ZoneStor[,'nh4'] = ZoneStor[,'nh4'] - ZoneBen$Pnitri*(1.0-ZoneBen$PnitriFrac); #*ZoneBen[,'colonizedBenthicArea_1'];
                        ZoneStor[ZoneStor[,'nh4']<0,'nh4']=0
                    
                    ZoneBen$algPUpNO3[] = ZoneBen$algPUpNO3 * NO3limit_scalar *ZoneBen[,'colonizedBenthicArea_1']# unit area
                    ZoneBen$algPUpNH4[] = ZoneBen$algPUpNH4 * NH4limit_scalar *ZoneBen[,'colonizedBenthicArea_1']# unit area
                    ZoneBen$algPUpPO4[] = ZoneBen$algPUpPO4 * Plimit_scalar *ZoneBen[,'colonizedBenthicArea_1']# unit area
                
                	ZoneBen[,'wN'] = ZoneBen[,'wN'] + (ZoneBen[,'netmineral']>0)*ZoneBen$im_PUpNO3*NO3limit_scalar*ZoneBen[,'colonizedBenthicArea_1']
                	ZoneBen[,'wN'] = ZoneBen[,'wN'] + (ZoneBen[,'netmineral']>0)*ZoneBen$im_PUpNH4*NH4limit_scalar*ZoneBen[,'colonizedBenthicArea_1']
                	ZoneBen[,'wP'] = ZoneBen[,'wP'] + (ZoneBen[,'netmineralP']>0)*ZoneBen$im_PUpPO4*Plimit_scalar*ZoneBen[,'colonizedBenthicArea_1']
                	
				ZoneBen$Pdenitri[] = ZoneBen$Pdenitri * (1.0 - ZoneBen$PdenitriFrac + NO3limit_scalar*ZoneBen$PdenitriFrac) *ZoneHyd$wba_1# unit area
				ZoneBen$Pnitri[] = ZoneBen$Pnitri * (1.0 - ZoneBen$PnitriFrac + NH4limit_scalar*ZoneBen$PnitriFrac) *ZoneHyd$wba_1# unit area
				
				ZoneBen[,'algalupNO3'] <- ZoneBen[,'algalupNO3'] + ZoneBen$algPUpNO3 # mgN/m2 # unit area
				ZoneBen[,'algalupNH4'] <- ZoneBen[,'algalupNH4'] + ZoneBen$algPUpNH4 # mgN/m2 # unit area
				ZoneBen[,'algalupP'] <- ZoneBen[,'algalupP'] + ZoneBen$algPUpPO4 # mgP/m2 # unit area
				ZoneBen[,'nitrification'] <- ZoneBen[,'nitrification'] + ZoneBen$Pnitri # unit area
				ZoneBen[,'denitrification'] <- ZoneBen[,'denitrification'] + ZoneBen$Pdenitri # unit area
                #print(paste(STS, ZoneBen$Pdenitri,sep=":"))
                
                
                
                
                
                
                ## hyporheic exchange and process *******
                	# zoneEnv_matrix[,BTS,'exchangeProp2SRIP'] m3/m2/s
                	# zoneEnv_matrix[,BTS,'exchangeProp2STR'] m3/m2/s
                #hold = 1.5*zoneEnv_matrix[,BTS,'exchangeProp2SRIP']*ZoneHyd$wba/ZoneCol$q * STSstep; hold[hold>1]=1; ## from str 2 rip
				#zoneEnv_matrix[,BTS,'exchange'] = 1.5*zoneEnv_matrix[,BTS,'exchangeProp2STR']*ZoneHyd$wba/ZoneStor$q * STSstep ## from rip 2 str
				
				TMP$h[] = zoneEnv_matrix_exchange * TMP$STSstep; TMP$h[TMP$h>1]=1;
				TMP$h2[] = TMP$h[] * ZoneHyd[,'crossarea'] / ZoneHyd[,'As']
				
				
				concMassChange[] = 0
				concChange[] = (TMP$h3[] - TMP$concNO3[]) #*zoneEnv_matrix[TMP$cond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
					TMP$cond = concChange>0 # positive "concChange" mean increase channel water[NO3]
					if(sum(TMP$cond)>0) concMassChange[TMP$cond] = vectorMin(
						concChange[TMP$cond] * TMP$h[TMP$cond] * ZoneCol[TMP$cond,'q'],# NO3 mass in channel increase due to exchange
						ZoneStor[TMP$cond,'no3'] ) ## <<---- should time barea? it's total mass, not areal in the new model
					TMP$cond = concChange<0 # negative "concChange" mean increase riparian water[NO3]
					if(sum(TMP$cond)>0) concMassChange[TMP$cond] = -vectorMin(
						-concChange[TMP$cond] * TMP$h2[TMP$cond] *ZoneStor[TMP$cond,'q'], # NO3 mass in ztorage increase due to exchange
						ZoneCol[TMP$cond,'no3']  ) ##<<-----
                        ZoneStor[,'no3'] = ZoneStor[,'no3'] - concMassChange[] ##<<---- correct with barea? it's total mass, not areal in the new model
					ZoneCol[,'no3'] = ZoneCol[,'no3'] + concMassChange[]
					ZoneBen[,'exchangeNO3'] = ZoneBen[,'exchangeNO3'] + concMassChange[]*ZoneHyd$wba_1 # areal
					
									
				concMassChange[] = 0	
				concChange[] = (TMP$h4[] - TMP$concNH4[]) #*zoneEnv_matrix[TMP$cond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
				TMP$cond = concChange>0 # positive "concChange" mean increase channel water[NO3]
					if(sum(TMP$cond)>0) concMassChange[TMP$cond] = vectorMin(
						concChange[TMP$cond] * TMP$h[TMP$cond] * ZoneCol[TMP$cond,'q'],
						ZoneStor[TMP$cond,'nh4'] )
					TMP$cond = concChange<0 # negative "concChange" mean increase riparian water[NO3]
					if(sum(TMP$cond)>0) concMassChange[TMP$cond] = -vectorMin(
						-concChange[TMP$cond] * TMP$h2[TMP$cond] *ZoneStor[TMP$cond,'q'],
						ZoneCol[TMP$cond,'nh4'] )
					ZoneStor[,'nh4'] = ZoneStor[,'nh4'] - concMassChange[]
					ZoneCol[,'nh4'] = ZoneCol[,'nh4'] + concMassChange[]
					ZoneBen[,'exchangeNH4'] = ZoneBen[,'exchangeNH4'] + concMassChange[]*ZoneHyd$wba_1
					
				concMassChange[] = 0	
				concChange[] = (TMP$h5[] - TMP$concPO4[]) #*zoneEnv_matrix[TMP$cond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
				TMP$cond = concChange>0 # positive "concChange" mean increase channel water[NO3]
					if(sum(TMP$cond)>0) concMassChange[TMP$cond] = vectorMin(
						concChange[TMP$cond] * TMP$h[TMP$cond] * ZoneCol[TMP$cond,'q'],
						ZoneStor[TMP$cond,'po4'] )
					TMP$cond = concChange<0 # negative "concChange" mean increase riparian water[NO3]
					if(sum(TMP$cond)>0) concMassChange[TMP$cond] = -vectorMin(
						-concChange[TMP$cond] * TMP$h2[TMP$cond] *ZoneStor[TMP$cond,'q'],
						ZoneCol[TMP$cond,'po4'] )
					ZoneStor[,'po4'] = ZoneStor[,'po4'] - concMassChange[]
					ZoneCol[,'po4'] = ZoneCol[,'po4'] + concMassChange[]
                	ZoneBen[,'exchangePO4'] = ZoneBen[,'exchangePO4'] + concMassChange[]*ZoneHyd$wba_1
                	
  					# ZoneCol
  					# ZoneStor
  					# TMP$concNO3 = ZoneCol[,'no3']/ZoneCol$q; TMP$concNO3[is.na(TMP$concNO3)|is.infinite(TMP$concNO3)]=0
  					# TMP$h3 = ZoneStor[,'no3']/ZoneStor$q
  					# cbind(TMP$concNO3, TMP$h3)
  					# plot(TMP$concNO3, type='b'); abline(v=8,col='red')
  					# plot(ZoneBen[,'algalupNO3'], type='b')
  					
  				ZoneStor[,'baseNH4balance'] = ZoneStor[,'nh4']-ZoneStor[,'baseNH4']
                ZoneStor[,'baseNO3balance'] = ZoneStor[,'no3']-ZoneStor[,'baseNO3']
          		TMP$cond = ZoneStor[,'baseNO3balance']>0
          		TMP$h[] = 0
                if(sum(TMP$cond)>0){
	                # hyporheicDenitrification[TMP$cond] = hyporheicDenitrification[TMP$cond] + vectorMin(
	                	# zoneEnv_matrix[TMP$cond,BTS,'Nuptake']*zoneInfo[TMP$cond,'len'],
	                	# ZoneStor[TMP$cond,'baseNO3balance'])
	                # ZoneStor[TMP$cond,'baseNO3balance'] = ZoneStor[TMP$cond,'baseNO3balance'] - vectorMin(
	                	# zoneEnv_matrix[TMP$cond,BTS,'Nuptake']*zoneInfo[TMP$cond,'len'],
	                	# ZoneStor[TMP$cond,'baseNO3balance'])
	                	
	                TMP$h[TMP$cond] = vectorMin(
		                # zoneEnv_matrix[,,'Nuptake']
	                	#zoneEnv_matrix$landuptake[BTS, TMP$cond]*zoneInfo[TMP$cond,'len'],
                        zoneEnv_matrix[TMP$cond, BTS,'Nuptake']*zoneInfo[TMP$cond,'len'], # unit for this uptake is mgN/m2/m, right?
	                	ZoneStor[TMP$cond,'baseNO3balance']*ZoneHyd[TMP$cond,'wba_1'] )
                	
	                hyporheicDenitrification[TMP$cond] = hyporheicDenitrification[TMP$cond] + TMP$h[TMP$cond]
                    ZoneStor[TMP$cond,'baseNO3balance'] = ZoneStor[TMP$cond,'baseNO3balance'] - TMP$h[TMP$cond]*ZoneHyd[TMP$cond,'wba'] # consistent to the total mass assumption
                }# if
                	
                	TMP$h = ZoneBen$algPUpNO3*ZoneBen[,'colonizedBenthicArea']*ZoneHyd$wba_1 + ZoneBen$Pdenitri + ZoneBen$im_PUpNO3 - ZoneBen$Pnitri # <<----------this is all NO3
                    ## here ZoneBen$algPUpNO3 was corrected for colonizedBenthicArea, not all BenthicArea!
				ZoneBen[,'totaluptake'] = ZoneBen[,'totaluptake'] + TMP$h #areal uptake
					TMP$h = ZoneCol[,'no3']/ZoneCol$q*ZoneHyd[,'depth']*ZoneHyd[,'velocity']/TMP$h
                ZoneBen[,'Sw'] = ZoneBen[,'Sw'] + TMP$h*TMP$STSstep # uptake length
                
                ###################################################################################### update hydraulic parameters
                crossSection = crossSection + ZoneHyd[,'crossarea']
                
                hydraulicHOLD=sapply(numZoneIndex,function(i){ c(
					Morphology_matrix[[i]]@Vol2Z(ZoneCol[i,'q']), #1
		            Morphology_matrix[[i]]@Vol2V(ZoneCol[i,'q']), #2
		            Morphology_matrix[[i]]@Vol2WP(ZoneCol[i,'q']), #3
		            Morphology_matrix[[i]]@Vol2Q(ZoneCol[i,'q']), #4
                    Morphology_matrix[[i]]@Vol2CA(ZoneCol[i,'q']), # 5
                    Morphology_matrix[[i]]@maxQ, #5 -> 6
                    Morphology_matrix[[i]]@maxVol ) }) #6 -> 7
				ZoneHyd[,'depth'] = hydraulicHOLD[1,]
				ZoneHyd[,'velocity'] = hydraulicHOLD[2,] #m/s
				ZoneHyd$wba = zoneInfo[,'len']*hydraulicHOLD[3,]
				ZoneHyd$wba_1 = 1.0 / ZoneHyd$wba
                ZoneHyd[,'crossarea'] = hydraulicHOLD[5,]
				ZoneHyd$q = hydraulicHOLD[4,]

                #print(paste(c(BTS,STS, ZoneCol$q),collapse=' '))
                if(is.na(sum(ZoneCol$q))|sum(ZoneCol$q<0)>0)
                	print(paste(c(BTS,STS, TMP$STSstep,round(abs(ZoneHyd$q)/ZoneCol$q*TMP$STSstep,6)),collapse=' '))
                	
                if(sum(ZoneStor[,'nh4']<0)>0)
                	print(paste(c(BTS,STS, TMP$STSstep, ZoneStor[,'nh4']),collapse=' '))
                # debug
                # print(	paste(
                	# STS, '|',
                	# paste(round(ZoneCol[1:3,'q'],3),collapse=','),'|',
                	# paste(round(ZoneExport[1:3,'q'],6),collapse=','),'|',
                	# paste(round(ZoneInput[1:3,'q'],6),collapse=',')
                	# ))		
                ZoneInput[ ] = 0
             
             	STS = STS + TMP$STSstep
             	TMP$h = ZoneHyd$q/ZoneCol$q;
             	TMP$STSstep = min(timemax_STS-STS,max(1,apply(TMP$h %o% STSstepOption<0.25,2,prod)* STSstepOption))
		}#STS
		
		ZoneBen[,'Sw'] = ZoneBen[,'Sw']* 0.0002777778 ## correcting aggregating
			
		## debugging
		# ZoneStor
		# ZoneCol
		# ZoneBen
		# ZoneBen[,'requiredN']
		# ZoneBen[,'dN']+ZoneBen[,'wN']
		# mean(ZoneBen[,'denitrification']*24) # 17.34134
		
		
		# end.time <- Sys.time()
		# time.taken <- end.time - start.time
		# time.taken	
		# ZoneBank
		# ZoneCol
		# ZoneCol[,'no3']/ZoneCol$q
		# ZoneStor
		# STSaggregateVar/3600
		# ZoneHyd
		# cbind(TMP$concNO3, TMP$concNH4, TMP$concPO4, ZoneStor)
		# cbind(TMP$concNO3,ZoneStor[,'no3'],ZoneStor[,'no3'] + hyporheic_nitrification)
		# 
			
		## algal process (hourly) # unit area
		
			## NPP
			ZoneBen[,'algalNF'] = 1-(ZoneBen[,'algalCN']-algPar['minCN'])/(algPar['maxCN']-algPar['minCN']);
				ZoneBen[ZoneBen[,'algalNF']<0,'algalNF']=0
				ZoneBen[ZoneBen[,'algalNF']>1,'algalNF']=1
				
			ZoneBen[,'algalPF'] = 1-(ZoneBen[,'algalCP']-algPar['minCP'])/(algPar['maxCP']-algPar['minCP'])
				ZoneBen[ZoneBen[,'algalPF']<0,'algalPF']=0
				ZoneBen[ZoneBen[,'algalPF']>1,'algalPF']=1
			ZoneBen[,'algalnutrientF'] <- apply(ZoneBen[,c('algalNF','algalPF')],1,min)
			
			ZoneBen[,'algalselfF'] <- 1 / (1 + algPar['selflimitcoef']* ZoneBen[,'algalc']) # algalc mgC/m2; selflimitcoef m2/mgC
			
			# exp(-(zoneEnv_matrix[,'PAR',BTS]-b)/a)
			# month = as.numeric(format(SimulationTime[as.integer(BTS/24)],'%m'))
            # if( as.numeric(format(SimulationTime[ceiling(BTS/24)],'%m')) %in% c(6,7,8,9,10)){
            #    ZoneBen[,'algallightF'] <- zoneEnv_matrix[,BTS,'PARf'] #* lightAdjust # 6,7,8,9,10 month
            # }else{
            # 	ZoneBen[,'algallightF'] <- zoneEnv_matrix[,BTS,'PARf'] #* lightAdjustII
            # }#
            ZoneBen[,'algallightF'] <- zoneEnv_matrix[,BTS,'PARf'] #*(0.9+0.1*exp(-0.1666667*ZoneCol[,'POC']/ZoneCol$q));
            
            
			ZoneBen[,'algalNPP'] <- ZoneBen[,'algalc']* algPar['maxgrowthRate'] * ZoneBen[,'algalnutrientF'] * ZoneBen[,'algalselfF'] * ZoneBen[,'algallightF'] * algalQ10f * timemax_STS * (1-alghaDmgIndex)
					
	
		ZoneBen[,'resp'] <- ZoneBen[,'algalNPP']*(1-algPar['gcoef']) + ZoneBen[,'algalmialC']
		ZoneBen[,'algalc'] <- ZoneBen[,'algalc'] + algPar['gcoef']*ZoneBen[,'algalNPP']
		ZoneBen[,'algaln'] <- ZoneBen[,'algaln'] + ZoneBen[,'algalupNO3'] + ZoneBen[,'algalupNH4']
		ZoneBen[,'algalp'] <- ZoneBen[,'algalp'] + ZoneBen[,'algalupP']
		
		## convert areal flux to total flux to exchange with water column resources
        #ZoneBen[,'algalmialC'] = ZoneBen[,'algalmialC'] * ZoneBen[,'colonizedBenthicArea'] # from unit area to full area
        #ZoneBen[,'algalmialN'] = ZoneBen[,'algalmialN'] * ZoneBen[,'colonizedBenthicArea']
        #ZoneBen[,'algalmialP'] = ZoneBen[,'algalmialP'] * ZoneBen[,'colonizedBenthicArea']
        #ZoneBen$algEntC = ZoneBen$algEntC * ZoneBen[,'colonizedBenthicArea'] # from unit area to full area
        #ZoneBen[,'algEntN'] = ZoneBen[,'algEntN'] * ZoneBen[,'colonizedBenthicArea']
        #ZoneBen[,'algEntP'] = ZoneBen[,'algEntP'] * ZoneBen[,'colonizedBenthicArea']
		
		## update C:N and C:P
		ZoneBen[,'algalCN'] <- ZoneBen[,'algalc']/ZoneBen[,'algaln']
		ZoneBen[,'algalCP'] <- ZoneBen[,'algalc']/ZoneBen[,'algalp']
        ZoneBen[,'algUpNF'] <- ZoneBen[,'algalc']*(1-ZoneBen[,'algalNF'])*ZoneBen[,'colonizedBenthicArea']
        ZoneBen[,'algUpPF'] <- ZoneBen[,'algalc']*(1-ZoneBen[,'algalPF'])*ZoneBen[,'colonizedBenthicArea']
        
		dailyAveragewba[] = dailyAveragewba  + ZoneHyd$wba* 0.04166667
		
	
		# decomposition process (hourly ... add later)	
		TMP$cond[] = ZoneBen[,'netmineral']>0 # 0.5/cn - 1/CN (uptake)
        if(sum(TMP$cond)>0) ZoneBen[TMP$cond,'decay'] = vectorMin(ZoneBen[TMP$cond,'decay'], ZoneBen[TMP$cond,'wN']/ZoneBen[TMP$cond,'Ncoef']) ## N can limit decay
        TMP$cond[] = ZoneBen[,'netmineralP']>0 # 0.5/cn - 1/CN (uptake)
        if(sum(TMP$cond)>0) ZoneBen[TMP$cond,'decay'] = vectorMin(ZoneBen[TMP$cond,'decay'], ZoneBen[TMP$cond,'wP']/ZoneBen[TMP$cond,'Pcoef']) ## P can limit decay
        
		TMP$h[] = ZoneBen[,'decay']/(ZoneBen$LabC+ZoneBen$IntC* immo_class2decayRatio +ZoneBen$RecC* immo_class3decayRatio) # decay by immo and mi
		TMP$h2[] = ZoneBen[,'mi_decayC']/(ZoneBen$LabC+ZoneBen$IntC* mi_class2decayRatio +ZoneBen$RecC* mi_class3decayRatio)
		TMP$h2[TMP$h2>1] = 1.0
		TMP$h[TMP$h>1]=1.0
			# cbind(
				# ZoneBen$LabC+ZoneBen$IntC* mi_class2decayRatio +ZoneBen$RecC* mi_class3decayRatio,
				# ZoneBen[,'mi_decayC']
			# )
		ZoneBen$LabC = ZoneBen$LabC * (1 - TMP$h) * (1 - TMP$h2) # mgC/m2 unit area
		ZoneBen$LabN = ZoneBen$LabN * (1 - TMP$h) * (1 - TMP$h2)# mgN/m2 unit area
		ZoneBen$LabP = ZoneBen$LabP * (1 - TMP$h) * (1 - TMP$h2)
		ZoneBen$IntC = ZoneBen$IntC * (1 - TMP$h* immo_class2decayRatio) * (1 - TMP$h2* mi_class2decayRatio)
		ZoneBen$IntN = ZoneBen$IntN * (1 - TMP$h* immo_class2decayRatio) * (1 - TMP$h2* mi_class2decayRatio)
		ZoneBen$IntP = ZoneBen$IntP * (1 - TMP$h* immo_class2decayRatio) * (1 - TMP$h2* mi_class2decayRatio)
		ZoneBen$RecC = ZoneBen$RecC * (1 - TMP$h* immo_class3decayRatio) * (1 - TMP$h2* mi_class3decayRatio)
		ZoneBen$RecN = ZoneBen$RecN * (1 - TMP$h* immo_class3decayRatio) * (1 - TMP$h2* mi_class3decayRatio)
		ZoneBen$RecP = ZoneBen$RecP * (1 - TMP$h* immo_class3decayRatio) * (1 - TMP$h2* mi_class3decayRatio)
		
		ZoneBen[,'mi'] = ZoneBen[,'mi'] + ZoneBen[,'decay']*micPar['gcoef']
		ZoneBen[,'mi'] = ZoneBen[,'mi'] + ZoneBen[,'decay']*micPar['gcoef']
        ZoneBen[,'resp'] = ZoneBen[,'resp'] + ZoneBen[,'decay']*(1-micPar['gcoef']) + ZoneBen[,'mi_decayC']-ZoneBen[,'mi_grow']
		ZoneBen[,'mial'] = 0
		TMP$cond = ZoneBen[,'netmineral']<0
		if(sum(TMP$cond)>0) ZoneBen[TMP$cond,'mial'] = -ZoneBen[TMP$cond,'netmineral']
		ZoneBen[,'mial'] = ZoneBen[,'mial'] + ZoneBen[,'mi_decayN'] ## for mi to release
		
		# hourly aggregate
        BTSaggregateVar[,'velocity'] = BTSaggregateVar[,'velocity'] + ZoneHyd[,'velocity'] ## need to be /24
        BTSaggregateVar[,'depth'] = BTSaggregateVar[,'depth'] + ZoneHyd[,'depth'] ## need to be /24
        BTSaggregateVar[,'benthicArea'] = BTSaggregateVar[,'benthicArea'] + ZoneHyd$wba ## need to be /24
        BTSaggregateVar[,'denitrification'] = BTSaggregateVar[,'denitrification'] + ZoneBen[,'denitrification'] # unit area
        BTSaggregateVar[,'nitrification'] = BTSaggregateVar[,'nitrification'] + ZoneBen[,'nitrification'] # unit area
        BTSaggregateVar[,'algalNO3uptake'] = BTSaggregateVar[,'algalNO3uptake'] + ZoneBen[,'algalupNO3'] # unit area
        BTSaggregateVar[,'algalNH4uptake'] = BTSaggregateVar[,'algalNH4uptake'] + ZoneBen[,'algalupNH4'] # unit area
        BTSaggregateVar[,'algalNPP'] = BTSaggregateVar[,'algalNPP'] + ZoneBen[,'algalNPP'] # unit area
        BTSaggregateVar[,'mialN'] = BTSaggregateVar[,'mialN'] + ZoneBen[,'algalmialN']
        BTSaggregateVar[,'algalCN'] = BTSaggregateVar[,'algalCN'] + ZoneBen[,'algalCN'] ## need to be /24
        BTSaggregateVar[,'algalCP'] = BTSaggregateVar[,'algalCP'] + ZoneBen[,'algalCP'] ## need to be /24
        BTSaggregateVar[,'algalC'] = BTSaggregateVar[,'algalC'] + ZoneBen[,'algalc'] ## need to be /24
        BTSaggregateVar[,'algalCloss'] = BTSaggregateVar[,'algalCloss'] + ZoneBen$algEntC + ZoneBen[,'algaldeadC']
        BTSaggregateVar[,'algalR'] = BTSaggregateVar[,'algalR'] + ZoneBen[,'resp']
        BTSaggregateVar[,'detritusR'] = BTSaggregateVar[,'detritusR'] + ZoneBen[,'resp']
        BTSaggregateVar[,'Sw'] = BTSaggregateVar[,'Sw'] + ZoneBen[,'Sw']## need to be /24
        
		BTSaggregateVar[,'netmineral'] = BTSaggregateVar[,'netmineral'] + ZoneBen[,'mial'] - ZoneBen[,'wN'] # unit area
		BTSaggregateVar[,'detric'] = BTSaggregateVar[,'detric'] + (ZoneBen$LabC+ ZoneBen$IntC+ ZoneBen$RecC)  ## need to be /24
		
        BTSaggregateVar[,'detriCN'] = BTSaggregateVar[,'detriCN'] + 1/ZoneBen[,'foodNC'] ## need to be /24 <<----------------------- debugging
        BTSaggregateVar[,'detriCP'] = BTSaggregateVar[,'detriCP'] + 1/ZoneBen[,'foodPC'] ## need to be /24 <<----------------------- debugging
		
		BTSaggregateVar[,'microbC'] = BTSaggregateVar[,'microbC'] + ZoneBen[,'mi']  ## need to be /24 <<----- ***
		BTSaggregateVar[,'miC'] = BTSaggregateVar[,'miC'] + ZoneBen[,'mi']
        
		BTSaggregateVar[,'tmperature'] = BTSaggregateVar[,'tmperature'] + zoneEnv_matrix[,BTS,'temperature'] ## need to be /24
		BTSaggregateVar[,'PARf'] = BTSaggregateVar[,'PARf'] + zoneEnv_matrix[,BTS,'PARf'] ## need to be /24
		BTSaggregateVar[,'hdmg'] = BTSaggregateVar[,'hdmg'] + alghaDmgIndex
        
        BTSaggregateVar[,'algalNlimit'] = BTSaggregateVar[,'algalNlimit'] + ZoneBen[,'algalNF']
        BTSaggregateVar[,'algalPlimit'] = BTSaggregateVar[,'algalPlimit'] + ZoneBen[,'algalPF']
        
		## reset
		ZoneBen[,'algalupP']=0;
		ZoneBen[,'algalupNO3']=0;
		ZoneBen[,'algalupNH4']=0;
		# ZoneBen[,'resp'] # calculated hourly
		# ZoneBen[,'algalmialC'] # calculated hourly
		# ZoneBen[,'algalmialN'] # calculated hourly
		# ZoneBen[,'algalmialP'] # calculated hourly
		# ZoneBen$algEntC # calculated hourly
		# ZoneBen[,'algEntN'] # calculated hourly
		# ZoneBen[,'algEntP'] # calculated hourly
		ZoneBen[,'nitrification']=0
		ZoneBen[,'denitrification']=0
		ZoneBen[,'wN'] = 0 ## N imilization
        ZoneBen[,'wP'] = 0
		# ZoneBen[,'resp'] # calculated hourly
		# ZoneBen[,'Ncoef'] # calculated hourly
		# ZoneBen[,'netmineral'] # calculated hourly
		# ZoneBen[,'decay'] # calculated hourly
		# ZoneBen[,'mial'] # calculated hourly
		
		if(BTS%%24==0){

			## adjust algal colonized area
				# cond = ZoneBen[,'colonizedBenthicArea'] < dailyAveragewba
				# # crowded, need expansion; ZoneBen[cond,'algalselfF']=low
				# ZoneBen[cond,'colonizedBenthicArea'] =
					# (0.9666667 - 0.03333333*(1-ZoneBen[cond,'algalselfF']) )* ZoneBen[cond,'colonizedBenthicArea'] +
					# (0.03333333 + 0.03333333*(1-ZoneBen[cond,'algalselfF']) )* dailyAveragewba[cond] #m2
					
				# cond = ZoneBen[,'colonizedBenthicArea'] > 1.5*dailyAveragewba & dailyAveragewba>0
				# # no density effect
				# ZoneBen[cond,'colonizedBenthicArea'] =
					# 0.99666667*ZoneBen[cond,'colonizedBenthicArea'] +
					# 0.003333333* dailyAveragewba[cond] #m2
					
				# cond = ZoneBen[,'colonizedBenthicArea'] > maxChannelBenthicArea
				# ZoneBen[cond,'colonizedBenthicArea'] = maxChannelBenthicArea[cond]
				# ZoneBen[,'colonizedBenthicArea_1'] = 1.0/ZoneBen[,'colonizedBenthicArea']
			
			# ZoneBen[,'algalc'] = ZoneBen[,'algalc'] * ZoneBen[,'colonizedBenthicArea_1']
			# ZoneBen[,'algaln'] = ZoneBen[,'algaln'] * ZoneBen[,'colonizedBenthicArea_1']
			# ZoneBen[,'algalp'] = ZoneBen[,'algalp'] * ZoneBen[,'colonizedBenthicArea_1']
			
			
            
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
            BTSaggregateVar[,'miC'] = BTSaggregateVar[,'miC']* 0.04166667
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
                ZoneInputTrack[upZs[1],] + colSums(STSLateralinputTrack[upZs[2]:upZs[length(upZs)],]), # inputs: discharge waterNO3 waterNH4 waterPO4 POC PON POP (flux/load) from upsteram & lateral IN
                STSaggregateVar[upZs[length(upZs)],], # outputs: discharge waterNO3 waterNH4 waterPO4 POC PON POP (flux/load)
                sum(BTSaggregateVar[upZs,'benthicArea']),
                (zoneInfo[upZs,'len']/sum(zoneInfo[upZs,'len'])) %*% BTSaggregateVar[upZs,c('velocity','depth')],
                (BTSaggregateVar[upZs,'benthicArea']/sum(BTSaggregateVar[upZs,'benthicArea'])) %*% BTSaggregateVar[upZs, BTSaggregateVarPrint],
               	100*(1-(STSaggregateVar[upZs[length(upZs)],'nitrate'])/(ZoneInputTrack[upZs[1],'no3']+sum(STSLateralinputTrack[upZs[2]:upZs[length(upZs)],'no3'])) )/sum(zoneInfo[upZs,'len']),
                
                #downstream (daily) full 11 - 34; odd = pool
                ZoneInputTrack[downZs[1],] + colSums(STSLateralinputTrack[downZs[2]:downZs[length(downZs)],]), # inputs: discharge waterNO3 waterNH4 waterPO4 POC PON POP (flux/load) from upsteram & lateral IN
                STSaggregateVar[downZs[length(downZs)],], # outputs: discharge waterNO3 waterNH4 waterPO4 POC PON POP (flux/load)
                sum(BTSaggregateVar[downZs,'benthicArea']),
                (zoneInfo[downZs,'len']/sum(zoneInfo[downZs,'len'])) %*% BTSaggregateVar[downZs,c('velocity','depth')],
                (BTSaggregateVar[downZs,'benthicArea']/sum(BTSaggregateVar[downZs,'benthicArea'])) %*% BTSaggregateVar[downZs, BTSaggregateVarPrint],
                100*(1-(STSaggregateVar[downZs[length(downZs)],'nitrate'])/(ZoneInputTrack[downZs[1],'no3']+sum(STSLateralinputTrack[downZs[2]:downZs[length(downZs)],'no3'])))/sum(zoneInfo[downZs,'len']),
                
                sum(ZoneBen[upZs,'exchangeNO3']* (BTSaggregateVar[upZs,'benthicArea']/sum(BTSaggregateVar[upZs,'benthicArea'])) ), # it's areal
                sum(ZoneBen[downZs,'exchangeNO3']* (BTSaggregateVar[downZs,'benthicArea']/sum(BTSaggregateVar[downZs,'benthicArea'])) ),
                
                ## input nitrate
                #ZoneCol[,'no3'] /ZoneCol$q
                #STSaggregateVar[,'nitrate'] /STSaggregateVar$q ## export
                #STSaggregateVar # tracking local output
                #ZoneInputTrack # tracking lcoal lateral and upstream inputs <-----
                #STSLateralinputTrack # tracking lcoal lateral only  <-----
                #downZs = 11:60 
				#upZs = 1:8
                sapply(upZs,function(i){
                    totalinput = ZoneInputTrack[upZs[1],'no3']
                    if(i>upZs[1]) totalinput = totalinput + sum(STSLateralinputTrack[upZs[2]:i,'no3'])
                    #return <- 100*(1-STSaggregateVar[i,'nitrate']/totalinput)
                    return <- totalinput
                }), # upstream
                sapply(9:10,function(i){
                    totalinput = ZoneInputTrack[9,'no3']
                    if(i>9) totalinput = totalinput + sum(STSLateralinputTrack[10:i,'no3'])
                    #return <- 100*(1-STSaggregateVar[i,'nitrate']/totalinput)
                    return <- totalinput
                }), # mid transition
                sapply(downZs[1]:max(zoneInfo[,'zoneID']),function(i){
                    totalinput = ZoneInputTrack[downZs[1],'no3']
                    if(i>downZs[1]) totalinput = totalinput + sum(STSLateralinputTrack[downZs[2]:i,'no3'])
                    #return <- 100*(1-STSaggregateVar[i,'nitrate']/totalinput)
                    return <- totalinput
                }), # downstream
                
                ## output nitrate
                STSaggregateVar[,'nitrate']
                ), '\n', file=outputFile_buff,sep=',')

                
                # downstream (daily) -- last 60 m [27,29,31,33]
                # ZoneInputTrack[27,] + colSums(STSLateralinputTrack[28:34,]), # inputs: discharge waterNO3 waterNH4 waterPO4 POC PON POP (flux/load) from upsteram & lateral IN
                # STSaggregateVar[34,], # outputs: discharge waterNO3 waterNH4 waterPO4 POC PON POP (flux/load)
                # sum(BTSaggregateVar[c(27,29,31,33),'benthicArea']),
                # (zoneInfo[c(27,29,31,33),'len']/sum(zoneInfo[c(27,29,31,33),'len'])) %*% BTSaggregateVar[c(27,29,31,33),c('velocity','depth')],
                # (BTSaggregateVar[c(27,29,31,33),'benthicArea']/sum(BTSaggregateVar[c(27,29,31,33),'benthicArea'])) %*% BTSaggregateVar[c(27,29,31,33),c(
                	# 'nitrification','denitrification','algalNO3uptake','algalNH4uptake','algalNPP','mialN','algalCN','algalCP','algalC',
                	# 'algalCloss','netmineral','detric','detriCN','microbC','algalR','detritusR','Sw')],
                # 100*(1-(STSaggregateVar[34,'nitrate'])/(ZoneInputTrack[27,'no3']+sum(STSLateralinputTrack[c(27,29,31,33),'no3'])))/sum(zoneInfo[27:34,'len'])
                # ), '\n', file=outputFile_buff,sep=',')
				
			BTSaggregateVar[ ] = 0
			STSaggregateVar[ ] = 0
			dailyAveragewba[ ] = 0
			ZoneBen[,'totaluptake'] = 0
			ZoneBen[,'Sw'] = 0
			ZoneBen[,'exchangeNO3'] = 0
			ZoneBen[,'exchangeNH4'] = 0
			ZoneBen[,'exchangePO4'] = 0
			hyporheicDenitrification[] =0
			ZoneInputTrack[] = 0
			STSLateralinputTrack[] = 0
		}## daily output -- BTS%%24==0	
		# ZoneBank
		# cbind(zoneInfo,ZoneCol)
		# STSaggregateVar/3600
		# ZoneHyd
		# check = data.frame(
			# detritus1 = ZoneBen$LabC,
			# detritus2 = ZoneBen$IntC,
			# detritus3 = ZoneBen$RecC,
			# detritus1cn = ZoneBen$LabC/ZoneBen$LabN,
			# detritus2cn = ZoneBen$IntC/ZoneBen$IntN,
			# detritus3cn = ZoneBen$RecC/ZoneBen$RecN,
			# microbes = ZoneBen[,'mi'],
			# mi = ZoneBen[,'mi'],
			# mi_grow = ZoneBen[,'mi_grow'],
			# mi_decayC = ZoneBen[,'mi_decayC'],
			# mi_decayN = ZoneBen[,'mi_decayN'],
			# mi_decayP = ZoneBen[,'mi_decayP'],
			# netmineral = ZoneBen[,'netmineral'],
			# netmineralP = ZoneBen[,'netmineralP'],
			# decay = ZoneBen[,'decay'],
			# foodNC = ZoneBen[,'foodNC'],
			# foodPC = ZoneBen[,'foodPC'],
			# mi_resNC = ZoneBen[,'mi_resNC'],
			# mi_resPC = ZoneBen[,'mi_resPC'],
			# resp = ZoneBen[,'resp'],
			# algaldeadc = ZoneBen[,'algaldeadC']
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
# result$inputQ = daily$q
# result$inputNO3 = daily[,'no3']
# result$inputNO3c = daily[,'no3c']
# result$inputNH4c = daily[,'nh4c']
# result$outputNO3c = check[,'nitrate']/(check$q*1000) # (mgN/s) / (L/s) = mgN/L
# result$outputNH4c = check[,'ammonium']/(check$q*1000) # (mgN/s) / (L/s) = mgN/L

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


