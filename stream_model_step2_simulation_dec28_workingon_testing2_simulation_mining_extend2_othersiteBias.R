source('stream_model_step0_predefine_CBT_extend2.R') 
arg=commandArgs(T)
#arg=c(172, 100, 100, 60, 201, 123, 96, 96, 20,'test11c2.csv')
# Rscript stream_model_step2_simulation_dec28_workingon_testing2_simulation_mining_extend2_othersiteBias.R 172 100 100 60 201 123 96 96 20 6
source("https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_misc.r")	
detritusCinput = 0 #219650.7/24 # per hour input <<----- what is it?
detritusCN = 40 # <<----- what is it?
detritusCP = 444.5
p0 = 0.485 # soil 8 silt_loam
p_d = 1/4000

#####################################sindex = list() # sensitivity index [1-200]
site = c('GFGL','GFGB','GFVN','GFCP','POBR','BARN','MCDN','DRKR')
sindex = list() # sensitivity index [1-200]
sindex$flow = as.numeric(arg[1]) #1 [used]
sindex$ninput = as.numeric(arg[2]) #2 [used]
sindex$PAR_temp = as.numeric(arg[3]) #3 [used]
sindex$exchange = as.numeric(arg[4]) #4 [used]
sindex$velocity = as.numeric(arg[5]) #5 [used]
sindex$denitri = as.numeric(arg[6]) #6 [used]
sindex$algal = as.numeric(arg[7]) #7 [used]
sindex$decomp = as.numeric(arg[8]) #8 [used]
sindex$nitri = as.numeric(arg[9]) #9 [used]
sindex$outputFile = paste(site[as.numeric(arg[10])],'_SLB_11_bias.csv',sep='')
print(arg[10])
print(sindex$outputFile)

#yy=c(seq(0,30,1),40,seq(70,450,1)); zz=1/(1+exp( -(yy-50)/(10*(1+4*(yy>50))))); # plot(yy, zz,pch=21 ) # old time
yy=seq(0,3000); zz = 1/ ( 1+exp(  -(yy-300)/50 ) );
PAR2PARF = splinefun(yy,zz)

# need uptake metric as % removal
	sensitivity = list()
	sensitivity$zoneinfo = 'channel/GIS_channel_lateral_full.csv'	
	sensitivity$pathway = 'sensitivity_inputs/'
	sensitivity$AjQ = paste('other_BES_sites/',site[as.numeric(arg[10])],'2SLB_qnp_bias.csv',sep='') # <-- * change 2SLB_qnp_bias.csv
    sensitivity$AjN = paste('other_BES_sites/',site[as.numeric(arg[10])],'2SLB_qnp_bias.csv',sep='') # <-- * change 2SLB_qnp_bias.csv
	sensitivity$AjP = paste('other_BES_sites/',site[as.numeric(arg[10])],'2SLB_qnp_bias.csv',sep='') # <-- * change 2SLB_qnp_bias.csv
	sensitivity$rPAR = paste(sensitivity$pathway, 'restored_Parps.csv', sep='')
	sensitivity$uPAR = paste(sensitivity$pathway, 'unrestored_Parps.csv', sep='')
	sensitivity$rwater = paste(sensitivity$pathway, 'restored_water_temperature.csv', sep='')
	sensitivity$uwater = paste(sensitivity$pathway, 'unrestored_water_temperature.csv', sep='')
	sensitivity$vPathway = paste(sensitivity$pathway, 'sensitivity_channel/', sep='')
	sensitivity$denitrification_coef = paste(sensitivity$pathway, 'denitrification_ceof.csv', sep='')
	sensitivity$totaluptake_coef = paste(sensitivity$pathway, 'totaluptake_ceof.csv', sep='')
	sensitivity$algaluptake_coef = paste(sensitivity$pathway, 'algaluptake_ceof.csv', sep='')
	sensitivity$AjS365Q = paste(sensitivity$pathway, 'AJ_S365Q.csv', sep='')
	sensitivity$uPARS365 = paste(sensitivity$pathway, 'unrestored_Parps_S365.csv', sep='')
	sensitivity$rPARS365 = paste(sensitivity$pathway, 'restored_Parps_S365.csv', sep='')
	sensitivity$rwaterS365 = paste(sensitivity$pathway, 'restored_water_temperature_S365.csv', sep='')
	sensitivity$uwaterS365 = paste(sensitivity$pathway, 'unrestored_water_temperature_S365.csv', sep='')
	sensitivity$vOrder = paste(sensitivity$pathway, 'velocity_order.csv', sep='')
	sensitivity$otherAnnualPattern = paste(sensitivity$pathway, 'annualpatterns_long2.csv', sep='') # annualpatterns.csv
	
	# format of the files
	sensit_velocity = read.csv(sensitivity$vOrder) ## velocity combination, upstream,rift,pool (keep)
	sensit_denitrification = read.csv(sensitivity$denitrification_coef) # --- (keep)
	sensit_algaluptake = read.csv(sensitivity$algaluptake_coef) # --- (keep)
	
    ## ---- here is the new inputs
    # sensitivity$AjQ = 'AJ_synthicFlow.csv' is a spreadsheet: year month day flow (input as L/s -> m3/s)
    # sensitivity$AjN = 'usgs01589290_WRTDS_GFGB.csv' is a spreadsheet: year month day conc(mgN/L); same time series as Q
    # sensitivity$AjP = 'usgs01589290_WRTDS_GFGB_po4_yr1.csv' (for now): same as N but conc(µgP/L); same time series as Q
	sensit_lateralQ = read.csv(sensitivity$AjQ); sensit_lateralQ$date = as.Date(paste(sensit_lateralQ[,'day'], sensit_lateralQ[,'month'], sensit_lateralQ[,'year'],sep="-"),format="%d-%m-%Y"); # shift col by 3
	sensit_lateralN = read.csv(sensitivity$AjN); sensit_lateralN$date = as.Date(paste(sensit_lateralN[,'day'], sensit_lateralN[,'month'], sensit_lateralN[,'year'],sep="-"),format="%d-%m-%Y"); # shift col by 3
    sensit_lateralP = read.csv(sensitivity$AjP); sensit_lateralP$date = as.Date(paste(sensit_lateralP[,'day'], sensit_lateralP[,'month'], sensit_lateralP[,'year'],sep="-"),format="%d-%m-%Y"); # shift col by 3
	## ----------------------------



	sensit_rePAR = read.csv(sensitivity$rPAR) # from 2012/1/1 to 2017/12/31 hourly; no col shift
	sensit_retemp = read.csv(sensitivity$rwater) # from 2012/1/1 to 2017/12/31 hourly; no col shift
	sensit_unPAR = read.csv(sensitivity$uPAR)
	sensit_untemp = read.csv(sensitivity$uwater)
	sensit_otherpattern = read.csv(sensitivity$otherAnnualPattern)
    # sensit_PO4 = read.csv(sensitivity$AjP)
	
	# sensit_denitrification[sindex,]
		# nitrification 0.00002339181 s-1 in my manuscript
		# nitrification 0 - 0.0002314815 mgN/m2/s in Bernhardt et al. 2002 
	    # median(c(1.13,0.67,11.95,8.85,4.55,10.14,0.23,17.76,13.69,0,2.35,2.18,5.16,1.29,7.85,2.20,4.19,1.26,1.29)) # 2.35 mgN/m2/day
		# mean(c(1.13,0.67,11.95,8.85,4.55,10.14,0.23,17.76,13.69,0,2.35,2.18,5.16,1.29,7.85,2.20,4.19,1.26,1.29)) # 5.09 mgN/m2/day
	    biochem_parameter_max_nitrif = seq(0, 2.444973e-05, length.out=200)[sindex$nitri]/23 #17 # exactly reasonable
        if(biochem_parameter_max_nitrif<=0) biochem_parameter_max_nitrif = 5.341868e-10;
	    biochem_parameter_Mulholland_max_denitr = 10^sensit_denitrification[sindex$denitri,1]*0.01 #1.059254e-05 # cm/s -> m/s
		biochem_parameter_Mulholland_pow_denitr = sensit_denitrification[sindex$denitri,2] #-0.493 # per second
	
	# sensit_algaluptake[sindex,]
		algal_parameter_Mulholland_max = 10^sensit_algaluptake[sindex$algal,1]*0.01 #0.0598444/2/25000 #1.49611e-05 #unit: mgN/mgC/s 
		algal_parameter_Mulholland_pow = sensit_algaluptake[sindex$algal,2] #10*1000#100 # mgN/m3 AJ's data --> 10mg/L
		algal_parameter_1ms_left = 0.15

		detritus_parameter_Mulholland_max = 10^sensit_algaluptake[sindex$decomp,1]*0.01 #6.134259e-07 # mgN/mgC/s in Jack's (chapter)
		detritus_parameter_Mulholland_pow = sensit_algaluptake[sindex$decomp,2] #100 # mgN/m3
	    detritus_parameter_1ms_left = 0.15

	SimulationTime = seq.Date(from=as.Date('2012-1-1'), to=as.Date('2017-12-31') ,by="day")    ## second round
	timemax_BTS=24 * length(SimulationTime);
	timemax_STS=3600
	commonDates = intersectDate(list(SimulationTime, sensit_lateralQ$date, sensit_lateralN$date, sensit_lateralP$date))
    sensit_lateralQ.dtsMatch = match(commonDates, sensit_lateralQ$date)
    sensit_lateralN.dtsMatch = match(commonDates, sensit_lateralN$date)
	sensit_lateralP.dtsMatch = match(commonDates, sensit_lateralP$date)

####################################### reading inputs ###########################################	
	
	## -- building channel zone -- traveling time distribution	
	Morphology_matrix=list()	
	maxChannelBenthicArea = rep(NA, numZone)
	for(i in numZoneIndex){
		
		speedindex = sensit_velocity[sindex$velocity,c('poolV_order','riftV_order','upV_order')] 
		ratingTable = read.csv(paste(sensitivity$vPathway,'channelZone', zoneInfo[i,'zoneClass'],'_ratingTable', speedindex[zoneInfo[i,'zoneClass']],'.csv',sep=''))
		# ratingTable = read.csv(paste('channel/whole_channel_new/channelZone', zoneInfo[i,'zoneID'],'_ratingTable1.csv',sep=''))
		maxChannelBenthicArea[i] = max(ratingTable[,'wettedBperimeter'])*zoneInfo[i,'len']
		ratingTable$Corrected_volume = ratingTable$volume*0.1*zoneInfo[i,'len']
			#cbind(ratingTable$volume, ratingTable$Corrected_volume, ratingTable$crossarea*zoneInfo$len[i])
		
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
	}#i
	
	
	
	
	
	
	zoneEnv_matrix = array(0, dim=c(numZone, timemax_BTS,length(zoneEnv_variable)), dimnames=list(numZoneIndex, 1:timemax_BTS, zoneEnv_variable) ); 
	zoneEnv_matrix[,,'temperature'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_untemp[,sindex$PAR_temp] # upstream
		}else{
			return <- sensit_retemp[,sindex$PAR_temp] # downstream
		}
	}))	
	zoneEnv_matrix[,,'PARf']= t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- PAR2PARF(sensit_unPAR[,sindex$PAR_temp]) # upstream
		}else{
			return <- PAR2PARF(sensit_rePAR[,sindex$PAR_temp]) # downstream
		}
	}))

    blocking = seq(0.0003,0.003, length.out=11)
    upperblock = blocking[-1];
    lowerblock = blocking[-length(blocking)]
    exchangeRateList = sapply(seq_len(length(blocking)-1),function(i){ return <- seq(lowerblock[i],upperblock[i],length.out=20) })
    #zoneEnv_matrix_exchange = seq(0.0002,0.003, length.out=200)[sindex$exchange] # just a value with A/As (later)
    zoneEnv_matrix_exchange = exchangeRateList[sindex$exchange] # just a value with A/As (later)
    zoneEnv_matrix_exchangeAreaScale = c(seq(1.0/3,1, length.out=100),seq(1,3, length.out=101)[-1])[sindex$exchange]



	zoneEnv_matrix[,,'Nuptake']= t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_uptake # upstream	<<----?
		}else{
			return <- sensit_otherpattern$Down_Nuptake # downstream <<---?
		}
	}))







	indexlen = unlist(lapply(table(format(SimulationTime, '%Y')),function(g){seq_len(g)})); indexlen[indexlen==366]=365
	lateralResource_matrix = array(0, dim=c(numZone, timemax_BTS,length(lateralResource_variable)), dimnames=list(numZoneIndex,1:timemax_BTS, lateralResource_variable) ); 
    
    totalflow_m3s = sensit_lateralQ[sensit_lateralQ.dtsMatch,'Qlps1']*0.001 #<<------ here (new inputs)
	baseflow = fivedayblockbaseflow(totalflow_m3s)
	stormQ = sapply(totalflow_m3s-baseflow,function(x){return <- ifelse(x>0,x,0)})
		# plot(totalflow_m3s, type='l', log='y'); lines(baseflow, col='blue')
		
	lateralResource_matrix[,,'discharge'] = t(sapply(numZoneIndex, function(kk){ zoneInfo[kk,'lateralFrac']*rep(baseflow,each=24) }))
	lateralResource_matrix[,,'stormQ'] = t(sapply(numZoneIndex, function(kk){ zoneInfo[kk,'stormFrac']*rep(stormQ,each=24) }))
	
	lateralResource_matrix[,,'waterNO3'] = t(1000*sapply(numZoneIndex, function(kk){ 
		lateralResource_matrix[kk, ,'discharge']*rep(sensit_lateralN[sensit_lateralN.dtsMatch,'Nmgpl'],each=24)
	}))
	lateralResource_matrix[,,'waterNH4'] = 0*lateralResource_matrix[, ,'discharge']
	lateralResource_matrix[,,'waterPO4'] = t(1000*sapply(numZoneIndex, function(kk){
		lateralResource_matrix[kk, ,'discharge']*rep(sensit_lateralP[sensit_lateralP.dtsMatch,'Pmgpl'],each=24)
		# 1.936937*0.001*
	}))
	
	
	lateralResource_matrix[,,'stormNO3'] = t(1000*sapply(numZoneIndex, function(kk){ 
		lateralResource_matrix[kk,, 'stormQ']*rep(sensit_lateralN[sensit_lateralN.dtsMatch,'Nmgpl'],each=24)
	}))
	lateralResource_matrix[,,'stormNH4'] = 0*lateralResource_matrix[,, 'stormQ']
	lateralResource_matrix[,,'stormPO4'] = t(1000*sapply(numZoneIndex, function(kk){
		lateralResource_matrix[kk,, 'stormQ']*rep(sensit_lateralP[sensit_lateralN.dtsMatch,'Pmgpl'],each=24)
	}))#200 *lateralResource_matrix[,, 'stormQ'] 1.936937*0.001*
	
	lateralResource_matrix[,,'detritus1c'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_detritus1c # upstream
		}else{
			return <- sensit_otherpattern$Down_detritus1c # downstream
		}
	}))
	lateralResource_matrix[,,'detritus1n'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_detritus1n # upstream
		}else{
			return <- sensit_otherpattern$Down_detritus1n # downstream
		}
	}))
	
	lateralResource_matrix[,,'detritus2c'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_detritus2c # upstream
		}else{
			return <- sensit_otherpattern$Down_detritus2c # downstream
		}
	}))
	lateralResource_matrix[,,'detritus2n'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_detritus2n # upstream
		}else{
			return <- sensit_otherpattern$Down_detritus2n # downstream
		}
	}))
	
	lateralResource_matrix[,,'detritus3c'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_detritus3c # upstream
		}else{
			return <- sensit_otherpattern$Down_detritus3c # downstream
		}
	}))
	lateralResource_matrix[,,'detritus3n'] = t(sapply(numZoneIndex, function(kk){
		if(zoneInfo[kk,'zoneClass']==3){
			return <- sensit_otherpattern$Up_detritus3n # upstream
		}else{
			return <- sensit_otherpattern$Down_detritus3n # downstream
		}
	}))
			
	## here
rm(sensit_velocity)
rm(sensit_denitrification)
rm(sensit_algaluptake)
rm(sensit_lateralQ)
rm(sensit_lateralN)
rm(sensit_rePAR)
rm(sensit_retemp)
rm(sensit_unPAR)
rm(sensit_untemp)
rm(sensit_otherpattern)
		
		
	
		
####################################### initial ###########################################			
		iniDischarge = cumsum(lateralResource_matrix[,1,'discharge'] + lateralResource_matrix[, 1,'stormQ'] )
		iniWater = sapply(numZoneIndex,function(i){ min(Morphology_matrix[[i]]@Q2Vol(iniDischarge[i]),Morphology_matrix[[i]]@maxVol,na.rm=T) })	
		meanWater = sapply(numZoneIndex,function(i){ min(Morphology_matrix[[i]]@Q2Vol(2*0.001), Morphology_matrix[[i]]@maxVol,na.rm=T) }) # assume baseflow is 2 L/s
	
	hydraulic_matrix = matrix(0, numZone,length(hydraulic_variable)); colnames(hydraulic_matrix)= hydraulic_variable
		hydraulicHOLD=sapply(numZoneIndex,function(i){ c(
            Morphology_matrix[[i]]@Vol2Z(iniWater[i]), #1
            Morphology_matrix[[i]]@Vol2V(iniWater[i]), #2
            Morphology_matrix[[i]]@Vol2WP(iniWater[i]), #3
            Morphology_matrix[[i]]@Vol2Q(iniWater[i]), #4
            Morphology_matrix[[i]]@Vol2CA(iniWater[i]), # 5
            Morphology_matrix[[i]]@maxQ, #5 -> 6
            Morphology_matrix[[i]]@maxVol, # 7
            	# Morphology_matrix[[i]]@Vol2CA(meanWater[i])/Morphology_matrix[[i]]@Vol2Z(meanWater[i])+3 = width
            	# (1-exp(-p_d*1.5))*p0/p_d space below the 1.5 m
            (Morphology_matrix[[i]]@Vol2CA(meanWater[i])/Morphology_matrix[[i]]@Vol2Z(meanWater[i])+3) * (1-exp(-p_d*1.5))*p0/p_d,  # porosity
            Morphology_matrix[[i]]@Vol2Z(meanWater[i])
            ) })
		hydraulic_matrix[,'depth'] = hydraulicHOLD[1,]
		hydraulic_matrix[,'velocity'] = hydraulicHOLD[2,] #m/s
		hydraulic_matrix[,'wetBenthicArea'] = zoneInfo[,'len']*hydraulicHOLD[3,]
		hydraulic_matrix[,'wetBenthicArea_1'] = 1.0 / hydraulic_matrix[,'wetBenthicArea']
        hydraulic_matrix[,'crossarea'] = hydraulicHOLD[5,]
        hydraulic_matrix[,'maxVol'] = hydraulicHOLD[7,]*0.99
		hydraulic_matrix[,'discharge'] = hydraulicHOLD[4,]
		hydraulic_matrix[,'As'] = hydraulicHOLD[8,] * zoneEnv_matrix_exchangeAreaScale
		#cbind(hydraulic_matrix[,c('crossarea','As')], iniWater, meanWater)
		
####################################### in-channel process parameters ###########################################
	
	bioprocess_matrix = matrix(0, numZone,length(bioprocess_variable)); colnames(bioprocess_matrix)= bioprocess_variable
			zoneBioAssign = (zoneInfo[,'zoneClass']%%3); zoneBioAssign[zoneBioAssign==2]=0.0
			
		bioprocess_matrix[,'algalc'] = 10*1000*zoneBioAssign + 0.01; #50 gC/m2 as inintial; model unit is mg/m2
		bioprocess_matrix[,'algaln'] = bioprocess_matrix[,'algalc']*0.5*(algal_parameter['maxNC']+algal_parameter['minNC'])
		bioprocess_matrix[,'algalp'] = bioprocess_matrix[,'algalc']*0.5*(algal_parameter['maxPC']+algal_parameter['minPC'])
		bioprocess_matrix[,'algalCN'] = bioprocess_matrix[,'algalc']/bioprocess_matrix[,'algaln']
		bioprocess_matrix[,'algalCP'] = bioprocess_matrix[,'algalc']/bioprocess_matrix[,'algalp']
		bioprocess_matrix[,'algalselfF'] = 1
		bioprocess_matrix[,'algalnutrientF'] = 0
		bioprocess_matrix[,'algallightF'] = 1
		bioprocess_matrix[,'algalNuptakeF'] = 1
		bioprocess_matrix[,'algalPuptakeF'] = 1
		bioprocess_matrix[,'colonizedBenthicArea']=hydraulic_matrix[,'wetBenthicArea']*0.9*sapply((zoneInfo[,'zoneClass']%%3),function(x){ifelse(x>0,1,0.1)}) #m2
		bioprocess_matrix[,'colonizedBenthicArea_1'] = 1.0/bioprocess_matrix[,'colonizedBenthicArea']
	hyportheicResource_matrix = matrix(0, numZone,length(hyportheicResource_variable)); colnames(hyportheicResource_matrix)= hyportheicResource_variable
	benthicResource_matrix = matrix(0, numZone,length(benthicResource_variable)); colnames(benthicResource_matrix)= benthicResource_variable	
		
		
####################################### model computational settings	

	columnResource_matrix = matrix(0,numZone,length(columnResource_variable)); colnames(columnResource_matrix)= columnResource_variable
		columnResource_matrix[,'discharge'] = iniWater
		columnResource_matrix[,'waterNO3'] = (lateralResource_matrix[,1,'waterNO3']+lateralResource_matrix[, 1,'stormNO3'])/(lateralResource_matrix[,1,'discharge'] + lateralResource_matrix[, 1,'stormQ'])*iniWater
		columnResource_matrix[,'waterNO3'] = columnResource_matrix[,'waterNO3']
		columnResource_matrix[,'waterNH4'] = (lateralResource_matrix[,1,'waterNH4']+lateralResource_matrix[, 1,'stormNH4'])/(lateralResource_matrix[,1,'discharge'] + lateralResource_matrix[, 1,'stormQ'])*iniWater
		columnResource_matrix[,'waterPO4'] = (lateralResource_matrix[,1,'waterPO4']+lateralResource_matrix[, 1,'stormPO4'])/(lateralResource_matrix[,1,'discharge'] + lateralResource_matrix[, 1,'stormQ'])*iniWater
	hyporheicResource_matrix = matrix(0,numZone,length(columnResource_variable)+5); 
		colnames(hyporheicResource_matrix)=c(columnResource_variable,'baseNO3','baseNH4','baseQ','baseNO3balance','baseNH4balance')
		hyporheicResource_matrix[,'discharge'] = hydraulic_matrix[,'As'] * zoneInfo[,'len']
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
	
	
	## -------- holders / containers
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
        hold[] =  1 - (1-(1-algal_parameter_1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(1/24) # removal rate
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
                # algal_parameter_1ms_left is the x% of material remain within a day under velocity V -> removal rate = f(V)
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
                holdIII[] = ((1-detritus_parameter_1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(STSstep/86400)
                hold[] =  1 - (1-(1-detritus_parameter_1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(STSstep/86400) # removal rate
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
				Pot_algalupNO3[] =  algal_parameter_Mulholland_max*(concNO3 ^(1+ algal_parameter_Mulholland_pow)) * bioprocess_matrix[,'algalNuptakeF'] * algalQ10f * STSstep *4e-05 # correct by 25gC/m2
					tmpCond[] = Pot_algalupNO3>columnResource_matrix[,'waterNO3']
					Pot_algalupNO3[tmpCond] = columnResource_matrix[tmpCond,'waterNO3']
					# mgN/mgC/s * mgC/m2 * s  * m2 = mgN; algalNuptakeF = mgC/m2 * m2 * fraction
				       
				#Pot_algalupNH4[] = algal_parameter['umaxN']*(concNH4+concNO3)/(algal_parameter['Nhalf']+concNH4+concNO3)*bioprocess_matrix[,'algalNuptakeF']*algalQ10f * STSstep # mgN
				Pot_algalupNH4[] = algal_parameter_Mulholland_max*((concNH4+concNO3)^(1+ algal_parameter_Mulholland_pow)) *bioprocess_matrix[,'algalNuptakeF']*algalQ10f * STSstep *4e-05 
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
				Pot_immobNO3[] = detritus_parameter_Mulholland_max*(concNO3 ^(1+ detritus_parameter_Mulholland_pow)) * hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep *detritusMatrix[,'microbes'] *4e-05;
					tmpCond[] = Pot_immobNO3>columnResource_matrix[,'waterNO3']
					Pot_immobNO3[tmpCond] = columnResource_matrix[tmpCond,'waterNO3']
				Pot_immobNH4[] = detritus_parameter_Mulholland_max*((concNH4+concNO3)^(1+ detritus_parameter_Mulholland_pow)) * hydraulic_matrix[,'wetBenthicArea']*bioQ10f*STSstep * detritusMatrix[,'microbes'] *4e-05 - Pot_immobNO3;
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


