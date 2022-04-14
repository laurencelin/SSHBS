################################################################################################################
# SSHBS_library.R contains all library parameters
# SSHBS_initiation.R  reads in channel morphology & hydraulic inputs and constructs model channels and hydraulic profile
# SSHBS_simulation.R (this file) starts stream ecosystem and flow simulations in the constructed and initiated model channels 
################################################################################################################
source('SSHBS_initiation.R')
#
################################################################################################################
outputFile = './test.csv'


outputTitleDownstreamExport = c('discharge','nitrate','ammonium','phosphate','FOC','FON','FOP')
outputTitleInChannel = c('storedWater','storedNitrate','storedAmmonium','storedPhosphate')
outputTitleWaterCond = c('velocity','depth','wetbenthicarea')#
outputTitleBio = c('colonizedBenthicArea','algalc','algalCN','algalCP','algalnutrientF','algalselfF',
                   'algalNPP', 'algalupNO3','algalupNH4','algalupP','algalMineralN','algalMineralP',
                   'nitrification','denitrification')

BTSaggregateVarPrint = c(
	'nitrification','denitrification','algalNO3uptake',
	'algalNH4uptake','algalNPP','mineralN',
	'algalCN','algalCP','algalC','algalCloss',
	'netmineral','detric','detriCN','detriCP','microbC','minerC',
	'algalR','detritusR','Sw','tmperature','PARf','hdmg',
	'algalNlimit','algalPlimit')
propfluxTitle = c(
  'benthicArea','velocity','depth', 
  'nitrification','denitrification',
  'algalNO3uptake','algalNH4uptake','algalNPP','mineralN','algalCN','algalCP','algalC','algalCloss',
  'netmineral','detric','detriCN','detriCP','microbC','minerC','algalR','detritusR',
  'Sw','tmperature','PARf','hdmg','algalNlimit','algalPlimit','rpm')


## ----- model time aggregators 
BTSaggregateVar = matrix(0,nrow=numZone, ncol=length(propfluxTitle)); colnames(BTSaggregateVar)=propfluxTitle
STSaggregateVar = matrix(0,nrow=numZone, ncol=length(outputTitleDownstreamExport)); colnames(STSaggregateVar)= outputTitleDownstreamExport
crossSection = rep(0, numZone)

## ----- output zone aggregators
zoneOutputAgg = list(
	# in this example, we aggregate zones into two output units
	outputUnit1 = c(1,2,3,4), 
	outputUnit2 = c(5,6,7,8)
) 

zoneOutputAggColname = lapply(seq_along(zoneOutputAgg), function(ii){
	xx = names(zoneOutputAgg)[ii]
	c(
		sprintf('%s_In_%s',xx,columnResource_variable),
		sprintf('%s_Out_%s',xx,columnResource_variable),
		sprintf('%s_%s',xx,propfluxTitle),
		sprintf('%s_%s',xx,'exchange'),
		sprintf('%s_z%s_%s',xx,zoneOutputAgg[[ii]],'NO3in'),
		sprintf('%s_z%s_%s',xx,zoneOutputAgg[[ii]],'NO3out')
	)
})


## ----- set output file and write out column names
outputFile_buff <- file(outputFile,'w')
cat( c('year','month','day', do.call(c,zoneOutputAggColname)), '\n', file=outputFile_buff,sep=',')
	
################################################################################################################
system.time({

	for(BTS in 1:timemax_BTS ){	 
		
		## system check
		if(sum(detritusMatrix[,'detritus1c']<0)>0){print(BTS);break}
		
		exportProp = hydraulic_matrix[,'discharge']/columnResource_matrix[,'discharge'];
		STSstep = max(1, apply(exportProp %o% STSstepOption<0.25,2,prod)* STSstepOption)	
		
    	algalQ10f = (algal_parameter['Q10']^(0.1*zoneRes_matrix[BTS,'temp',]-0.2) ) ## Jack's equation for algae # version 8 testing
		bioQ10f = (biochem_parameter['Q10']^(0.1*zoneRes_matrix[BTS,'temp',]-0.2) ) ## Jack's equation for microbes
		
        ##-------------------------------------------------------------------------------------------------------------------------------------------------------------------##
		## algal death (biological with some velocity related)
        hold[] = algalQ10f * algal_parameter['moralityRate'] * timemax_STS # moralityRate = death rate inthe modle now
		bioprocess_matrix[,'algaldeadC'] <- bioprocess_matrix[,'algalc'] * hold[]  
		bioprocess_matrix[,'algaldeadN'] <- bioprocess_matrix[,'algaln'] * hold[] 
		bioprocess_matrix[,'algaldeadP'] <- bioprocess_matrix[,'algalp'] * hold[]
       
        ## wash out; 
        hold[] =  1 - (1-(1-algal_parameter_1ms_left)/(1+exp(-(hydraulic_matrix[,'velocity']-0.4)/0.06)))^(1/24) # removal rate
        hold[hydraulic_matrix[,'velocity']<0.1]=0
		bioprocess_matrix[,'algalentC'] <- bioprocess_matrix[,'algalc'] * hold[]  
		bioprocess_matrix[,'algalentN'] <- bioprocess_matrix[,'algaln'] * hold[] 
		bioprocess_matrix[,'algalentP'] <- bioprocess_matrix[,'algalp'] * hold[] 

		## damaged algal recovery (numerical approximation; rationale: 
		## a) algal community reestablishment; b) recolonization on benthic surface; c) regrowth phase under stress
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
		detritusMatrix[,'detritus1c'] = detritusMatrix[,'detritus1c'] + 
            zoneRes_matrix[BTS,'detritus1c',]*0.04166667 +
            hold*0.5 + bioprocess_matrix[,'algaldeadC']*0.5 # mgN/m2 unit area
		detritusMatrix[,'detritus1n'] = detritusMatrix[,'detritus1n'] +
            zoneRes_matrix[BTS,'detritus1n',]*0.04166667 +
			holdIV*0.5 + bioprocess_matrix[,'algaldeadN']*0.5 # mgN/m2 unit area
		detritusMatrix[,'detritus1p'] = detritusMatrix[,'detritus1p'] + 
            zoneRes_matrix[BTS,'detritus1c',]*0.04166667/detritusCP +
			holdV*0.5 + bioprocess_matrix[,'algaldeadP']*0.5 # mgN/m2 unit area	
            
		detritusMatrix[,'detritus2c'] = detritusMatrix[,'detritus2c'] + 
            zoneRes_matrix[BTS,'detritus2c',]*0.04166667 +
			hold*0.5 + bioprocess_matrix[,'algaldeadC']*0.5
		detritusMatrix[,'detritus2n'] = detritusMatrix[,'detritus2n'] + 
            zoneRes_matrix[BTS,'detritus2n',]*0.04166667 +
			holdIV*0.5 + bioprocess_matrix[,'algaldeadN']*0.5
		detritusMatrix[,'detritus2p'] = detritusMatrix[,'detritus2p'] + 
            zoneRes_matrix[BTS,'detritus2c',]*0.04166667/detritusCP +
			holdV*0.5 + bioprocess_matrix[,'algaldeadP']*0.5

		detritusMatrix[,'detritus3c'] = detritusMatrix[,'detritus3c'] +
            zoneRes_matrix[BTS,'detritus3c',]*0.04166667
		detritusMatrix[,'detritus3n'] = detritusMatrix[,'detritus3n'] + 
            zoneRes_matrix[BTS,'detritus3n',]*0.04166667
		detritusMatrix[,'detritus3p'] = detritusMatrix[,'detritus3p'] + 
            zoneRes_matrix[BTS,'detritus3c',]*0.04166667/detritusCP
        
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
				))
			);
				
			
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
            
            
      
    ##---------------------------------------------------------------------##
		STS = 0
		while(STS < timemax_STS){ #1:timemax_STS
			
			tmpCond[] = zoneRes_matrix[BTS,'temp',]>0 ## freezing issue
			## export
			exportProp = hydraulic_matrix[,'discharge']/columnResource_matrix[,'discharge']; exportProp[!tmpCond]=0
			STSexport = columnResource_matrix * exportProp * STSstep
			STSaggregateVar = STSaggregateVar + STSexport
			columnResource_matrix[] = columnResource_matrix[] * (1-exportProp * STSstep)
			
								
			## input (per second)
			STSinput[tmpCond, interested_variable] = #STSinput[tmpCond, interested_variable] + 
				t(zoneRes_matrix[BTS,c('baseflow','baseflowNO3','baseflowNH4','baseflowPO4'),tmpCond] * STSstep + 
				zoneRes_matrix[BTS,c('stormQ','stormQNO3','stormQNH4','stormQPO4'),tmpCond] * STSstep) #lateral/spring
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
			
					
			## overflow to riparian	
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
			hyporheicResource_matrix[,'waterNO3'] = hyporheicResource_matrix[,'baseNO3'] + hyporheicResource_matrix[,'baseNO3balance']
			hyporheicResource_matrix[,'waterNH4'] = hyporheicResource_matrix[,'baseNH4'] + hyporheicResource_matrix[,'baseNH4balance']
			hyporheicResource_matrix[,'baseNO3'] = hyporheicResource_matrix[,'waterNO3']
			hyporheicResource_matrix[,'baseNH4'] = hyporheicResource_matrix[,'waterNH4'] 
			
			hypoconcNO3 = hyporheicResource_matrix[,'waterNO3']/hyporheicResource_matrix[,'discharge'] ##<<----- come back here
			hypoconcNH4 = hyporheicResource_matrix[,'waterNH4']/hyporheicResource_matrix[,'discharge'] ##<<----- come back here
			hypoconcPO4 = hyporheicResource_matrix[,'waterPO4']/hyporheicResource_matrix[,'discharge'] ##<<----- come back here
			tmpCond[] = concNO3 > hypoconcNO3;
			denitrify_concNO3[tmpCond] = concNO3[tmpCond]; denitrify_concNO3[!tmpCond]=hypoconcNO3[!tmpCond]
			

			## potential uptake in hyporheic zone
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
			bioQ10f = (biochem_parameter['Q10']^(0.1*zoneRes_matrix[BTS,'temp',]-1.2) ) ## Jack's equation
			
			# ## mgN *1/s *s = mgN
			#Pot_nitrification[] = biochem_parameter['max_nitrif'] * bioQ10f * 0.5*(columnResource_matrix[,'waterNH4']+hydraulic_matrix[,'wetBenthicArea']*hyporheicResource_matrix[,'waterNH4']) * STSstep 
			Pot_nitrification[] = biochem_parameter_max_nitrif * bioQ10f * 0.5*(columnResource_matrix[,'waterNH4']+hyporheicResource_matrix[,'waterNH4']) * STSstep
			Pot_nitrificationFrac = (0.5*Pot_nitrification - hyporheicResource_matrix[,'waterNH4'])/Pot_nitrification
			Pot_nitrificationFrac[Pot_nitrificationFrac<0] = 0
			Pot_nitrificationFrac = Pot_nitrificationFrac+0.5 
					 
			## immobilization --  mgN/m2/s* m2 * s = mgN
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
				 
			## denitrification
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
 
            ## hyporheic exchange and process ******* 
			hold[] = alpha_exchange * STSstep; hold[hold>1]=1; 
			holdII[] = hold[] * hydraulic_matrix[,'crossarea'] / hydraulic_matrix[,'As'] 
			
			concMassChange[] = 0
			concChange[] = (hypoconcNO3[] - concNO3[]) #*zoneRes_matrix[tmpCond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
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
			concChange[] = (hypoconcNH4[] - concNH4[]) #*zoneRes_matrix[tmpCond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
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
			concChange[] = (hypoconcPO4[] - concPO4[]) #*zoneRes_matrix[tmpCond,BTS,'exchangeRate']/hydraulic_parameter #row=zone; col=solute
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
				
			hyporheicResource_matrix[,'baseNH4balance'] = hyporheicResource_matrix[,'waterNH4']-hyporheicResource_matrix[,'baseNH4']	
			hyporheicResource_matrix[,'baseNO3balance'] = hyporheicResource_matrix[,'waterNO3']-hyporheicResource_matrix[,'baseNO3']  
			tmpCond = hyporheicResource_matrix[,'baseNO3balance']>0
			hold[] = 0
			if(sum(tmpCond)>0){	
				# uptake from terrestrial tree root in hyporheic zone
				hold[tmpCond] = vectorMin(
					zoneRes_matrix[BTS,'uptake',tmpCond]*zoneInfo$len[tmpCond], 
					hyporheicResource_matrix[tmpCond,'baseNO3balance']*hydraulic_matrix[tmpCond,'wetBenthicArea_1'] 
				)#min
				hyporheicDenitrification[tmpCond] = hyporheicDenitrification[tmpCond] + hold[tmpCond]
				hyporheicResource_matrix[tmpCond,'baseNO3balance'] = hyporheicResource_matrix[tmpCond,'baseNO3balance'] - hold[tmpCond]*hydraulic_matrix[tmpCond,'wetBenthicArea'] # consistent to the total mass assumption
			}# if
                	
            hold = Pot_algalupNO3*bioprocess_matrix[,'colonizedBenthicArea']*hydraulic_matrix[,'wetBenthicArea_1'] + Pot_denitrification + Pot_immobNO3 - Pot_nitrification # <<----------this is all NO3
			bioprocess_matrix[,'totaluptake'] = bioprocess_matrix[,'totaluptake'] + hold #areal uptake
			hold = columnResource_matrix[,'waterNO3']/columnResource_matrix[,'discharge']*hydraulic_matrix[,'depth']*hydraulic_matrix[,'velocity']/hold
			bioprocess_matrix[,'Sw'] = bioprocess_matrix[,'Sw'] + hold*STSstep # uptake length
			
			###################################################################################### update hydraulic parameters
			crossSection = crossSection + hydraulic_matrix[,'crossarea']
			
			hydraulicHOLD=sapply(numZoneIndex,function(i){ c(
				depth = Morphology_matrix[[i]]@Vol2Z(columnResource_matrix[i,'discharge']), #1
				velocity = Morphology_matrix[[i]]@Vol2V(columnResource_matrix[i,'discharge']), #2
				waterPerimeter = Morphology_matrix[[i]]@Vol2WP(columnResource_matrix[i,'discharge']), #3
				discharge = Morphology_matrix[[i]]@Vol2Q(columnResource_matrix[i,'discharge']), #4
				crossArea = Morphology_matrix[[i]]@Vol2CA(columnResource_matrix[i,'discharge']), # 5
				maxDischarge = Morphology_matrix[[i]]@maxQ, #5 -> 6
				maxWater = Morphology_matrix[[i]]@maxVol ) }) #6 -> 7
			hydraulic_matrix[,'depth'] = hydraulicHOLD[1,]
			hydraulic_matrix[,'velocity'] = hydraulicHOLD[2,] #m/s
			hydraulic_matrix[,'wetBenthicArea'] = zoneInfo$len*hydraulicHOLD[3,]
			hydraulic_matrix[,'wetBenthicArea_1'] = 1.0 / hydraulic_matrix[,'wetBenthicArea']
			hydraulic_matrix[,'crossarea'] = hydraulicHOLD[5,]
			hydraulic_matrix[,'discharge'] = hydraulicHOLD[4,]

			#print(paste(c(BTS,STS, columnResource_matrix[,'discharge']),collapse=' '))
			if(is.na(sum(columnResource_matrix[,'discharge']))|sum(columnResource_matrix[,'discharge']<0)>0) 
				print(paste(c(BTS,STS, STSstep,round(abs(hydraulic_matrix[,'discharge'])/columnResource_matrix[,'discharge']*STSstep,6)),collapse=' '))
				
			if(sum(hyporheicResource_matrix[,'waterNH4']<0)>0)
				print(paste(c(BTS,STS, STSstep, hyporheicResource_matrix[,'waterNH4']),collapse=' ')) 

			## reset STS	
			STSinput[ ] = 0	
			STS = STS + STSstep
			exportProp = hydraulic_matrix[,'discharge']/columnResource_matrix[,'discharge'];
			STSstep = min(timemax_STS-STS,max(1,apply(exportProp %o% STSstepOption<0.25,2,prod)* STSstepOption))
		}# end of STS loop
		
		## back to BTS loop
		bioprocess_matrix[,'Sw'] = bioprocess_matrix[,'Sw']* 0.0002777778 ## correcting aggregating
				
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
		bioprocess_matrix[,'algallightF'] <- zoneRes_matrix[BTS,'PARf',] #*(0.9+0.1*exp(-0.1666667*columnResource_matrix[,'FOC']/columnResource_matrix[,'discharge']));
        bioprocess_matrix[,'algalNPP'] <- bioprocess_matrix[,'algalc']* algal_parameter['maxgrowthRate'] * bioprocess_matrix[,'algalnutrientF'] * bioprocess_matrix[,'algalselfF'] * bioprocess_matrix[,'algallightF'] * algalQ10f * timemax_STS * (1-algalhabitatDmgIndex)
		bioprocess_matrix[,'resp'] <- bioprocess_matrix[,'algalNPP']*(1-algal_parameter['gcoef']) + bioprocess_matrix[,'algalMineralC']
		bioprocess_matrix[,'algalc'] <- bioprocess_matrix[,'algalc'] + algal_parameter['gcoef']*bioprocess_matrix[,'algalNPP'] 
		bioprocess_matrix[,'algaln'] <- bioprocess_matrix[,'algaln'] + bioprocess_matrix[,'algalupNO3'] + bioprocess_matrix[,'algalupNH4']
		bioprocess_matrix[,'algalp'] <- bioprocess_matrix[,'algalp'] + bioprocess_matrix[,'algalupP']
		
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
        
		BTSaggregateVar[,'tmperature'] = BTSaggregateVar[,'tmperature'] + zoneRes_matrix[BTS,'temp',] ## need to be /24
		BTSaggregateVar[,'PARf'] = BTSaggregateVar[,'PARf'] + zoneRes_matrix[BTS,'PARf',] ## need to be /24
		BTSaggregateVar[,'hdmg'] = BTSaggregateVar[,'hdmg'] + algalhabitatDmgIndex
        
        BTSaggregateVar[,'algalNlimit'] = BTSaggregateVar[,'algalNlimit'] + bioprocess_matrix[,'algalNF']
        BTSaggregateVar[,'algalPlimit'] = BTSaggregateVar[,'algalPlimit'] + bioprocess_matrix[,'algalPF']
        
		## reset
		bioprocess_matrix[,'algalupP']=0;
		bioprocess_matrix[,'algalupNO3']=0;
		bioprocess_matrix[,'algalupNH4']=0;
		bioprocess_matrix[,'nitrification']=0
		bioprocess_matrix[,'denitrification']=0
		detritusMatrix[,'wN'] = 0 ## N immobilization
        detritusMatrix[,'wP'] = 0
		
		## forming output
		if(BTS%%24==0){

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
				zoneRes_matrix[BTS,'year',1],
				zoneRes_matrix[BTS,'month',1],
				zoneRes_matrix[BTS,'day',1], #3
               
				do.call(c,lapply(zoneOutputAgg,function(ii){
					return <- c(
						# input (to the whole stream)
						#STSinputTrack # from upstream (accumulative) [zone, variable]
						#STSLateralinputTrack # from lateral sides [zone, variable]
						STSinputTrack[ii[1],] + colSums(STSLateralinputTrack[ii[-1],]),
						# output to downstream
						STSaggregateVar[ii[length(ii)], ],
						# channel physical
						sum(BTSaggregateVar[ii,'benthicArea']),
						(zoneInfo$len[ii]/sum(zoneInfo$len[ii]) ) %*% BTSaggregateVar[ii,c('velocity','depth')],
						# ecology
						(BTSaggregateVar[ii,'benthicArea']/sum(BTSaggregateVar[ii,'benthicArea'])) %*% BTSaggregateVar[ii, BTSaggregateVarPrint],
						# nitrate reduction
						100*(1-STSaggregateVar[ii[length(ii)],'nitrate']/(STSinputTrack[ii[1],'waterNO3']+sum(STSLateralinputTrack[ii[-1],'waterNO3'])))/sum(zoneInfo$len[ii]),
						# nitrate exchange to hyporheic 
						sum(bioprocess_matrix[ii,'exchangeNO3']* (BTSaggregateVar[ii,'benthicArea']/sum(BTSaggregateVar[ii,'benthicArea'])) ),
						# nitrate input (up and side) by each zone within the output aggregated unit
						STSinputTrack[ii[1],'waterNO3']+c(0,STSLateralinputTrack[ii[-1],'waterNO3']),
						# nitrate output by each zone within the output aggregated unit
						STSaggregateVar[ii,'nitrate']
					)# return 42 variables
				}))
                ), '\n', file=outputFile_buff,sep=',')

            
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
			
			# progresss message
			print(sprintf('%s: %s-%s-%s %02d:00 is done',
				BTS, 
				zoneRes_matrix[BTS,'year',1],
				zoneRes_matrix[BTS,'month',1],
				zoneRes_matrix[BTS,'day',1],
				zoneRes_matrix[BTS,'hour',1]))# print
		}		
		#if(BTS%%720==0) print(paste(BTS,"is done"))#monthly

	}# end of BTS loop
	close(outputFile_buff)
})  

