library(methods)
## all mass unit in the model is "mgC" or "mgN"; all distance unit is metric, "m"
####################################### cross-section morphology
	# http://www.fsl.orst.edu/geowater/FX3/help/8_Hydraulic_Reference/Mannings_n_Tables.htm
	# https://www.engineeringtoolbox.com/mannings-roughness-d_799.html
	# http://xpsolutions.com/webhelp/Manning_s_N_Roughness_Coefficients_Open_Channel.htm
	# https://www.researchgate.net/publication/245287209_Estimation_of_Mean_Velocity_in_Natural_Channels_Based_on_Chiu%27s_Velocity_Distribution_Equation?enrichId=rgreq-0bfd7c5348e6a22d1dec8f49a33fb6a3-XXX&enrichSource=Y292ZXJQYWdlOzI0NTI4NzIwOTtBUzoxMDQxMjU4MTk3ODUyMTlAMTQwMTgzNjkzODM3MA%3D%3D&el=1_x_2&_esc=publicationCoverPdf
	# phi(M) = 1/0.6668 = 1.4997, --> M = 2.1475
	setClass('CrossSectionQrelation',slots=list(
		Q2Vol='function',
		Vol2Z='function', 
		Vol2V='function', 
		Vol2WP='function', 
		Vol2Q='function',
		Vol2CA='function',
		maxQ='numeric',
		maxVol='numeric') )
	ma <- function(x,n=5){filter(x,rep(1/n,n), sides=2)}
	

####################################### define functions and constants
	
	
	hydraulic_variable=c('velocity','wetBenthicArea','wetBenthicArea_1','len','depth','maxVol','discharge','crossarea','As')
	
	interested_variable =c('discharge','waterNO3','waterNH4','waterPO4') 
	columnResource_variable=c(interested_variable,'FOC','FON','FOP')
	zoneEnv_variable=c('temperature','PAR','PARf','exchangeProp2STR','exchangeProp2SRIP','Nuptake','baseNO3','baseQ','exchange')
	lateralResource_variable = c(interested_variable,c('stormQ','stormNO3','stormNH4','stormPO4'),c('detritus1c','detritus1n','detritus2c','detritus2n','detritus3c','detritus3n'))
	hyportheicResource_variable = c(interested_variable)
	
	benthicResource_variable = c('cbomC_fast','cbomN_fast','cbomP_fast','cbomC_med','cbomN_med','cbomP_med','cbomC_slow','cbomN_slow','cbomP_slow','FOC','FON','FOP')
		algal_variable = c(
			'algalc','algaln','algalp','algalCN','algalCP',
			'algalNPP','algalupP', 'algalupNO3','algalupNH4',
			'algalMineralC','algalMineralN','algalMineralP',
			'algalnutrientF','algalPF', 'algalNF','algalselfF','algallightF','algalNuptakeF', 'algalPuptakeF',
			'algalentC','algalentN','algalentP', 'colonizedBenthicArea','colonizedBenthicArea_1','resp','algaldeadC','algaldeadN','algaldeadP')
		biochem_variable = c('nitrification','denitrification')
		detritusName = c(
			'detritus1c','detritus1n','detritus1p',
			'detritus2c','detritus2n','detritus2p',
			'detritus3c','detritus3n','detritus3p', 'microbes', ## tracking variable
			'miner','growth_miner','decay_miner','decay_minerN','decay_minerP', # adding miner
			'Ncoef','Pcoef','netmineral','netmineralP','mineral','decay','wN','wP', 'foodNC','foodPC','foodNC_miner','foodPC_miner', ## temporaly; wN is uptake result
			'resp','respN', 'respP','age')
			
	bioprocess_variable = c(algal_variable, biochem_variable,'totaluptake','Sw','exchangeNO3','exchangeNH4','exchangePO4')
	
	
	outputTitleDownstreamExport = c('discharge','nitrate','ammonium','phosphate','FOC','FON','FOP')
	outputTitleInChannel = c('storedWater','storedNitrate','storedAmmonium','storedPhosphate')
	outputTitleWaterCond = c('velocity','depth','wetbenthicarea')#
	outputTitleBio = c('colonizedBenthicArea','algalc','algalCN','algalCP','algalnutrientF','algalselfF',
			'algalNPP', 'algalupNO3','algalupNH4','algalupP','algalMineralN','algalMineralP',
			'nitrification','denitrification')
	






	# Nconc = c(0.75,1.51,1.73,1.52,12.02,11.07) #mgN/L
    # uptake = c(0.013105182,0.005403471,0.016731141,0.052403754,0.35492126,0.822029703) #mgN/m2/s
	# dataset = as.data.frame(cbind(
		# 1/uptake,1/Nconc
	# )); colnames(dataset)=c('y','x')
    # # (1/U) = k/Vmax * (1/[c]) + 1/Vmax
	# result = lm(y~x,dataset);summary(result)
	# plot(dataset$x, dataset$y); abline(result)
	# Coefficients:
    # 				Estimate Std. Error t value Pr(>|t|)
	# (Intercept)    16.71      46.93   0.356    0.740
	# x              71.65      66.39   1.079    0.341 
	# --> Vmax = 1/16.71 = 0.0598444 mgN/m2/s; with 25000 mgC/m2  algal carbon standing stock =>  2.393776e-06 mgN/mgC/s
	# --> K = 71.65* 0.0598444 = 4.287851 mg/L = 4287.851 mg/m3
	algal_parameterName = c(
		'umaxN','Nhalf','umaxP','Phalf',
		'maxCN','maxNC','maxCP','maxPC','minCN','minNC','minCP','minPC','gcoef',
		'selflimitcoef','maxgrowthRate','quotaN','quotaP','mineralRate','Q10','moralityRate')
	algal_parameter = rep(NA,length(algal_parameterName)); names(algal_parameter) = algal_parameterName
	
		# 2.99222e-05 (using AJ's total uptake data; much higher than in Jack's chapter; half of it) --> 1.49611e-05(above);
		# 0.0530 mgN/mgC/d = 6.134259e-07 mgN/mgC/s in Jack's
	algal_parameter['umaxN'] = 1.49611e-05 #0.0598444/2/25000 #1.49611e-05 #unit: mgN/mgC/s 
	algal_parameter['Nhalf'] = 100 #10*1000 #100 # mgN/m3 AJ's data --> 10mg/L
		# 0.3/86400 #mgP/mgC/s [from daily] # 0.00731 mgP/mgC/d = 8.460648e-08 mgP/mgC/s in Jack's [no P data, so make it not limiting factor]
	algal_parameter['umaxP'] = 8.460648e-08 * 2
	algal_parameter['Phalf'] = 2 # mgP/m3 --> 2µgP/L
	
		## double check here to correct notes!!
		algal_parameter['maxCN'] = (1/0.0606) *12/14 # from molar ratio to mass ratio ~ 10.28 -- 14.14427 -- 24 cross et al. 2005; min 3.6857
		algal_parameter['maxNC'] = 1.0 / algal_parameter['maxCN']
		algal_parameter['maxCP'] = (1/0.00377) *12/31 # from molar ratio to mass ratio ~ 102.6782
		algal_parameter['maxPC'] = 1.0 / algal_parameter['maxCP']
		
		algal_parameter['minCN'] = (10.404) *12/14 # cross et al. 2005; 4.308, 10.404; min ~> 3.692571, 8.917714 #<<------
		algal_parameter['minNC'] = 1.0 / algal_parameter['minCN']
		algal_parameter['minCP'] = (123.817) *12/31 # cross et al. 2005; 25.536, 123.817; min ~ 9.884903, 47.92916
		algal_parameter['minPC'] = 1.0 / algal_parameter['minCP']
	
	algal_parameter['gcoef'] = 0.8
	algal_parameter['selflimitcoef'] = 0.0015 # m2/mgC from Jack's 667mg/m2 ;  1/(1+#*c); so # is about 1% of maxC --> 1/101 -> 0.00990099
	algal_parameter['maxgrowthRate'] = 2*2*1.0/86400 # 1/s [from Jack's]; [count for daily to hourly; half day is in dark; convert NPP to GPP] 
	algal_parameter['mineralRate'] = 0.01/86400 # 1/s from Jack's
	algal_parameter['Q10'] = 2 # from Jack's
	algal_parameter['moralityRate'] = 0.004*algal_parameter['maxgrowthRate'] #0.02/86400 = 0.005*algal_parameter['maxgrowthRate'] (1/s) from Jack's -->> turn into mortality rate; this is very sensitivity parameter to control algal mass
	 
	#biochem_parameter = c(0.0002339181, 1.059254e-05, -0.493, 2)# per second 
	biochem_parameter = c(0.00002339181/17, 1.059254e-05, -0.493, 2)# per second 
	names(biochem_parameter)=c('max_nitrif','Mulholland_max_denitr','Mulholland_pow_denitr',"Q10")	
		# Mulholland et al. 2008
		# log10(vf_cmps) = -0.493 * log10(ugNpL) - 2.975; cmps = cm per second
		# vf_cmps = (ugNpL)^-0.493 * 10^(-2.975) = 0.001059254 * ugNpL^-0.493
		# cm/s = m/s * 1/100
		# vf_cmps = vf_mps * 1/100
		# vf_mps = 0.00001059254 * (ugNpL)^-0.493 (convert unit from cmps to mps)
		# ugNpL = mgN/m3
		# vf_mps * mgN/m3 = mgN/m2/s = areal_den
		# areal_den = mgN/m2/s = mgN/m3 * 0.00001059254 * (mgN/m3)^-0.493 = 0.00001059254 * (mgN/m3)^(1-0.493)
		# e.g., 80 mgN/m2/day -> 80/3600/24 mg/m2/s -> 80/3600/24/(2000) m/s = 4.62963e-07 = 4.62963e-05 cm/s

	detritus_parameter = list()
	#detritus_parameter$growthRate = as.numeric(algal_parameter['maxgrowthRate']) #1.925427e-05 # unit = 1/s; double by 10 hours (1+g)^(10*3600) A = 2A -> g = exp(log(2)/(10*3600))-1
	detritus_parameter$growthRate = 4.456975e-06   	   # exp(log(2)/(1.8*24*3600))-1				#immob_growth=log(2)/(1.8*3600*24) # 4.456965e-06 1/s
	detritus_parameter$growthRate_miner = 3.342729e-06 # exp(log(2)/(2.4*24*3600))-1				#miner_growth=log(2)/(2.4*3600*24) # 3.342724e-06 1/s
	#detritus_parameter$deadRate = 0.02*detritus_parameter$growthRate
	detritus_parameter$baseR = 1.16e-7 # my manuscript #<---------------- maintain respiration
	detritus_parameter$cn = 7.0 #<-------------------------- cross et al. 2005
	detritus_parameter$nc = 1.0 / detritus_parameter$cn
	detritus_parameter$cp = 188
	detritus_parameter$pc = 1.0 / detritus_parameter$cp
		detritus_parameter$miner_cn = 5.0 #<-------------------------- cross et al. 2005
		detritus_parameter$miner_nc = 1.0 / detritus_parameter$miner_cn
		detritus_parameter$miner_cp = 20
		detritus_parameter$miner_pc = 1.0 / detritus_parameter$miner_cp
	detritus_parameter$gcoef = 0.5
	detritus_parameter$umaxN = 6.134259e-07 # mgN/mgC/s in Jack's (chapter)
	detritus_parameter$Nhalf = 100 # mgN/m3
	detritus_parameter = unlist(detritus_parameter)


	source("https://raw.githubusercontent.com/laurencelin/Date_analysis/master/LIB_dailytimeseries3.R")
	#SimulationTime = seq.Date(from=as.Date('2009-1-1'), to=as.Date('2017-12-31') ,by="day")	
	SimulationTime = seq.Date(from=as.Date('2012-1-1'), to=as.Date('2017-12-31') ,by="day")	## second round
	timemax_BTS=24 * length(SimulationTime); 
	timemax_STS=3600
	
	## -- building channel zone -- traveling time distribution
	channelInfo = read.csv('channel/GIS_channel_lateral_full_extend2.csv') # 1.01484728
	zoneInfo = t(simplify2array(tapply(seq_len(dim(channelInfo)[1]), INDEX=channelInfo[,'zoneID'],function(ii){
		return <- c(
			mean(channelInfo[ii,'zoneID']),
			mean(channelInfo[ii,'zoneClass']),
			#max(channelInfo[ii,'accDist'])-min(channelInfo[ii,'accDist']), # zone length
			sum(channelInfo[ii,'dist']), # zone length
			sum(channelInfo[ii,'laterFrac']), # zone lateral input fraction
			0
		)
	}))); colnames(zoneInfo)=c('zoneID','zoneClass','len','lateralFrac','stormFrac')
	numZone = dim(zoneInfo)[1]
	numZoneIndex = seq_len(numZone)
	
		classcond = zoneInfo[,'zoneClass']==3
	zoneInfo[classcond,'stormFrac'] = 0.9* zoneInfo[classcond,'len'] / sum(zoneInfo[classcond,'len'])
		classcond = ! classcond; classcond[9]=F
	zoneInfo[classcond,'stormFrac'] = 0.05* zoneInfo[classcond,'len'] / sum(zoneInfo[classcond,'len'])	
		classcond = rep(F,numZone); classcond[zoneInfo[,'zoneID']==9]=T
	zoneInfo[classcond,'stormFrac'] = 0.05
	
	
	#### previous before AGU
		# outputTitle = c(
		# 'year','month','day',
		# outputTitleDownstreamExport,
		# outputTitleWaterCond,
		# outputTitleBio,
		# paste('sw',zoneInfo[,1],sep=''),
		# paste('exNO3_',zoneInfo[,1],sep=''),
		# paste('no3_',zoneInfo[,1],sep='')
		# )
		
	#### after AGU
		# outputTitle = c(
		# 'year','month','day',
		# outputTitleDownstreamExport,
		# outputTitleWaterCond,
		# outputTitleBio,
		# paste('ca',zoneInfo[,1],sep=''),
		# paste('zoneQin_',zoneInfo[,1],sep=''),
		# paste('zoneNO3in_',zoneInfo[,1],sep=''),
		# paste('zoneNH4in_',zoneInfo[,1],sep=''),
		# paste('zoneFONin_',zoneInfo[,1],sep=''),
		# paste('zoneNO3out_',zoneInfo[,1],sep=''),
		# paste('zoneNH4out_',zoneInfo[,1],sep=''),
		# paste('zoneFONout_',zoneInfo[,1],sep='')
		# )
		
		propfluxTitle = c(
			'benthicArea','velocity','depth', 
			'nitrification','denitrification',
			'algalNO3uptake','algalNH4uptake','algalNPP','mineralN','algalCN','algalCP','algalC','algalCloss',
			'netmineral','detric','detriCN','detriCP','microbC','minerC','algalR','detritusR',
			'Sw','tmperature','PARf','hdmg','algalNlimit','algalPlimit','rpm')
		
		outputTitle = c(
		'year','month','day',
		paste('UpIn',columnResource_variable,sep='_'), 
		paste('UpOut',columnResource_variable, sep='_'),
        paste('Up',propfluxTitle,sep='_'), #<<----------------------
		paste('DownIn',columnResource_variable, sep='_'),
		paste('DownOut',columnResource_variable, sep='_'),
		paste('Down',propfluxTitle,sep='_'), #<<----------------------
		'upexchange','downexchange',
		paste('no3',1:numZone,sep='_')
		)
		
		

		
		
		
	
