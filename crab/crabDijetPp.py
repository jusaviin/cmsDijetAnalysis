from WMCore.Configuration import Configuration
config = Configuration()

card = 'cardDijetPp.input'
infoString = 'dijet_pp_highForest_part3_2018-07-27'
output = infoString + '.root'

config.section_("General")
config.General.requestName = infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 400
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.userInputFiles = open('pp_5TeVdata_forest.txt').readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 10
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetTestPpHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

