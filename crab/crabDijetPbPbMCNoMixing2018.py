from WMCore.Configuration import Configuration
config = Configuration()

correlationType = 'RecoGen'
filledHistograms = 'noUncorr'
infoString = 'akPu4CaloJets_' + filledHistograms + '_noMixing_subeNon0_wtaAxis_JECv5b_2019-09-08'
card='cardDijetPbPbMCNoMixing.input'
output='PbPbMC_' + correlationType + '_' + infoString + '.root'
inputFile='PbPbMC2018_commonForest.txt'
fileLocation='1'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

config.section_("General")
config.General.requestName = 'PbPbMC_' + correlationType + '_' + infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 600
config.JobType.maxMemoryMB = 2200

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetPbPbMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_*']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

