from WMCore.Configuration import Configuration
config = Configuration()

system = 'RecoGen'
infoString = 'pfCsJets_noUncOrInc_5eveStrictMix_xj_2019-06-12_lastFailed'
card='cardDijetPbPbMC.input'
output='PbPbMC_' + system + '_' + infoString + '.root'
inputFile='lastFailedRecoGen.txt'
fileLocation='0'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

config.section_("General")
config.General.requestName = 'PbPbMC_' + system + '_' + infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 500
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetPbPbMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

