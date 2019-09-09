from WMCore.Configuration import Configuration
config = Configuration()

system = 'RecoGen'
infoString = 'skims_pfJets_noUncorr_5moreMixedEvents_xj_2019-03-22_combine4_failed'
card='cardDijetPbPbMCMoreMixing.input'
output='PbPbMC_' + system + '_' + infoString + '.root'
inputFile='PbPbMC_Pythia6HydjetCymbal_list_combine04_toResubmit.txt'
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
config.JobType.maxJobRuntimeMin = 1440
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
config.Site.whitelist = ['T2_US_Purdue','T2_US_Wisconsin','T2_US_Nebraska']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

