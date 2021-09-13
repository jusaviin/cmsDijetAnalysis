from WMCore.Configuration import Configuration
config = Configuration()

correlationType = 'RecoGen'
filledHistograms = 'onlyRegular'
infoString = 'akCaloJet_' + filledHistograms + '_4pCentShift_subeNon0_noMixing_jet80Trigger_2021-09-08'
card='cardDijetPbPbMCNoMixing.input'
#card='cardDijetPbPbMCNoMixingClosure2018.input'
output='PbPbMC2018_' + correlationType + '_' + infoString + '.root'
#inputFile='PbPbMC2018_commonForestSkim.txt'
#inputFile='newPbPbMcForest.txt'
inputFile='newPbPbMcSkim.txt'
#inputFile='PbPbMC2018_sampleWithEventPlane.txt'  # No flow subtraction
#inputFile='PbPbMC2018_onlyJetsAndEventPlane_2021-01-26.txt' # Flow subtraction
fileLocation='1'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

config.section_("General")
config.General.requestName = 'PbPbMC2018_' + correlationType + '_' + infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 200
config.JobType.maxMemoryMB = 2000

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetPbPbMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
#config.Site.whitelist = ['T2_US_Purdue','T2_US_Wisconsin','T2_US_Nebraska']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

