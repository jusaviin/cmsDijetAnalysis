from WMCore.Configuration import Configuration
config = Configuration()

correlationType = 'GenGen'
infoString = 'Pythia8_pfJets_wtaAxis_jetsNtracks_jetClosures_JECv3_2019-09-03'
inputFile='ppMC2017_mergedFiles.txt'
card='cardDijetPpMC2017.input'
output='ppMC2017_'+correlationType+'_'+infoString+'.root'
fileLocation='0'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

config.section_("General")
config.General.requestName = 'ppMC2017_'+correlationType+'_'+infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 180
config.JobType.maxMemoryMB = 1500

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 20
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetTestPpMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue','T2_US_Wisconsin','T2_US_Nebraska']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

