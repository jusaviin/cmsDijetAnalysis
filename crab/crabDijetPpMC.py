from WMCore.Configuration import Configuration
config = Configuration()

correlationType = 'GenGen'
infoString = 'mergedSkims_Pythia6_2018-07-27_part1'
inputFile='mergedSkimPpPythia5TeV.txt'
card='cardDijetPpMC.input'
output='dijet_ppMC_'+correlationType+'_'+infoString+'.root'

config.section_("General")
config.General.requestName = 'ppMC_'+correlationType+'_'+infoString
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [output]
config.JobType.maxJobRuntimeMin = 800
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.userInputFiles = open(inputFile).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetTestPpMCHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T3_US_FNALLPC']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

