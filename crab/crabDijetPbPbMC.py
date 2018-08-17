from WMCore.Configuration import Configuration
config = Configuration()

system = 'RecoReco'
infoString = 'skims_allButInclusive_2018-08-15_part6'
card='cardDijetPbPbMC.input'
output='PbPbMC_' + system + '_' + infoString + '.root'
inputFile='PbPbMC_Pythia6HydjetCymbal_list_part05.txt'

config.section_("General")
config.General.requestName = 'PbPbMC_' + system + '_' + infoString
config.General.workArea = config.General.requestName 
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+output]
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
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

