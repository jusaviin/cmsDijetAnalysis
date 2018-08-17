from WMCore.Configuration import Configuration
config = Configuration()

config.section_("General")
config.General.requestName = 'dijetTestPp_2019-04-18'
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
#config.JobType.scriptArgs = ['ppMC_Pythia6_forest_5TeV.txt','cardDijet.input','dijetSpectraTest.root','2']
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeVpp.tar.gz']
config.JobType.outputFiles = ['dijetSpectraTestPp.root']
config.JobType.maxJobRuntimeMin = 600

config.section_("Data")
config.Data.userInputFiles = open('ppMC_Pythia6_forest_5TeV.txt').readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 100
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetTestHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

