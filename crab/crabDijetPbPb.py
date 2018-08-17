from WMCore.Configuration import Configuration
config = Configuration()

card='cardDijetPbPb.input'
output='dijetTestPbPb_noMixing_2018-06-13.root'

config.section_("General")
config.General.requestName = 'dijetTestPbPb_noMixing_2018-06-13'
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
config.Data.userInputFiles = open('PbPb2015_data_Marta_csidfix_full.txt').readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 50
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetTestPbPbHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

