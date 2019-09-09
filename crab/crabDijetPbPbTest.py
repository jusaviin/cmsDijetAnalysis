from WMCore.Configuration import Configuration
config = Configuration()

card='cardDijetPbPbTest.input'
jobTag='dijetPbPb_skims_pfJets_quickTest_2019-01-24'
inputList='PbPb_5TeV_skim_test.txt'
outputFile=jobTag+'.root'
fileLocation='2'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

config.section_("General")
config.General.requestName = jobTag
config.General.workArea = config.General.requestName 

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'PSet.py'
config.JobType.scriptExe = 'compileAndRun.sh'
config.JobType.scriptArgs = ['card='+card,'output='+outputFile,'location='+fileLocation]
config.JobType.inputFiles = ['FrameworkJobReport.xml','dijet5TeV.tar.gz',card]
config.JobType.outputFiles = [outputFile]
config.JobType.maxJobRuntimeMin = 1440
config.JobType.maxMemoryMB = 2500

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetPbPbHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_*']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

