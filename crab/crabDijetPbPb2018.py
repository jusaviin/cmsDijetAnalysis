from WMCore.Configuration import Configuration
config = Configuration()

card='cardDijetPbPb2018.input'
jobTag='dijetPbPb2018_highForest_akFlowPuCs4PfJets_5eveMix_calo80Trigger_wtaAxis_JECv5b_2019-09-05_part5'
inputList='PbPbData2018_CaloJet80_newJets_part5.txt'
#inputList='PbPbData2018_newJets_combine0.txt'
outputFile=jobTag+'.root'
fileLocation='1'  # Locations: 0 = Purdue, 1 = CERN, 2 = Search with xrootd

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
config.JobType.maxJobRuntimeMin = 400
config.JobType.maxMemoryMB = 2000

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetPbPbHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_US_Purdue','T2_US_Wisconsin','T2_US_Nebraska']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

