from WMCore.Configuration import Configuration
config = Configuration()

card='cardDijetPbPbNoMixing2018.input'
jobTag='dijetPbPb2018_akCaloJet_onlyJets_jet100TriggerEta1v3_2021-09-09'
#inputList='PbPbData2018_MinBiasFiles.txt'
inputList='triggerEfficiencyFiles2018.txt'
#inputList='PbPbData2018_flowJetSkims_80and100triggers.txt'
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
config.JobType.maxJobRuntimeMin = 120
config.JobType.maxMemoryMB = 300

config.section_("Data")
config.Data.userInputFiles = open(inputList).readlines() 
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1000
config.Data.totalUnits = len(config.Data.userInputFiles)
config.Data.outputPrimaryDataset = 'diJetPbPbHistograms'
config.Data.outLFNDirBase = '/store/user/jviinika/'+config.General.requestName
config.Data.publication = False

config.section_("Site")
config.Site.whitelist = ['T2_CH_CERN']
#config.Site.whitelist = ['T2_US_Purdue','T2_US_Wisconsin','T2_US_Nebraska']
config.Site.storageSite = 'T3_US_FNALLPC'

#"really" force crab to only run at whitelisted sites
config.section_("Debug")
config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']

