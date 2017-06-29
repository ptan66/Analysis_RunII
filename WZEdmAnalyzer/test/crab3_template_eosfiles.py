from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

#config.General.requestName = 'testMC'
config.General.workArea = '?subdir'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
#config.JobType.psetName = 'had_ntuplizer.py'
config.JobType.psetName = '?config_py'





config.Data.outputPrimaryDataset = '?output_tag'
config.Data.userInputFiles = open('?eos_files').readlines()
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob =1
config.Data.outLFNDirBase = '/store/user/ptan/noreplica/Private_Ntuple/'
config.Data.publication = False
config.Data.outputDatasetTag = '?dataset_tag'

config.Site.storageSite = 'T3_US_FNALLPC'
#config.Site.whitelist   = ['T3_US_FNALLPC', 'T1_US_FNAL']
#config.Debug.extraJDL = ['+CMS_ALLOW_OVERFLOW=False']
