from CRABClient.UserUtilities import config
import datetime
import glob
now = datetime.datetime.now()

config = config()

## era 16_UL
dataset_in = '/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD'

config.General.requestName = 'test_data_'+now.strftime("%Y_%m_%d")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'construct_secondary_vertices_cfg.py'

## era 16_UL
#config.JobType.scriptArgs = ['tag=crab,genmatchalltracks,era16_UL_APV','dataset='+dataset_in]
#config.JobType.scriptArgs = ['tag=crab,genmatchalltracks,era16_UL_APV,cleanleptons', 'dataset='+dataset_in]
config.JobType.scriptArgs = ['tag=crab,era16_UL,data']

## era 16_07Aug17
#config.JobType.scriptArgs = ['tag=crab,era16_07Aug17,data,cleanleptons']

 
config.JobType.inputFiles = ['/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/plantTrees.py', 
'/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/commons.py',
'/nfs/dust/cms/user/wolfmor/NTupleStuff/']

config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 3500
config.JobType.scriptExe = 'allinonejob.sh'

config.JobType.outputFiles = ['crab_NTuple.root']

## era 16_UL
config.Data.inputDataset = dataset_in
config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
NJOBS = 1
config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

config.Data.outLFNDirBase = '/store/user/altews/NTuples/NTuplesVXXX/16UL_preAPV/' 
config.Data.publication = False
config.Data.outputDatasetTag = now.strftime("%Y_%m_%d")

#config.Site.blacklist = ['T1_RU_JINR']
config.Site.storageSite = "T2_DE_DESY" 
