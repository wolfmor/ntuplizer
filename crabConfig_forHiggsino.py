from CRABClient.UserUtilities import config
import datetime
import glob
now = datetime.datetime.now()

config = config()
#tag = 'test'

## era 16_UL
#tag = 'DYJetsToLL_M-50_Zpt-200toInf_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'TTJets_DiLept_RunIISummer20UL16RECOAPV_preVFP'
tag = 'ZJetsToNuNu_Zpt-200toInf_BPSFilter_RunIISummer20UL16RECOAPV_preVFP'

#tag = 'WJetsToLNu_HT-200To400_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-70To100_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-100to200_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-400to600_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-600to800_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-800to1200_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-1200to2500_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'WJetsToLNu_HT-2500toInf_RunIISummer20UL16RECOAPV_preVFP'

#tag = 'ZJetsTNuNu_HT-200To400_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'ZJetsTNuNu_HT-800To1200_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'ZJetsTNuNu_HT-1200To2500_RunIISummer20UL16RECOAPV_preVFP'
#tag = 'ZJetsTNuNu_HT-2500ToInf_RunIISummer20UL16RECOAPV_preVFP'

#tag = 'SingleMuon_Run2016B_ver1_UL2016_HIPM'
#tag = 'SingleMuon_Run2016B_ver2_UL2016_HIPM'
#tag = 'Test_MET_Run2016E-21Feb2020_UL2016_HIPM'


#tag = 'MET_Run2016B-21Feb2020_ver1_UL2016_HIPM'
#tag = 'MET_Run2016B-21Feb2020_ver2_UL2016_HIPM'
#tag = 'MET_Run2016C-21Feb2020_UL2016_HIPM'
#tag = 'MET_Run2016D-21Feb2020_UL2016_HIPM'
#tag = 'MET_Run2016E-21Feb2020_UL2016_HIPM'
#tag = 'MET_Run2016F-21Feb2020_UL2016_HIPM'
#tag = 'MET_Run2016F-21Feb2020_UL2016'
#tag = 'MET_Run2016G-21Feb2020_UL2016'
#tag = 'MET_Run2016H-21Feb2020_UL2016'

#tag = 'ZJetsToNuNu_Zpt-200toInf_BPSFilter_RunIISummer20UL16RECO'
#

## era 16_07Aug17
#tag = 'ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#tag = 'DYJetsToLL_M-50_Zpt-150toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8'
#tag = 'SingleMuon_Run2016B-07Aug17_ver2-v1'

config.General.requestName = tag+'_'+now.strftime("%Y_%m_%d")
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'construct_secondary_vertices_cfg.py'

## era 16_UL
#config.JobType.scriptArgs = ['tag=crab,era16_UL_APV']
config.JobType.scriptArgs = ['tag=crab,genmatchalltracks,era16_UL_APV','dataset=/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM']
#config.JobType.scriptArgs = ['tag=crab,era16_UL_APV,TTJets_DiLept']
#config.JobType.scriptArgs = ['tag=crab,genmatchalltracks,era16_UL_APV,ZJetsToNuNu_Zpt-200toInf_UL16_preVFP']
#config.JobType.scriptArgs = ['tag=crab,genmatchalltracks,era16_UL,ZJetsToNuNu_Zpt-200toInf_UL16']
#config.JobType.scriptArgs = ['tag=crab,era16_UL_APV,cleanleptons,DYJetsToLL_M-50_Zpt-200toInf_UL16_preVFP']
#config.JobType.scriptArgs = ['tag=crab,era16_UL,data']


## era 16_07Aug17
#config.JobType.scriptArgs = ['tag=crab,era16_07Aug17,data,cleanleptons']

 
config.JobType.inputFiles = ['/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/plantTrees.py', # if use absolute paths here, only the filename is needed in NTupleizer; not the full path again, max. 120 MB to sandbox
'/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/commons.py',
'/nfs/dust/cms/user/wolfmor/NTupleStuff/']

config.JobType.allowUndistributedCMSSW = True
#config.JobType.maxMemoryMB = 3500
config.JobType.scriptExe = 'allinonejob.sh'

#config.JobType.outputFiles = ['simpleoutput.txt',', 'crab_NTuple.root']
config.JobType.outputFiles = ['crab_NTuple.root']

## era 16_UL

#config.Data.inputDataset = '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8_ext1-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'

#config.Data.inputDataset = '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM'
#config.Data.inputDataset = '/DYJetsToLL_M-50_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'

config.Data.inputDataset = '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'

#config.Data.inputDataset = '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM'

#config.Data.inputDataset = '/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD'

#config.Data.inputDataset = '/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016F-21Feb2020_UL2016-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016G-21Feb2020_UL2016-v1/AOD'
#config.Data.inputDataset = '/MET/Run2016H-21Feb2020_UL2016-v2/AOD'

#config.Data.inputDataset = '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM'

## era 16_07Aug17
#config.Data.inputDataset = '/SingleMuon/Run2016B-07Aug17_ver2-v1/AOD'
#config.Data.inputDataset = '/DYJetsToLL_M-50_Zpt-150toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'
#config.Data.inputDataset = '/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIISummer16DR80Premix-PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/AODSIM'

config.Data.inputDBS = 'global'

config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#NJOBS = 1
#NJOBS = 1000
#config.Data.totalUnits = config.Data.unitsPerJob * NJOBS

#config.Data.outLFNDirBase = '/store/user/altews/' 
config.Data.publication = False # this must be false for additional non-root outputfiles
config.Data.outputDatasetTag = tag + '_'+now.strftime("%Y_%m_%d")

#config.Site.blacklist = ['T1_RU_JINR']
config.Site.storageSite = "T2_DE_DESY" # uncomment or put 
