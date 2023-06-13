#! /usr/bin/env python

import os, sys, shutil
from glob import glob
from random import shuffle 
from time import sleep
import datetime
now = datetime.datetime.now()


datasets_era16_UL_APV_MC = [
    # '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',

    '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',

    '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    #'/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8_ext1-v1/AODSIM',

    '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',

    '/DYJetsToLL_M-50_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
]

datasets_era16_UL_MC = [
    # '/ZJetsToNuNu_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',

    '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',

    '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v1/AODSIM',

    '/DYJetsToLL_M-50_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',

    '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
]

datasets_era16_UL_MET = [
    '/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD',
    '/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016F-21Feb2020_UL2016-v1/AOD',
    # '/MET/Run2016G-21Feb2020_UL2016-v1/AOD',
    # '/MET/Run2016H-21Feb2020_UL2016-v2/AOD',
]

datasets_era16_UL_SingleMuon = [
    '/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016F-21Feb2020_UL2016-v1/AOD',
    '/SingleMuon/Run2016G-21Feb2020_UL2016-v1/AOD',
    '/SingleMuon/Run2016H-21Feb2020_UL2016-v1/AOD',
]

# TODO: add low-mass DY?
# TODO: add ttW?

test = False
limitNJobs = False
moritz = True

data = True
sameSign = False
cleanDY = True

datasets = datasets_era16_UL_MC
if data:
    datasets = datasets_era16_UL_SingleMuon  # datasets_era16_UL_MET  #

ntuplizer_tag = "crab,genmatchalltracks,era16_UL"  # era16_UL_APV  # TODO: adapt
if data:
    ntuplizer_tag = ntuplizer_tag.replace('genmatchalltracks', 'data')
if sameSign:
    ntuplizer_tag += ",samesign"
if cleanDY:
    ntuplizer_tag += ",cleanleptons"

if moritz:
    crabprojectpath = "/nfs/dust/cms/user/wolfmor/"
    inputpath = "/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_6_34/src/SoftDisplacedPion/ntuplizer/"
    outputpath = "/store/user/mowolf/NTuples/NTuplesV12/16UL/"
else:
    crabprojectpath = ""
    inputpath = "/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/"
    outputpath = "/store/user/altews/NTuples/NTuplesV12/16UL_preAPV/"


def submitToCrab():
    nJobs = 0
    for ijob, dataset in enumerate(datasets):

        if data:
            task_tag = ((dataset.split('/')[1]).split("/")[0])+"_Run"+(dataset.split('/Run')[1]).split("/")[0]
        else:
            task_tag = ((dataset.split('/')[1]).split("_")[0])+((dataset.split('/')[1]).split("_")[1])+(dataset.split('/')[2]).split("-")[0]

        if sameSign:
            task_tag += "_sameSign"
        if cleanDY:
            task_tag += "_cleaned"

        print "use tag:", task_tag

        if os.path.exists(crabprojectpath + "crab_projects/crab_" + task_tag):
            shutil.rmtree(crabprojectpath + "crab_projects/crab_" + task_tag)
        if os.path.exists("crabSubmissionScript_"+str(ijob)+".py"): os.remove("crabSubmissionScript_"+str(ijob)+".py")

        f = open("crabSubmissionScript_"+str(ijob)+".py", "w")

        f.write("from CRABClient.UserUtilities import config\n")
        f.write("import datetime\n")
        f.write("import glob\n")
        f.write("now = datetime.datetime.now()\n")
        f.write("config = config()\n")
        f.write("config.General.requestName = \"" + task_tag +"\"\n")
        f.write("config.General.workArea = \"" + crabprojectpath + "crab_projects\"\n")
        f.write("config.General.transferOutputs = True\n")
        f.write("config.General.transferLogs = True\n")
        f.write("config.JobType.pluginName = \"Analysis\"\n")
        if sameSign:
            f.write("config.JobType.psetName = \"construct_secondary_vertices_ss_cfg.py\"\n")
        else:
            f.write("config.JobType.psetName = \"construct_secondary_vertices_cfg.py\"\n")
        f.write("config.JobType.scriptArgs = [\"tag=" + ntuplizer_tag + "\" , \"dataset=" + dataset + "\"]\n")
        f.write("config.JobType.inputFiles = [\"" + inputpath + "plantTrees.py\",\n")
        f.write("\"" + inputpath + "commons.py\",\n")
        f.write("\"/nfs/dust/cms/user/wolfmor/NTupleStuff/\"]\n")
        f.write("config.JobType.allowUndistributedCMSSW = True\n")
        f.write("config.JobType.scriptExe = \"allinonejob.sh\"\n")
        if sameSign:
            if data: f.write("config.JobType.outputFiles = [\"crab_ss_NTuple.root\", \"crab_ss_NTuple.json\"]\n")
            else: f.write("config.JobType.outputFiles = [\"crab_ss_NTuple.root\"]\n")
        elif data: f.write("config.JobType.outputFiles = [\"crab_NTuple.root\", \"crab_NTuple.json\"]\n")
        else: f.write("config.JobType.outputFiles = [\"crab_NTuple.root\"]\n")
        f.write("config.Data.inputDataset = \"" + dataset + "\"\n")
        f.write("config.Data.inputDBS = \"global\"\n")
        f.write("config.Data.splitting = \"FileBased\"\n")
        f.write("config.Data.unitsPerJob = 1\n")
        f.write("config.Data.outLFNDirBase = \"" + outputpath + "\"\n")
        if sameSign:
            f.write("config.Data.outLFNDirBase = \"" + outputpath + "samesign/\"\n")
        if limitNJobs:
            f.write("NJOBS = 1 \n")
            f.write("config.Data.totalUnits = config.Data.unitsPerJob * NJOBS \n")
        f.write("config.Data.publication = False \n")
        f.write("config.Data.outputDatasetTag = now.strftime(\"%Y_%m_%d\")\n")
        f.write("config.Site.storageSite = \"T2_DE_DESY\"\n")
        f.close()

        f = open("crabSubmissionScript_" +str(ijob)+ ".py", "r")
        command = "crab submit -c crabSubmissionScript_" +str(ijob)+ ".py & "
        print "command", command
        if not test: os.system(command)

        sleep(120)
        nJobs += 1

    if not test: print "submitted " + str(nJobs) + " tasks"

submitToCrab()

