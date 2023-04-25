#! /usr/bin/env python

import os, sys,shutil
from glob import glob
from random import shuffle 
from time import sleep
import datetime
now = datetime.datetime.now()


limitNJobs = True
data = False
tag = ""
test = False
sameSign = False


datasets = [
			'/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM'
			,'/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM'
			,'/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM'
			,'/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM'
			
			,'/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			#,'/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
			,'/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8_ext1-v1/AODSIM'
			,'/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM'
			,'/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM'
			,'/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM'
			]
datasets = [
			'/DYJetsToLL_M-50_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
]
datasets = [
			'/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
]
if data:
	datasets = [
				'/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD'
				,'/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD'
				,'/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD'
				,'/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD'
				,'/MET/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD'
				,'/MET/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD'
				,'/MET/Run2016F-21Feb2020_UL2016-v1/AOD'
				,'/MET/Run2016G-21Feb2020_UL2016-v1/AOD'
				,'/MET/Run2016H-21Feb2020_UL2016-v2/AOD'
				]

def submitToCrab():
	nJobs = 0
	for ijob, dataset in enumerate(datasets):

		
		if data:
			tag = ((dataset.split('/')[1]).split("/")[0])+"_Run"+(dataset.split('/Run')[1]).split("/")[0]
		else:
			tag = ((dataset.split('/')[1]).split("_")[0])+((dataset.split('/')[1]).split("_")[1])+(dataset.split('/')[2]).split("-")[0]
		print "use tag:", tag

		if os.path.exists("crab_projects/crab_" + tag): 
			shutil.rmtree("crab_projects/crab_" + tag)
		if os.path.exists("crabSubmissionScript_"+str(ijob)+".py"): os.remove("crabSubmissionScript_"+str(ijob)+".py")
		
		f = open("crabSubmissionScript_"+str(ijob)+".py", "w")
		
		f.write("from CRABClient.UserUtilities import config\n")
		f.write("import datetime\n")
		f.write("import glob\n")
		f.write("now = datetime.datetime.now()\n")
		f.write("config = config()\n")
		f.write("config.General.requestName = \"" + tag +"\"\n")
		f.write("config.General.workArea = \"crab_projects\"\n")
		f.write("config.General.transferOutputs = True\n")
		f.write("config.General.transferLogs = True\n")
		f.write("config.JobType.pluginName = \"Analysis\"\n")
		if sameSign:
			f.write("config.JobType.psetName = \"construct_secondary_vertices_ss_cfg.py\"\n")
		else:
			f.write("config.JobType.psetName = \"construct_secondary_vertices_cfg.py\"\n")
		
		if sameSign and not data:
			f.write("config.JobType.scriptArgs = [\"tag=crab,genmatchalltracks,era16_UL_APV,samesign\" , \"dataset="+ dataset + "\"]\n")
		#f.write("config.JobType.scriptArgs = [\"tag=crab,genmatchalltracks,era16_UL_APV,cleanleptons\" , \"dataset="+ dataset + "\"]\n")
		elif sameSign and data:
			f.write("config.JobType.scriptArgs = [\"tag=crab,data,era16_UL_APV,samesign\"]\n")
		elif data:
			f.write("config.JobType.scriptArgs = [\"tag=crab,data,era16_UL_APV\"]\n")
		else:
			f.write("config.JobType.scriptArgs = [\"tag=crab,genmatchalltracks,era16_UL_APV\" , \"dataset="+ dataset + "\"]\n")
			
		#f.write("config.JobType.inputFiles = [\"/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/plantTrees.py\",\n") #toDo fix once pushed ntupelizer
		f.write("config.JobType.inputFiles = [\"/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/plantTrees_cp.py\",\n") #toDo fix once pushed ntupelizer
		f.write("\"/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/commons.py\",\n")
		f.write("\"/nfs/dust/cms/user/wolfmor/NTupleStuff/\"]\n")
		f.write("config.JobType.allowUndistributedCMSSW = True\n")
		f.write("config.JobType.scriptExe = \"allinonejob.sh\"\n")
		f.write("config.JobType.outputFiles = [\"crab_NTuple.root\"]\n")
		if data:
			f.write("config.JobType.outputFiles = [\"crab_NTuple.root\", \"crab_NTuple.json\"]\n")
		f.write("config.Data.inputDataset = \"" + dataset +"\"\n")
		f.write("config.Data.inputDBS = \"global\"\n")
		f.write("config.Data.splitting = \"FileBased\"\n")
		f.write("config.Data.unitsPerJob = 1\n")
		f.write("config.Data.outLFNDirBase = \"/store/user/altews/NTuples/NTuplesV12/16UL_preAPV/\"\n")
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
		nJobs +=1

	print "submitted " + str(nJobs) + " jobs"

submitToCrab()

