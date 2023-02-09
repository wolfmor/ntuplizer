#! /usr/bin/env python

import os, sys,shutil
from glob import glob
from random import shuffle 
from time import sleep
import datetime
now = datetime.datetime.now()


limitNJobs = False
tag = ""
test = False

datasets = [
			"/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM"
			,"/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM"
			# ,"/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM"
			# ,"/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD"
			# ,"/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD"
			# ,"/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD"
			# ,"/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD"
			]

### for back. files on pnfs	
def submitToCrab():
	nJobs = 0
	for ijob, dataset in enumerate(datasets):

		tag = ((dataset.split('/')[1]).split("_")[0])+((dataset.split('/')[1]).split("_")[1])+(dataset.split('/')[2]).split("-")[0]
		#tag = dataset.split('/')[1]+(dataset.split('/')[2]).split("-")[0]
		print "use tag:", tag

		if os.path.exists("crab_projects/crab_" + tag): 
			shutil.rmtree("crab_projects/crab_" + tag)
			#os.rmdir("crab_projects/crab_" + tag)
		if os.path.exists("crabSubmissionScript_"+str(ijob)+".py"): os.remove("crabSubmissionScript_"+str(ijob)+".py")
		
		f = open("crabSubmissionScript_"+str(ijob)+".py", "w")
		
		f.write("from CRABClient.UserUtilities import config\n")
		f.write("import datetime\n")
		f.write("import glob\n")
		f.write("now = datetime.datetime.now()\n")
		f.write("config = config()\n")
		#f.write("config.General.requestName = \"" + tag +"_ \" +now.strftime(\"%Y_%m_%d\")\n")
		f.write("config.General.requestName = \"" + tag +"\"\n")
		f.write("config.General.workArea = \"crab_projects\"\n")
		f.write("config.General.transferOutputs = True\n")
		f.write("config.General.transferLogs = True\n")
		f.write("config.JobType.pluginName = \"Analysis\"\n")
		f.write("config.JobType.psetName = \"construct_secondary_vertices_cfg.py\"\n")
		f.write("config.JobType.scriptArgs = [\"tag=crab,genmatchalltracks,era16_UL_APV\" , \"dataset="+ dataset + "\"]\n")
		f.write("config.JobType.inputFiles = [\"/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/plantTrees.py\",\n")
		f.write("\"/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/commons.py\",\n")
		f.write("\"/nfs/dust/cms/user/wolfmor/NTupleStuff/\"]\n")
		f.write("config.JobType.allowUndistributedCMSSW = True\n")
		f.write("config.JobType.scriptExe = \"allinonejob.sh\"\n")
		f.write("config.JobType.outputFiles = [\"crab_NTuple.root\"]\n")
		f.write("config.Data.inputDataset = \"" + dataset +"\"\n")
		f.write("config.Data.inputDBS = \"global\"\n")
		f.write("config.Data.splitting = \"FileBased\"\n")
		f.write("config.Data.unitsPerJob = 1\n")
		if limitNJobs:
			f.write("NJOBS = 1 \n")
			f.write("config.Data.totalUnits = config.Data.unitsPerJob * NJOBS \n")
		f.write("config.Data.publication = False \n")
		f.write("config.Data.outputDatasetTag = \"" + tag +"_\" +now.strftime(\"%Y_%m_%d\")\n")
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
