#! /usr/bin/env python

import os, sys
from glob import glob
from random import shuffle 
from time import sleep

test = False

### for signal files (local)
def submitNTuples():
	#totalJobs = 124
	totalJobs = 64
	numPerJob = 1 #set to 1 for one file per job
	startAt = 0
	nJobs = 0
	for ifile in range(startAt, startAt+totalJobs):

		if os.path.exists("submissionScriptO_"+str(ifile+1)+".sh"): os.remove("submissionScriptO_"+str(ifile+1)+".sh")
			
		f = open("submissionScriptO_"+str(ifile+1)+".sh", "w")
		f.write("source /etc/profile.d/modules.sh\n")
		f.write("source /afs/desy.de/user/t/tewsalex/.bash_profile\n")
		f.write("module use -a /afs/desy.de/group/cms/modulefiles/\n")
		f.write("module load cmssw\n")
		f.write("cd /nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src\n")
		f.write("export SCRAM_ARCH=slc7_amd64_gcc820\n")
		f.write("cmsenv\n")
		### ToDo cange signal/background in plant_nTuplez.py and ctau
		#f.write("python plant_nTuplez_signal.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/edmfiles/testDeltaM1p14_MChi115_ctau0p005_"+str(ifile+1)+"of124_v2.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_testDeltaM1p14_MChi115_ctau0p005_MuTrack_"+str(ifile+1)+"_v3.root\" maxEvents=1000 startEvent=0 \n")
		#f.write("python plant_nTuplez_signal.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/edmfiles/testDeltaM1p14_MChi115_ctau0p005_"+str(ifile+1)+"of124_v2.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_testDeltaM1p14_MChi115_ctau0p005_MuTrack_"+str(ifile+1)+"_v3.root\" maxEvents=1000 startEvent=0 datatier=\"signal\" ctau=0.05 deltaM=\"1p14\" \n")
		f.write("python plant_nTuplez_signal.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/edmfiles/FullSim_mChipm115GeV_dm0p768GeV_ctauAll_"+str(ifile+1)+"of64.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_FullSim_testDeltaM1p54_MChi115_ctauAll_"+str(ifile+1)+".root\" maxEvents=1000 startEvent=0 datatier=\"signal\" ctau=-1 deltaM=\"0p786\" \n")
		#f.write("python plant_nTuplez_signal.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/edmfiles/testDeltaM1p54_MChi115_ctauAll_"+str(ifile+1)+"of44.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_testDeltaM1p54_MChi115_ctauAll_"+str(ifile+1)+".root\" maxEvents=1000 startEvent=0 datatier=\"signal\" ctau=-1 deltaM=\"0p786\" \n")
		#f.write("python plant_nTuplez_signal.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/edmfiles/testDeltaM1p54_MChi115_ctau0p005_"+str(ifile+1)+"of124_v2.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_testDeltaM1p54_MChi115_ctau0p005_"+str(ifile+1)+".root\" maxEvents=1000 startEvent=0 datatier=\"signal\" ctau=0.005 deltaM=\"1p54\" \n")
		#f.write("python plant_nTuplez_signal.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/edmfiles/testDeltaM1p14_MChi115_ctau0p001_MuTrack_*of2.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_testDeltaM1p14_MChi115_ctau0p001_MuTrack_"+str(nJobs+1)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((nJobs)*(numPerJob))))+"\n")
		#f.write("python plant_nTuplez.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/testDeltaM0p94_MChi115_ctau0p001_v2.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_testDeltaM0p94_MChi115_ctau0p001_v2_"+str(nJobs+1)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((nJobs)*(numPerJob))))+"\n")
		#f.write("python plant_nTuplez.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/testDeltaM1p14_MChi115_ctau0p005_MuTrack.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_pdgIDs_"+str(nJobs+1)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((nJobs)*(numPerJob))))+"\n")
		#f.write("python plant_nTuplez.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/test_SingleMuonRunC_large.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_ZLlDataCleaned_"+str(nJobs+1)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((nJobs)*(numPerJob))))+"\n")

		f.close()

		f = open("submissionScriptO_"+str(ifile+1)+".sh", "r")
		#print "Read file " +str(ifile+1)
		#print(f.read()) 
			
		command = 'condor_qsub -l h_vmem=2G -l h_rt=02:59:00 -cwd submissionScriptO_'+str(ifile+1)+'.sh &'

		print 'command', command
		if not test: os.system(command)

		sleep(1)
		nJobs +=1

	print 'submitted ' + str(nJobs) + ' jobs'
	
### for back. files on pnfs	
def submitNTuples_back():
	totalJobs = 514 #754 
	#numPerJob = 514 #set to 1 for one file per job !!!!!!
	startAt = 0
	nJobs = 0
	for ifile in range(startAt, startAt+totalJobs):

		if os.path.exists("submissionScriptA_"+str(ifile+1)+".sh"): os.remove("submissionScriptA_"+str(ifile+1)+".sh")
			
		f = open("submissionScriptA_"+str(ifile+1)+".sh", "w")
		f.write("source /etc/profile.d/modules.sh\n")
		f.write("source /afs/desy.de/user/t/tewsalex/.bash_profile\n")
		f.write("module use -a /afs/desy.de/group/cms/modulefiles/\n")
		f.write("module load cmssw\n")
		f.write("cd /nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src\n")
		f.write("export SCRAM_ARCH=slc7_amd64_gcc820\n")
		f.write("cmsenv\n")
		#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/sameSign_19/210524_091327/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_test_ZJetsToNuNu_SC_"+str(ifile+1)+"of"+str(totalJobs)+".root\" maxEvents=1000 startEvent=0 \n")		
		#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/WJetsToLNu_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/OC_W_10/210611_114514/0002/testSubmission2_"+str(ifile+2001)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_test_WJetsToNLNu_"+str(ifile+2001)+"of"+str(totalJobs)+"_v5.root\" maxEvents=2000 startEvent=8000 \n")		
		f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/DYJetsToLL_Zpt-200toInf_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_Zpt-200toInf/210709_113425/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_test_DYJetsMC_"+str(ifile+1)+"of"+str(totalJobs)+".root\" maxEvents=8000 startEvent=0 \n")		
		#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/WJetsToLNu_HT-400To600_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/OC_W_HT400to600_2/210611_114552/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_test_WJetsToNLNu_HT400to600_"+str(ifile+1)+"of"+str(totalJobs)+"_v4.root\" maxEvents=2000 startEvent=6000 \n")		
		#f.write("python plant_nTuplez.py inputFiles=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/test_ZJetsToNuNu.root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_test_ZJetsToNuNu_OC_"+str(ifile+1)+"of"+str(totalJobs)+".root\" maxEvents=1000 startEvent="+str(ifile*1000)+" \n")		
		
		### !!!
		### ToDo change signal/background/clean DY in plant_nTuplez.py
		### !!!
		
		f.close()

		f = open("submissionScriptA_"+str(ifile+1)+".sh", "r")
		#print "Read file " +str(ifile+1)
		#print(f.read()) 
			
		command = 'condor_qsub -l h_vmem=2G -l h_rt=06:00:00 -cwd submissionScriptA_'+str(ifile+1)+'.sh &'

		print 'command', command
		if not test: os.system(command)

		sleep(1)
		nJobs +=1

	print 'submitted ' + str(nJobs) + ' jobs'
	
### for back. files on pnfs	
def submitNTuples_back_2():
	nFiles = 1000
	#nFiles = 320
	#nEventsTotal = 8000
	#nEventsTotal = 20000
	#nEventsTotal = 16000 #for DY
	nEventsTotal =  9000 #for MET 1600, 4000
	numPerJob = 4500
	perFile = nEventsTotal/numPerJob
	startAt = 0
	nJobs = 0
	for ifile in range(nFiles):
		for ijob in range(perFile):

			if os.path.exists("submissionScriptB_"+str(ifile+1)+"_"+str(ijob)+".sh"): os.remove("submissionScriptB_"+str(ifile+1)+"_"+str(ijob)+".sh")
				
			f = open("submissionScriptB_"+str(ifile+1)+"_"+str(ijob)+".sh", "w")
			f.write("source /etc/profile.d/modules.sh\n")
			f.write("source /afs/desy.de/user/t/tewsalex/.bash_profile\n")
			f.write("module use -a /afs/desy.de/group/cms/modulefiles/\n")
			f.write("module load cmssw\n")
			f.write("cd /nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src\n")
			f.write("export SCRAM_ARCH=slc7_amd64_gcc820\n")
			f.write("cmsenv\n")
			#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/MET/data_MET_Run2016D/210625_091322/0001/testSubmission2_"+str(ifile+1373)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_MET_Run2016D_"+str(ifile+1373)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/MET/data_MET_Run2016H/210625_093941/0004/testSubmission2_"+str(ifile+4001)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_MET_Run2016H_"+str(ifile+4001)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/MET/data_MET_Run2016E/210625_091935/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_MET_Run2016E_"+str(ifile+1)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/MET/data_MET_Run2016G/210625_093334/0001/testSubmission2_"+str(ifile+1001)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_MET_Run2016G_"+str(ifile+1001)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/MET/data_MET_Run2016B/210624_143014/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_MET_Run2016C_"+str(ifile+1)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			#f.write("python plant_nTuplez.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/SingleMuon/DY_9/210715_094044/0002/testSubmission2_"+str(ifile+2751)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_DYNoCleaning_Run2016E_"+str(ifile+2751)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			#f.write("python plant_nTuplez_DYdata.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/SingleMuon/DY_8/210621_084356/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_DY_Run2016B_"+str(ifile+1)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	
			#f.write("python plant_nTuplez_DY.py inputFiles=\"/pnfs/desy.de/cms/tier2/store/user/altews/DYJetsToLL_Zpt-200toInf_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/DYJetsToLL_Zpt-200toInf/210709_113425/0000/testSubmission2_"+str(ifile+1)+".root\" outputFile=\"/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/NTuple_test_DYJetsMC_"+str(ifile+1)+"of"+str(nFiles)+"_v"+str(ijob)+".root\" maxEvents="+str(numPerJob)+" startEvent="+str((startAt+((ijob)*(numPerJob))))+"\n")	

			### !!!
			### ToDo change signal/background/clean DY in plant_nTuplez.py
			### !!!
			
			f.close()

			f = open("submissionScriptB_"+str(ifile+1)+"_"+str(ijob)+".sh", "r")
			#print "Read file " +str(ifile+1)
			#print(f.read()) 
				
			#command = "condor_qsub -l h_vmem=1G -l h_rt=00:30:00 -cwd submissionScriptB_"+str(ifile+1)+"_"+str(ijob)+".sh"
			command = "condor_qsub -l mem=2Gb -l h_rt=05:59:00 -cwd submissionScriptB_"+str(ifile+1)+"_"+str(ijob)+".sh"

			print 'command', command
			if not test: os.system(command)

			sleep(1)
			nJobs +=1

	print 'submitted ' + str(nJobs) + ' jobs'

### to produce slimmed edms for signal 
def submitEDMFiles():
	nJobs = 0
	#for ifile in range(totalJobs):
	myfiles = glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v*/higgsino94x_susyall_mChipm115GeV_dm0p77GeV_pu35_part*of2?.root")
	#myfiles = glob("/nfs/dust/cms/user/beinsam/CommonSamples/MC_BSM/CompressedHiggsino/RadiativeMu_2016Full/v2/higgsino94xfull_susyall_mChipm115GeV_dm0p768GeV_pu35_part*of100.root")


	for ifile in range(len(myfiles)):	
		if os.path.exists("submissionScriptT"+str(ifile+1)+".sh"): os.remove("submissionScriptT"+str(ifile+1)+".sh")
		
		f = open("submissionScriptT"+str(ifile+1)+".sh", "w")
		f.write("source /etc/profile.d/modules.sh\n")
		f.write("source /afs/desy.de/user/t/tewsalex/.bash_profile\n")
		f.write("module use -a /afs/desy.de/group/cms/modulefiles/\n")
		f.write("module load cmssw\n")
		f.write("cd /nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src\n")
		f.write("export SCRAM_ARCH=slc7_amd64_gcc820\n")
		f.write("scram b -j12\n")
		f.write("cmsenv\n")
		f.write("cmsRun construct_secondary_vertices_cfg.py inputFiles=\"file:"+myfiles[nJobs]+"\" maxEvents=-1 outputFile=\"/nfs/dust/cms/user/tewsalex/rootfiles/edmfiles/signal/higgsino94xfast_slimmededm_mChipm115GeV_dm0p768GeV_pu35_"+str(nJobs)+"of"+str(len(myfiles))+".root\"\n")
		#f.write("cmsRun construct_secondary_vertices_cfg.py inputFiles=\"file:"+myfiles[nJobs]+"\" maxEvents=-1 outputFile=\"rootfiles/testDeltaM1p14_MChi115_ctauAll_"+str(nJobs)+"of"+str(len(myfiles))+".root\"\n")
		#print "cmsRun construct_secondary_vertices_cfg.py inputFiles=\"file:"+myfiles[nJobs]+"\" maxEvents=-1 outputFile=\"rootfiles/testDeltaM1p94_MChi115_ctau0p005_"+str(nJobs)+"of"+str(len(myfiles))+"_v2.root\"\n"

		f.close()
		
		f = open("submissionScriptT"+str(ifile+1)+".sh", "r")

		command = 'condor_qsub -l h_vmem=4G -l h_rt=06:00:00 -cwd submissionScriptT'+str(ifile+1)+'.sh &'
		
		print 'command', command
		if not test: os.system(command)

		sleep(1)
		nJobs +=1

	print 'submitted 1 jobs'
	
def main():
	#submitNTuples()
	#submitNTuples_back()
	#submitNTuples_back_2()
	submitEDMFiles()

	
main()

