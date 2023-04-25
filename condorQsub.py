#! /usr/bin/env python

import os, sys
from glob import glob
from random import shuffle 
from time import sleep

test = False
sameSign = False 
	
### to produce slimmed edms for signal 
def submit_allInOneJob():
	nJobs = 0
	#myfiles = glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v4/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p7*GeV_Chi20ctau5MM_part11of25.root")
	#myfiles = glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v4/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm*GeV_Chi20ctau*MM_part*of100.root")
	myfiles = glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v4/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm*GeV_dm*GeV_part*of*.root")

	for ifile in range(len(myfiles)):	
		if os.path.exists("submissionScriptU"+str(ifile+1)+".sh"): os.remove("submissionScriptU"+str(ifile+1)+".sh")
		nameout = myfiles[nJobs].split('/')[-1]
		
		f = open("submissionScriptU"+str(ifile+1)+".sh", "w")
		f.write("source /etc/profile.d/modules.sh\n")
		f.write("source /afs/desy.de/user/t/tewsalex/.bash_profile\n")
		f.write("module use -a /afs/desy.de/group/cms/modulefiles/\n")
		f.write("module load cmssw\n")
		f.write("cd /nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src\n")
		f.write("export SCRAM_ARCH=slc7_amd64_gcc820\n")
		f.write("scram b -j12 RecoVertex\n")
		f.write("cmsenv\n")
		if sameSign: #toDo test!!
			f.write("cmsRun construct_secondary_vertices_cfg.py inputFiles=\"file://"+myfiles[nJobs]+"\" maxEvents=-1 crab=False sameSign=True outputFile=\"/nfs/dust/cms/user/tewsalex/rootfiles/edmfiles/svfile_ss_"+nameout+"\"\n")
			f.write("python /nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/plantTrees.py inputFiles=\"file:"+myfiles[nJobs]+"\" tag=\"era16_UL, fastsim, local, signal, samesign\"\n")	
		else:
			f.write("cmsRun construct_secondary_vertices_cfg.py inputFiles=\"file://"+myfiles[nJobs]+"\" maxEvents=-1 crab=False outputFile=\"/nfs/dust/cms/user/tewsalex/rootfiles/edmfiles/svfile_"+nameout+"\"\n")
			f.write("python /nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/plantTrees_cp.py inputFiles=\"file:"+myfiles[nJobs]+"\" tag=\"era16_UL, fastsim, local, signal\"\n")	
			
		#f.write("python /nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/plantTrees.py inputFiles=\"file:"+myfiles[nJobs]+"\" tag=\"era16_UL, local, signal, fastsim, samesign\"\n")	
		#ToDo: add hadd of NTuplefiles to this workflow 
			
		f.close()
		
		f = open("submissionScriptU"+str(ifile+1)+".sh", "r")

		command = 'condor_qsub -l h_vmem=4G -l h_rt=06:00:00 -cwd submissionScriptU'+str(ifile+1)+'.sh &'
		
		print 'command', command
		if not test: os.system(command)

		sleep(1)
		nJobs +=1

	print 'submitted 1 jobs'
	
	
def haddNTuple():
	
	myfiles = glob("/nfs/dust/cms/user/tewsalex/CRAB3-tutorial/CMSSW_10_6_18/src/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm*GeV_dm*GeV_part*of*.root")
	samefiles ={}
	for afile in myfiles:	
		#mass = ((ifile.split('/')[-1])
		mass = ((afile.split('/')[-1]).split('mChipm')[-1]).split('GeV_dm')[0]
		#print mass
		dm = afile.split('/')[-1].split('GeV_dm')[-1].split('GeV_part')[0]

		
		if not mass in samefiles.keys(): 
			samefiles[mass] = {}
			samefiles[mass][dm] = []
			if afile not in samefiles[mass][dm]:samefiles[mass][dm].append(afile)
			
		else: 
			if not dm in samefiles[mass].keys():
				samefiles[mass][dm] = []
				samefiles[mass][dm].append(afile)
			else: 
				if afile not in samefiles[mass][dm]:samefiles[mass][dm].append(afile)
				
	for mass in samefiles.keys():
		print 'mass ', mass 
		for dm in samefiles[mass].keys():
			print '        dm ', dm
			for idx, afile in enumerate(samefiles[mass][dm]):
				#print '              file ', afile
				if idx == 0:
					command = 'hadd -f /nfs/dust/cms/user/tewsalex/rootfiles/ntuple/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm' + mass +'GeV_dm' + dm + '_NTuple.root'+ ' ' + afile
					print command
					os.system(command)
				else: 
					command = 'hadd -a /nfs/dust/cms/user/tewsalex/rootfiles/ntuple/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm' + mass +'GeV_dm' + dm + '_NTuple.root'+ ' ' + afile
					print command
					os.system(command)
			

	
def main():
	submit_allInOneJob()
	#haddNTuple()

	
main()

