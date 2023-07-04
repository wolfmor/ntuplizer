
from ROOT import *
from array import array
from math import sqrt, pi
import re, array
import glob
import cms_figure
import numpy as np
from numpy import linspace
from decimal import *
from glob import glob
from tmva_include import book_tmva, fill_tmva
import uproot
import gfal2
import os 

gSystem.Load("libFWCoreFWLite.so")
gSystem.Load("libDataFormatsFWLite.so")
FWLiteEnabler.enable()
from DataFormats.FWLite import Events, Handle

gROOT.SetBatch(True)
gStyle.SetOptStat(111111) 

import uproot
	
### , 
### , 
### , 
def addBDTScores(
									fin ='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/test.root', 
									MVATree='Data/TestTree', 
									MVAFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/simpleTMVA.root', 
									weightFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/Data/weights/TMVAClassification_DNN.weights.xml',  
									
									):
	fout = TFile(fin)
	#fout = TFile(fin, "update")
	#fout = TFile("root://dcache-cms-xrootd.desy.de:1094/"+fin, "update")
	tree = fout.Get('tEvent') #todo: i think the problem lies here: I should make this tree is still bound th the fut, but gets a new branch: that gives an error
	#tree = TTree()
	#fout.Get('tEvent').Copy(tree)
	
	print "Retrieve tree" 

	max_length = -1
	
	#bdt_scores = array.array("f", 1000*[-999]) #too Fix
	bdt_scores = std.vector('double')()
	#bdt_branch = tree.Branch("bdt_score", bdt_scores, "bdt_score[n_sv]/F")
	bdt_branch = tree.Branch("bdt_score", bdt_scores)
	
	print "Retrieve MVA values "
	print "								"

	chain_BDT = TChain('Data/TestTree')
	chain_BDT.Add(MVAFile)
	
	### Booking the reader for all variables in the BDT test tree
	reader = TMVA.Reader()
	print chain_BDT
	print chain_BDT.GetListOfBranches()
	var = {}
	BDTvalues = {}
	for b, branch in enumerate(chain_BDT.GetListOfBranches()):
		branch_name = branch.GetName()
		if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
		var[branch_name] = array.array("f",[0])
		reader.AddVariable(branch_name,var[branch_name])
		
	reader.BookMVA("BDT", weightFile)

	nevents = tree.GetEntries()
	for ievent in range (nevents):
		tree.GetEntry(ievent)
		nsv = tree.n_sv
		BDTvalues = [None]*nsv
		bdt_scores.clear()
		
		#if nsv==0: continue
		print ievent, nsv
		for isv in range(nsv):
			for b, branch in enumerate(chain_BDT.GetListOfBranches()):
				branch_name = branch.GetName()
				if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
				var[branch_name][0] = getattr(tree, branch_name)[isv]
				
			mybdt_for_this_sv = reader.EvaluateMVA("BDT")
			print ievent, isv, mybdt_for_this_sv
			#bdt_scores[isv]= mybdt_for_this_sv
			bdt_scores.push_back(mybdt_for_this_sv)
			
		bdt_branch.Fill()
		if nsv>0:
			tree.Show(ievent)
		#print " single val"
		#print bdt_scores
			for val in tree.bdt_score: print val
			#exit()
	print "end of the script" 
	
	fout.cd()
	tree.Write("tEvent", TObject.kOverwrite)

	return tree

def joinTreeByFriend(tree, friendFile, friendTree):
	
	try:
		ff = ROOT.TFile(friendFile)
		tf = ff.Get(friendTree)
		tree.AddFriend(tf)
	except:
		print "failed to add friend", friendFile, "to tree"

### , 
### , 
### , 
# # # >>> import ROOT
# # # >>> f = ROOT.TFile("/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm1p168GeV_Chi20ctau5MM_part5of100_NTuple.root")
# # # >>> t = f.Get("tEvent")
# # # >>> ff = ROOT.TFile("/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm1p168GeV_Chi20ctau5MM_part5of100_NTuple_friend.root")
# # # >>> tf = ff.Get("tFriend")
# # # >>> t.AddFriend(tf)
# # # <ROOT.TFriendElement object ("tFriend") at 0x54748f0>
# # # >>> t.Draw("deltaR")
# # # Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
# # # >>> t.Draw("bdt_score")

def makeFriendTree(
									fin ='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm1p168GeV_Chi20ctau5MM_part5of100_NTuple.root', 
									outpath = '/nfs/dust/cms/user/tewsalex/rootfiles/friendTrees_V12/',
									MVATree='Data/TestTree', 
									MVAFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/simpleTMVA.root', 
									weightFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/Data/weights/TMVAClassification_DNN.weights.xml',  
									
									):

	
	finTree = TFile(fin)
	tree = finTree.Get('tEvent') 
	
	outname = fin.split('/')[-1].split('.root')[0]
	fout = TFile(outpath+outname+'_friend.root', 'recreate')	
	
	friendTree = TTree('tFriend', 'tFriend')


	max_length = -1
	
	bdt_scores = std.vector('double')()
	is_maxscore = std.vector('int')()
	bdt_branch = friendTree.Branch("bdt_score", bdt_scores)
	bdt_maxbranch = friendTree.Branch("is_maxscoring_sv", is_maxscore)
	
	print "Retrieve MVA values "
	print "								"

	chain_BDT = TChain('Data/TestTree')
	chain_BDT.Add(MVAFile)
	
	### Booking the reader for all variables in the BDT test tree
	reader = TMVA.Reader()
	print chain_BDT
	print chain_BDT.GetListOfBranches()
	var = {}
	BDTvalues = {}
	for b, branch in enumerate(chain_BDT.GetListOfBranches()):
		branch_name = branch.GetName()
		if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
		var[branch_name] = array.array("f",[0])
		reader.AddVariable(branch_name,var[branch_name])
		
	reader.BookMVA("BDT", weightFile)

	nevents = tree.GetEntries()
	for ievent in range (nevents):
		tree.GetEntry(ievent)
		nsv = tree.n_sv
		BDTvalues = [None]*nsv
		bdt_scores.clear()
		is_maxscore.clear()
		
		maxvalue = double(-999)
		#print ievent, nsv
		for isv in range(nsv):
			for b, branch in enumerate(chain_BDT.GetListOfBranches()):
				branch_name = branch.GetName()
				if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
				var[branch_name][0] = getattr(tree, branch_name)[isv]
				
			mybdt_for_this_sv = reader.EvaluateMVA("BDT")
			#print ievent, isv, mybdt_for_this_sv
			if mybdt_for_this_sv > maxvalue: maxvalue = mybdt_for_this_sv

			bdt_scores.push_back(mybdt_for_this_sv)

		if nsv>0:
			for val in friendTree.bdt_score: 
				#print val
				if val == maxvalue:
					is_maxscore.push_back(1)
				else: 
					is_maxscore.push_back(0)
		else: is_maxscore.push_back(-1)
			
		
		friendTree.Fill()
	fout.Write('', TObject.kWriteDelete)
	#fout.cd()
	#friendTree.Write("tFriend", TObject.kOverwrite)

	print "added a tFriend to", fin, fout.GetName()


### give the path to the folder that contains the ntuples for datasets on pnfs 
### or give a path to files with wildcards for signal files on dust / local
list_of_datasets = [
	#'/pnfs/desy.de/cms/tier2/store/user/altews/NTuples/NTuplesV12/16UL_preAPV/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/2023_04_28/230428_123855/0000/'
	'/nfs/dust/cms/user/tewsalex/rootfiles/ntuple_V12/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm*GeV_dm*GeV_Chi20ctau5MM_part*of*_NTuple.root'
	
	]
outpath = '/nfs/dust/cms/user/tewsalex/rootfiles/friendTrees_V12/'

for dataset in list_of_datasets:
	
	print "----"
	print "adding friends to", dataset
	print "----"
	if 'pnfs' in dataset: 
		inpath = dataset+'crab_NTuple_*.root'
		infiles = glob(inpath)
		fout_folder = inpath.split('16UL_preAPV/')[-1].split('/crab_')[0].replace( '/', '_')
	else: 
		inpath = dataset
		infiles = glob(inpath)
		fout_folder = ""

	try:
		# Create the folder
		os.mkdir(outpath+fout_folder)
		print("Folder created successfully.")
	except:
		print("Folder already exists.")

	for afile in infiles: makeFriendTree(afile, outpath+fout_folder+"/")
