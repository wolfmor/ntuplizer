import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import glob

process = cms.Process("SVS")
options = VarParsing.VarParsing("analysis")


deltaMChi0LSP = '0p56'
deltaM = '1p12'
runOnData = False
### bad files list


fnames = []
#fnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm*GeV_pu35_part*of25.root") #signal take awaz 1
fnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_20ctau0p005GeV_pu35_part*of500.root") 
#fnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p97GeV_20ctau0p005GeV_pu35_part*.root") # ZtoLL Dataset
#fnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_20ctau0p005m_pu35_part*of100.root")
#fnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/altews/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/test_11/210504_092021/0000/testSumbmission_67.root")
#fnames3GeV = glob.glob("/nfs/dust/cms/user/wolfmor/testsamples/WJetsToLNu/80A3D525-0FBC-E611-AF19-549F35AC7EA4.root") #signal take awaz 10 dm0p57GeV_20ctau0p005GeV_pu35_part
#sfnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_pu35_part*of25.root") #signal take awaz 10
#fnames3GeV = glob.glob("/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_pu35_part*of25.root")
#fnames3GeV = glob.glob("root://cmsxrootd.fnal.gov//store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/002F49A1-2AC9-E611-8D3F-0025905A612C.root") #signal take awaz 10

for name in fnames3GeV:
	fnames.append('file:'+name)
	
if True: 
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p77GeV_20ctau0p001GeV_pu35_part459of500.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p77GeV_20ctau0p001GeV_pu35_part459of500.root")
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_20ctau0p001GeV_pu35_part171of500.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_20ctau0p001GeV_pu35_part171of500.root")
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_pu35_part13of25.root")
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm0p98GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm0p98GeV_pu35_part13of25.root")
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm100GeV_dm0p26GeV_pu35_part5of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm100GeV_dm0p26GeV_pu35_part5of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm3p3GeV_pu35_part15of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm3p3GeV_pu35_part15of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm100GeV_dm4p26GeV_pu35_part25of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm100GeV_dm4p26GeV_pu35_part25of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm4p3GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm4p3GeV_pu35_part13of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm100GeV_dm4p3GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm100GeV_dm4p3GeV_pu35_part13of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm400GeV_dm0p224GeV_pu35_part18of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm400GeV_dm0p224GeV_pu35_part18of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm2p31GeV_pu35_part19of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm2p31GeV_pu35_part19of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm1p8GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm1p8GeV_pu35_part13of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm2p31GeV_pu35_part19of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm2p31GeV_pu35_part19of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm225GeV_dm1p0GeV_pu35_part10of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm225GeV_dm1p0GeV_pu35_part10of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p03GeV_pu35_part10of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p03GeV_pu35_part10of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm6p31GeV_pu35_part17of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm6p31GeV_pu35_part17of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm0p81GeV_pu35_part20of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm0p81GeV_pu35_part20of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm0p51GeV_pu35_part14of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm0p51GeV_pu35_part14of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm0p41GeV_pu35_part11of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm300GeV_dm0p41GeV_pu35_part11of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm0p43GeV_pu35_part9of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm0p43GeV_pu35_part9of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p83GeV_pu35_part9of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p83GeV_pu35_part9of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p03GeV_pu35_part23of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p03GeV_pu35_part23of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p33GeV_pu35_part8of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p33GeV_pu35_part8of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p03GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm500GeV_dm1p03GeV_pu35_part13of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm275GeV_dm0p512GeV_pu35_part5of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm275GeV_dm0p512GeV_pu35_part5of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm1p81GeV_pu35_part11of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm1p81GeV_pu35_part11of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm0p51GeV_pu35_part6of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm0p51GeV_pu35_part6of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm6p31GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm6p31GeV_pu35_part13of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm2p31GeV_pu35_part8of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm2p31GeV_pu35_part8of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm2p31GeV_pu35_part11of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm250GeV_dm2p31GeV_pu35_part11of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p8GeV_pu35_part25of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p8GeV_pu35_part25of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p8GeV_pu35_part23of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p8GeV_pu35_part23of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p6GeV_pu35_part22of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p6GeV_pu35_part22of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p8GeV_pu35_part18of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm0p8GeV_pu35_part18of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm1p8GeV_pu35_part21of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm1p8GeV_pu35_part21of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm3p3GeV_pu35_part16of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm200GeV_dm3p3GeV_pu35_part16of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm0p99GeV_pu35_part14of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm0p99GeV_pu35_part14of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm2p29GeV_pu35_part2of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm2p29GeV_pu35_part2of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm0p39GeV_pu35_part13of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm0p39GeV_pu35_part13of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm3p29GeV_pu35_part8of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm3p29GeV_pu35_part8of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm1p79GeV_pu35_part20of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm1p79GeV_pu35_part20of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm0p99GeV_pu35_part23of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm180GeV_dm0p99GeV_pu35_part23of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm0p79GeV_pu35_part20of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm0p79GeV_pu35_part20of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm1p79GeV_pu35_part18of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm1p79GeV_pu35_part18of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm5p29GeV_pu35_part25of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm5p29GeV_pu35_part25of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm0p29GeV_pu35_part4of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm160GeV_dm0p29GeV_pu35_part4of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm3p28GeV_pu35_part18of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm3p28GeV_pu35_part18of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm3p28GeV_pu35_part3of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm3p28GeV_pu35_part3of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm0p18GeV_pu35_part19of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm140GeV_dm0p18GeV_pu35_part19of25.root")	
	while "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm1p77GeV_pu35_part1of25.root" in fnames: 
		fnames.remove("file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm1p77GeV_pu35_part1of25.root")	


# any default options
options.outputFile = "rootfiles/testDeltaM1p14_MChi115_test.root"
options.inputFiles = "file:/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMuLongerCtau_2016Fast/AODSIM/higgsino94x_susyall_mChipm115GeV_dm0p57GeV_20ctau0p005GeV_pu35_part396of500.root"
options.maxEvents = -1

options.parseArguments()


process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(options.inputFiles),
    #fileNames = cms.untracked.vstring(fnames),
    #fileNames = cms.untracked.vstring(
    #'file:/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/WToLNu/0651F0F0-3209-E911-B7D5-1866DAEA6E1C.root', 
    ##'file:/nfs/dust/cms/user/wolfmor/testsamples/WJetsToLNu/80A3D525-0FBC-E611-AF19-549F35AC7EA4.root',
    #'file:/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/WToLNu/08A8411E-8309-E911-9A59-0025905B85A2.root'),
    #fileNames = cms.untracked.vstring(
    #'file:/nfs/dust/cms/user/wolfmor/testsamples/ZJetsToNuNu_Zpt-200toInf/B044CEA0-F8C9-E611-8F67-0CC47AD990C4.root',
    #'file:/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/ZToNuNu/002F49A1-2AC9-E611-8D3F-0025905A612C.root', 
    #'file:/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/ZToNuNu/0280D3B7-27C9-E611-A435-0025905A605E.root'),
    # duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
    #firstEvent = cms.untracked.uint32(142)
 )
    
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(14000))
#process.firstEvent = cms.untracked.PSet(input = cms.untracked.int32(102000))
#process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(20))

# Suppress messages that are less important than ERRORs.
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

# Load part of the CMSSW reconstruction sequence to make vertexing possible.
# We'll need the CMS geometry and magnetic field to follow the true, non-helical
# shapes of tracks through the detector.
process.load("Configuration/StandardSequences/FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
if (runOnData):
	process.GlobalTag =  GlobalTag(process.GlobalTag, "auto:run2_data")
else:
	process.GlobalTag =  GlobalTag(process.GlobalTag, "100X_upgrade2018_realistic_v10")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")

# Copy most of the vertex producer's parameters, but accept tracks with
# progressively more strict quality.
process.load("RecoVertex.V0Producer.generalV0Candidates_cfi")

process.SecondaryVerticesFromLooseTracks = process.generalV0Candidates.clone(
	#trackRecoAlgorithm = cms.InputTag("MuonTag","selectedTracks"),   ### ToDO: adjust for tt    
	trackRecoAlgorithm = cms.InputTag("TrackTag1","selectedTracks"),       
	muons = cms.InputTag("muons"),       
	#trackRecoAlgorithm = cms.InputTag("generalTracks"),       
	trackQualities = cms.string("loose"),
	vertexFitter = cms.bool(True),  
	#tkChi2Cut = cms.double(10.),
	tkChi2Cut = cms.double(999.),
	#tkNHitsCut = cms.int32(3), 
	tkNHitsCut = cms.int32(0), 
	tkPtMax = cms.double(15.0),
	tkPtCut = cms.double(0.180),
	tkIPSigXYCut = cms.double(0.),
	useVertex = cms.bool(True),	
	sameSign = cms.bool(False),	
	#vtxChi2Cut = cms.double(2.),
	vtxChi2Cut = cms.double(999.),
	vtxDecaySigXYCut = cms.double(0.),
	vtxDecaySigXYZCut = cms.double(-1.),

	#tkDCACut = cms.double(0.004),
	tkDCACut = cms.double(0.5),
	mPiPiCut = cms.double(3), 
	innerHitPosCut = cms.double(-1),
	cosThetaXYCut = cms.double(-1),
	cosThetaXYZCut = cms.double(-1.),

	doKShorts = cms.bool(True),
	doLambdas = cms.bool(False),
	kShortMassCut = cms.double(5),
	lambdaMassCut = cms.double(999)
    )
                                     
process.EventSelection1 = cms.EDFilter("EventSelectionFilter",
                                      inputTracks = cms.InputTag("generalTracks"),
                                      muons = cms.InputTag("muons"), 
                                      selectedJets = cms.InputTag("ak4PFJetsCHS"),
									  PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
									  genParticles = cms.InputTag("genParticles"),
									  selectedMET = cms.InputTag("pfMet"),
                                      minJetPt = cms.double(100),
                                      minMET = cms.double(250),
                                      minDeltaPhi = cms.double(0.5),
                                      maxJetEta = cms.double(2.5),
                                      debug =cms.bool(True),
                                      runOnData =cms.bool(True),
                                     )
process.TrackTag1 = cms.EDProducer("TrackTag",
                                      inputTracks = cms.InputTag("generalTracks"),
                                      useMuonTag = cms.bool(True),
                                      muons = cms.InputTag("muons"),  
									  genParticles = cms.InputTag("genParticles"),
									  particleFlowCandidates = cms.InputTag("particleFlow"),
									  PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
									  selectedJets = cms.InputTag("ak4PFJetsCHS"),
									  dEdxData = cms.InputTag("dedxHarmonic2"),
									  matchToPFs =cms.bool(False),
									  runOnData =cms.bool(True), #false means run on signal, uses signal matched tracks (gen dependence)
									  coneSizeJets = cms.double(0.4), #0.4
									  evaluatedDeltaM = cms.double(1.12), 
									  minPtJets = cms.double(15),
									  minTrackPt = cms.double(0.180),
									  minTrackDxy = cms.double(0.),
									  maxTrackDz = cms.double(3), 
									  #minMVAScoreSingle = cms.double(-0.2),
									  minMVAScoreSingle = cms.double(-1.),
									  maxDeDx = cms.double(999),
									  weightFileSingleTrack = cms.string("/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/BDT_singleBDT/weights/TMVAClassification_BDT70.weights.xml"),
									  debug =cms.bool(False),
                                     )
process.MuonTag = cms.EDProducer("MuonTag",
                                      inputTracks = cms.InputTag("generalTracks"),
									  genParticles = cms.InputTag("genParticles"),
									  muons = cms.InputTag("muons"),
									  PrimaryVertex = cms.InputTag("offlinePrimaryVertices"),
									  minTrackPt = cms.double(0),
									  minTrackDxy = cms.double(0),
									  maxTrackDz = cms.double(1.),
                                      debug =cms.bool(True),
                                     )
# Run the algorithm.
#process.path = cms.Path(process.EventSelection1*process.TrackTag1*process.SecondaryVerticesFromLooseTracks)
#process.path = cms.Path(process.EventSelection1*process.TrackTag1*process.SecondaryVerticesFromLooseTracks)
#process.path = cms.Path(process.TrackTag1*process.SecondaryVerticesFromLooseTracks)
process.path = cms.Path(process.TrackTag1*process.SecondaryVerticesFromLooseTracks)
#process.path = cms.Path(process.EventSelection1*process.SecondaryVerticesFromLooseTracks)
#process.path = cms.Path(process.EventSelection1*process.MuonTag*process.SecondaryVerticesFromLooseTracks)

# Writer to a new file called output.root.  Save only the new K-shorts and the
# primary vertices (for later exercises).
process.output = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_*_selectedTracks",
        "keep *_generalTracks_*_*",
        "keep *_gedGsfElectrons_*_*",
        "keep *_*Gsf*_*_*",
        "keep *_electronGsfTracks_*_*",
        "keep *_hpsPFT*_*_*",
        "keep *_muons_*_*",
        "keep *_*_*_bkgTracks",
        "keep *_*_*_KSHORTS",
        "keep *_*_*_SVS",
        "keep *_offlineBeamSpot_*_*",
        "keep *_offlinePrimaryVertices_*_*",
        "keep *_offlinePrimaryVerticesWithBS_*_*",
        "keep *_genParticles_*_*",
        "keep *_ak4PFJetsCHS_*_*",
        "keep *_ak4PFJets_*_*",
        "keep *_pfCand_*_*",
        "keep *_gedPhotons_*_*",
        "keep *_pfMet_*_*",
        "keep *_particleFlow_*_*",
        "keep *_fixedGridRho*_*_*",
        "keep *_*JetTags*_*_*",
        "keep *_*Cluster*_*_*",
        "keep *Cluster*_*_*_*",
        "keep *_*_*Cluster*_*",
        "keep *_*Track*_*_*",
        "keep *_genMetTrue_*_*",
        "keep *_generator_*_*",
        "keep *_ak4GenJets_*_*",
        #"keep *_TriggerResults_*_*",
        ),
    #fileName = cms.untracked.string("output.root")
    #fileName = cms.untracked.string("2020_05_only_signal_tracks_dca0p5_large.root") #This file is when only passing signal tracks to the SV builder. For this the TrackTag is used
    #fileName = cms.untracked.string("2020_05_all_tracks_dca0p5.root") #This is when building SV out of all sorts of tracks. For this the TrackTag is skiped. 
    #fileName = cms.untracked.string("rootfiles/2021_01_all_tracks_dca0p5_dm"+deltaM+"_m0.2BDT.root")  
    #fileName = cms.untracked.string("rootfiles/2021_01_all_tracks_dca0p5_mChipm200.root")  
    #fileName = cms.untracked.string("rootfiles/testDeltaM1p2.root")  
    #fileName = cms.untracked.string("rootfiles/testDeltaM0p92_ZJetsToNuNu.root")  
    #fileName = cms.untracked.string("rootfiles/test_ZJetsToNuNu.root")  
    #fileName = cms.untracked.string("rootfiles/test_WJetsToLNu.root")  
    #fileName = cms.untracked.string("rootfiles/testDeltaM0p92_WJetsToLNu.root")  
    #fileName = cms.untracked.string("rootfiles/testDeltaM1p14_MChi115_ctauAll.root")  
    fileName = cms.untracked.string(options.outputFile)  
    )
process.endpath = cms.EndPath(process.output)
