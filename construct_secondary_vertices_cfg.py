"""

Runs over AOD files and builds SVs for pairs of OS tracks.

----------------------------------------------------------------------
cmsRun construct_secondary_vertices_cfg.py inputFiles="file1, file2,..." options="..."
minimal example: 

#cmsRun construct_secondary_vertices_cfg.py inputFiles="file:///nfs/dust/cms/user/wolfmor/testsamples/ZJetsToNuNu_Zpt-200toInf/B044CEA0-F8C9-E611-8F67-0CC47AD990C4.root"  maxEvents=2
#cmsRun construct_secondary_vertices_cfg.py inputFiles="file:///pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p77GeV_pu35_part22of25.root" maxEvents=-1
#cmsRun construct_secondary_vertices_cfg.py inputFiles="file:///nfs/dust/cms/user/beinsam/CommonSamples/MC_BSM/CompressedHiggsino/RadiativeMu_2016Full/v2/higgsino94xfull_susyall_mChipm115GeV_dm0p768GeV_pu35_part22of100.root" maxEvents=-1
---------------------------------------------------------------------

options
-----
crab -> is this script run via CRAB or local / on HTCondor; default True
outputFile -> needed if submission NOT via CRAB; if use crab static name for outputs is needed; crab adds job ID to each job; if local / HTCondor is used name of outputfile is given here
maxEvents -> only look at X events; default no limit (-1)
inputFiles -> list of inputfiles; takes wildcards

"""

import FWCore.ParameterSet.Config as cms
#import FWCore.ParameterSet.VarParsing as VarParsing
from FWCore.ParameterSet.VarParsing import VarParsing
import glob

process = cms.Process("SVS")
options = VarParsing("analysis")

options.register('crab', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool)

# any default options
options.maxEvents = -1
options.outputFile = "test.root"

runOnData = False

options.parseArguments()

# output name
nameout = "svfiles_forCrab.root" ### needed for CRAB submission!
if not (options.crab): nameout = options.outputFile

#nameout = options.outputFile

process.source = cms.Source("PoolSource",

    fileNames = cms.untracked.vstring(options.inputFiles),
 )

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))

# Suppress messages that are less important than ERRORs.
process.MessageLogger = cms.Service("MessageLogger",
    destinations = cms.untracked.vstring("cout"),
    cout = cms.untracked.PSet(threshold = cms.untracked.string("ERROR")))

SimpleMemoryCheck = cms.Service("SimpleMemoryCheck",
	ignoreTotal = cms.untracked.int32(1)) 

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
	#trackRecoAlgorithm = cms.InputTag("TrackTag1","selectedTracks"),       
	muons = cms.InputTag("muons"),       
	trackRecoAlgorithm = cms.InputTag("generalTracks"),       
	trackQualities = cms.string("loose"),
	vertexFitter = cms.bool(True),  
	#tkChi2Cut = cms.double(10.),
	tkChi2Cut = cms.double(999.),
	#tkNHitsCut = cms.int32(3), 
	tkNHitsCut = cms.int32(0), 
	tkPtMax = cms.double(999.0),
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

# Run the algorithm.

#process.path = cms.Path(process.TrackTag1*process.SecondaryVerticesFromLooseTracks)
process.path = cms.Path(process.SecondaryVerticesFromLooseTracks)


# Writer to a new file called output.root.  Save only the new K-shorts and the
# primary vertices (for later exercises).
process.output = cms.OutputModule(
    "PoolOutputModule",
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("path")),
    outputCommands = cms.untracked.vstring(
        "drop *",
        "keep *_*_mvaScore_SVS",
        "keep *_*_DcaKshort_SVS",
        "keep *_*_Kshort_SVS",
        "keep *_*_selectedTrackIDs_SVS",
        ),
 
    fileName = cms.untracked.string(nameout)  
    )
process.endpath = cms.EndPath(process.output)
