import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import glob

process = cms.Process("SVS")
options = VarParsing.VarParsing("analysis")
#  cmsRun construct_secondary_vertices_cfg.py inputFiles="file:///nfs/dust/cms/user/wolfmor/testsamples/ZJetsToNuNu_Zpt-200toInf/B044CEA0-F8C9-E611-8F67-0CC47AD990C4.root"  maxEvents=2


# any default options
options.maxEvents = -1
runOnData = False

options.parseArguments()

# output name
nameout = cms.untracked.vstring(options.inputFiles).value()[0].split("/")[-2]+"_"+cms.untracked.vstring(options.inputFiles).value()[0].split("/")[-1]
#nameout = "/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/test_forCrab.root"
#print cms.untracked.vstring(options.inputFiles).value()[0], nameout

process.source = cms.Source("PoolSource",

    #fileNames = cms.untracked.vstring(options.inputFiles),
    fileNames = cms.untracked.vstring(options.inputFiles),
 )
    
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(options.maxEvents))


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

process.path = cms.Path(process.TrackTag1*process.SecondaryVerticesFromLooseTracks)


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
