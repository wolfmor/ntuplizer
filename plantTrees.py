#! /usr/bin/env python2


"""

Runs over AOD files and writes file with histos and tree.

----------------------------------------------------------------------
python plantTrees.py inputFiles="file1, file2,..." tag="tag1 tag2 ..."

minimal example: 
python plantTrees.py inputFiles="/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p77GeV_pu35_part22of25.root" tag="test, local, era16_07Aug17, signal, skipSVs, fastsim"
python plantTrees.py inputFiles="/nfs/dust/cms/user/wolfmor/testsamples/ZJetsToNuNu_Zpt-200toInf/B044CEA0-F8C9-E611-8F67-0CC47AD990C4.root" tag="test, local, era16_07Aug17"
----------------------------------------------------------------------

tags:
-----

local -> don't open file via xrootd
redirinfn -> use infn redirector for xrootd
redirfnal -> use fnal redirector for xrootd

debug -> Do some helpful printouts for debugging

fastsim -> FastSim correction for MET/JEC
data -> check trigger flags and runnum/lumisec, save json file
cleanleptons -> perform DY cleaning
pmssm -> save pMSSM IDs
signal -> save signal information, if signal files is in inputFiles

skipSVs -> skip collections from SV building, do not save sv information 
max50SVs -> save only 50 svs per event

DATASETNAME -> save GEN info accordingly (e.g. Signal...)

veto_jet100 -> veto events with no jet with pT > 100 GeV
veto_dphimetjets -> veto events with dphi(MET, jets) < 0.5
veto_isolepton -> veto events with isolated leptons

genmatchtracks -> try to find GEN match to every 10th track
genmatchalltracks -> try to find GEN match to every track

era16_07Aug17 -> use corresponding golden json, JECs, jetID, working points,...
era16_UL -> ...
era16_UL_APV -> ...
era17_17Nov2017 -> ...
era18_17Sep2018 -> ...

"""


import sys
import json
import re
import random
import numpy as np
from array import array

import ROOT

ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background

ROOT.gSystem.Load('libFWCoreFWLite.so')
ROOT.gSystem.Load('libDataFormatsFWLite.so')
ROOT.FWLiteEnabler.enable()

from ROOT import gROOT, gSystem, FWLiteEnabler, TFile, TH1F, TMath, TLorentzVector, TH2F, TTree, TVector3, TRandom, TCanvas, TLegend, TColor, TEfficiency, Math, TMVA
from math import ceil,  fabs, sqrt, acos, cos, asin, degrees, sin, pi
import random
import re
import sys

from DataFormats.FWLite import Events, Lumis, Handle
from FWCore.ParameterSet.VarParsing import VarParsing

from DataFormats.Candidate import *
#from DataFormats.Candidate import VertexCompositeCandidate
from commons import *
from copy import copy, deepcopy

def cleanZllEvent(zl1Idx, zl2Idx, collection, tracks, pfcands, jets, met, hZllLeptonPt, hZllDrTrack, hZllDrPfc, hZllDrJet):
    """DY cleaning: replace leptons with "neutrinos" and update collections: clean jets, tracks, pfcands and adapt MET.
    """

    leptonscleaning = [collection[zl1Idx], collection[zl2Idx]]

    badtracks = []
    badpfcands = []
    badjets = []
    for l in leptonscleaning:

        hZllLeptonPt.Fill(l.pt())


        idx, drmin = findMatch_track_old(l, tracks)

        hZllDrTrack.Fill(drmin)

        if idx != -1 and drmin < 0.05:

            badtracks.append(idx)


        idx, drmin = findMatch_pfc_old(l, pfcands)

        hZllDrPfc.Fill(drmin)

        if idx != -1 and drmin < 0.05:

            badpfcands.append(idx)


        idx, drmin = findMatch_jet_old(l, jets)

        hZllDrJet.Fill(drmin)

        if idx != -1 and drmin < 0.2:

            badjets.append(idx)


        met.setP4(met.p4() + ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(l.px(), l.py(), 0, l.energy()))

    jets = [j for (ij, j) in enumerate(jets) if ij not in badjets]
    tracks = [t for (it, t) in enumerate(tracks) if it not in badtracks]
    pfcands = [p for (ip, p) in enumerate(pfcands) if ip not in badpfcands]

    collection = [c for (ic, c) in enumerate(collection) if ic not in [zl1Idx, zl2Idx]]

    if not len(badtracks) == 2: tracks = None

    return collection, tracks, pfcands, jets, met


'''
###############################################################################################
# define functions/class for jet energy corrections
# from: https://github.com/cmsb2g/B2GDAS/blob/master/test/b2gdas_fwlite.py
###############################################################################################
'''


def createJEC(jecSrc, jecLevelList, jetAlgo):

    jecParameterList = ROOT.vector('JetCorrectorParameters')()

    # Load the different JEC levels (the order matters!)
    for jecLevel in jecLevelList:
        jecParameter = ROOT.JetCorrectorParameters('%s_%s_%s.txt' % (jecSrc, jecLevel, jetAlgo))
        jecParameterList.push_back(jecParameter)

    # Chain the JEC levels together
    return ROOT.FactorizedJetCorrector(jecParameterList)


def getJEC(jecSrc, jet, area, rho, nPV):

    jecSrc.setJetEta(jet.Eta())
    jecSrc.setJetPt(jet.Perp())
    jecSrc.setJetE(jet.E())
    jecSrc.setJetA(area)
    jecSrc.setRho(rho)
    jecSrc.setNPV(nPV)

    return jecSrc.getCorrection()


class DataJEC:

    JECList = []

    def __init__(self, inputmap, jettype):
        for minrun, maxrun, version in inputmap:
            JECMap = {}
            JECMap['jecAK4'] = createJEC('/nfs/dust/cms/user/wolfmor/JECs/'+version+'/'+version, ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
            self.JECList.append([minrun, maxrun, JECMap])

    def GetJECMap(self, run):
        for minrun,maxrun,returnmap in self.JECList:
            if run >= minrun and run <= maxrun:
                return returnmap
        raise Exception('Error! Run '+str(run)+' not found in run ranges')

    def jecAK4(self, run):
        JECMap = self.GetJECMap(run)
        return JECMap['jecAK4']



'''
###############################################################################################
# get input arguments
###############################################################################################
'''

options = VarParsing('python')
options.parseArguments()



nMaxTracksPerEvent = 10000

# TODO: check if saveoutputfile
saveOutputFile = True

if 'test' in options.tag: isTest = True
else: isTest = False
nEventsTest = 100 # number of events that are analyzed in case of test
printevery = 1

# TODO: check thresholds for "new" matching
matchingDrThreshold = 0.05
matchingDxyzThreshold = 0.2

# TODO: check met threshold for event selection
metthreshold = 200
metthresholdtrackgenmatch = 200  # only relevant if "genmatch(all)tracks" in options.tag


'''
###############################################################################################
# define output file with histos and variables for tree
###############################################################################################
'''

if True:

    if saveOutputFile:
        nameout = 'NTuple'
        if len(options.inputFiles) == 1: nameout = options.inputFiles[0].split('/')[-1].strip().replace('.root', '') + '_' + nameout
        if isTest: nameout += '_test'
        if 'skipSVs' in options.tag: nameout += '_noSVs'
        fout = ROOT.TFile(nameout + '.root', 'recreate')

        fout.cd()

    tCounter = ROOT.TTree('tCounter', 'tCounter')

    tEvent = ROOT.TTree('tEvent', 'tEvent')


    event_level_var_names = []

    if 'pmssm' in options.tag:
        event_level_var_names += [('pMSSMid1', 'F'), ('pMSSMid2', 'F')]

    var_names_gen_signal = [
        ('deltamFile', 'F'), ('mchipmFile', 'F'), ('mstopFile', 'F')

        , ('chiC1_deltamN1', 'F'), ('chiC1_m', 'F')
        , ('chiN2_deltamN1', 'F'), ('chiN2_m', 'F')

        , ('n_chiC1', 'I'), ('n_chiN2', 'I'), ('n_chiN1', 'I')
        , ('n_chiC1Daughter', 'I'), ('n_chiN2Daughter', 'I')
        , ('n_chiDaughter', 'I')

        , ('stop_mass', 'F'), ('antistop_mass', 'F')
        , ('stop_pt', 'F'), ('antistop_pt', 'F')
        , ('stop_eta', 'F'), ('antistop_eta', 'F')
        , ('stop_phi', 'F'), ('antistop_phi', 'F')
        , ('stop_decay', 'F'), ('antistop_decay', 'F')
        ]
    event_level_var_names += var_names_gen_signal

    var_names_gen_background = [
        ('gen_met_pt', 'F'), ('gen_met_phi', 'F')
        , ('gen_ht', 'F'), ('gen_ht5', 'F'), ('gen_htMiss', 'F')
        , ('gen_neutrinoSumPt', 'F')

        , ('n_zGamma', 'I'), ('zGamma_pdgId', 'F')
        , ('zGamma_pt', 'F'), ('zGamma_eta', 'F'), ('zGamma_phi', 'F')
        , ('zGamma_neutrinoSumPt', 'F'), ('zGamma_tauDecayMode', 'F')
        , ('n_zDaughter', 'I')

        , ('n_wBoson', 'I'), ('wBoson_pdgId', 'F')
        , ('wBoson_pt', 'F'), ('wBoson_eta', 'F'), ('wBoson_phi', 'F')
        , ('wBoson_neutrinoPt', 'F')
        , ('wBoson_tauDecayMode', 'F')
        , ('wBoson_tauDecaylengthXYZ', 'F'), ('wBoson_tauDecaylengthXY', 'F'), ('wBoson_tauDecaylengthZ', 'F')
        , ('wBoson_tauPt', 'F'), ('wBoson_tauEta', 'F'), ('wBoson_tauPhi', 'F')
        , ('wBoson_tauPtVis', 'F'), ('wBoson_tauEtaVis', 'F'), ('wBoson_tauPhiVis', 'F')
        , ('n_wDaughter', 'I')

        , ('n_genParticle', 'I')
        ]
    event_level_var_names += var_names_gen_background

    var_names_cleaning = [
        ('cleaning_electronsCleaned', 'I'), ('cleaning_muonsCleaned', 'I')
        , ('cleaning_invm', 'F'), ('cleaning_zPt', 'F')
        , ('cleaning_l1Pt', 'F'), ('cleaning_l2Pt', 'F')
        , ('cleaning_l1Eta', 'F'), ('cleaning_l2Eta', 'F')
        , ('cleaning_l1Phi', 'F'), ('cleaning_l2Phi', 'F')
        , ('cleaning_l1dBetaAbsIso', 'F'), ('cleaning_l2dBetaAbsIso', 'F')
        , ('cleaning_l1dBetaRelIso', 'F'), ('cleaning_l2dBetaRelIso', 'F')
        , ('cleaning_metPtBeforeCleaning', 'F'), ('cleaning_metPhiBeforeCleaning', 'F')
    ]
    event_level_var_names += var_names_cleaning

    var_names_event = [
        ('cutflow', 'I'), ('random', 'I')

        , ('crossSection', 'F'), ('numSimEvents', 'F')

        , ('weight_fastSimBug', 'F')
        , ('weight_PU_FastFull', 'F'), ('weight_PU_FastFull_rebin', 'F')
        , ('weight_PU_SigBkg', 'F'), ('weight_PU_SigBkg_rebin', 'F')

        , ('n_pv', 'I'), ('rho', 'F')

        , ('met_pt', 'F'), ('met_phi', 'F')
        , ('met_ptNoFastSimCorr', 'F'), ('met_phiNoFastSimCorr', 'F')

        , ('ht', 'F'), ('ht5', 'F'), ('htMiss', 'F')

        , ('n_genJet', 'I')

        , ('badJets_n', 'I'), ('badJets_minEta', 'F'), ('badJets_nForEventVeto', 'I')
        , ('badJets_lepVeto_n', 'I'), ('badJets_lepVeto_minEta', 'F'), ('badJets_lepVeto_nForEventVeto', 'I')

        , ('hasISRJet', 'I'), ('leadingJet_pt', 'F'), ('leadingJet_eta', 'F'), ('leadingJet_phi', 'F')
        ,('JetMetdeltaPhi1', 'F'),('JetMetdeltaPhi2', 'F'),('JetMetdeltaPhi3', 'F'),('JetMetdeltaPhi4', 'F')
        ,('JetPt1', 'F'),('JetPt2', 'F'),('JetPt3', 'F'),('JetPt4', 'F')
        ,('JetEta1', 'F'),('JetEta2', 'F'),('JetEta3', 'F'),('JetEta4', 'F')
        , ('n_jet', 'I')
        , ('n_jet_15', 'I'), ('n_jet_30', 'I'), ('n_jet_50', 'I'), ('n_jet_100', 'I'), ('n_jet_200', 'I')

        , ('n_jet_30_btagloose', 'I'), ('n_jet_15_btagloose', 'I')
        , ('n_jet_30_btagmedium', 'I'), ('n_jet_15_btagmedium', 'I')
        , ('n_jet_30_btagtight', 'I'), ('n_jet_15_btagtight', 'I')

        , ('mtMetLeadingJet', 'F')
        , ('dphiminMetJets', 'F')

        , ('n_photon', 'I'), ('n_photon_iso', 'I')
        , ('n_pfLepton', 'I'), ('n_pfLepton_iso', 'I')
        , ('n_electron', 'I'), ('n_electron_iso', 'I')
        , ('n_muon', 'I'), ('n_muon_iso', 'I')
        , ('n_lepton', 'I'), ('n_lepton_iso', 'I')
        , ('n_tau', 'I'), ('n_tau_vloose', 'I'), ('n_tau_loose', 'I'), ('n_tau_medium', 'I'), ('n_tau_tight', 'I'), ('n_tau_vtight', 'I'), ('n_tau_vvtight', 'I')
        , ('n_tau_20', 'I'), ('n_tau_20_vloose', 'I'), ('n_tau_20_loose', 'I'), ('n_tau_20_medium', 'I'), ('n_tau_20_tight', 'I'), ('n_tau_20_vtight', 'I'), ('n_tau_20_vvtight', 'I')
        , ('n_track_total', 'I'), ('n_track_basic', 'I'), ('n_track', 'I')
        , ('numSVs', 'I'), ('n_sv_total', 'I'), ('n_sv', 'I'),
        ]
    event_level_var_names += var_names_event

    var_names_data = [
        ('runNum', 'F'), ('lumiSec', 'F'), ('eventNum', 'F')
        ]
    event_level_var_names += var_names_data

    event_level_var_array = {}
    for n in event_level_var_names:
        event_level_var_array[n[0]] = array(n[1].lower(), [0])
        tEvent.Branch(nice_string(n[0]), event_level_var_array[n[0]], nice_string(n[0]) + '/' + n[1])


    # TODO: implement the correct MET filters:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    trigger_flags = []
    if 'era16_07Aug17' in options.tag:
        trigger_flags = [
            'goodVertices'
            , 'globalSuperTightHalo2016Filter'
            , 'HBHENoiseFilter'
            , 'HBHENoiseIsoFilter'
            , 'EcalDeadCellTriggerPrimitiveFilter'
            , 'BadPFMuonFilter'
            # , 'BadChargedCandidateFilter'
            , 'eeBadScFilter'
        ]
    elif 'era16_UL' in options.tag:
        pass
    elif 'era17_17Nov2017' in options.tag:
        pass
    elif 'era18_17Sep2018' in options.tag:
        pass
    else:
        raise NotImplementedError('MET filters: era unknown or not specified')

    """
    flags available in /nfs/dust/cms/user/wolfmor/testsamples/DataMET_16C/8AFF8257-2BAF-E711-BBAC-0CC47A4D7614.root
    
    raw2digi_step
    L1Reco_step
    reconstruction_step
    EXOMONOPOLEPath
    condPath
    pathALCARECOHcalCalNoise
    eventinterpretaion_step
    HBHENoiseFilter
    HBHENoiseIsoFilter
    CSCTightHaloFilter
    CSCTightHaloTrkMuUnvetoFilter
    CSCTightHalo2015Filter
    globalTightHalo2016Filter
    globalSuperTightHalo2016Filter
    HcalStripHaloFilter
    hcalLaserEventFilter
    EcalDeadCellTriggerPrimitiveFilter
    EcalDeadCellBoundaryEnergyFilter
    goodVertices
    eeBadScFilter
    ecalLaserCorrFilter
    trkPOGFilters
    chargedHadronTrackResolutionFilter
    muonBadTrackFilter
    BadChargedCandidateFilter
    BadPFMuonFilter
    BadChargedCandidateSummer16Filter
    BadPFMuonSummer16Filter
    trkPOG_manystripclus53X
    trkPOG_toomanystripclus53X
    trkPOG_logErrorTooManyClusters
    METFilters
    dqmoffline_step
    dqmoffline_1_step
    dqmoffline_2_step
    dqmoffline_3_step
    dqmoffline_4_step
    """

    # trigger_flags = [
    #     'globalSuperTightHalo2016Filter'
    #     , 'HBHENoiseFilter'
    #     , 'HBHEIsoNoiseFilter'
    #     , 'eeBadScFilter'
    #     , 'EcalDeadCellTriggerPrimitiveFilter'
    #     , 'ecalBadCalibFilter'
    #     , 'BadChargedHadronFilter'
    #     , 'BadPFMuonFilter'
    #     , 'globalTightHalo2016Filter'
    #     , 'PrimaryVertexFilter'
    #     , 'CSCTightHaloFilter'
    #     ]

    for tf in trigger_flags:
        event_level_var_array[tf] = array('i', [0])
        tEvent.Branch(tf, event_level_var_array[tf], tf + '/I')

    # TODO: add triggers for SingleMuon datastream
    trigger_hlt = [
        'HLT_PFMET90_PFMHT90_IDTight_v'
        , 'HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned_v'
        , 'HLT_PFMET100_PFMHT100_IDTight_v'
        , 'HLT_PFMET110_PFMHT110_IDTight_v'
        , 'HLT_PFMET120_PFMHT120_IDTight_v'
        , 'triggerfired'
    ]

    if 'SingleMuon' in options.tag:
        trigger_hlt = [
            'HLT_IsoMu24_v'
            , 'HLT_IsoMu27_v'
            , 'HLT_Mu50_v'
            , 'triggerfired'
        ]

    for t_hlt in trigger_hlt:
        event_level_var_array[t_hlt] = array('i', [0])
        tEvent.Branch(t_hlt, event_level_var_array[t_hlt], t_hlt + '/I')

    var_names_chiC1 = [
        ('chiC1_pt', 'F'), ('chiC1_eta', 'F'), ('chiC1_phi', 'F')
        , ('chiC1_decaylengthXYZ', 'F'), ('chiC1_decaylengthXY', 'F'), ('chiC1_chidecaylengthZ', 'F')
        , ('chiC1_log10(decaylengthXYZ)', 'F'), ('chiC1_log10(decaylengthXY)', 'F'), ('chiC1_log10(chidecaylengthZ)', 'F')

        , ('chiC1_hasPion', 'I'), ('chiC1_hasMatchedTrackPion', 'I')

        , ('chiC1_pionPt', 'F'), ('chiC1_pionEta', 'F'), ('chiC1_pionPhi', 'F'), ('chiC1_pionCharge', 'F')
        , ('chiC1_pionMatching_tmin', 'F'), ('chiC1_pionMatching_dxyzmin', 'F'), ('chiC1_pionMatching_drmin', 'F')
        , ('chiC1_pionMatching_dxyzminrandom', 'F'), ('chiC1_pionMatching_drminrandom', 'F')
        , ('chiC1_pionMatching_drminold', 'F'), ('chiC1_pionMatching_drminoldrandom', 'F')
        ]
    
    chiC1_var_array = {}
    for n in var_names_chiC1:
        chiC1_var_array[n[0]] = array('f', 10*[0.])
        tEvent.Branch(nice_string(n[0]), chiC1_var_array[n[0]], nice_string(n[0]) + '[n_chiC1]/F')


    var_names_chiN2 = [
        ('chiN2_pt', 'F'), ('chiN2_eta', 'F'), ('chiN2_phi', 'F')
        ,('chiN2_mass', 'F')
        , ('chiN2_decaylengthXYZ', 'F'), ('chiN2_decaylengthXY', 'F'), ('chiN2_chidecaylengthZ', 'F')
        , ('chiN2_log10(decaylengthXYZ)', 'F'), ('chiN2_log10(decaylengthXY)', 'F'), ('chiN2_log10(chidecaylengthZ)', 'F')

		,('chi02_pz', 'F')
        ]
    
    chiN2_var_array = {}
    for n in var_names_chiN2:
        chiN2_var_array[n[0]] = array('f', 10*[0.])
        tEvent.Branch(nice_string(n[0]), chiN2_var_array[n[0]], nice_string(n[0]) + '[n_chiN2]/F')
        
    var_names_chiN2ToObjects = [
		('deltaEtaChi0sToLeptons','F'), ('absdeltaEtaChi0sToLeptons','F'), ('deltaPhiChi0sToLeptons','F'), 
		('Chi0sToPV_phi','F'), ('Chi0sToPV_eta','F'),
	
		('mZstar', 'F'),('abspZstar', 'F'),('absnormalVector', 'F'),
		('beta', 'F'), # angle from normal Vector to pZstar
		('ZstarBoost', 'F'),('Chi01Boost', 'F'),('Chi02Boost', 'F'),
		('num', 'F'),('hasSignalSV', 'F'),
		('mtransverse2', 'F'),('mtransverse2_paper', 'F'),('ptZstar', 'F')
        ]

    for n in var_names_chiN2ToObjects:
        chiN2_var_array[n[0]] = array('f', 10*[0.])
        tEvent.Branch(nice_string(n[0]), chiN2_var_array[n[0]], nice_string(n[0]) + '[n_chiN2]/F')
        
        
    var_names_chiN2Leptons = [
		('leptonID', 'I'),
		('leptonBoost', 'F'),
		('lepton_Low_eta', 'F'),('lepton_High_eta', 'F'),
		('lepton_Low_pt', 'F'),('lepton_High_pt', 'F') 
		]

    for n in var_names_chiN2Leptons:
        chiN2_var_array[n[0]] = array('f', 10*[0.])
        tEvent.Branch(nice_string(n[0]), chiN2_var_array[n[0]], nice_string(n[0]) + '[n_chiN2]/F')
        
        
    var_names_chiN1 = [
        ('chiN1_pt', 'F'), ('chiN1_eta', 'F'), ('chiN1_phi', 'F')
       ,('chiN1_mass', 'F')
               
        ,('chi01_vx', 'F'),('chi01_vy', 'F'),('chi01_vz', 'F')
		,('chi01_dx', 'F'),('chi01_dy', 'F'),('chi01_dz', 'F')
        ]
    
    chiN1_var_array = {}
    for n in var_names_chiN1:
        chiN1_var_array[n[0]] = array('f', 10*[0.])
        tEvent.Branch(nice_string(n[0]), chiN1_var_array[n[0]], nice_string(n[0]) + '[n_chiN1]/F')


    var_names_chidaughter = [
        ('chiDaughter_pdgIdMother', 'I'), ('chiDaughter_pdgId', 'F')
        , ('chiDaughter_pt', 'F'), ('chiDaughter_eta', 'F'), ('chiDaughter_phi', 'F')
        , ('chiDaughter_hasMatchedTrack', 'I')
        ]

    chidaughter_var_array = {}
    for n in var_names_chidaughter:
        chidaughter_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), chidaughter_var_array[n[0]], nice_string(n[0]) + '[n_chiDaughter]/F')


    var_names_zdaughter = [
        ('zDaughter_pdgId', 'F')
        , ('zDaughter_pt', 'F'), ('zDaughter_eta', 'F'), ('zDaughter_phi', 'F')
        ]

    zdaughter_var_array = {}
    for n in var_names_zdaughter:
        zdaughter_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), zdaughter_var_array[n[0]], nice_string(n[0]) + '[n_zDaughter]/F')


    var_names_wdaughter = [
        ('wDaughter_pdgId', 'F')
        , ('wDaughter_pt', 'F'), ('wDaughter_eta', 'F'), ('wDaughter_phi', 'F')
        , ('wDaughter_ptVis', 'F'), ('wDaughter_etaVis', 'F'), ('wDaughter_phiVis', 'F')
        ]

    wdaughter_var_array = {}
    for n in var_names_wdaughter:
        wdaughter_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), wdaughter_var_array[n[0]], nice_string(n[0]) + '[n_wDaughter]/F')
        
    
    var_names_genparticle = [
        ('genParticle_pdgId', 'F'), ('genParticle_motherPdgId', 'F'), ('genParticle_status', 'F')
        , ('genParticle_mass', 'F'), ('genParticle_energy', 'F')
        , ('genParticle_pt', 'F'), ('genParticle_eta', 'F'), ('genParticle_phi', 'F')
        ]

    genparticle_var_array = {}
    for n in var_names_genparticle:
        genparticle_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), genparticle_var_array[n[0]], nice_string(n[0]) + '[n_genParticle]/F')
        

    var_names_pv = [
        ('pv_idx', 'I'), ('pv_numTracks', 'F')
        , ('pv_x', 'F'), ('pv_y', 'F'), ('pv_z', 'F')
        ]

    pv_var_array = {}
    for n in var_names_pv:
        pv_var_array[n[0]] = array('f', 500*[0.])
        tEvent.Branch(nice_string(n[0]), pv_var_array[n[0]], nice_string(n[0]) + '[n_pv]/F')


    var_names_genjet = [
        ('genJet_pt', 'F'), ('genJet_eta', 'F'), ('genJet_phi', 'F'), ('genJet_mass', 'F')
        ]

    genjet_var_array = {}
    for n in var_names_genjet:
        genjet_var_array[n[0]] = array('f', 1000*[0.])
        tEvent.Branch(nice_string(n[0]), genjet_var_array[n[0]], nice_string(n[0]) + '[n_genJet]/F')


    var_names_jet = [
        ('jet_px', 'F'), ('jet_py', 'F'), ('jet_pz', 'F')
        , ('jet_pt', 'F'), ('jet_energy', 'F')
        , ('jet_eta', 'F'), ('jet_phi', 'F')
        , ('jet_numConstituents', 'I')
        , ('jet_btag', 'F')
        , ('jet_drminLepton', 'F'), ('jet_ptClosestLepton', 'F'), ('jet_isLepton', 'F')
        , ('jet_drminGenJet', 'F'), ('jet_ptClosestGenJet', 'F'), ('jet_isGenJet', 'F')
        ]

    jet_var_array = {}
    for n in var_names_jet:
        jet_var_array[n[0]] = array('f', 1000*[0.])
        tEvent.Branch(nice_string(n[0]), jet_var_array[n[0]], nice_string(n[0]) + '[n_jet]/F')


    var_names_photon = [
        ('photon_px', 'F'), ('photon_py', 'F'), ('photon_pz', 'F')
        , ('photon_pt', 'F'), ('photon_energy', 'F')
        , ('photon_eta', 'F'), ('photon_phi', 'F')
        , ('photon_pfAbsIso', 'F'), ('photon_pfAbsIsoMini', 'F')
        , ('photon_chPfAbsIso', 'F'), ('photon_chPfAbsIsoMini', 'F')
        , ('photon_jetIso', 'F'), ('photon_jetIsoMulti', 'F'), ('photon_drminJet', 'F'), ('photon_minvJet', 'F')
        , ('photon_chHadIso', 'F'), ('photon_neHadIso', 'F'), ('photon_photIso', 'F')
        , ('photon_absIso', 'F'), ('photon_relIso', 'F')
        ]

    photon_var_array = {}
    for n in var_names_photon:
        photon_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), photon_var_array[n[0]], nice_string(n[0]) + '[n_photon]/F')


    var_names_electron = [
        ('electron_charge', 'I')
        , ('electron_px', 'F'), ('electron_py', 'F'), ('electron_pz', 'F')
        , ('electron_pt', 'F'), ('electron_energy', 'F')
        , ('electron_eta', 'F'), ('electron_phi', 'F')
        , ('electron_dz', 'F'), ('electron_dxy', 'F')
        , ('electron_pfAbsIso', 'F'), ('electron_pfAbsIsoMini', 'F')
        , ('electron_chPfAbsIso', 'F'), ('electron_chPfAbsIsoMini', 'F')
        , ('electron_jetIso', 'F'), ('electron_jetIsoMulti', 'F'), ('electron_drminJet', 'F'), ('electron_minvJet', 'F')
        , ('electron_chHadIso', 'F'), ('electron_chAllIso', 'F')
        , ('electron_neHadIso', 'F'), ('electron_photIso', 'F')
        , ('electron_puChHadIso', 'F')
        , ('electron_dBetaAbsIso', 'F'), ('electron_dBetaRelIso', 'F')
        ]

    electron_var_array = {}
    for n in var_names_electron:
        electron_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), electron_var_array[n[0]], nice_string(n[0]) + '[n_electron]/F')


    var_names_muon = [
        ('muon_charge', 'I')
        , ('muon_px', 'F'), ('muon_py', 'F'), ('muon_pz', 'F')
        , ('muon_pt', 'F'), ('muon_energy', 'F')
        , ('muon_eta', 'F'), ('muon_phi', 'F')
        , ('muon_dz', 'F'), ('muon_dxy', 'F')
        , ('muon_pfAbsIso', 'F'), ('muon_pfAbsIsoMini', 'F')
        , ('muon_chPfAbsIso', 'F'), ('muon_chPfAbsIsoMini', 'F')
        , ('muon_jetIso', 'F'), ('muon_jetIsoMulti', 'F'), ('muon_drminJet', 'F'), ('muon_minvJet', 'F')
        , ('muon_chHadIso', 'F'), ('muon_chAllIso', 'F')
        , ('muon_neHadIso', 'F'), ('muon_photIso', 'F')
        , ('muon_puChHadIso', 'F')
        , ('muon_dBetaAbsIso', 'F'), ('muon_dBetaRelIso', 'F')
        ]

    muon_var_array = {}
    for n in var_names_muon:
        muon_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), muon_var_array[n[0]], nice_string(n[0]) + '[n_muon]/F')


    var_names_tau = [
        ('tau_charge', 'I')
        , ('tau_px', 'F'), ('tau_py', 'F'), ('tau_pz', 'F')
        , ('tau_pt', 'F'), ('tau_energy', 'F'), ('tau_mass', 'F')
        , ('tau_ptScaled', 'F'), ('tau_energyScaled', 'F'), ('tau_massScaled', 'F')
        , ('tau_eta', 'F'), ('tau_phi', 'F')
        , ('tau_dz', 'F'), ('tau_dxy', 'F')
        , ('tau_chHadIso', 'F'), ('tau_photIso', 'F')
        , ('tau_decayMode', 'F'), ('tau_decayModeFinding', 'F'), ('tau_mvaDiscr', 'F')
        , ('tau_isvloose', 'I'), ('tau_isloose', 'I'), ('tau_ismedium', 'I'), ('tau_istight', 'I'), ('tau_isvtight', 'I'), ('tau_isvvtight', 'I')
        , ('tau_elRejection', 'F'), ('tau_muRejection', 'F')
        , ('tau_leadPfChHadCandPt', 'F'), ('tau_leadPfChHadCandEta', 'F'), ('tau_leadPfChHadCandPhi', 'F')
        , ('tau_genMatch', 'F'), ('tau_genMatchPdgId', 'F'), ('tau_genMatchDr', 'F')
        , ('tau_sanityDm', 'F'), ('tau_sanityRaw', 'F'), ('tau_sanityVloose', 'F')
        ]

    tau_var_array = {}
    for n in var_names_tau:
        tau_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), tau_var_array[n[0]], nice_string(n[0]) + '[n_tau]/F')


    var_names_pflepton = [
        ('pfLepton_charge', 'I')
        , ('pfLepton_pdgId', 'F')
        , ('pfLepton_px', 'F'), ('pfLepton_py', 'F'), ('pfLepton_pz', 'F')
        , ('pfLepton_pt', 'F'), ('pfLepton_energy', 'F')
        , ('pfLepton_eta', 'F'), ('pfLepton_phi', 'F')
        , ('pfLepton_dz', 'F'), ('pfLepton_dxy', 'F')
        , ('pfLepton_pfRelIso', 'F'), ('pfLepton_pfRelIsoMini', 'F')
        , ('pfLepton_chPfRelIso', 'F'), ('pfLepton_chPfRelIsoMini', 'F')
        , ('pfLepton_jetIso', 'F'), ('pfLepton_jetIsoMulti', 'F'), ('pfLepton_drminJet', 'F'), ('pfLepton_minvJet', 'F')
        ]

    pflepton_var_array = {}
    for n in var_names_pflepton:
        pflepton_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), pflepton_var_array[n[0]], nice_string(n[0]) + '[n_pfLepton]/F')


    track_level_var_names = [
        ('track_random', 'I')

        , ('track_charge', 'F')
        , ('track_px', 'F'), ('track_py', 'F'), ('track_pz', 'F')
        , ('track_pt', 'F')
        , ('track_ptError', 'F'), ('track_log10(ptError)', 'F')
        , ('track_ptError/pt', 'F'), ('track_log10(ptError/pt)', 'F')
        , ('track_eta', 'F'), ('track_phi', 'F')
        , ('track_etaError', 'F'), ('track_phiError', 'F')

        , ('track_isPfCand', 'I')
        , ('track_pfCandPt', 'F'), ('track_pfCandEta', 'F'), ('track_pfCandPhi', 'F'), ('track_pfCandPdgId', 'F'), ('track_pfCandParticleId', 'F')
        , ('track_pfCandEnergy', 'F'), ('track_pfCandEcalEnergy', 'F'), ('track_pfCandHcalEnergy', 'F')

        , ('track_associatedPV', 'I')
        , ('track_associatedPU', 'I')
        , ('track_associatedPUAssPV', 'I')

        , ('track_distPVAssPVxy', 'F'), ('track_distPVAssPVz', 'F')

        , ('track_IPsig', 'F'), ('track_IPxyz', 'F'), ('track_IPxy', 'F'), ('track_IPz', 'F')
        , ('track_log10(IPsig)', 'F'), ('track_log10(IPxyz)', 'F'), ('track_log10(IPxy)', 'F'), ('track_log10(IPz)', 'F')

        , ('track_IPsigPU', 'F'), ('track_IPxyzPU', 'F'), ('track_IPxyPU', 'F'), ('track_IPzPU', 'F')
        , ('track_log10(IPsigPU)', 'F'), ('track_log10(IPxyzPU)', 'F'), ('track_log10(IPxyPU)', 'F'), ('track_log10(IPzPU)', 'F')

        , ('track_dxy0', 'F'), ('track_dz0', 'F')
        , ('track_log10(dxy0)', 'F'), ('track_log10(dz0)', 'F')

        , ('track_dxyNoAbs', 'F'), ('track_dzNoAbs', 'F')
        , ('track_dxySign', 'F'), ('track_dzSign', 'F')

        , ('track_dxy', 'F'), ('track_dxyHandmade', 'F'), ('track_dxyPU', 'F')
        , ('track_dz', 'F'), ('track_dzHandmade', 'F'), ('track_dzPU', 'F')
        , ('track_log10(dxy)', 'F'), ('track_log10(dxyHandmade)', 'F'), ('track_log10(dxyPU)', 'F')
        , ('track_log10(dz)', 'F'), ('track_log10(dzHandmade)', 'F'), ('track_log10(dzPU)', 'F')

        , ('track_IPsigAssPV', 'F'), ('track_IPxyzAssPV', 'F'), ('track_IPxyAssPV', 'F'), ('track_IPzAssPV', 'F')
        , ('track_log10(IPsigAssPV)', 'F'), ('track_log10(IPxyzAssPV)', 'F'), ('track_log10(IPxyAssPV)', 'F'), ('track_log10(IPzAssPV)', 'F')

        , ('track_IPsigPUAssPV', 'F'), ('track_IPxyzPUAssPV', 'F'), ('track_IPxyPUAssPV', 'F'), ('track_IPzPUAssPV', 'F')
        , ('track_log10(IPsigPUAssPV)', 'F'), ('track_log10(IPxyzPUAssPV)', 'F'), ('track_log10(IPxyPUAssPV)', 'F'), ('track_log10(IPzPUAssPV)', 'F')

        , ('track_dxyNoAbsAssPV', 'F'), ('track_dzNoAbsAssPV', 'F')
        , ('track_dxySignAssPV', 'F'), ('track_dzSignAssPV', 'F')

        , ('track_dxyAssPV', 'F'), ('track_dxyHandmadeAssPV', 'F'), ('track_dxyPUAssPV', 'F')
        , ('track_dzAssPV', 'F'), ('track_dzHandmadeAssPV', 'F'), ('track_dzPUAssPV', 'F')
        , ('track_log10(dxyAssPV)', 'F'), ('track_log10(dxyHandmadeAssPV)', 'F'), ('track_log10(dxyPUAssPV)', 'F')
        , ('track_log10(dzAssPV)', 'F'), ('track_log10(dzHandmadeAssPV)', 'F'), ('track_log10(dzPUAssPV)', 'F')

        , ('track_dxyError', 'F'), ('track_dzError', 'F')
        , ('track_log10(dxyError)', 'F'), ('track_log10(dzError)', 'F')

        , ('track_pfAbsIso', 'F'), ('track_pfRelIso', 'F'), ('track_drminPf', 'F'), ('track_numneighboursPf', 'I')
        , ('track_chPfAbsIso', 'F'), ('track_chPfRelIso', 'F'), ('track_drminChPf', 'F'), ('track_numneighboursChPf', 'I')
        
        , ('track_tkAbsIso0', 'F'), ('track_tkRelIso0', 'F'), ('track_drminTrack0', 'F'), ('track_numneighboursTrack0', 'I')
        , ('track_tkAbsIso1', 'F'), ('track_tkRelIso1', 'F'), ('track_drminTrack1', 'F'), ('track_numneighboursTrack1', 'I')
        , ('track_tkAbsIso5', 'F'), ('track_tkRelIso5', 'F'), ('track_drminTrack5', 'F'), ('track_numneighboursTrack5', 'I')
        , ('track_tkAbsIso10', 'F'), ('track_tkRelIso10', 'F'), ('track_drminTrack10', 'F'), ('track_numneighboursTrack10', 'I')
        
        , ('track_jetIso0', 'F'), ('track_jetIsoMulti0', 'F'), ('track_drminJet0', 'F'), ('track_btagJet0', 'F'), ('track_minvJet0', 'F')
        , ('track_jetIso10', 'F'), ('track_jetIsoMulti10', 'F'), ('track_drminJet10', 'F'), ('track_btagJet10', 'F'), ('track_minvJet10', 'F')
        , ('track_jetIso15', 'F'), ('track_jetIsoMulti15', 'F'), ('track_drminJet15', 'F'), ('track_btagJet15', 'F'), ('track_minvJet15', 'F')
        , ('track_jetIso20', 'F'), ('track_jetIsoMulti20', 'F'), ('track_drminJet20', 'F'), ('track_btagJet20', 'F'), ('track_minvJet20', 'F')
        , ('track_jetIso30', 'F'), ('track_jetIsoMulti30', 'F'), ('track_drminJet30', 'F'), ('track_btagJet30', 'F'), ('track_minvJet30', 'F')
        
        , ('track_jetIsoNoLepton15', 'F'), ('track_jetIsoMultiNoLepton15', 'F'), ('track_drminJetNoLepton15', 'F'), ('track_btagJetNoLepton15', 'F'), ('track_minvJetNoLepton15', 'F')

        , ('track_bjetLooseIso15', 'F'), ('track_bjetLooseIsoMulti15', 'F'), ('track_drminBjetLoose15', 'F'), ('track_btagBjetLoose15', 'F'), ('track_minvBjetLoose15', 'F')
        , ('track_bjetLooseIso30', 'F'), ('track_bjetLooseIsoMulti30', 'F'), ('track_drminBjetLoose30', 'F'), ('track_btagBjetLoose30', 'F'), ('track_minvBjetLoose30', 'F')
        , ('track_bjetMediumIso15', 'F'), ('track_bjetMediumIsoMulti15', 'F'), ('track_drminBjetMedium15', 'F'), ('track_btagBjetMedium15', 'F'), ('track_minvBjetMedium15', 'F')
        , ('track_bjetMediumIso30', 'F'), ('track_bjetMediumIsoMulti30', 'F'), ('track_drminBjetMedium30', 'F'), ('track_btagBjetMedium30', 'F'), ('track_minvBjetMedium30', 'F')
        , ('track_bjetTightIso15', 'F'), ('track_bjetTightIsoMulti15', 'F'), ('track_drminBjetTight15', 'F'), ('track_btagBjetTight15', 'F'), ('track_minvBjetTight15', 'F')
        , ('track_bjetTightIso30', 'F'), ('track_bjetTightIsoMulti30', 'F'), ('track_drminBjetTight30', 'F'), ('track_btagBjetTight30', 'F'), ('track_minvBjetTight30', 'F')
        
        , ('track_neHadAbsIso0', 'F'), ('track_drminNeHad0', 'F'), ('track_invmNeHad0', 'F')
        , ('track_neHadAbsIso1', 'F'), ('track_drminNeHad1', 'F'), ('track_invmNeHad1', 'F')
        , ('track_neHadAbsIso5', 'F'), ('track_drminNeHad5', 'F'), ('track_invmNeHad5', 'F')
        , ('track_neHadAbsIso10', 'F'), ('track_drminNeHad10', 'F'), ('track_invmNeHad10', 'F')
        
        , ('track_drminPhoton', 'F'), ('track_drminElectron', 'F'), ('track_drminMuon', 'F')

        , ('track_isTauLeadPfChHadCand0', 'I'), ('track_drminTau0', 'F')
        , ('track_mvaDiscrTau0', 'F'), ('track_decayModeTau0', 'F')
        , ('track_dr3highestWpTau0', 'I'), ('track_dr4highestWpTau0', 'I'), ('track_dr5highestWpTau0', 'I')
        
        , ('track_isTauLeadPfChHadCand20', 'I'), ('track_drminTau20', 'F')
        , ('track_mvaDiscrTau20', 'F'), ('track_decayModeTau20', 'F')
        , ('track_dr3highestWpTau20', 'I'), ('track_dr4highestWpTau20', 'I'), ('track_dr5highestWpTau20', 'I')

        , ('track_detaLeadingJet', 'F'), ('track_dphiLeadingJet', 'F')
        , ('track_dphiMet', 'F'), ('track_dphiMetPca', 'F')

        , ('track_chi2', 'F')
        , ('track_quality', 'I')
        , ('track_numValidHits', 'I'), ('track_numLostHits', 'I')

        , ('track_isSignalTrack', 'I'), ('track_isSusyTrack', 'I'), ('track_susyTrackPdgIdMother', 'I'), ('track_susyTrackPdgId', 'I')

        , ('track_hasGenMatch', 'I'), ('track_genMatchTmin', 'F')
        , ('track_genMatchPdgId', 'F'), ('track_genMatchPt', 'F'), ('track_genMatchStatus', 'F')
        , ('track_genMatchIsHardProcess', 'F'), ('track_genMatchIsFromHardProcess', 'F')
        , ('track_genMatchIsPrompt', 'F'), ('track_genMatchIsDirectHadronDecayProduct', 'F'), ('track_genMatchIsDirectTauDecayProduct', 'F')
        , ('track_genMatchMotherPdgId', 'F'), ('track_genMatchMotherPt', 'F'), ('track_genMatchMotherStatus', 'F')
        , ('track_genMatchMotherIsHardProcess', 'F')
        , ('track_genMatchMotherIsTheTau', 'I'), ('track_genMatchMotherTauDecay', 'F')

        , ('track_drminGenTauJet', 'F'), ('track_genTauJetPt', 'F')
        ]

    track_level_var_array = {}
    for n in track_level_var_names:
        track_level_var_array[n[0]] = array('f', nMaxTracksPerEvent*[0.])
        tEvent.Branch(nice_string(n[0]), track_level_var_array[n[0]], nice_string(n[0]) + '[n_track]/F')

        
    SV_level_var_names = [
        ('isSignal','I'),
        ('vtxDCA','I'),
        ('numberofdaughters','I'),
        ('log10vtxChi2','F'), ('vtxChi2Ndof','F'), ('vtxVx','F'), ('vtxVy','F'), ('vtxVz','F'), 
        ('log10vtxdxy','F'), ('log10vtxdz','F'), ('vtxdx','F'), ('vtxdy','F'), ('vtxdz','F'), 
        ('vtxiso','F'), ('vtxdrmin','F'), ('vtxnumneighbours','I'),
        ('PVVtxEta','F'), ('PVVtxPhi','F'),	
        ('deltaEtaPVVtxToTrack_Low','F'), ('deltaEtaPVVtxToTrack_High','F'), 
        ('deltaEtaPVVtxToTrackSum','F'),('absdeltaEtaPVVtxToTrackSum','F'),
        ('deltaPhiPVVtxToTrack_High','F'),('deltaPhiPVVtxToTrack_Low','F'), ('deltaPhiPVVtxToTrackSum','F'),
        ('deltaEtaPVVtxToMET','F'),('absdeltaEtaPVVtxToMET','F'),
        ('deltaPhiPVVtxToMET','F'),('absdeltaPhiPVVtxToMET','F'),
        ('deltaM','F'), ('deltaPhi','F'), ('deltaR','F'), ('deltaEta','F'), ('invMass','F'), 
        ('deltaInnerPos','F'), ('sumCharge','I'), ('vectorSumPt','D'),
        ('vectorSumPxy','D'), ('sumPt','F'),


        ### muon related stuff
        ('muonMatched_Low','F'), ('muonMatched_High','F'),	
        ('numberOfChambers_Low','F'), ('numberOfChambers_High','F'),
        ('numberOfMatchedStations_Low','F'), ('numberOfMatchedStations_High','F'),
        ('numberOfSegments_Low','F'), ('numberOfSegments_High','F'),
        ('isGlobalMuon_Low','F'), ('isGlobalMuon_High','F'),	
        ('isGoodMuon_Low','F'), ('isGoodMuon_High','F'),
        ('isTrackerMuon_Low','F'), ('isTrackerMuon_High','F'),
        ('normalizedChi2Muon_Low','F'), ('normalizedChi2Muon_High','F'),
        ('trackerLayersWithMeasurementMuon_Low','F'), ('trackerLayersWithMeasurementMuon_High','F'),
        ('pixelLayersWithMeasurementMuon_Low','F'), ('pixelLayersWithMeasurementMuon_High','F'),
        ('isSoftMuon_Low','F'), ('isSoftMuon_High','F'),
        ('dxyPVMuon_Low','F'), ('dxyPVMuon_High','F'),
        ('dzPVMuon_Low','F'), ('dzPVMuon_High','F'),


        ('pt_Low','F'), ('eta_Low','F'),
        ('log10PttrackerrorPttrack_Low','F'), ('log10dxy_Low','F'), 
        ('log10dxyerrorDxy_Low','F'), ('log10dzerrorDz_Low','F'),
        ('log10dz_Low','F'), ('nvalidhits_Low','I'), ('absChi2_Low','F'),
        ('mvaSingle_Low','F'), ('quality_Low','I'), ('trackiso_Low','F'), 
        ('trackdrmin_Low','F'), ('tracknumneighbours_Low','I'),('trackisoLoose_Low','F'), 
        ('trackdrminLoose_Low','F'), ('tracknumneighboursLoose_Low','I'),
        ('pt_High','F'), ('eta_High','F'), ('isLow_High','I'), ('isHigh_High','I'), 
        ('log10PttrackerrorPttrack_High','F'), ('log10dxy_High','F'), 
        ('log10dxyerrorDxy_High','F'), ('log10dzerrorDz_High','F'),
        ('log10dz_High','F'), ('nvalidhits_High','I'), ('absChi2_High','F'),
        ('mvaSingle_High','F'), ('quality_High','I'), ('trackiso_High','F'), 
        ('trackdrmin_High','F'), ('tracknumneighbours_High','I'),('trackisoLoose_High','F'), 
        ('trackdrminLoose_High','F'), ('tracknumneighboursLoose_High','I'),
        ('jetdrmin_Low','F'),('jetrelpt_Low','F'),('jetnum_Low','F'), 
        ('jetdrmin_High','F'), ('jetrelpt_High','F'), ('jetnum_High','F'), 

        ('IPsignificance_Low','F'),('IPxyz_Low','F'), ('IPxy_Low','F'),  ('IPz_Low','F'), 
        ('log10IPsignificance_Low','F'),('log10IPxyz_Low','F'), ('log10IPxy_Low','F'),  ('log10IPz_Low','F'), 

        ('IPsignificance_High','F'),('IPxyz_High','F'), ('IPxy_High','F'),  ('IPz_High','F'), 
        ('log10IPsignificance_High','F'),('log10IPxyz_High','F'), ('log10IPxy_High','F'),  ('log10IPz_High','F'), 

        ('IPsignificancePU_Low','F'),('IPxyzPU_Low','F'), ('IPxyPU_Low','F'),  ('IPzPU_Low','F'), 
        ('log10IPsignificancePU_Low','F'),('log10IPxyzPU_Low','F'), ('log10IPxyPU_Low','F'),  ('log10IPzPU_Low','F'), 

        ('IPsignificancePU_High','F'),('IPxyzPU_High','F'), ('IPxyPU_High','F'),  ('IPzPU_High','F'), 
        ('log10IPsignificancePU_High','F'),('log10IPxyzPU_High','F'), ('log10IPxyPU_High','F'),  ('log10IPzPU_High','F'), 

        ('hasTrackMatch_Low', 'F'),('hasTrackMatch_High', 'F'), #ToDo move to gen SV set of variables
        ('hasGenMatch_Low', 'F'),('hasGenMatch_High', 'F'), 
        ('hasGenMatchWithSameMother', 'F'),
        ('events_hasGenMatchWithSameMother1', 'F'),
        ('events_hasGenMatchWithSameMotherm1', 'F'),
        ('pdgID_Low', 'F'), ('pdgID_High', 'F'),
        ('pdgIDMother_Low', 'F'), ('pdgIDMother_High', 'F'),
        ('numDaughtersOfMother_Low', 'F'), ('numDaughtersOfMother_High', 'F'),
        ('hasEWancestor_Low', 'F'),('hasEWancestor_High', 'F'),

        ('mZstar_reco_muon', 'F'),('mZstar_reco_electron', 'F'),('abspZstar_reco', 'F'),('absnormalVector_reco', 'F'),
        ('beta_reco', 'F'), # angle from normal Vector to pZstar
        ('theta_reco', 'F'), # angle from normal Vector to pZstar
        ('error_mtransverse2_reco_muon', 'F'),
        ('mtransverse2_reco_muon', 'F'),('mtransverse2_reco_electron', 'F'),
        ('mtransverse2_hybrid', 'F'),('ptZstar_reco', 'F'),

        # BDT scores
        ('score', 'F'),
        ('deltaPhiMetSumPt', 'F'),('deltaEtaLeadingJetSumPt', 'F'),('deltaPhiLeadingJetSumPt', 'F'),
        #('deltaPhiMetSumPt', 'F'),('Met', 'F'), 
    ]	

    SV_level_var_array = {}
    for n in SV_level_var_names:
        if 'max50SVs' in options.tag: SV_level_var_array[n[0]] = array('f', 51*[0.])
        else: SV_level_var_array[n[0]] = array('f', 500*[0.])
        tEvent.Branch(nice_string(n[0]), SV_level_var_array[n[0]], nice_string(n[0])+ '[n_sv]/F')

    sv_vars = {}
    for n in SV_level_var_names:
        if 'max50SVs' in options.tag: sv_vars[n[0]] = array('f', 51*[0.])
        else: sv_vars[n[0]] = array('f', 500*[0.])    
    
    
    var_names_sv_resolution = [
        ('res_vx', 'F'),('res_vy', 'F'),('res_vz', 'F'),
        ('res_dx', 'F'),('res_dy', 'F'),('res_dz', 'F'),
        ('res_PVx', 'F'),('res_PVy', 'F'),('res_PVz', 'F'),
        ('res_deltaEta','F'), ('res_deltaPhi','F'), 
        ('res_phi','F'), ('res_eta','F'), 
        ('res_ncrossn','F'),('res_alphan','F'),('res_mtransverse2','F'),
        ('res_theta','F'),('gen_theta','F'),('reco_theta','F'),('res_mtransverse','F')
    ]
    for n in var_names_sv_resolution:
        if 'max50SVs' in options.tag: SV_level_var_array[n[0]] = array('f', 51*[0.])
        else: SV_level_var_array[n[0]] = array('f', 500*[0.])
        tEvent.Branch(nice_string(n[0]), SV_level_var_array[n[0]], nice_string(n[0])+ '[n_sv]/F')


    # cutflow histos

    hCutflow = ROOT.TH1F('hCutflow', 'hCutflow', 10, 0., 10.)

    hMetptRaw = ROOT.TH1F('hMetptRaw', 'hMetptRaw', 5000, 0., 5000.)
    hMetptBeforeLeptonCleaning = ROOT.TH1F('hMetptBeforeLeptonCleaning', 'hMetptBeforeLeptonCleaning', 5000, 0., 5000.)
    hMetpt = ROOT.TH1F('hMetpt', 'hMetpt', 5000, 0., 5000.)

    hMindphimetjets = ROOT.TH1F('hMindphimetjets', 'hMindphimetjets', 1000, -1., 4.)
    hMtmetleadingjet = ROOT.TH1F('hMtmetleadingjet', 'hMtmetleadingjet', 5000, 0., 5000.)
    hNjetsbtagmedium = ROOT.TH1F('hNjetsbtagmedium', 'hNjetsbtagmedium', 10, 0., 10.)
    hNumphotons = ROOT.TH1F('hNumphotons', 'hNumphotons', 10, 0., 10.)
    hNumleptons = ROOT.TH1F('hNumleptons', 'hNumleptons', 10, 0., 10.)
    hNumtaus = ROOT.TH1F('hNumtaus', 'hNumtaus', 10, 0., 10.)

    hNumjets = ROOT.TH1F('hNumjets', 'hNumjets', 100, 0., 100.)
    hNumjets30 = ROOT.TH1F('hNumjets30', 'hNumjets30', 20, 0., 20.)
    hNumjets50 = ROOT.TH1F('hNumjets50', 'hNumjets50', 20, 0., 20.)
    hNumjets100 = ROOT.TH1F('hNumjets100', 'hNumjets100', 20, 0., 20.)
    hNumjets200 = ROOT.TH1F('hNumjets200', 'hNumjets200', 20, 0., 20.)

    hBtagjets = ROOT.TH1F('hBtagjets', 'hBtagjets', 400, -2., 2.)


    # Zll cleaning histos

    hZllLeptonPt = ROOT.TH1F('hZllLeptonPt', 'hZllLeptonPt', 1000, 0., 1000.)

    hZllDrTrack = ROOT.TH1F('hZllDrTrack', 'hZllDrTrack', 100, 0., 1.)
    hZllDrPfc = ROOT.TH1F('hZllDrPfc', 'hZllDrPfc', 100, 0., 1.)
    hZllDrJet = ROOT.TH1F('hZllDrJet', 'hZllDrJet', 100, 0., 1.)


    # PV histos

    hNPVsPerEvent = ROOT.TH1F('hNPVsPerEvent', 'hNPVsPerEvent', 100, 0., 100.)

    hPV0x = ROOT.TH1F('hPV0x', 'hPV0x', 500, -0.5, 0.5)
    hPV0y = ROOT.TH1F('hPV0y', 'hPV0y', 500, -0.5, 0.5)
    hPV0z = ROOT.TH1F('hPV0z', 'hPV0z', 500, -25, 25)

    hPVsx = ROOT.TH1F('hPVsx', 'hPVsx', 500, -0.5, 0.5)
    hPVsy = ROOT.TH1F('hPVsy', 'hPVsy', 500, -0.5, 0.5)
    hPVsz = ROOT.TH1F('hPVsz', 'hPVsz', 500, -25, 25)


    # jetID histos

    jetIDvars = [('nhfjet', (50, 0., 1.)), ('nefjet', (50, 0., 1.)), ('chfjet', (50, 0., 1.)), ('cefjet', (50, 0., 1.)), ('mefjet', (50, 0., 1.)), ('nconstituentsjet', (50, 0., 50.)), ('cmjet', (50, 0., 50.)), ('nmjet', (50, 0., 50.))]
    jetIDhistos = {}
    for v in jetIDvars:
        jetIDhistos[v[0] + 'all'] = ROOT.TH1F('h' + v[0] + 'all', 'h' + v[0] + 'all', *v[1])
        jetIDhistos[v[0] + 'pass'] = ROOT.TH1F('h' + v[0] + 'pass', 'h' + v[0] + 'pass', *v[1])


'''
###############################################################################################
# define... stuff and local variables
###############################################################################################
'''

if True:



    jettype = 'AK4PFchs'

    # this discriminator (MVArun2v1DBnewDM) seems to be the one referred to as "MVA2017v2"???
    # see https://cms-nanoaod-integration.web.cern.ch/integration/master/mc102X_doc.html#Tau
    # better use oldDM, see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2
    tauIDdecaymode = 'OldDMs'
    tauIDalgo = 'MVArun2v1DBoldDMwLT'
    # tauIDdecaymode = 'NewDMs'
    # tauIDalgo = 'MVArun2v1DBnewDMwLT'

    if 'era16_07Aug17' in options.tag:

        # https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2016Analysis#Re_reco_datasets_07Aug17
        with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/goldenjson_era16_07Aug17.json') as goldenjsonfile:
            goldenjson = json.load(goldenjsonfile)

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        if 'data' in options.tag:  # data
            jet_energy_corrections = [
                [1, 276811, 'Summer16_07Aug2017BCD_V11_DATA'],
                [276831, 278801, 'Summer16_07Aug2017EF_V11_DATA'],
                [278802, float('inf'), 'Summer16_07Aug2017GH_V11_DATA']]
            DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'fastsim' in options.tag:
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer16_FastSimV1_MC/Summer16_FastSimV1_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
        else:  # FullSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

        # tau energy scale (TES)
        # from https://github.com/cms-tau-pog/TauIDSFs#dm-dependent-tau-energy-scale
        if tauIDalgo == 'MVArun2v1DBoldDMwLT':
            tesfile = ROOT.TFile('/nfs/dust/cms/user/wolfmor/TES/TauES_dm_MVAoldDM2017v2_2016Legacy.root')
            teshist = tesfile.Get('tes')
        else:
            raise NotImplementedError('tauIDalgo unknown or not specified')


    elif 'era16_UL' in options.tag:

        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage
        with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/goldenjson_era16_UL.json') as goldenjsonfile:
            goldenjson = json.load(goldenjsonfile)

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # and https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDataReprocessingUL2016
        if 'data' in options.tag:  # data
            jet_energy_corrections = [
                [1, 276811, 'Summer19UL16APV_RunBCD_V7_DATA'],
                [276831, 278807, 'Summer19UL16APV_RunEF_V7_DATA'],
                [278769, float('inf'), 'Summer19UL16_RunFGH_V7_DATA']]
            DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'fastsim' in options.tag:
            raise NotImplementedError('no JECs for UL FastSim')
        else:  # FullSim
            if 'era16_UL_APV' in options.tag:
                jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer19UL16_V7_MC/Summer19UL16_V7_MC',
                                   ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
            else:
                jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC',
                                   ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

        # tau energy scale (TES)
        # from https://github.com/cms-tau-pog/TauIDSFs#dm-dependent-tau-energy-scale
        if tauIDalgo == 'MVArun2v1DBoldDMwLT':
            tesfile = ROOT.TFile('/nfs/dust/cms/user/wolfmor/TES/TauES_dm_MVAoldDM2017v2_2016Legacy.root')  # TODO: this is not UL but ok...
            teshist = tesfile.Get('tes')
        else:
            raise NotImplementedError('tauIDalgo unknown or not specified')

    elif 'era17_17Nov2017' in options.tag:

        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage
        with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/goldenjson_era17_17Nov2017.json') as goldenjsonfile:
            goldenjson = json.load(goldenjsonfile)

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        if 'data' in options.tag:  # data
            raise NotImplementedError('no JECs yet for 2017 data')
            # TODO: implement JECs for 2017 data
            # jet_energy_corrections = []
            # DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'fastsim' in options.tag:
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
        else:  # FullSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

        # tau energy scale (TES)
        # from https://github.com/cms-tau-pog/TauIDSFs#dm-dependent-tau-energy-scale
        if tauIDalgo == 'MVArun2v1DBoldDMwLT':
            tesfile = ROOT.TFile('/nfs/dust/cms/user/wolfmor/TES/TauES_dm_MVAoldDM2017v2_2017ReReco.root')
            teshist = tesfile.Get('tes')
        else:
            raise NotImplementedError('tauIDalgo unknown or not specified')

    elif 'era18_17Sep2018' in options.tag:

        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage
        with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/goldenjson_era18_17Sep2018.json') as goldenjsonfile:
            goldenjson = json.load(goldenjsonfile)

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        if 'data' in options.tag:  # data
            raise NotImplementedError('no JECs yet for 2018 data')
            # TODO: implement JECs for 2018 data
            # jet_energy_corrections = []
            # DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'fastsim' in options.tag:
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
        else:  # FullSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Autumn18_V19_MC/Autumn18_V19_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

        # tau energy scale (TES)
        # from https://github.com/cms-tau-pog/TauIDSFs#dm-dependent-tau-energy-scale
        if tauIDalgo == 'MVArun2v1DBoldDMwLT':
            tesfile = ROOT.TFile('/nfs/dust/cms/user/wolfmor/TES/TauES_dm_MVAoldDM2017v2_2018ReReco.root')
            teshist = tesfile.Get('tes')
        else:
            raise NotImplementedError('tauIDalgo unknown or not specified')

    else:

        raise NotImplementedError('era unknown or not specified')


'''
###############################################################################################
# define handles and labels
###############################################################################################
'''

if 'pmssm' in options.tag:
    handle_lumis = Handle('GenLumiInfoHeader')
    
elif 'fastsim' not in options.tag:

    handle_trigger_hlt = Handle('edm::TriggerResults')
    label_trigger_hlt = ('TriggerResults', '', 'HLT')

if 'data' in options.tag:

    handle_trigger_flags = Handle('edm::TriggerResults')
    label_trigger_flags = ('TriggerResults', '', 'RECO')

else:

    handle_genparticles = Handle('std::vector<reco::GenParticle>')
    label_genparticles = ('genParticles')

    handle_genmet = Handle('std::vector<reco::GenMET>')
    label_genmet = ('genMetTrue')

    handle_genjets = Handle('std::vector<reco::GenJet>')
    label_genjets = ('ak4GenJets')

handle_tracks = Handle('std::vector<reco::Track>')
label_tracks = ('generalTracks')

handle_pfcands = Handle('std::vector<reco::PFCandidate>')
label_pfcands = ('particleFlow')

handle_photons = Handle('std::vector<reco::Photon>')
label_photons = ('gedPhotons')

handle_electrons = Handle('std::vector<reco::GsfElectron>')
label_electrons = ('gedGsfElectrons')

handle_muons = Handle('std::vector<reco::Muon>')
label_muons = ('muons')

handle_taus = Handle('std::vector<reco::PFTau>')
label_taus = ('hpsPFTauProducer')

handle_taudiscriminatorDM = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorDM = ('hpsPFTauDiscriminationByDecayModeFinding' + tauIDdecaymode)

handle_taudiscriminatorMVA_VLoose = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_VLoose = ('hpsPFTauDiscriminationByVLooseIsolation' + tauIDalgo)

handle_taudiscriminatorMVA_Loose = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_Loose = ('hpsPFTauDiscriminationByLooseIsolation' + tauIDalgo)

handle_taudiscriminatorMVA_Medium = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_Medium = ('hpsPFTauDiscriminationByMediumIsolation' + tauIDalgo)

handle_taudiscriminatorMVA_Tight = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_Tight = ('hpsPFTauDiscriminationByTightIsolation' + tauIDalgo)

handle_taudiscriminatorMVA_VTight = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_VTight = ('hpsPFTauDiscriminationByVTightIsolation' + tauIDalgo)

handle_taudiscriminatorMVA_VVTight = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_VVTight = ('hpsPFTauDiscriminationByVVTightIsolation' + tauIDalgo)

handle_taudiscriminatorMVAraw = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVAraw = ('hpsPFTauDiscriminationByIsolation' + tauIDalgo + 'raw')

handle_taudiscriminatorElectronRej = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorElectronRej = ('hpsPFTauDiscriminationByMVA6VLooseElectronRejection')

handle_taudiscriminatorMuonRej = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMuonRej = ('hpsPFTauDiscriminationByLooseMuonRejection3')

handle_rhos = Handle('double')
label_rhos = ('fixedGridRhoFastjetAll')

handle_pv = Handle('std::vector<reco::Vertex>')
label_pv = ('offlinePrimaryVertices')

handle_met = Handle('std::vector<reco::PFMET>')
label_met = ('pfMet')

handle_btag = Handle('edm::AssociationVector<edm::RefToBaseProd<reco::Jet>,vector<float>,edm::RefToBase<reco::Jet>,unsigned int,edm::helper::AssociationIdenticalKeyReference>')
label_btag = ('pfCombinedSecondaryVertexV2BJetTags')

handle_jets = Handle('std::vector<reco::PFJet>')
label_jets = ('ak4PFJetsCHS')

handle_sv = Handle("std::vector<reco::VertexCompositeCandidate>")
label_sv = ('SecondaryVerticesFromLooseTracks','Kshort', 'SVS')

handle_dca = Handle("std::vector<float>")
label_dca = ('SecondaryVerticesFromLooseTracks','DcaKshort', 'SVS')

handle_mva = Handle("std::vector<float>")
label_mva = ('TrackTag1', 'mvaScore', 'SVS')

handle_selected_tracks = Handle("std::vector<int>")
label_selected_tracks = ('TrackTag1' , 'selectedTrackIDs', 'SVS')


runs = {}
lastlumi = -1
lastrun = -1


'''
###############################################################################################
# start with SV files
###############################################################################################
'''
    
if not 'skipSVs' in options.tag:
    
    localpath = '/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/'

    filesWithSV = np.array([None]*len(options.inputFiles))
    filesWithTrackID = np.array([None]*len(options.inputFiles))
    filesWithDCA = np.array([None]*len(options.inputFiles))
    filesWithMVA = np.array([None]*len(options.inputFiles))
    
    for ifile, f in enumerate(options.inputFiles):
        
        vertexfile = localpath+f.split("/")[-2]+"_"+f.split("/")[-1]

        '''
        ###############################################################################################
        # get sv-level info
        ###############################################################################################
        '''
        if 'local' in options.tag:

            v_fname = vertexfile.strip()
            v_fin = ROOT.TFile.Open(v_fname)
            print vertexfile

        else:

            redir = 'cms-xrd-global.cern.ch'
            if 'redirinfn' in options.tag: redir = 'xrootd-cms.infn.it'
            if 'redirfnal' in options.tag: redir = 'cmsxrootd.fnal.gov'

            v_fname = 'root://' + redir + '/' + vertexfile.strip()
            v_fin = ROOT.TFile.Open(v_fname)

            retry = 1
            while not v_fin or v_fin.IsZombie() or not v_fin.IsOpen():

                print 'retry open file ', retry

                retry += 1
                if retry > 5: break

                v_fin = ROOT.TFile.Open('root://' + redir + '/' + vertexfile.strip())

        v_events = Events(v_fin)
        v_nevents = v_events.size()
        filesWithSV[ifile] = np.array([None]*v_events.size())
        filesWithTrackID[ifile] = np.array([None]*v_events.size())
        filesWithDCA[ifile] = np.array([None]*v_events.size())
        filesWithMVA[ifile] = np.array([None]*v_events.size())

        for v_ievent, v_event in enumerate(v_events):
            
            if isTest and v_ievent >= nEventsTest:
                print 'nEventsTest boundary'
                break
                               
            ### get collections for same event
            v_event.getByLabel(label_sv, handle_sv)
            v_event.getByLabel(label_dca, handle_dca)
            v_event.getByLabel(label_mva, handle_mva)
            v_event.getByLabel(label_selected_tracks, handle_selected_tracks)

            secondaryVertices = handle_sv.product()
            dcas = handle_dca.product()
            mvaScores = handle_mva.product()
            selectedTrackIDs = handle_selected_tracks.product()

            filesWithSV[ifile][v_ievent] = np.array([None]*secondaryVertices.size())
            filesWithTrackID[ifile][v_ievent] = np.array([-1.]*selectedTrackIDs.size())
            filesWithDCA[ifile][v_ievent] = np.array([-1.]*dcas.size())
            filesWithMVA[ifile][v_ievent] = np.array([-1.]*mvaScores.size())
            
            dca = 0.0
            
            for nSV, secondary in enumerate(secondaryVertices):
                
                filesWithSV[ifile][v_ievent][nSV] = deepcopy(secondary)
                dca = float(dcas[nSV])
                filesWithDCA[ifile][v_ievent][nSV] = dca
                #print "file ", ifile, "event ", v_ievent, "SV", nSV,  "dca", dca, filesWithDCA[ifile][v_ievent][nSV]
                
            for idx, index in enumerate(selectedTrackIDs):
                
                filesWithTrackID[ifile][v_ievent][idx] = index
                filesWithMVA[ifile][v_ievent][idx] = mvaScores[idx]
                #print "file ", ifile, "event ", v_ievent, "selected track", idx, "track ID ", index, filesWithTrackID[ifile][v_ievent][idx]

    print "----------Finished loop over SV files----------------"


if 'debug' in options.tag:
    for ifile , afile in enumerate(filesWithSV):
        #print 'File loop ', ifile, afile
        for ievent, event in enumerate(filesWithSV[ifile]):
            
            if isTest and ievent >= nEventsTest:
                print 'nEventsTest boundary'
                break
            #print 'Event loop', ievent, event
            
            #for isv, sv in enumerate(filesWithSV[ifile][ievent]):
                #print 'file', ifile, 'event', ievent, 'SV loop', isv, sv.vx()
                
            for isv, sv in enumerate(filesWithDCA[ifile][ievent]):
                print 'file', ifile, 'event', ievent, 'DCA loop', isv, ' DCA', sv
#sys.exit()
'''
###############################################################################################
# start with AOD files 
###############################################################################################
'''

print "----------Start loop over AOD files----------------"
print ''
print 'n input files: ' + str(len(options.inputFiles))

for ifile, f in enumerate(options.inputFiles):

    print ''
    print f

    '''
    ###############################################################################################
    # open file
    ###############################################################################################
    '''

    if 'local' in options.tag:

        fname = f.strip()
        fin = ROOT.TFile.Open(fname)

    else:

        redir = 'cms-xrd-global.cern.ch'
        if 'redirinfn' in options.tag: redir = 'xrootd-cms.infn.it'
        if 'redirfnal' in options.tag: redir = 'cmsxrootd.fnal.gov'

        fname = 'root://' + redir + '/' + f.strip()
        fin = ROOT.TFile.Open(fname)

        retry = 1
        while not fin or fin.IsZombie() or not fin.IsOpen():

            print 'retry open file ', retry

            retry += 1
            if retry > 5: break

            fin = ROOT.TFile.Open('root://' + redir + '/' + f.strip())

    events = Events(fin)

    nevents = events.size()

    # don't just silently skip a file
    # try:
    #     nevents = events.size()
    # except:
    #     print 'skipping file ' + f
    #     continue

    print '### with ' + str(nevents) + ' events'
    print '### printing every ' + str(printevery) + '. event'

    if saveOutputFile: fout.cd()

    nEventsPerFile = 0

    phifirsttrack = 0
    etafirsttrack = 0

    if 'pmssm' in options.tag:
        lumis = Lumis(fname)
        pMSSMname = ''
        pMSSMid1 = -1
        pMSSMid2 = -1
        

    '''
    ###############################################################################################
    # event loop
    ###############################################################################################
    '''

    for ievent, event in enumerate(events):
            
        if isTest and ievent >= nEventsTest:
            print 'nEventsTest boundary'
            break

        if 'local' not in options.tag:
            if fin.IsZombie() or not fin.IsOpen():
                print 'file not usable'
                sys.exit(1)

        if saveOutputFile and ievent % 100 == 0: fout.Write('', ROOT.TObject.kWriteDelete)

        if ievent % printevery == 0: print 'analyzing event %d of %d' % (ievent, nevents)

        tCounter.Fill()

        cutflow = 0
        hCutflow.Fill(cutflow)

        random.seed()
        event_level_var_array['random'][0] = random.randrange(10)

        '''
        ###############################################################################################
        # check runnum/lumisec and add to list for json
        ###############################################################################################
        '''

        runnum = event.eventAuxiliary().run()
        lumisec = event.eventAuxiliary().luminosityBlock()
        eventnum = event.eventAuxiliary().event()

        if 'pmssm' in options.tag:

            if lumisec != lastlumi:
                lastlumi = lumisec
                for lumi in lumis:
                    lumi.getByLabel('generator', handle_lumis)
                    if lumi.luminosityBlockAuxiliary().luminosityBlock() != lumisec: continue
                    pMSSMname = handle_lumis.product().configDescription()
                    break
                pMSSMid1 = float(pMSSMname.replace('.slha', '').split('_')[-2])
                pMSSMid2 = float(pMSSMname.replace('.slha', '').split('_')[-1])

            event_level_var_array['pMSSMid1'][0] = pMSSMid1
            event_level_var_array['pMSSMid2'][0] = pMSSMid2

        if 'data' in options.tag:

            if runnum != lastrun or lumisec != lastlumi:
                if str(runnum) in goldenjson: goodlumisecs = goldenjson[str(runnum)]
                else: goodlumisecs = []
                isgood = False
                for gls in goodlumisecs:
                    if lumisec in range(gls[0], gls[1]+1): isgood = True
                    if isgood: break

            if not isgood:
                lastrun = runnum
                lastlumi = lumisec
                continue
        # ########################################################################################### veto

            if runnum != lastrun:
                lastrun = runnum
                if runnum not in runs:
                    runs[runnum] = []
            if lumisec != lastlumi:
                lastlumi = lumisec
                if runnum not in runs:
                    runs[runnum] = []
                if lumisec not in runs[runnum]:
                    runs[runnum].append(lumisec)

        cutflow = 1
        hCutflow.Fill(cutflow)

        event_level_var_array['runNum'][0] = runnum
        event_level_var_array['lumiSec'][0] = lumisec
        event_level_var_array['eventNum'][0] = eventnum

        numC1 = 0
        numN2 = 0
        numN1 = 0

        chipmmFILE = -1
        deltamFILE = -1
        mstopFILE = -1

        chipmmGEN = -1
        deltamGEN = -1

        chiN2mGEN = -1
        deltamN2GEN = -1

        chipmnumdaughters = 0
        chiN2numdaughters = 0
        numchidaughters = 0
        
        theChi01 = None
        theChi02 = None

        '''
        ###############################################################################################
        # get products
        ###############################################################################################
        '''

        event.getByLabel(label_rhos, handle_rhos)
        rhos = handle_rhos.product()
        if not len(rhos) > 0: continue
        rho = rhos[0]

        event.getByLabel(label_pv, handle_pv)
        primaryvertices = handle_pv.product()
        primaryvertices = [pv for pv in primaryvertices if pv.isValid()]
        n_pv = len(primaryvertices)
        if not n_pv > 0: continue
        pv_pos = primaryvertices[0].position()

        event.getByLabel(label_tracks, handle_tracks)
        tracks = handle_tracks.product()
        if not len(tracks) > 0: continue

        # ########################################################################################### veto

        cutflow = 2
        hCutflow.Fill(cutflow)

        '''
        ###############################################################################################
        # trigger
        ###############################################################################################
        '''

        trigger_flags_accept = {}
        for tf in trigger_flags:
            trigger_flags_accept[tf] = -1

        allfine = True

        # TODO: do this also for MC??
        if 'data' in options.tag:

            event.getByLabel(label_trigger_flags, handle_trigger_flags)
            triggerresults = handle_trigger_flags.product()

            triggernames = event.object().triggerNames(triggerresults)

            for i in range(triggerresults.size()):

                tn = triggernames.triggerName(i).replace('Flag_', '')
                if tn in trigger_flags:
                    if triggerresults.accept(i):
                        trigger_flags_accept[tn] = 1
                    else:
                        trigger_flags_accept[tn] = 0
                        allfine = False

        for tf in trigger_flags:
            event_level_var_array[tf][0] = trigger_flags_accept[tf]

        if not allfine: continue

        # ########################################################################################### veto

        cutflow = 3
        hCutflow.Fill(cutflow)

        trigger_hlt_accept = {}
        for t_hlt in trigger_hlt:
            if t_hlt == 'triggerfired': continue
            trigger_hlt_accept[t_hlt] = -1

        triggerfired = 0

        if 'fastsim' not in options.tag:

            event.getByLabel(label_trigger_hlt, handle_trigger_hlt)
            triggerresults_hlt = handle_trigger_hlt.product()

            triggernames_hlt = event.object().triggerNames(triggerresults_hlt)

            for i in range(triggerresults_hlt.size()):

                tn_hlt = triggernames_hlt.triggerName(i)

                for t_hlt in trigger_hlt:

                    if re.match(r'' + t_hlt + '.*', tn_hlt):
                        if triggerresults_hlt.accept(i):
                            trigger_hlt_accept[t_hlt] = 1
                            triggerfired = 1
                        else:
                            trigger_hlt_accept[t_hlt] = 0

        for t_hlt in trigger_hlt:
            if t_hlt == 'triggerfired': continue
            event_level_var_array[t_hlt][0] = trigger_hlt_accept[t_hlt]
        event_level_var_array['triggerfired'][0] = triggerfired


        if 'data' not in options.tag:
            event.getByLabel(label_genparticles, handle_genparticles)
            event.getByLabel(label_genmet, handle_genmet)
            event.getByLabel(label_genjets, handle_genjets)
        event.getByLabel(label_pfcands, handle_pfcands)
        event.getByLabel(label_jets, handle_jets)
        event.getByLabel(label_met, handle_met)
        event.getByLabel(label_btag, handle_btag)
        event.getByLabel(label_photons, handle_photons)
        event.getByLabel(label_electrons, handle_electrons)
        event.getByLabel(label_muons, handle_muons)
        event.getByLabel(label_taus, handle_taus)

        if 'data' not in options.tag:
            genparticles = handle_genparticles.product()
            genmet = handle_genmet.product().front()
            genjets = handle_genjets.product()
        pfcands = handle_pfcands.product()
        jets = handle_jets.product()
        met = handle_met.product().front()
        nbtags = len(handle_btag.product())
        photons = handle_photons.product()
        electrons = handle_electrons.product()
        muons = handle_muons.product()
        taus = handle_taus.product()
        

        event.getByLabel(label_taudiscriminatorDM, handle_taudiscriminatorDM)
        taudiscriminatorDM = handle_taudiscriminatorDM.product()

        event.getByLabel(label_taudiscriminatorMVA_VLoose, handle_taudiscriminatorMVA_VLoose)
        taudiscriminatorMVA_VLoose = handle_taudiscriminatorMVA_VLoose.product()

        event.getByLabel(label_taudiscriminatorMVA_Loose, handle_taudiscriminatorMVA_Loose)
        taudiscriminatorMVA_Loose = handle_taudiscriminatorMVA_Loose.product()

        event.getByLabel(label_taudiscriminatorMVA_Medium, handle_taudiscriminatorMVA_Medium)
        taudiscriminatorMVA_Medium = handle_taudiscriminatorMVA_Medium.product()

        event.getByLabel(label_taudiscriminatorMVA_Tight, handle_taudiscriminatorMVA_Tight)
        taudiscriminatorMVA_Tight = handle_taudiscriminatorMVA_Tight.product()

        event.getByLabel(label_taudiscriminatorMVA_VTight, handle_taudiscriminatorMVA_VTight)
        taudiscriminatorMVA_VTight = handle_taudiscriminatorMVA_VTight.product()

        event.getByLabel(label_taudiscriminatorMVA_VVTight, handle_taudiscriminatorMVA_VVTight)
        taudiscriminatorMVA_VVTight = handle_taudiscriminatorMVA_VVTight.product()

        event.getByLabel(label_taudiscriminatorMVAraw, handle_taudiscriminatorMVAraw)
        taudiscriminatorMVAraw = handle_taudiscriminatorMVAraw.product()

        event.getByLabel(label_taudiscriminatorElectronRej, handle_taudiscriminatorElectronRej)
        taudiscriminatorElectronRej = handle_taudiscriminatorElectronRej.product()
        
        event.getByLabel(label_taudiscriminatorMuonRej, handle_taudiscriminatorMuonRej)
        taudiscriminatorMuonRej = handle_taudiscriminatorMuonRej.product()

### ask Moritz
        tauswithdiscriminators = []
        # tauswithdiscriminators = [
            # (tau, taudiscriminatorDM.value(itau)
                # , taudiscriminatorMVAraw.value(itau)
                # , taudiscriminatorMVA_VLoose.value(itau)
                # , taudiscriminatorMVA_Loose.value(itau)
                # , taudiscriminatorMVA_Medium.value(itau)
                # , taudiscriminatorMVA_Tight.value(itau)
                # , taudiscriminatorMVA_VTight.value(itau)
                # , taudiscriminatorMVA_VVTight.value(itau)
                # , taudiscriminatorElectronRej.value(itau)
                # , taudiscriminatorMuonRej.value(itau)
                # , taudiscriminatorDM.key(itau).get().pt()
                # , taudiscriminatorMVAraw.key(itau).get().pt()
                # , taudiscriminatorMVA_VLoose.key(itau).get().pt())
            # for itau, tau in enumerate(taus)
        # ]

        # first dummy entry needed for jitted jet iso calculation function
        btagvalues = [(-2., 0., 0.)]


        '''
        ###############################################################################################
        # lepton and photon selection
        ###############################################################################################
        '''

        photons = [p for p in photons
                   if passesPhotID(p, wp='loose')
                   and p.pt() > 15
                   and abs(p.eta()) < 2.5]

        electrons = [e for e in electrons
                     if passesEleID(e, pv_pos, rho, wp='veto')
                     and e.pt() > 10
                     and abs(e.eta()) < 2.5]

        muons = [m for m in muons
                 if passesMuonID(m, wp='loose')
                 and m.pt() > 10
                 and abs(m.eta()) < 2.4]

        tauswithdiscriminators = [t for t in tauswithdiscriminators
                                  if t[1] > 0.5  # decaymodefinding
                                  and t[3] > 0.5  # vloose working point
                                  and t[9] > 0.5  # electron rejection
                                  and t[10] > 0.5  # muon rejection
                                  # and t[0].pt() > 20
                                  and abs(t[0].eta()) < 2.3]

        '''
        ###############################################################################################
        # jetID and JECs
        ###############################################################################################
        '''

        jetsP4Raw = []
        jetsP4Corr = []
        jetsIdxGood = []
        numBadJets = 0
        minetaabsbadjets = 9
        numBadJetsEventVeto = 0
        numBadJetsLepVeto = 0
        minetaabsbadjetsLepVeto = 9
        numBadJetsLepVetoEventVeto = 0
        for ijet, jet in enumerate(jets):

            jetP4Raw = ROOT.TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy())
            jetsP4Raw.append(jetP4Raw)

            if 'data' in options.tag:
                correction = getJEC(DataJECs.jecAK4(runnum), jetP4Raw, jet.jetArea(), rho, n_pv)
            else:
                correction = getJEC(jecAK4, jetP4Raw, jet.jetArea(), rho, n_pv)

            jetP4Corr = jetP4Raw * correction
            jetsP4Corr.append(jetP4Corr)

            nhf = jet.neutralHadronEnergy() / jetP4Raw.E()
            nef = jet.neutralEmEnergy() / jetP4Raw.E()
            chf = jet.chargedHadronEnergy() / jetP4Raw.E()
            cef = jet.chargedEmEnergy() / jetP4Raw.E()
            mef = jet.muonEnergy() / jetP4Raw.E()
            nconstituents = jet.numberOfDaughters()
            cm = jet.chargedMultiplicity()
            nm = jet.neutralMultiplicity()

            goodJet = jetID(options.tag, jet.eta(), nhf, nef, chf, cef, mef, nconstituents, cm, nm, lepveto=False)
            goodJetLepVeto = jetID(options.tag, jet.eta(), nhf, nef, chf, cef, mef, nconstituents, cm, nm, lepveto=True)

            if jet.pt() > 30 and abs(jet.eta()) < 5.0:
                for v in jetIDvars:
                    jetIDhistos[v[0] + 'all'].Fill(globals()[v[0].replace('jet', '')])
                    if goodJet: jetIDhistos[v[0] + 'pass'].Fill(globals()[v[0].replace('jet', '')])

            if goodJetLepVeto:
                jetsIdxGood.append(ijet)

            if not goodJet:
                numBadJets += 1
                if abs(jet.eta()) < minetaabsbadjets: minetaabsbadjets = abs(jet.eta())
                if jet.pt() > 30 and abs(jet.eta()) < 5.0:
                    numBadJetsEventVeto += 1

            if not goodJetLepVeto:
                numBadJetsLepVeto += 1
                if abs(jet.eta()) < minetaabsbadjetsLepVeto: minetaabsbadjetsLepVeto = abs(jet.eta())
                if jet.pt() > 30 and abs(jet.eta()) < 5.0:
                    numBadJetsLepVetoEventVeto += 1

        if numBadJetsEventVeto > 0: continue

        # ########################################################################################### veto

        cutflow = 4
        hCutflow.Fill(cutflow)

        event_level_var_array['badJets_n'][0] = numBadJets
        event_level_var_array['badJets_minEta'][0] = minetaabsbadjets
        event_level_var_array['badJets_nForEventVeto'][0] = numBadJetsEventVeto

        event_level_var_array['badJets_lepVeto_n'][0] = numBadJetsLepVeto
        event_level_var_array['badJets_lepVeto_minEta'][0] = minetaabsbadjetsLepVeto
        event_level_var_array['badJets_lepVeto_nForEventVeto'][0] = numBadJetsLepVetoEventVeto

        jets = [j for ij, j in enumerate(jets) if ij in jetsIdxGood and abs(j.eta()) < 5.]

        if not len(jets) > 0: continue

        # ########################################################################################### veto

        '''
        ###############################################################################################
        # Type-I MET correction
        ###############################################################################################
        '''

        hMetptRaw.Fill(met.pt())

        for j in jetsP4Raw:
            met.setP4(met.p4() + ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(j.Px(), j.Py(), 0, j.Energy()))

        for j in jetsP4Corr:
            met.setP4(met.p4() - ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(j.Px(), j.Py(), 0, j.Energy()))

        '''
        ###############################################################################################
        # DY lepton cleaning
        ###############################################################################################
        '''

        hMetptBeforeLeptonCleaning.Fill(met.pt())

        metptBeforeCleaning = met.pt()
        metphiBeforeCleaning = met.phi()

        electronsCleaned = 0
        muonsCleaned = 0
        invm = -1
        zpt = -1
        l1pt = -1
        l2pt = -1
        l1eta = -1
        l2eta = -1
        l1phi = -1
        l2phi = -1
        l1absisodbeta, l1relisodbeta = -1, -1
        l2absisodbeta, l2relisodbeta = -1, -1

        if 'cleanleptons' in options.tag:

            l1Idx = -1
            l2Idx = -1

            if len(electrons) > 1:

                for ie1, e1 in enumerate(electrons):

                    e1Tlv = ROOT.TLorentzVector()
                    e1Tlv.SetPxPyPzE(e1.px(), e1.py(), e1.pz(), e1.energy())

                    for ie2, e2 in enumerate(electrons):

                        if ie2 == ie1: continue

                        if e1.charge() * e2.charge() > 0: continue

                        e2Tlv = ROOT.TLorentzVector()
                        e2Tlv.SetPxPyPzE(e2.px(), e2.py(), e2.pz(), e2.energy())

                        invm = (e1Tlv + e2Tlv).M()
                        zpt = (e1Tlv + e2Tlv).Pt()

                        if not invm > 75: continue
                        if not invm < 105: continue

                        l1Idx = ie1
                        l2Idx = ie2

                        electronsCleaned = 1

                        break

                    break

            if len(muons) > 1:

                for im1, m1 in enumerate(muons):

                    m1Tlv = ROOT.TLorentzVector()
                    m1Tlv.SetPxPyPzE(m1.px(), m1.py(), m1.pz(), m1.energy())

                    for im2, m2 in enumerate(muons):

                        if im2 == im1: continue

                        if m1.charge() * m2.charge() > 0: continue

                        m2Tlv = ROOT.TLorentzVector()
                        m2Tlv.SetPxPyPzE(m2.px(), m2.py(), m2.pz(), m2.energy())

                        invm = (m1Tlv + m2Tlv).M()
                        zpt = (m1Tlv + m2Tlv).Pt()

                        if not invm > 75: continue
                        if not invm < 105: continue

                        l1Idx = im1
                        l2Idx = im2

                        muonsCleaned = 1
                        electronsCleaned = 0

                        break

                    break

            if l1Idx == -1 or l2Idx == -1: continue

            if electronsCleaned: collection = electrons
            elif muonsCleaned: collection = muons
            else:
                print 'not tidy...'
                sys.exit(1)

            l1pt = collection[l1Idx].pt()
            l2pt = collection[l2Idx].pt()
            l1eta = collection[l1Idx].eta()
            l2eta = collection[l2Idx].eta()
            l1phi = collection[l1Idx].phi()
            l2phi = collection[l2Idx].phi()

            if electronsCleaned:
                l1absisodbeta = calcIso_dBeta(collection[l1Idx].pfIsolationVariables())
                l1relisodbeta = calcIso_dBeta(collection[l1Idx].pfIsolationVariables()) / collection[l1Idx].pt()
                l2absisodbeta = calcIso_dBeta(collection[l2Idx].pfIsolationVariables())
                l2relisodbeta = calcIso_dBeta(collection[l2Idx].pfIsolationVariables()) / collection[l2Idx].pt()
            else:
                l1absisodbeta = calcIso_dBeta(collection[l1Idx].pfIsolationR03())
                l1relisodbeta = calcIso_dBeta(collection[l1Idx].pfIsolationR03()) / collection[l1Idx].pt()
                l2absisodbeta = calcIso_dBeta(collection[l2Idx].pfIsolationR03())
                l2relisodbeta = calcIso_dBeta(collection[l2Idx].pfIsolationR03()) / collection[l2Idx].pt()

            collection, tracks, pfcands, jets, met = cleanZllEvent(l1Idx, l2Idx, collection, tracks, pfcands, jets, met, hZllLeptonPt, hZllDrTrack, hZllDrPfc, hZllDrJet)

            if tracks is None: continue

            if electronsCleaned: electrons = collection
            else: muons = collection

        # ########################################################################################### veto

        cutflow = 5
        hCutflow.Fill(cutflow)

        event_level_var_array['cleaning_electronsCleaned'][0] = electronsCleaned
        event_level_var_array['cleaning_muonsCleaned'][0] = muonsCleaned
        event_level_var_array['cleaning_invm'][0] = invm
        event_level_var_array['cleaning_zPt'][0] = zpt
        event_level_var_array['cleaning_l1Pt'][0] = l1pt
        event_level_var_array['cleaning_l2Pt'][0] = l2pt
        event_level_var_array['cleaning_l1Eta'][0] = l1eta
        event_level_var_array['cleaning_l2Eta'][0] = l2eta
        event_level_var_array['cleaning_l1Phi'][0] = l1phi
        event_level_var_array['cleaning_l2Phi'][0] = l2phi
        event_level_var_array['cleaning_l1dBetaAbsIso'][0] = l1absisodbeta
        event_level_var_array['cleaning_l1dBetaRelIso'][0] = l1relisodbeta
        event_level_var_array['cleaning_l2dBetaAbsIso'][0] = l2absisodbeta
        event_level_var_array['cleaning_l2dBetaRelIso'][0] = l2relisodbeta
        event_level_var_array['cleaning_metPtBeforeCleaning'][0] = metptBeforeCleaning
        event_level_var_array['cleaning_metPhiBeforeCleaning'][0] = metphiBeforeCleaning

        '''
        ###############################################################################################
        # GEN MET and HT(miss) and FastSim MET correction
        ###############################################################################################
        '''

        pTneutrinosum = -1
        genmetpt = -1
        genmetphi = -1
        genht = -1
        genht5 = -1
        genhtmiss = -1
        if 'data' not in options.tag:

            genht = 0
            genht5 = 0

            allneutrinos = [gp for gp in genparticles if gp.status() == 1
                            and (abs(gp.pdgId()) == 12 or abs(gp.pdgId()) == 14 or abs(gp.pdgId()) == 16)]

            sumneutrinos = ROOT.TLorentzVector()
            for n in allneutrinos:
                nTlv = ROOT.TLorentzVector(n.px(), n.py(), n.pz(), n.energy())
                sumneutrinos += nTlv

            pTneutrinosum = sumneutrinos.Pt()

            genmetpt = genmet.pt()
            genmetphi = genmet.phi()
            genhtmissTlv = ROOT.TLorentzVector()
            for genjet in genjets:
                if abs(genjet.eta()) < 2.4 and genjet.pt() > 30:
                    genht += genjet.pt()
                if abs(genjet.eta()) < 5 and genjet.pt() > 30:
                    genht5 += genjet.pt()
                    genjetTlv = ROOT.TLorentzVector(genjet.px(), genjet.py(), genjet.pz(), genjet.energy())
                    genhtmissTlv -= genjetTlv
            genhtmiss = genhtmissTlv.Pt()

        nofastsimcorrmetpt = met.pt()
        nofastsimcorrmetphi = met.phi()

        if 'fastsim' in options.tag:

            met.setP4(ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.5 * (genmet.px() + met.px())
                                                                               , 0.5 * (genmet.py() + met.py())
                                                                               , 0
                                                                               , 0.5 * (genmet.energy() + met.energy())))

        event_level_var_array['gen_neutrinoSumPt'][0] = pTneutrinosum
        event_level_var_array['gen_met_pt'][0] = genmetpt
        event_level_var_array['gen_met_phi'][0] = genmetphi
        event_level_var_array['gen_ht'][0] = genht
        event_level_var_array['gen_ht5'][0] = genht5
        event_level_var_array['gen_htMiss'][0] = genhtmiss
        event_level_var_array['met_ptNoFastSimCorr'][0] = nofastsimcorrmetpt
        event_level_var_array['met_phiNoFastSimCorr'][0] = nofastsimcorrmetphi
        event_level_var_array['met_phi'][0] = met.phi()
        event_level_var_array['met_pt'][0] = met.pt()

        hMetpt.Fill(met.pt())


        '''
        ###############################################################################################
        # get event-level info
        ###############################################################################################
        '''

        event_level_var_array['n_pv'][0] = n_pv
        event_level_var_array['rho'][0] = rho

        hasISRJet = False
        passedMet = False

        deltaPhi1 = 999  
        deltaPhi2 = 999
        deltaPhi3 = 999
        deltaPhi4 = 999

        minPt1 = -1
        minPt2 = -1
        minPt3 = -1
        minPt4 = -1

        eta1 = 999
        eta2 = 999
        eta3 = 999
        eta4 = 999

        lJet1 = None
        lJet2 = None
        lJet3 = None
        lJet4 = None

        numjets = len(jets)
        numjets15 = 0
        numjets30 = 0
        numjets50 = 0
        numjets100 = 0
        numjets200 = 0
        ht = 0
        ht5 = 0
        htmissTlv = ROOT.TLorentzVector()
        idxhighestptjet = 0
        for ijet, jet in enumerate(jets):

            pt = jet.pt()
            eta = fabs(jet.eta())
            if ((eta < 2.5) and (pt >100) ): hasISRJet = True
            
            if (pt > 30): 
                if (pt > minPt1):
                    minPt1 = pt
                    lJet1 = jet
                    eta1 = eta
                elif (pt > minPt2):
                    minPt2 = pt
                    lJet2 = jet
                elif (pt > minPt3):
                    minPt3 = pt
                    lJet3 = jet
                elif (pt > minPt4):
                    minPt4 = pt
                    lJet4 = jet
                    
            if abs(jet.eta()) < 2.4 and jet.pt() > 30: ht += jet.pt()
            if abs(jet.eta()) < 5. and jet.pt() > 30:
                ht5 += jet.pt()
                jetTlv = ROOT.TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy())
                htmissTlv -= jetTlv

            jetpt = jet.pt()
            if jetpt > 15: numjets15 += 1
            if jetpt > 30: numjets30 += 1
            if jetpt > 50: numjets50 += 1
            if jetpt > 100: numjets100 += 1
            if jetpt > 200: numjets200 += 1

            if jetpt > jets[idxhighestptjet].pt(): idxhighestptjet = ijet

        event_level_var_array['hasISRJet'][0] = hasISRJet
        if len(jets) > 0:
            event_level_var_array['leadingJet_pt'][0] = jets[idxhighestptjet].pt()
            event_level_var_array['leadingJet_eta'][0] = jets[idxhighestptjet].eta()
            event_level_var_array['leadingJet_phi'][0] = jets[idxhighestptjet].phi()
        else:
            event_level_var_array['leadingJet_pt'][0] = -1
            event_level_var_array['leadingJet_eta'][0] = -1
            event_level_var_array['leadingJet_phi'][0] = -1
                        
            
        if (not(minPt1 == -1)): deltaPhi1 = deltaPhi(met.phi(),lJet1.phi())
        if (not(minPt2 == -1)): deltaPhi2 = deltaPhi(met.phi(),lJet2.phi())
        if (not(minPt3 == -1)): deltaPhi3 = deltaPhi(met.phi(),lJet3.phi())
        if (not(minPt4 == -1)): deltaPhi4 = deltaPhi(met.phi(),lJet4.phi())

        event_level_var_array['JetMetdeltaPhi1'][0] = deltaPhi1
        event_level_var_array['JetMetdeltaPhi2'][0] = deltaPhi2
        event_level_var_array['JetMetdeltaPhi3'][0] = deltaPhi3
        event_level_var_array['JetMetdeltaPhi4'][0] = deltaPhi4

        event_level_var_array['JetPt1'][0] = minPt1
        event_level_var_array['JetPt2'][0] = minPt2
        event_level_var_array['JetPt3'][0] = minPt3
        event_level_var_array['JetPt4'][0] = minPt4

        event_level_var_array['JetEta1'][0] = eta1
        event_level_var_array['JetEta2'][0] = eta1
        event_level_var_array['JetEta3'][0] = eta1
        event_level_var_array['JetEta4'][0] = eta1

        event_level_var_array['n_jet'][0] = numjets
        event_level_var_array['n_jet_15'][0] = numjets15
        event_level_var_array['n_jet_30'][0] = numjets30
        event_level_var_array['n_jet_50'][0] = numjets50
        event_level_var_array['n_jet_100'][0] = numjets100
        event_level_var_array['n_jet_200'][0] = numjets200
        event_level_var_array['ht'][0] = ht
        event_level_var_array['ht5'][0] = ht5
        event_level_var_array['htMiss'][0] = htmissTlv.Pt()

        hNumjets.Fill(numjets)
        hNumjets30.Fill(numjets30)
        hNumjets50.Fill(numjets50)
        hNumjets100.Fill(numjets100)
        hNumjets200.Fill(numjets200)


        dphimetjets = []
        for jet in jets:
            if jet.pt() > 30 and abs(jet.eta()) < 2.4:

                dphimetjet = abs(deltaPhi(met.phi(), jet.phi()))
                dphimetjets.append(dphimetjet)

            if len(dphimetjets) > 3: break

        if len(dphimetjets) > 0:
            event_level_var_array['dphiminMetJets'][0] = min(dphimetjets)
            hMindphimetjets.Fill(min(dphimetjets))
        else:
            event_level_var_array['dphiminMetJets'][0] = -1
            hMindphimetjets.Fill(-1)


        hNPVsPerEvent.Fill(n_pv)

        hPV0x.Fill(primaryvertices[0].x())
        hPV0y.Fill(primaryvertices[0].y())
        hPV0z.Fill(primaryvertices[0].z())

        for i in range(1, n_pv):

            hPVsx.Fill(primaryvertices[i].x())
            hPVsy.Fill(primaryvertices[i].y())
            hPVsz.Fill(primaryvertices[i].z())


        '''
        ###############################################################################################
        # event selection
        ###############################################################################################
        '''

        if not met.pt() > metthreshold: continue

        # ########################################################################################### veto

        cutflow = 6
        hCutflow.Fill(cutflow)

        if 'veto_jet100' in options.tag:
            if not numjets100 > 0: continue

        # ########################################################################################### veto

        cutflow = 7
        hCutflow.Fill(cutflow)

        if 'veto_dphimetjets' in options.tag:
            if len(dphimetjets) > 0:
                if not min(dphimetjets) > 0.5: continue
            else:
                continue

        # ########################################################################################### veto

        cutflow = 8
        hCutflow.Fill(cutflow)

        ## ask Moritz
        #chpfcandsforiso0 = np.array([(p.pt(), p.eta(), p.phi()) for p in pfcands if passesPreselection_iso_chpf(p, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=0.)])
        chpfcandsforiso0 = np.array([(p.pt(), p.eta(), p.phi()) for p in pfcands if passesPreselection_iso_pf(p, pt_threshold=0.)])
        pfcandsforiso0 = np.array([(p.pt(), p.eta(), p.phi()) for p in pfcands if passesPreselection_iso_pf(p, pt_threshold=0.)])
        jetsforiso15 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets if passesPreselection_iso_jet(j, pt_threshold=15.)])


        event_level_var_array['n_electron'][0] = len(electrons)

        numelectronsiso = 0
        for ie, e in enumerate(electrons):

            electron_var_array['electron_charge'][ie] = e.charge()
            electron_var_array['electron_px'][ie] = e.px()
            electron_var_array['electron_py'][ie] = e.py()
            electron_var_array['electron_pz'][ie] = e.pz()
            electron_var_array['electron_pt'][ie] = e.pt()
            electron_var_array['electron_energy'][ie] = e.energy()
            electron_var_array['electron_eta'][ie] = e.eta()
            electron_var_array['electron_phi'][ie] = e.phi()
            electron_var_array['electron_dz'][ie] = abs(e.gsfTrack().dz(pv_pos))
            electron_var_array['electron_dxy'][ie] = abs(e.gsfTrack().dxy(pv_pos))
            electron_var_array['electron_pfAbsIso'][ie], _, _, _ = calcIso_pf_or_track_new(e, pfcandsforiso0, subtractObject=True)
            electron_var_array['electron_pfAbsIsoMini'][ie], _, _, _ = calcIso_pf_or_track_new(e, pfcandsforiso0, isMini=True, subtractObject=True)
            electron_var_array['electron_chPfAbsIso'][ie], _, _, _ = calcIso_pf_or_track_new(e, chpfcandsforiso0, subtractObject=(abs(e.gsfTrack().dz(pv_pos)) < 0.1 and abs(e.gsfTrack().dxy(pv_pos)) < 0.1))
            electron_var_array['electron_chPfAbsIsoMini'][ie], _, _, _ = calcIso_pf_or_track_new(e, chpfcandsforiso0, isMini=True, subtractObject=(abs(e.gsfTrack().dz(pv_pos)) < 0.1 and abs(e.gsfTrack().dxy(pv_pos)) < 0.1))
            electron_var_array['electron_jetIso'][ie], electron_var_array['electron_jetIsoMulti'][ie], electron_var_array['electron_drminJet'][ie], _, electron_var_array['electron_minvJet'][ie] = calcIso_jet_new(e, jetsforiso15, isTrack=False, btagvalues=btagvalues)
            electron_var_array['electron_chHadIso'][ie] = e.pfIsolationVariables().sumChargedHadronPt
            electron_var_array['electron_chAllIso'][ie] = e.pfIsolationVariables().sumChargedParticlePt
            electron_var_array['electron_neHadIso'][ie] = e.pfIsolationVariables().sumNeutralHadronEt
            electron_var_array['electron_photIso'][ie] = e.pfIsolationVariables().sumPhotonEt
            electron_var_array['electron_puChHadIso'][ie] = e.pfIsolationVariables().sumPUPt
            electron_var_array['electron_dBetaAbsIso'][ie] = calcIso_dBeta(e.pfIsolationVariables())
            electron_var_array['electron_dBetaRelIso'][ie] = calcIso_dBeta(e.pfIsolationVariables()) / e.pt()

            if calcIso_dBeta(e.pfIsolationVariables()) / e.pt() < 0.2: numelectronsiso += 1

        event_level_var_array['n_electron_iso'][0] = numelectronsiso


        event_level_var_array['n_muon'][0] = len(muons)

        nummuonsiso = 0
        for im, m in enumerate(muons):

            muon_var_array['muon_charge'][im] = m.charge()
            muon_var_array['muon_px'][im] = m.px()
            muon_var_array['muon_py'][im] = m.py()
            muon_var_array['muon_pz'][im] = m.pz()
            muon_var_array['muon_pt'][im] = m.pt()
            muon_var_array['muon_energy'][im] = m.energy()
            muon_var_array['muon_eta'][im] = m.eta()
            muon_var_array['muon_phi'][im] = m.phi()
            muon_var_array['muon_dz'][im] = abs(m.muonBestTrack().dz(pv_pos))
            muon_var_array['muon_dxy'][im] = abs(m.muonBestTrack().dxy(pv_pos))
            muon_var_array['muon_pfAbsIso'][im], _, _, _ = calcIso_pf_or_track_new(m, pfcandsforiso0, subtractObject=True)
            muon_var_array['muon_pfAbsIsoMini'][im], _, _, _ = calcIso_pf_or_track_new(m, pfcandsforiso0, isMini=True, subtractObject=True)
            muon_var_array['muon_chPfAbsIso'][im], _, _, _ = calcIso_pf_or_track_new(m, chpfcandsforiso0, subtractObject=(abs(m.muonBestTrack().dz(pv_pos)) < 0.1 and abs(m.muonBestTrack().dxy(pv_pos)) < 0.1))
            muon_var_array['muon_chPfAbsIsoMini'][im], _, _, _ = calcIso_pf_or_track_new(m, chpfcandsforiso0, isMini=True, subtractObject=(abs(m.muonBestTrack().dz(pv_pos)) < 0.1 and abs(m.muonBestTrack().dxy(pv_pos)) < 0.1))
            muon_var_array['muon_jetIso'][im], muon_var_array['muon_jetIsoMulti'][im], muon_var_array['muon_drminJet'][im], _, muon_var_array['muon_minvJet'][im] = calcIso_jet_new(m, jetsforiso15, isTrack=False, btagvalues=btagvalues)
            muon_var_array['muon_chHadIso'][im] = m.pfIsolationR03().sumChargedHadronPt
            muon_var_array['muon_chAllIso'][im] = m.pfIsolationR03().sumChargedParticlePt
            muon_var_array['muon_neHadIso'][im] = m.pfIsolationR03().sumNeutralHadronEt
            muon_var_array['muon_photIso'][im] = m.pfIsolationR03().sumPhotonEt
            muon_var_array['muon_puChHadIso'][im] = m.pfIsolationR03().sumPUPt
            muon_var_array['muon_dBetaAbsIso'][im] = calcIso_dBeta(m.pfIsolationR03())
            muon_var_array['muon_dBetaRelIso'][im] = calcIso_dBeta(m.pfIsolationR03()) / m.pt()

            if calcIso_dBeta(m.pfIsolationR03()) / m.pt() < 0.2: nummuonsiso += 1

        event_level_var_array['n_muon_iso'][0] = nummuonsiso

        event_level_var_array['n_lepton'][0] = len(electrons) + len(muons)
        event_level_var_array['n_lepton_iso'][0] = numelectronsiso + nummuonsiso

        hNumleptons.Fill(numelectronsiso+nummuonsiso)

        if 'veto_isolepton' in options.tag:
            if not (numelectronsiso + nummuonsiso) == 0: continue

        # ########################################################################################### veto

        cutflow = 9
        hCutflow.Fill(cutflow)


        '''
        ###############################################################################################
        # GEN info background
        ###############################################################################################
        '''

        numZgamma = 0
        pdgidZgamma = -1
        ptZgamma = -1
        etaZgamma = -1
        phiZgamma = -1
        numZgammaDaughters = 0
        decayZtau = -1
        ptsumZgammaNeutrinos = -1
        if 'ZJetsToNuNu' in options.tag or 'DYJetsToLL' in options.tag:

            Zgammas = [gp for gp in genparticles if gp.isLastCopy() and gp.statusFlags().fromHardProcess()
                       and (abs(gp.pdgId()) == 22 or abs(gp.pdgId()) == 23)]

            numZgamma = len(Zgammas)

            if numZgamma == 1:  # else there's something funny going on...

                Zgamma = Zgammas[0]

                pdgidZgamma = Zgamma.pdgId()
                ptZgamma = Zgamma.pt()
                etaZgamma = Zgamma.eta()
                phiZgamma = Zgamma.phi()

                numZgammaDaughters = Zgamma.numberOfDaughters()

                ZgammaNeutrinos = ROOT.TLorentzVector()

                firstZtau = True
                for i in range(numZgammaDaughters):

                    daughter = Zgamma.daughter(i)
                    if not abs(daughter.pdgId()) == 15:
                        daughter = getLastCopyStatusOne(daughter)
                    else:
                        daughter = getLastCopy(daughter)

                    if daughter is None: continue

                    zdaughter_var_array['zDaughter_pdgId'][i] = daughter.pdgId()
                    zdaughter_var_array['zDaughter_pt'][i] = daughter.pt()
                    zdaughter_var_array['zDaughter_eta'][i] = daughter.eta()
                    zdaughter_var_array['zDaughter_phi'][i] = daughter.phi()

                    if abs(daughter.pdgId()) == 15 and firstZtau:
                        decayZtau = getTauDecayMode(daughter, decayZtau)
                        firstZtau = False

                    if abs(daughter.pdgId()) == 12 or abs(daughter.pdgId()) == 14 or abs(daughter.pdgId()) == 16:
                        nTlv = ROOT.TLorentzVector(daughter.px(), daughter.py(), daughter.pz(), daughter.energy())
                        ZgammaNeutrinos += nTlv

                ptsumZgammaNeutrinos = ZgammaNeutrinos.Pt()

        event_level_var_array['n_zGamma'][0] = numZgamma
        event_level_var_array['zGamma_pdgId'][0] = pdgidZgamma
        event_level_var_array['zGamma_pt'][0] = ptZgamma
        event_level_var_array['zGamma_eta'][0] = etaZgamma
        event_level_var_array['zGamma_phi'][0] = phiZgamma
        event_level_var_array['zGamma_neutrinoSumPt'][0] = ptsumZgammaNeutrinos
        event_level_var_array['zGamma_tauDecayMode'][0] = decayZtau
        event_level_var_array['n_zDaughter'][0] = numZgammaDaughters


        numW = 0
        pdgidW = -1
        ptW = -1
        etaW = -1
        phiW = -1
        numWDaughters = 0
        ptWneutrino = -1
        decayWtau = -1
        decaylengthXYZWtau = -1
        decaylengthXYWtau = -1
        decaylengthZWtau = -1
        ptWtau = -1
        etaWtau = -1
        phiWtau = -1
        ptVisWtau = -1
        etaVisWtau = -1
        phiVisWtau = -1
        thetau = None
        if 'WJetsToLNu' in options.tag:

            Ws = [gp for gp in genparticles if gp.isLastCopy() and gp.statusFlags().fromHardProcess()
                  and abs(gp.pdgId()) == 24]

            numW = len(Ws)

            if numW == 1:  # else there's something funny going on...

                W = Ws[0]

                ptW = W.pt()
                etaW = W.eta()
                phiW = W.phi()
                pdgidW = W.pdgId()

                numWDaughters = W.numberOfDaughters()

                ptWneutrino = -1
                decayWtau = -1
                decaylengthXYZWtau = -1
                decaylengthXYWtau = -1
                decaylengthZWtau = -1
                ptWtau = -1
                etaWtau = -1
                phiWtau = -1
                ptVisWtau = -1
                etaVisWtau = -1
                phiVisWtau = -1

                for i in range(numWDaughters):

                    daughter = W.daughter(i)
                    if not abs(daughter.pdgId()) == 15:
                        daughter = getLastCopyStatusOne(daughter)
                    else:
                        daughter = getLastCopy(daughter)

                    if daughter is None: continue

                    viswdaughter = ROOT.TLorentzVector()
                    for igranddaughter in range(daughter.numberOfDaughters()):
                        if abs(daughter.daughter(igranddaughter).pdgId()) in [11, 12, 13, 14, 15, 16]: continue
                        viswdaughterTlv = ROOT.TLorentzVector(daughter.daughter(igranddaughter).px(), daughter.daughter(igranddaughter).py(), daughter.daughter(igranddaughter).pz(), daughter.daughter(igranddaughter).energy())
                        viswdaughter += viswdaughterTlv

                    wdaughter_var_array['wDaughter_pdgId'][i] = daughter.pdgId()

                    wdaughter_var_array['wDaughter_pt'][i] = daughter.pt()
                    wdaughter_var_array['wDaughter_eta'][i] = daughter.eta()
                    wdaughter_var_array['wDaughter_phi'][i] = daughter.phi()

                    wdaughter_var_array['wDaughter_ptVis'][i] = viswdaughter.Pt()
                    wdaughter_var_array['wDaughter_etaVis'][i] = viswdaughter.Eta()
                    wdaughter_var_array['wDaughter_phiVis'][i] = viswdaughter.Phi()

                    if abs(daughter.pdgId()) == 15:
                        thetau = daughter
                        decayWtau = getTauDecayMode(daughter, decayWtau)
                        decaylengthXYZWtau = ROOT.TMath.Sqrt(pow(daughter.vx() - daughter.daughter(0).vx(), 2)
                                                             + pow(daughter.vy() - daughter.daughter(0).vy(), 2)
                                                             + pow(daughter.vz() - daughter.daughter(0).vz(), 2))
                        decaylengthXYWtau = ROOT.TMath.Sqrt(pow(daughter.vx() - daughter.daughter(0).vx(), 2)
                                                            + pow(daughter.vy() - daughter.daughter(0).vy(), 2))
                        decaylengthZWtau = abs(daughter.vz() - daughter.daughter(0).vz())
                        ptWtau = daughter.pt()
                        etaWtau = daughter.eta()
                        phiWtau = daughter.phi()
                        ptVisWtau = viswdaughter.Pt()
                        etaVisWtau = viswdaughter.Eta()
                        phiVisWtau = viswdaughter.Phi()

                    if abs(daughter.pdgId()) == 12 or abs(daughter.pdgId()) == 14 or abs(daughter.pdgId()) == 16:
                        ptWneutrino = daughter.pt()

        event_level_var_array['n_wBoson'][0] = numW
        event_level_var_array['wBoson_pdgId'][0] = pdgidW
        event_level_var_array['wBoson_pt'][0] = ptW
        event_level_var_array['wBoson_eta'][0] = etaW
        event_level_var_array['wBoson_phi'][0] = phiW
        event_level_var_array['wBoson_neutrinoPt'][0] = ptWneutrino
        event_level_var_array['wBoson_tauDecayMode'][0] = decayWtau
        event_level_var_array['wBoson_tauDecaylengthXYZ'][0] = decaylengthXYZWtau
        event_level_var_array['wBoson_tauDecaylengthXY'][0] = decaylengthXYWtau
        event_level_var_array['wBoson_tauDecaylengthZ'][0] = decaylengthZWtau
        event_level_var_array['wBoson_tauPt'][0] = ptWtau
        event_level_var_array['wBoson_tauEta'][0] = etaWtau
        event_level_var_array['wBoson_tauPhi'][0] = phiWtau
        event_level_var_array['wBoson_tauPtVis'][0] = ptVisWtau
        event_level_var_array['wBoson_tauEtaVis'][0] = etaVisWtau
        event_level_var_array['wBoson_tauPhiVis'][0] = phiVisWtau
        event_level_var_array['n_wDaughter'][0] = numWDaughters


        '''
        ###############################################################################################
        # GEN info for signal and track matching
        ###############################################################################################
        '''

        matchedTrackIdxCharginoPion1 = -1
        matchedTrackIdxCharginoPion2 = -1

        stop_mass, antistop_mass = -1, -1
        stop_pt, antistop_pt = -1, -1
        stop_eta, antistop_eta = -1, -1
        stop_phi, antistop_phi = -1, -1
        stop_decay, antistop_decay = -1, -1

        susytracks = {}

        if 'signal' in options.tag:

            chipmmFILE = float(re.search(r'mChipm(.*?)GeV', f).group(1))
            deltamFILE = float(re.search(r'_dm(.*?)GeV', f).group(1).replace('p', '.'))
            try:
                mstopFILE = float(re.search(r'stopstop_(.*?)GeV', f).group(1))
            except AttributeError:
                pass

            C1s = [gp for gp in genparticles if gp.isLastCopy() and abs(gp.pdgId()) == 1000024]
            N2s = [gp for gp in genparticles if gp.isLastCopy() and abs(gp.pdgId()) == 1000023]
            N1s = [gp for gp in genparticles if gp.isLastCopy() and abs(gp.pdgId()) == 1000022]

            numC1 = len(C1s)
            numN2 = len(N2s)
            numN1 = len(N1s)

            if len(N2s) == 1: 
                theChi02 = N2s[0]
                if len(N1s) == 1: theChi01 = N1s[0]
                
            c1daughters = []
            n2daughters = []
            
            leptons = [None, None]
            signalIdx = -1

            for igp, gp in enumerate(C1s):

                chipmmGEN = round(gp.mass(), 2)

                deltamGEN = round((gp.mass() - gp.daughter(0).mass()), 2)

                chiC1_var_array['chiC1_pt'][igp] = gp.pt()
                chiC1_var_array['chiC1_eta'][igp] = gp.eta()
                chiC1_var_array['chiC1_phi'][igp] = gp.phi()

                chiC1_decaylengthXYZ = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                       + pow(gp.vy() - gp.daughter(0).vy(), 2)
                                                       + pow(gp.vz() - gp.daughter(0).vz(), 2))

                chiC1_decaylengthXY = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                      + pow(gp.vy() - gp.daughter(0).vy(), 2))

                chiC1_chidecaylengthZ = abs(gp.vz() - gp.daughter(0).vz())

                chiC1_var_array['chiC1_decaylengthXYZ'][igp] = chiC1_decaylengthXYZ
                chiC1_var_array['chiC1_decaylengthXY'][igp] = chiC1_decaylengthXY
                chiC1_var_array['chiC1_chidecaylengthZ'][igp] = chiC1_chidecaylengthZ

                if chiC1_decaylengthXYZ > 0:
                    chiC1_var_array['chiC1_log10(decaylengthXYZ)'][igp] = ROOT.TMath.Log10(chiC1_decaylengthXYZ)
                    chiC1_var_array['chiC1_log10(decaylengthXY)'][igp] = ROOT.TMath.Log10(chiC1_decaylengthXY)
                    chiC1_var_array['chiC1_log10(chidecaylengthZ)'][igp] = ROOT.TMath.Log10(chiC1_chidecaylengthZ)
                else:
                    chiC1_var_array['chiC1_log10(decaylengthXYZ)'][igp] = -10.
                    chiC1_var_array['chiC1_log10(decaylengthXY)'][igp] = -10.
                    chiC1_var_array['chiC1_log10(chidecaylengthZ)'][igp] = -10.

                thisc1daughters = findDaughters(gp)
                c1daughters += thisc1daughters

                hasPion = 0
                hasMatchedTrackPion = 0

                for c1d in thisc1daughters:

                    if abs(c1d.pdgId()) == 211:

                        pion = c1d

                        hasPion = 1

                        chiC1_var_array['chiC1_pionPt'][igp] = pion.pt()
                        chiC1_var_array['chiC1_pionEta'][igp] = pion.eta()
                        chiC1_var_array['chiC1_pionPhi'][igp] = pion.phi()
                        chiC1_var_array['chiC1_pionCharge'][igp] = pion.charge()

                        idxold, drminold = findMatch_track_old(pion, tracks)
                        _, drminoldrandom = findMatch_track_old_random(pion, tracks)

                        idx, dxyzmin, tminmatching, drmin = findMatch_track_new(pion, tracks)
                        _, dxyzminrandom, _, drminrandom = findMatch_track_new_random(pion, tracks)

                        chiC1_var_array['chiC1_pionMatching_tmin'][igp] = tminmatching
                        chiC1_var_array['chiC1_pionMatching_dxyzmin'][igp] = dxyzmin
                        chiC1_var_array['chiC1_pionMatching_drmin'][igp] = drmin
                        chiC1_var_array['chiC1_pionMatching_dxyzminrandom'][igp] = dxyzminrandom
                        chiC1_var_array['chiC1_pionMatching_drminrandom'][igp] = drminrandom
                        chiC1_var_array['chiC1_pionMatching_drminold'][igp] = drminold
                        chiC1_var_array['chiC1_pionMatching_drminoldrandom'][igp] = drminoldrandom

                        if not idx == -1:
                            if drmin < matchingDrThreshold and dxyzmin < matchingDxyzThreshold:

                                hasMatchedTrackPion = 1

                                if matchedTrackIdxCharginoPion1 == -1: matchedTrackIdxCharginoPion1 = idx
                                else: matchedTrackIdxCharginoPion2 = idx

                        break

                chiC1_var_array['chiC1_hasPion'][igp] = hasPion
                chiC1_var_array['chiC1_hasMatchedTrackPion'][igp] = hasMatchedTrackPion

            signalTrkIdx = [-1, -1]
            signalTrk = [None, None]
            
            for igp, gp in enumerate(N2s):

                chiN2mGEN = round(gp.mass(), 2)

                deltamN2GEN = round((gp.mass() - gp.daughter(0).mass()), 2)

                chiN2_var_array['chiN2_pt'][igp] = gp.pt()
                chiN2_var_array['chiN2_eta'][igp] = gp.eta()
                chiN2_var_array['chiN2_phi'][igp] = gp.phi()
                
                chiN2_var_array['chi02_pz'][igp] = gp.pz()


                chiN2_decaylengthXYZ = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                       + pow(gp.vy() - gp.daughter(0).vy(), 2)
                                                       + pow(gp.vz() - gp.daughter(0).vz(), 2))

                chiN2_decaylengthXY = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                      + pow(gp.vy() - gp.daughter(0).vy(), 2))

                chiN2_chidecaylengthZ = abs(gp.vz() - gp.daughter(0).vz())

                chiN2_var_array['chiN2_decaylengthXYZ'][igp] = chiN2_decaylengthXYZ
                chiN2_var_array['chiN2_decaylengthXY'][igp] = chiN2_decaylengthXY
                chiN2_var_array['chiN2_chidecaylengthZ'][igp] = chiN2_chidecaylengthZ

                if chiN2_decaylengthXYZ > 0:
                    chiN2_var_array['chiN2_log10(decaylengthXYZ)'][igp] = ROOT.TMath.Log10(chiN2_decaylengthXYZ)
                    chiN2_var_array['chiN2_log10(decaylengthXY)'][igp] = ROOT.TMath.Log10(chiN2_decaylengthXY)
                    chiN2_var_array['chiN2_log10(chidecaylengthZ)'][igp] = ROOT.TMath.Log10(chiN2_chidecaylengthZ)
                else:
                    chiN2_var_array['chiN2_log10(decaylengthXYZ)'][igp] = -10.
                    chiN2_var_array['chiN2_log10(decaylengthXY)'][igp] = -10.
                    chiN2_var_array['chiN2_log10(chidecaylengthZ)'][igp] = -10.

                n2daughters += findDaughters(gp)
                
                leptons[0],leptons[1] = findLeptons(gp)
                if not leptons[0] == None and not leptons[1]==None: 


                    chiN2_var_array['leptonID'][igp] = leptons[0].pdgId()
                    TLV_l1 = TLorentzVector()
                    TLV_l1.SetPxPyPzE(leptons[0].px(),leptons[0].py(),leptons[0].pz(),leptons[0].energy())
                    TLV_l2 = TLorentzVector()
                    TLV_l2.SetPxPyPzE(leptons[1].px(),leptons[1].py(),leptons[1].pz(),leptons[1].energy())
                
                    v_chi02 = TVector3(gp.vx(), gp.vy(), gp.vz())
                    v_chi01 = TVector3(gp.daughter(0).vx(), gp.daughter(0).vy(), gp.daughter(0).vz())
                    v_normal = v_chi01-v_chi02
                    normalvector = v_normal.Unit()
                    
            
                    ptZstar = leptons[0].pt()+leptons[1].pt()
                    mZstar = sqrt((TLV_l1+TLV_l2)*(TLV_l1+TLV_l2))
                    v_pZstar = TVector3(leptons[0].px()+leptons[1].px(), leptons[0].py()+leptons[1].py(), leptons[0].pz()+leptons[1].pz())
                    mtransverse2 = mZstar*mZstar+ ((v_pZstar.Cross(normalvector))*(v_pZstar.Cross(normalvector)))
                    mtransverse2_paper = 2*(leptons[0].pt())*(leptons[1].pt())*(1-cos(TLV_l1.Angle(TLV_l2.Vect())))
                    
                    chiN2_var_array['leptonBoost'][igp] = TLV_l1.Gamma()
                    chiN2_var_array['leptonBoost'][igp] = TLV_l2.Gamma()
                    if TLV_l1.Pt() > TLV_l2.Pt():
                        chiN2_var_array['lepton_High_pt'][igp] = TLV_l1.Pt()
                        chiN2_var_array['lepton_High_eta'][igp] = TLV_l1.Eta()
                        chiN2_var_array['lepton_Low_pt'][igp] = TLV_l2.Pt()
                        chiN2_var_array['lepton_Low_eta'][igp] = TLV_l2.Eta()
                    else:
                        chiN2_var_array['lepton_High_pt'][igp] = TLV_l2.Pt()
                        chiN2_var_array['lepton_High_eta'][igp] = TLV_l2.Eta()
                        chiN2_var_array['lepton_Low_pt'][igp] = TLV_l1.Pt()
                        chiN2_var_array['lepton_Low_eta'][igp] = TLV_l1.Eta()
                    
                    
                    chiN2_var_array['ZstarBoost'][igp] = (TLV_l1+TLV_l2).Gamma()

                    chiN2_var_array['mZstar'][igp] = mZstar
                    chiN2_var_array['ptZstar'][igp] = ptZstar
                    chiN2_var_array['abspZstar'][igp] = v_pZstar.Mag()
                    chiN2_var_array['absnormalVector'][igp] = normalvector.Mag()
                    chiN2_var_array['mtransverse2'][igp] = mtransverse2
                    chiN2_var_array['mtransverse2_paper'][igp] = mtransverse2_paper
                    chiN2_var_array['beta'][igp] = v_pZstar.Angle(normalvector)	

                    TLV_theChi01 = TLorentzVector()
                    TLV_theChi01.SetPxPyPzE(gp.daughter(0).px(), gp.daughter(0).py(), gp.daughter(0).pz(), gp.daughter(0).energy())
                    chiN2_var_array['Chi01Boost'][igp] = TLV_theChi01.Gamma()
                    TLV_theChi02 = TLorentzVector()
                    TLV_theChi02.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy())			
                    chiN2_var_array['Chi02Boost'][igp] = TLV_theChi02.Gamma()
                                                
                    PVSV = TVector3(gp.daughter(0).vx()-pv_pos.x(), gp.daughter(0).vy()-pv_pos.y(), gp.daughter(0).vz()-pv_pos.z())		
                    chiN2_var_array['Chi0sToPV_eta'][igp] = PVSV.Eta()
                    chiN2_var_array['Chi0sToPV_phi'][igp] = PVSV.Phi()					
                    summedLeptons = TLV_l1+TLV_l2
                    chiN2_var_array['deltaEtaChi0sToLeptons'][igp] = summedLeptons.Eta()-PVSV.Eta()
                    chiN2_var_array['absdeltaEtaChi0sToLeptons'][igp] = abs(summedLeptons.Eta()-PVSV.Eta()) #TODO> remove abs version, do that bz hnd later
                    chiN2_var_array['deltaPhiChi0sToLeptons'][igp] = deltaPhi(summedLeptons.Phi(), PVSV.Phi())

                    #####################
                    ### Match gen leptons to any tracks
                    #####################

                    idxTrk = [-1 ,-1]
                    drminTrk = [999,999]
                    dxyzminTrk = [999,999]
                    tminTrk = [-1, -1]

                    idxTrkClassic = [-1 ,-1]
                    drminTrkClassic = [999,999]	

                        
                    for ilepton, lepton in enumerate(leptons):
                        
                        idxTrk[ilepton], dxyzminTrk[ilepton], tminTrk[ilepton], drminTrk[ilepton] = findMatch_track_new(lepton, tracks)
                        _, idxTrkClassic[ilepton], drminTrkClassic[ilepton], _ = findMinDr(lepton, tracks, 20.)
                                
                    if  ((drminTrkClassic[0]<0.04) and (drminTrkClassic[1]<0.04)):	
                        for i, idx in enumerate(idxTrkClassic):
                            if 'debug' in options.tag : print i, idx, tracks.size()
                            signalTrk[i] = tracks[idx]
                            signalTrkIdx[i] = idx


                    elif  ((dxyzminTrk[0]<0.04 and drminTrk[0]<0.04) and (dxyzminTrk[1]<0.04 and drminTrk[1]<0.04)):
                        for i, idx in enumerate(idxTrk):
                            signalTrk[i] = tracks[idx]
                            signalTrkIdx[i] = idx

                    else: 
                        print "SV has no machting tracks"  ### todo: proper error handlling
                        
                        

            for igp, gp in enumerate(N1s):
				chiN1_var_array['chiN1_pt'][igp] = gp.pt()
				chiN1_var_array['chiN1_eta'][igp] = gp.eta()
				chiN1_var_array['chiN1_phi'][igp] = gp.phi()
				
				chiN1_var_array['chi01_vx'][igp] = gp.vx()
				chiN1_var_array['chi01_vy'][igp] = gp.vy()
				chiN1_var_array['chi01_vz'][igp] = gp.vz()
				chiN1_var_array['chi01_dx'][igp] = abs(gp.vx()-pv_pos.x())
				chiN1_var_array['chi01_dy'][igp] = abs(gp.vy()-pv_pos.y())
				chiN1_var_array['chi01_dz'][igp] = abs(gp.vz()-pv_pos.z())

            chipmnumdaughters = len(c1daughters)
            chiN2numdaughters = len(n2daughters)
            numchidaughters = chipmnumdaughters + chiN2numdaughters
            for ichid, chid in enumerate(c1daughters + n2daughters):
                ### ask Moritz
                #chidaughter_var_array['chiDaughter_pdgIdMother'][ichid] = chid.mother(0).pdgId()
                #chidaughter_var_array['chiDaughter_pdgId'][ichid] = chid.pdgId()
                #chidaughter_var_array['chiDaughter_pt'][ichid] = chid.pt()
                #chidaughter_var_array['chiDaughter_eta'][ichid] = chid.eta()
                #chidaughter_var_array['chiDaughter_phi'][ichid] = chid.phi()

                #chidaughter_var_array['chiDaughter_hasMatchedTrack'][ichid] = 0
                
                if chid.charge() != 0:

                    idxchid, dxyzminchid, _, drminchid = findMatch_track_new(chid, tracks)

                    if not idxchid == -1:
                        if drminchid < matchingDrThreshold and dxyzminchid < matchingDxyzThreshold:
                            susytracks[idxchid] = (chid.mother(0).pdgId(), chid.pdgId())
                            #chidaughter_var_array['chiDaughter_hasMatchedTrack'][ichid] = 1


            stops = [gp for gp in genparticles if gp.isLastCopy() and gp.pdgId() == 1000006]
            antistops = [gp for gp in genparticles if gp.isLastCopy() and gp.pdgId() == -1000006]

            if len(stops) == 1 and len(antistops) == 1:

                thestop = stops[0]
                theantistop = antistops[0]

                stop_mass = thestop.mass()
                stop_pt = thestop.pt()
                stop_eta = thestop.eta()
                stop_phi = thestop.phi()

                thestop_top = None
                for idx in range(thestop.numberOfDaughters()):
                    if abs(thestop.daughter(idx).pdgId()) == 6:
                        thestop_top = getLastCopy(thestop.daughter(idx))
                        break

                if thestop_top is not None:
                    thestop_w = None
                    for idx in range(thestop_top.numberOfDaughters()):
                        if abs(thestop_top.daughter(idx).pdgId()) == 24:
                            thestop_w = getLastCopy(thestop_top.daughter(idx))
                            break

                    if thestop_w is not None:
                        stop_decay = len([thestop_w.daughter(idx) for idx in range(thestop_w.numberOfDaughters()) if abs(thestop_w.daughter(idx).pdgId()) in [11, 13]])

                antistop_mass = theantistop.mass()
                antistop_pt = theantistop.pt()
                antistop_eta = theantistop.eta()
                antistop_phi = theantistop.phi()

                theantistop_top = None
                for idx in range(theantistop.numberOfDaughters()):
                    if abs(theantistop.daughter(idx).pdgId()) == 6:
                        theantistop_top = getLastCopy(theantistop.daughter(idx))
                        break

                if theantistop_top is not None:
                    theantistop_w = None
                    for idx in range(theantistop_top.numberOfDaughters()):
                        if abs(theantistop_top.daughter(idx).pdgId()) == 24:
                            theantistop_w = getLastCopy(theantistop_top.daughter(idx))
                            break

                    if theantistop_w is not None:
                        antistop_decay = len([theantistop_w.daughter(idx) for idx in range(theantistop_w.numberOfDaughters()) if abs(theantistop_w.daughter(idx).pdgId()) in [11, 13]])


        event_level_var_array['mchipmFile'][0] = chipmmFILE
        event_level_var_array['deltamFile'][0] = deltamFILE
        event_level_var_array['mstopFile'][0] = mstopFILE

        event_level_var_array['chiC1_m'][0] = chipmmGEN
        event_level_var_array['chiC1_deltamN1'][0] = deltamGEN

        event_level_var_array['chiN2_m'][0] = chiN2mGEN
        event_level_var_array['chiN2_deltamN1'][0] = deltamN2GEN

        event_level_var_array['n_chiC1'][0] = numC1
        event_level_var_array['n_chiN2'][0] = numN2
        event_level_var_array['n_chiN1'][0] = numN1

        event_level_var_array['n_chiC1Daughter'][0] = chipmnumdaughters
        event_level_var_array['n_chiN2Daughter'][0] = chiN2numdaughters
        event_level_var_array['n_chiDaughter'][0] = numchidaughters

        event_level_var_array['stop_mass'][0] = stop_mass
        event_level_var_array['stop_pt'][0] = stop_pt
        event_level_var_array['stop_eta'][0] = stop_eta
        event_level_var_array['stop_phi'][0] = stop_phi
        event_level_var_array['stop_decay'][0] = stop_decay
        
        event_level_var_array['antistop_mass'][0] = antistop_mass
        event_level_var_array['antistop_pt'][0] = antistop_pt
        event_level_var_array['antistop_eta'][0] = antistop_eta
        event_level_var_array['antistop_phi'][0] = antistop_phi
        event_level_var_array['antistop_decay'][0] = antistop_decay



        '''
        ###############################################################################################
        # Weights
        ###############################################################################################
        '''

        crossSection = 1.
        numSimEvents = 1.

        if 'SignalV1' in options.tag or 'SignalV2' in options.tag or 'SignalFullV2' in options.tag:

            # from https://twiki.cern.ch/twiki/bin/view/LHCPhysics/SUSYCrossSections13TeVhino
            higgsinoxsecfile = ROOT.TFile('/nfs/dust/cms/user/wolfmor/NTupleStuff/CN_hino_13TeV.root')

            if 100 <= chipmmFILE < 150:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_0')
            elif 150 <= chipmmFILE < 200:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_1')
            elif 200 <= chipmmFILE < 300:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_2')
            elif 300 <= chipmmFILE < 400:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_3')
            elif 400 <= chipmmFILE < 600:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_4')
            elif 600 <= chipmmFILE < 800:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_5')
            elif 800 <= chipmmFILE < 1000:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_6')
            elif 1000 <= chipmmFILE < 1200:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_7')
            elif 1200 <= chipmmFILE < 1500:
                higgsinoxsec = higgsinoxsecfile.Get('fit_nom_8')
            else:
                higgsinoxsec = None

            if higgsinoxsec is not None:
                crossSection = higgsinoxsec.Eval(chipmmFILE)

            higgsinoxsecfile.Close()

            fSimEventNumbers_Signal = None
            hSimEventNumbers_Signal = None
            if 'SignalV1' in options.tag:
                fSimEventNumbers_Signal = ROOT.TFile('/nfs/dust/cms/user/wolfmor/NTupleStuff/simEventNumbers_AOD_v1.root')
                hSimEventNumbers_Signal = fSimEventNumbers_Signal.Get('simEventNumbers_AOD_v1')  # is TH2F with sim. event numbers for each model point
            elif 'SignalV2' in options.tag:
                fSimEventNumbers_Signal = ROOT.TFile('/nfs/dust/cms/user/wolfmor/NTupleStuff/simEventNumbers_AOD_v2.root')
                hSimEventNumbers_Signal = fSimEventNumbers_Signal.Get('simEventNumbers_AOD_v2')  # is TH2F with sim. event numbers for each model point
            elif 'SignalFullV2' in options.tag:
                fSimEventNumbers_Signal = ROOT.TFile('/nfs/dust/cms/user/wolfmor/NTupleStuff/simEventNumbers_FullSim_AOD_v2.root')
                hSimEventNumbers_Signal = fSimEventNumbers_Signal.Get('simEventNumbers_FullSim_AOD_v2')  # is TH2F with sim. event numbers for each model point
                
            if fSimEventNumbers_Signal is not None and hSimEventNumbers_Signal is not None:

                binx = hSimEventNumbers_Signal.GetXaxis().FindBin(chipmmFILE)
                biny = hSimEventNumbers_Signal.GetYaxis().FindBin(deltamFILE)
                binglob = hSimEventNumbers_Signal.GetBin(binx, biny)
                numSimEvents = hSimEventNumbers_Signal.GetBinContent(binglob)

                fSimEventNumbers_Signal.Close()

        elif 'SignalStopV3' in options.tag:

            with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/SUSYCrossSections13TeVstopsbottom.json') as stopxsecfile:
                stopxsec = json.load(stopxsecfile)
                crossSection = stopxsec[str(int(mstopFILE))][0]

        else:


            # TODO: adapt x-sec and nsimevents for UL

            processid = None
            if 'ZJetsToNuNu_Zpt-200toInf' in options.tag:
                processid = 'ZJetsToNuNu_Zpt-200toInf'
            elif 'WJetsToLNu_HT-100To200' in options.tag:
                processid = 'WJetsToLNu_HT-100To200'
            elif 'WJetsToLNu_HT-200To400' in options.tag:
                processid = 'WJetsToLNu_HT-200To400'
            elif 'WJetsToLNu_HT-400To600' in options.tag:
                processid = 'WJetsToLNu_HT-400To600'
            elif 'WJetsToLNu_HT-600To800' in options.tag:
                processid = 'WJetsToLNu_HT-600To800'
            elif 'WJetsToLNu_HT-800To1200' in options.tag:
                processid = 'WJetsToLNu_HT-800To1200'
            elif 'WJetsToLNu_HT-1200To2500' in options.tag:
                processid = 'WJetsToLNu_HT-1200To2500'
            elif 'WJetsToLNu_HT-2500ToInf' in options.tag:
                processid = 'WJetsToLNu_HT-2500ToInf'
            elif 'DYJetsToLL_Zpt-150toInf' in options.tag:
                processid = 'DYJetsToLL_Zpt-150toInf'
            elif 'TTJets_DiLept' in options.tag:
                processid = 'TTJets_DiLept'
            elif 'TTJets_SingleLeptFromT' in options.tag:
                processid = 'TTJets_SingleLeptFromT'
            elif 'TTJets_SingleLeptFromTbar' in options.tag:
                processid = 'TTJets_SingleLeptFromTbar'
            else:
                pass

            if processid is not None:

                with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/BkgCrossSections.json') as bkgxsecfile:
                    bkgxsec = json.load(bkgxsecfile)
                    crossSection = bkgxsec[processid]

                with open('/nfs/dust/cms/user/wolfmor/NTupleStuff/simEventNumbers_Bkg.json') as bkgnsimfile:
                    bkgnsim = json.load(bkgnsimfile)
                    numSimEvents = bkgnsim[processid]

        event_level_var_array['crossSection'][0] = crossSection
        event_level_var_array['numSimEvents'][0] = numSimEvents


        eventWeightFastSimBug = 1.
        if 'fastsim' in options.tag:
            eventWeightFastSimBug = computeFastSimBugWeight(genjets, genparticles)
        event_level_var_array['weight_fastSimBug'][0] = eventWeightFastSimBug

        weight_PU_FastFull = 1.
        weight_PU_FastFull_rebin = 1.
        weight_PU_SigBkg = 1.
        weight_PU_SigBkg_rebin = 1.
        if 'fastsim' in options.tag and 'signal' in options.tag:

            fPUweights = ROOT.TFile('/nfs/dust/cms/user/wolfmor/NTupleStuff/PUweights.root')

            hPUweightFastFull = fPUweights.Get('NPVsPerEventSignalFull:SignalFast')
            binPUweightFastFull = hPUweightFastFull.GetXaxis().FindBin(n_pv)
            weight_PU_FastFull = hPUweightFastFull.GetBinContent(binPUweightFastFull)

            hPUweightFastFull_rebin = fPUweights.Get('NPVsPerEvent_rebinFastFullSignalFull:SignalFast')
            binPUweightFastFull_rebin = hPUweightFastFull_rebin.GetXaxis().FindBin(n_pv)
            weight_PU_FastFull_rebin = hPUweightFastFull_rebin.GetBinContent(binPUweightFastFull_rebin)
            
            hPUweightSigBkg = fPUweights.Get('NPVsPerEventZJetsToNuNu_Zpt-200toInf_16:SignalFast')
            binPUweightSigBkg = hPUweightSigBkg.GetXaxis().FindBin(n_pv)
            weight_PU_SigBkg = hPUweightSigBkg.GetBinContent(binPUweightSigBkg)

            hPUweightSigBkg_rebin = fPUweights.Get('NPVsPerEvent_rebinSigBkgZJetsToNuNu_Zpt-200toInf_16:SignalFast')
            binPUweightSigBkg_rebin = hPUweightSigBkg_rebin.GetXaxis().FindBin(n_pv)
            weight_PU_SigBkg_rebin = hPUweightSigBkg_rebin.GetBinContent(binPUweightSigBkg_rebin)

            fPUweights.Close()

        event_level_var_array['weight_PU_FastFull'][0] = weight_PU_FastFull
        event_level_var_array['weight_PU_FastFull_rebin'][0] = weight_PU_FastFull_rebin
        event_level_var_array['weight_PU_SigBkg'][0] = weight_PU_SigBkg
        event_level_var_array['weight_PU_SigBkg_rebin'][0] = weight_PU_SigBkg_rebin


        event_level_var_array['n_genParticle'][0] = 0
        event_level_var_array['n_genJet'][0] = 0

        if 'data' not in options.tag:

            genparticlesfinalstate = [genparticle for genparticle in genparticles if genparticle.fromHardProcessFinalState() == 1]
                            
            event_level_var_array['n_genParticle'][0] = len(genparticlesfinalstate)

            for igp, gp in enumerate(genparticlesfinalstate):


                genparticle_var_array['genParticle_pdgId'][igp] = gp.pdgId()
                genparticle_var_array['genParticle_status'][igp] = gp.status()
                genparticle_var_array['genParticle_mass'][igp] = gp.mass()
                genparticle_var_array['genParticle_energy'][igp] = gp.energy()
                genparticle_var_array['genParticle_pt'][igp] = gp.pt()
                genparticle_var_array['genParticle_eta'][igp] = gp.eta()
                genparticle_var_array['genParticle_phi'][igp] = gp.phi()

                mm = 0
                while not gp.statusFlags().isFirstCopy():
                    for idx in range(gp.numberOfMothers()):
                        if gp.pdgId() == gp.mother(idx).pdgId():
                            gp = gp.mother(idx)
                            break
                    mm += 1
                    if mm > 100: break

                genparticle_var_array['genParticle_motherPdgId'][igp] = gp.mother(0).pdgId()
            

            event_level_var_array['n_genJet'][0] = len(genjets)

            for igj, gj in enumerate(genjets):

                genjet_var_array['genJet_pt'][igj] = gj.pt()
                genjet_var_array['genJet_eta'][igj] = gj.eta()
                genjet_var_array['genJet_phi'][igj] = gj.phi()
                genjet_var_array['genJet_mass'][igj] = gj.mass()


        event_level_var_array['n_photon'][0] = len(photons)

        numphotonsiso = 0
        for ip, p in enumerate(photons):

            photon_var_array['photon_px'][ip] = p.px()
            photon_var_array['photon_py'][ip] = p.py()
            photon_var_array['photon_pz'][ip] = p.pz()
            photon_var_array['photon_pt'][ip] = p.pt()
            photon_var_array['photon_energy'][ip] = p.energy()
            photon_var_array['photon_eta'][ip] = p.eta()
            photon_var_array['photon_phi'][ip] = p.phi()
            photon_var_array['photon_pfAbsIso'][ip], _, _, _ = calcIso_pf_or_track_new(p, pfcandsforiso0, subtractObject=True)
            photon_var_array['photon_pfAbsIsoMini'][ip], _, _, _ = calcIso_pf_or_track_new(p, pfcandsforiso0, isMini=True, subtractObject=True)
            photon_var_array['photon_chPfAbsIso'][ip], _, _, _ = calcIso_pf_or_track_new(p, chpfcandsforiso0, subtractObject=False)
            photon_var_array['photon_chPfAbsIsoMini'][ip], _, _, _ = calcIso_pf_or_track_new(p, chpfcandsforiso0, isMini=True, subtractObject=False)
            photon_var_array['photon_jetIso'][ip], photon_var_array['photon_jetIsoMulti'][ip], photon_var_array['photon_drminJet'][ip], _, photon_var_array['photon_minvJet'][ip] = calcIso_jet_new(p, jetsforiso15, isTrack=False, btagvalues=btagvalues)
            photon_var_array['photon_chHadIso'][ip] = p.chargedHadronIso()
            photon_var_array['photon_neHadIso'][ip] = p.neutralHadronIso()
            photon_var_array['photon_photIso'][ip] = p.photonIso()
            absisophoton = p.chargedHadronIso() + p.neutralHadronIso() + p.photonIso()
            photon_var_array['photon_absIso'][ip] = absisophoton
            relisophoton = absisophoton / p.pt()
            photon_var_array['photon_relIso'][ip] = relisophoton

            if relisophoton < 0.2: numphotonsiso += 1

        event_level_var_array['n_photon_iso'][0] = numphotonsiso
        hNumphotons.Fill(numphotonsiso)


        pfleptons = [l for l in pfcands
                     if ((abs(l.pdgId()) == 11 or abs(l.pdgId()) == 13)
                         and l.pt() > 10
                         and abs(l.eta()) < 2.4
                         and not (1.4442 < abs(l.eta()) < 1.566))]

        event_level_var_array['n_pfLepton'][0] = len(pfleptons)

        numpfleptonsiso = 0
        for il, l in enumerate(pfleptons):

            pflepton_var_array['pfLepton_charge'][il] = l.charge()
            pflepton_var_array['pfLepton_pdgId'][il] = l.pdgId()
            pflepton_var_array['pfLepton_px'][il] = l.px()
            pflepton_var_array['pfLepton_py'][il] = l.py()
            pflepton_var_array['pfLepton_pz'][il] = l.pz()
            pflepton_var_array['pfLepton_pt'][il] = l.pt()
            pflepton_var_array['pfLepton_energy'][il] = l.energy()
            pflepton_var_array['pfLepton_eta'][il] = l.eta()
            pflepton_var_array['pfLepton_phi'][il] = l.phi()
            if l.trackRef().isNull():
                pflepton_var_array['pfLepton_dz'][il] = -1
                pflepton_var_array['pfLepton_dxy'][il] = -1
            else:
                pflepton_var_array['pfLepton_dz'][il] = abs(l.trackRef().get().dz(pv_pos))
                pflepton_var_array['pfLepton_dxy'][il] = abs(l.trackRef().get().dxy(pv_pos))
            _, pfrelisopflepton, _, _ = calcIso_pf_or_track_new(l, pfcandsforiso0)
            pflepton_var_array['pfLepton_pfRelIso'][il] = pfrelisopflepton
            _, pflepton_var_array['pfLepton_pfRelIsoMini'][il], _, _ = calcIso_pf_or_track_new(l, pfcandsforiso0, isMini=True)
            _, pflepton_var_array['pfLepton_chPfRelIso'][il], _, _ = calcIso_pf_or_track_new(l, chpfcandsforiso0)
            _, pflepton_var_array['pfLepton_chPfRelIsoMini'][il], _, _ = calcIso_pf_or_track_new(l, chpfcandsforiso0, isMini=True)
            pflepton_var_array['pfLepton_jetIso'][il], pflepton_var_array['pfLepton_jetIsoMulti'][il], pflepton_var_array['pfLepton_drminJet'][il], _, pflepton_var_array['pfLepton_minvJet'][il] = calcIso_jet_new(l, jetsforiso15, isTrack=False, btagvalues=btagvalues)

            if pfrelisopflepton < 0.2: numpfleptonsiso += 1

        event_level_var_array['n_pfLepton_iso'][0] = numpfleptonsiso


        event_level_var_array['n_tau'][0] = len(tauswithdiscriminators)

        # for info on MC matching and categories see
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/HiggsToTauTauWorking2016#MC_Matching
        # https://github.com/cms-tau-pog/TauIDSFs
        if 'data' not in options.tag and (len(tauswithdiscriminators) > 0 or 'genmatchtracks' in options.tag or 'genmatchalltracks' in options.tag):

            # key = category ID
            genparticlesfortaumatching = {}

            # prompt electrons
            genparticlesfortaumatching[1] = [gp for gp in genparticles if abs(gp.pdgId()) == 11
                                             and gp.pt() > 8 and gp.statusFlags().isPrompt()]

            # prompt muons
            genparticlesfortaumatching[2] = [gp for gp in genparticles if abs(gp.pdgId()) == 13
                                             and gp.pt() > 8 and gp.statusFlags().isPrompt()]

            # electrons from tau decay
            genparticlesfortaumatching[3] = [gp for gp in genparticles if abs(gp.pdgId()) == 11
                                             and gp.pt() > 8 and gp.statusFlags().isDirectPromptTauDecayProduct()]

            # muons from tau decay
            genparticlesfortaumatching[4] = [gp for gp in genparticles if abs(gp.pdgId()) == 13
                                             and gp.pt() > 8 and gp.statusFlags().isDirectPromptTauDecayProduct()]

            # real taus (visibile gen jet)
            genparticlesfortaumatching[5] = []
            gentaujets = []
            for gp in genparticles:
                if abs(gp.pdgId()) == 15 and gp.statusFlags().isPrompt():
                    gentaujet = ROOT.TLorentzVector()
                    for itaudaughter in range(gp.numberOfDaughters()):
                        if abs(gp.daughter(itaudaughter).pdgId()) in [11, 12, 13, 14, 15, 16]: continue
                        taudaughterTlv = ROOT.TLorentzVector(gp.daughter(itaudaughter).px(), gp.daughter(itaudaughter).py(), gp.daughter(itaudaughter).pz(), gp.daughter(itaudaughter).energy())
                        gentaujet += taudaughterTlv
                    if gentaujet.Pt() > 0:
                        gentaujets.append(Dummy(gentaujet))
                    if gentaujet.Pt() > 15:
                        genparticlesfortaumatching[5].append(Dummy(gentaujet))

        numtausvloose = 0
        numtausloose = 0
        numtausmedium = 0
        numtaustight = 0
        numtausvtight = 0
        numtausvvtight = 0
        numtausvloosePt20 = 0
        numtausloosePt20 = 0
        numtausmediumPt20 = 0
        numtaustightPt20 = 0
        numtausvtightPt20 = 0
        numtausvvtightPt20 = 0
        for it, t in enumerate(tauswithdiscriminators):

            genmatchtau = 6
            genmatchpdgidtau = -1.
            drmintaumatch = 9.

            if 'data' not in options.tag:
                for taumatchcategory in [1, 2, 3, 4, 5]:
                    idxtaumatch, drmintaumatch = findMatch_gen_old_easy(t[0], genparticlesfortaumatching[taumatchcategory])
                    if drmintaumatch < 0.2:
                        genmatchtau = taumatchcategory
                        genmatchpdgidtau = 15 if taumatchcategory == 5 else genparticlesfortaumatching[taumatchcategory][idxtaumatch].pdgId()

            tau_var_array['tau_genMatch'][it] = genmatchtau
            tau_var_array['tau_genMatchPdgId'][it] = genmatchpdgidtau
            tau_var_array['tau_genMatchDr'][it] = drmintaumatch

            tes = 1.
            if genmatchtau == 5:
                tes = teshist.GetBinContent(teshist.GetXaxis().FindBin(t[0].decayMode()))

            tau_var_array['tau_sanityDm'][it] = t[11] - t[0].pt()
            tau_var_array['tau_sanityRaw'][it] = t[12] - t[0].pt()
            tau_var_array['tau_sanityVloose'][it] = t[13] - t[0].pt()

            tau_var_array['tau_charge'][it] = t[0].charge()
            tau_var_array['tau_px'][it] = t[0].px()
            tau_var_array['tau_py'][it] = t[0].py()
            tau_var_array['tau_pz'][it] = t[0].pz()
            tau_var_array['tau_pt'][it] = t[0].pt()
            tau_var_array['tau_energy'][it] = t[0].energy()
            tau_var_array['tau_mass'][it] = t[0].mass()
            tau_var_array['tau_ptScaled'][it] = t[0].pt() * tes
            tau_var_array['tau_energyScaled'][it] = t[0].energy() * tes
            tau_var_array['tau_massScaled'][it] = t[0].mass() * tes
            tau_var_array['tau_eta'][it] = t[0].eta()
            tau_var_array['tau_phi'][it] = t[0].phi()
            try:
                if t[0].leadPFChargedHadrCand().trackRef().isNonnull() and t[0].leadPFChargedHadrCand().trackRef().get():
                    leadpfchhadcand = t[0].leadPFChargedHadrCand().trackRef().get()
                    tau_var_array['tau_dz'][it] = abs(leadpfchhadcand.dz(pv_pos))
                    tau_var_array['tau_dxy'][it] = abs(leadpfchhadcand.dxy(pv_pos))
                    tau_var_array['tau_leadPfChHadCandPt'][it] = leadpfchhadcand.pt()
                    tau_var_array['tau_leadPfChHadCandEta'][it] = leadpfchhadcand.eta()
                    tau_var_array['tau_leadPfChHadCandPhi'][it] = leadpfchhadcand.phi()
                else:
                    tau_var_array['tau_dz'][it] = -1.
                    tau_var_array['tau_dxy'][it] = -1.
                    tau_var_array['tau_leadPfChHadCandPt'][it] = -1.
                    tau_var_array['tau_leadPfChHadCandEta'][it] = -1.
                    tau_var_array['tau_leadPfChHadCandPhi'][it] = -1.
            except :
                tau_var_array['tau_dz'][it] = -1.
                tau_var_array['tau_dxy'][it] = -1.
                tau_var_array['tau_leadPfChHadCandPt'][it] = -1.
                tau_var_array['tau_leadPfChHadCandEta'][it] = -1.
                tau_var_array['tau_leadPfChHadCandPhi'][it] = -1.
            tau_var_array['tau_chHadIso'][it] = t[0].isolationPFChargedHadrCandsPtSum()
            tau_var_array['tau_photIso'][it] = t[0].isolationPFGammaCandsEtSum()

            tau_var_array['tau_decayMode'][it] = t[0].decayMode()
            tau_var_array['tau_decayModeFinding'][it] = t[1]
            tau_var_array['tau_mvaDiscr'][it] = t[2]

            tau_var_array['tau_isvloose'][it] = t[3]
            tau_var_array['tau_isloose'][it] = t[4]
            tau_var_array['tau_ismedium'][it] = t[5]
            tau_var_array['tau_istight'][it] = t[6]
            tau_var_array['tau_isvtight'][it] = t[7]
            tau_var_array['tau_isvvtight'][it] = t[8]

            if t[3] > 0.5:
                numtausvloose += 1
                if t[0].pt() > 20: numtausvloosePt20 += 1
            if t[4] > 0.5:
                numtausloose += 1
                if t[0].pt() > 20: numtausloosePt20 += 1
            if t[5] > 0.5:
                numtausmedium += 1
                if t[0].pt() > 20: numtausmediumPt20 += 1
            if t[6] > 0.5:
                numtaustight += 1
                if t[0].pt() > 20: numtaustightPt20 += 1
            if t[7] > 0.5:
                numtausvtight += 1
                if t[0].pt() > 20: numtausvtightPt20 += 1
            if t[8] > 0.5:
                numtausvvtight += 1
                if t[0].pt() > 20: numtausvvtightPt20 += 1

            tau_var_array['tau_elRejection'][it] = t[9]
            tau_var_array['tau_muRejection'][it] = t[10]

        hNumtaus.Fill(len(tauswithdiscriminators))
        event_level_var_array['n_tau_vloose'][0] = numtausvloose
        event_level_var_array['n_tau_loose'][0] = numtausloose
        event_level_var_array['n_tau_medium'][0] = numtausmedium
        event_level_var_array['n_tau_tight'][0] = numtaustight
        event_level_var_array['n_tau_vtight'][0] = numtausvtight
        event_level_var_array['n_tau_vvtight'][0] = numtausvvtight
        event_level_var_array['n_tau_20_vloose'][0] = numtausvloosePt20
        event_level_var_array['n_tau_20_loose'][0] = numtausloosePt20
        event_level_var_array['n_tau_20_medium'][0] = numtausmediumPt20
        event_level_var_array['n_tau_20_tight'][0] = numtaustightPt20
        event_level_var_array['n_tau_20_vtight'][0] = numtausvtightPt20
        event_level_var_array['n_tau_20_vvtight'][0] = numtausvvtightPt20


        tracksByPV = {}
        for ipv, pv in enumerate(primaryvertices):

            pv_var_array['pv_idx'][ipv] = ipv
            pv_var_array['pv_numTracks'][ipv] = pv.tracksSize()
            pv_var_array['pv_x'][ipv] = pv.position().x()
            pv_var_array['pv_y'][ipv] = pv.position().y()
            pv_var_array['pv_z'][ipv] = pv.position().z()
            tracksByPV[ipv] = [pv.trackRefAt(i).get() for i in range(pv.tracksSize())]


        if 'era16' in options.tag:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation80XReReco
            loosewp = 0.5426
            mediumwp = 0.8484
            tightwp = 0.9535
        elif 'era17' in options.tag:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            loosewp = 0.5803
            mediumwp = 0.8838
            tightwp = 0.9693
        elif 'era18' in options.tag:
            # same as for 17... because no info for 18
            loosewp = 0.5803
            mediumwp = 0.8838
            tightwp = 0.9693
        else:
            raise NotImplementedError('btag: era unknown or not specified')

        isleptonjetidxlist = []
        btaglist = []
        for ijet, jet in enumerate(jets):

            jet_var_array['jet_px'][ijet] = jet.px()
            jet_var_array['jet_py'][ijet] = jet.py()
            jet_var_array['jet_pz'][ijet] = jet.pz()
            jet_var_array['jet_pt'][ijet] = jet.pt()
            jet_var_array['jet_energy'][ijet] = jet.energy()
            jet_var_array['jet_eta'][ijet] = jet.eta()
            jet_var_array['jet_phi'][ijet] = jet.phi()
            jet_var_array['jet_numConstituents'][ijet] = jet.numberOfDaughters()

            dRmin = 0.1
            thisbtag = -2.
            for ib in range(nbtags):
                dR = deltaR(handle_btag.product().key(ib).get().eta(), jet.eta(), handle_btag.product().key(ib).get().phi(), jet.phi())
                if dR < dRmin:
                    dRmin = dR
                    thisbtag = max(-1.0, handle_btag.product().value(ib))

            jet_var_array['jet_btag'][ijet] = thisbtag
            btaglist.append(thisbtag)
            btagvalues.append((thisbtag, jet.pt(), jet.eta()))
            hBtagjets.Fill(thisbtag)

            drminleptonjet = 9.
            ptclosestleptonjet = -1.
            isleptonjet = 0.

            ptthreshold_leptonjet = 50

            for e in electrons:
                if e.pt() < ptthreshold_leptonjet: continue
                dR_jetlepton = deltaR(e.eta(), jet.eta(), e.phi(), jet.phi())
                if dR_jetlepton < drminleptonjet:
                    drminleptonjet = dR_jetlepton
                    ptclosestleptonjet = e.pt()

            for m in muons:
                if m.pt() < ptthreshold_leptonjet: continue
                dR_jetlepton = deltaR(m.eta(), jet.eta(), m.phi(), jet.phi())
                if dR_jetlepton < drminleptonjet:
                    drminleptonjet = dR_jetlepton
                    ptclosestleptonjet = m.pt()

            if drminleptonjet < 0.1 and ptclosestleptonjet / jet.pt() > 0.8:
                isleptonjet = 1.
                isleptonjetidxlist.append(ijet)

            jet_var_array['jet_drminLepton'][ijet] = drminleptonjet
            jet_var_array['jet_ptClosestLepton'][ijet] = ptclosestleptonjet
            jet_var_array['jet_isLepton'][ijet] = isleptonjet

            drmingenjetjet = 9.
            ptclosestgenjetjet = -1.
            isgenjetjet = 0.

            if 'data' not in options.tag:

                idx_genjet, drmingenjetjet = findMatch_gen_old_easy(jet, genjets)
                if not idx_genjet == -1:
                    ptclosestgenjetjet = genjets[idx_genjet].pt()

                if drmingenjetjet < 0.2:
                    isgenjetjet = 1.

            jet_var_array['jet_drminGenJet'][ijet] = drmingenjetjet
            jet_var_array['jet_ptClosestGenJet'][ijet] = ptclosestgenjetjet
            jet_var_array['jet_isGenJet'][ijet] = isgenjetjet


        event_level_var_array['n_jet_30_btagloose'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > loosewp and jetpt > 30 and abs(jeteta) < 2.4)])
        event_level_var_array['n_jet_15_btagloose'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > loosewp and jetpt > 15 and abs(jeteta) < 2.4)])

        njetsbtagmedium = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > mediumwp and jetpt > 30 and abs(jeteta) < 2.4)])
        hNjetsbtagmedium.Fill(njetsbtagmedium)
        event_level_var_array['n_jet_30_btagmedium'][0] = njetsbtagmedium
        event_level_var_array['n_jet_15_btagmedium'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > mediumwp and jetpt > 15 and abs(jeteta) < 2.4)])
        
        event_level_var_array['n_jet_30_btagtight'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > tightwp and jetpt > 30 and abs(jeteta) < 2.4)])
        event_level_var_array['n_jet_15_btagtight'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > tightwp and jetpt > 15 and abs(jeteta) < 2.4)])

        mtmetleadingjet = ROOT.TMath.Sqrt(2 * met.pt() * jets[idxhighestptjet].pt()
                                          * (1 - ROOT.TMath.Cos(deltaPhi(met.phi(), jets[idxhighestptjet].phi()))))
        hMtmetleadingjet.Fill(mtmetleadingjet)
        event_level_var_array['mtMetLeadingJet'][0] = mtmetleadingjet


        if phifirsttrack == tracks[0].phi() and etafirsttrack == tracks[0].eta():
            print 'suspicious... better get out of here!'
            sys.exit(1)
        phifirsttrack = tracks[0].phi()
        etafirsttrack = tracks[0].eta()


        '''
        ###############################################################################################
        # get sv-level info
        ###############################################################################################
        '''
        numsvsfinalpreselection = 0  
        n_sv_total = 0
        numSVs = 0

        if not 'skipSVs' in options.tag:
            n_sv_total = len(filesWithSV[ifile][ievent]) 
            numSVs = len(filesWithSV[ifile][ievent]) 
            n_selTracks = len(filesWithTrackID[ifile][ievent])
            selectedTracks = np.array([None]*len(filesWithTrackID[ifile][ievent]))
            
            for idx, index in enumerate(filesWithTrackID[ifile][ievent]):
                selectedTracks[idx]=deepcopy(tracks[int(index)])

            for nSV, secondary in enumerate(filesWithSV[ifile][ievent]):

                isSignal = 0
                proceed = False	
                
                
                idxHelical = [-1,-1]
                drmin = [999,999]
                dxyzmin = [999,999]
                tmin = [0,0]
                
                ClassicIdx = [-1,-1]
                ClassicDrmin = [999,999]
                dxyzmin = [999,999]
                tmin = [0,0]
                
                matchingTrk = [None, None]
                matchingTrkIdx = [None, None]

                SV_level_var_array['isSignal'][nSV] = isSignal
                SV_level_var_array['vtxDCA'][nSV] = filesWithDCA[ifile][ievent][nSV]
                SV_level_var_array['log10vtxChi2'][nSV] = ROOT.TMath.Log10(secondary.vertexChi2())
                SV_level_var_array['vtxChi2Ndof'][nSV] = ROOT.TMath.Log10(secondary.vertexNormalizedChi2())
                SV_level_var_array['vtxVx'][nSV] = secondary.vx()
                SV_level_var_array['vtxVy'][nSV] = secondary.vy()
                SV_level_var_array['vtxVz'][nSV] = secondary.vz()
                SV_level_var_array['vtxdx'][nSV] = abs(secondary.vx()-pv_pos.x())
                SV_level_var_array['vtxdy'][nSV] = abs(secondary.vy()-pv_pos.y())
                SV_level_var_array['vtxdz'][nSV] = abs(secondary.vz()-pv_pos.z())
                SV_level_var_array['log10vtxdxy'][nSV] = ROOT.TMath.Log10(sqrt(pow((secondary.vx()-pv_pos.x()),2)+pow((secondary.vy()-pv_pos.y()),2))+0.00001)
                SV_level_var_array['log10vtxdz'][nSV] = ROOT.TMath.Log10(sqrt(pow((secondary.vz()-pv_pos.z()),2))+0.00001)
                SV_level_var_array['numberofdaughters'][nSV] = secondary.numberOfDaughters()
                
                ### find matching tracks

                for k in range(secondary.numberOfDaughters()):
                    
                    if 'debug' in options.tag: print "SV no. ", nSV, "daughter no. ", k, " charge ", secondary.daughter(k).charge()
                    
                    idxHelical[k], dxyzmin[k], tmin[k], drmin[k] = findMatch_tracktrack_new(secondary.daughter(k), selectedTracks)
                    _, ClassicIdx[k], ClassicDrmin[k], _ = findMinDr_track(secondary.daughter(k), selectedTracks, 20.)
                    
                    if 'debug' in options.tag: print "matching", idxHelical[k], dxyzmin[k], tmin[k], drmin[k]
                    if 'debug' in options.tag: print "matching", ClassicIdx[k], ClassicDrmin[k]
                    
                            
                    if ((ClassicDrmin[0]<0.02) and (ClassicDrmin[1]<0.02)):
                        
                        
                        if 'debug' in options.tag: print "classic match"
                        
                        for i, idx in enumerate(ClassicIdx):
                            matchingTrk[i] = selectedTracks[idx]
                            matchingTrkIdx[i] = idx 
                            
                            
                        if 'signal' in options.tag:
                            if ClassicIdx[0] in signalTrkIdx and ClassicIdx[1] in signalTrkIdx:
                                if 'debug' in options.tag:print '---------------------'
                                if 'debug' in options.tag:print 'is Signal SV'
                                if 'debug' in options.tag:print '---------------------'
                                isSignal = 1
                                hasSignalSV = 1
                                signalIdx = nSV

                            else:
                                if 'debug' in options.tag:print '---------------------'
                                if 'debug' in options.tag:print 'is Back SV'
                                if 'debug' in options.tag:print '---------------------'


                        elif ((dxyzmin[0]<0.04 and drmin[0]<0.02) and (dxyzmin[1]<0.04 and drmin[1]<0.02)):
                            
                            if 'debug' in options.tag: print "dxyz match"
                            
                            for i, idx in enumerate(idxHelical):
                                matchingTrk[i] = selectedTracks[idx]
                                matchingTrkIdx[i] = idx 
                            
                            if 'signal' in options.tag:
                                if idxHelical[0] in signalTrkIdx and idxHelical[1] in signalTrkIdx:
                                    if 'debug' in options.tag:print '---------------------'
                                    if 'debug' in options.tag:print 'is Signal SV'
                                    if 'debug' in options.tag:print '---------------------'
                                    isSignal = 1
                                    hasSignalSV = 1
                                    signalIdx = nSV

                                else:
                                    if 'debug' in options.tag:print '---------------------'
                                    if 'debug' in options.tag:print 'is Back SV'
                                    if 'debug' in options.tag:print '---------------------'


                ######################################
                #### "filling tree on SV level"
                ######################################
                if None in matchingTrk: 
                    if matchingTrk[0] == matchingTrk[1]: 
                        SV_level_var_array['hasTrackMatch_Low'][nSV] = 0
                        SV_level_var_array['hasTrackMatch_High'][nSV] = 0
                    else: 
                        SV_level_var_array['hasTrackMatch_Low'][nSV] = 0
                        SV_level_var_array['hasTrackMatch_High'][nSV] = 1
                        
                    numsvsfinalpreselection += 1
                    continue
                    
                SV_level_var_array['hasTrackMatch_Low'][nSV] = 1
                SV_level_var_array['hasTrackMatch_High'][nSV] = 1
                
                if 'debug' in options.tag: print "filling tree on SV level, nSV", nSV
                TLV1 = TLorentzVector()
                TLV1.SetPxPyPzE(matchingTrk[0].px(),matchingTrk[0].py(),matchingTrk[0].pz(),matchingTrk[0].pt()*TMath.CosH(matchingTrk[0].eta()))
                TLV2 = TLorentzVector()
                TLV2.SetPxPyPzE(matchingTrk[1].px(),matchingTrk[1].py(),matchingTrk[1].pz(),matchingTrk[1].pt()*TMath.CosH(matchingTrk[1].eta()))
                
                PVVtx = TVector3(secondary.vx()-pv_pos.x(), secondary.vy()-pv_pos.y(), secondary.vz()-pv_pos.z())
                summedTracks = TLV1+TLV2
                
                if matchingTrk[0].pt() == min(matchingTrk[0].pt(), matchingTrk[1].pt()):
                    trackLow = matchingTrk[0]
                    idxLow = matchingTrkIdx[0]
                    trackHigh = matchingTrk[1]
                    idxHigh = matchingTrkIdx[1]
                else:
                    trackLow = matchingTrk[1]
                    idxLow = matchingTrkIdx[1]
                    trackHigh = matchingTrk[0]
                    idxHigh = matchingTrkIdx[0]
                                    
                SV_level_var_array['deltaPhiMetSumPt'][nSV] = deltaPhi(summedTracks.Phi(), met.phi())

                if lJet1 == None: SV_level_var_array['deltaPhiLeadingJetSumPt'][nSV] = -999
                else: 
                    SV_level_var_array['deltaPhiLeadingJetSumPt'][nSV] = deltaPhi(summedTracks.Phi(),lJet1.phi())
                if lJet1 == None: SV_level_var_array['deltaEtaLeadingJetSumPt'][nSV]= -999
                else:
                    SV_level_var_array['deltaEtaLeadingJetSumPt'][nSV]= summedTracks.Eta()-lJet1.eta()
                            

                SV_level_var_array['deltaEtaPVVtxToTrack_High'][nSV] = abs(trackHigh.eta()-PVVtx.Eta())
                SV_level_var_array['deltaEtaPVVtxToTrack_Low'][nSV] = abs(trackLow.eta()-PVVtx.Eta())
                SV_level_var_array['deltaEtaPVVtxToTrackSum'][nSV] = summedTracks.Eta()-PVVtx.Eta()
                SV_level_var_array['absdeltaEtaPVVtxToTrackSum'][nSV] = abs(summedTracks.Eta()-PVVtx.Eta())
                SV_level_var_array['deltaPhiPVVtxToTrack_High'][nSV] = deltaPhi(trackHigh.phi(), PVVtx.Phi())
                SV_level_var_array['deltaPhiPVVtxToTrack_Low'][nSV] = deltaPhi(trackLow.phi(), PVVtx.Phi())
                SV_level_var_array['deltaPhiPVVtxToTrackSum'][nSV] = deltaPhi(summedTracks.Phi(), PVVtx.Phi())
                SV_level_var_array['PVVtxEta'][nSV] = PVVtx.Eta()
                SV_level_var_array['PVVtxPhi'][nSV] = PVVtx.Phi()
                SV_level_var_array['vtxiso'][nSV] , SV_level_var_array['vtxdrmin'][nSV], SV_level_var_array['vtxnumneighbours'][nSV] = calcIso_vtx(secondary, secondaryVertices)	
                SV_level_var_array['deltaEtaPVVtxToMET'][nSV] = met.eta()-PVVtx.Eta()
                SV_level_var_array['absdeltaEtaPVVtxToMET'][nSV] = abs(met.eta()-PVVtx.Eta())
                SV_level_var_array['deltaPhiPVVtxToMET'][nSV] = met.phi()-PVVtx.Phi()
                SV_level_var_array['absdeltaPhiPVVtxToMET'][nSV] = abs(met.phi()-PVVtx.Phi())

                v_normal_reco = TVector3(secondary.vx()-pv_pos.x(),secondary.vy()-pv_pos.y(),secondary.vz()-pv_pos.z())
                normalvector_reco = v_normal_reco.Unit()
                
                ptZstar_reco = TLV1.Pt()+TLV2.Pt()
                v_pZstar_reco = TVector3(TLV1.Px()+TLV2.Px(), TLV1.Py()+TLV2.Py(), TLV1.Pz()+TLV2.Pz())

                SV_level_var_array['ptZstar_reco'][nSV] = ptZstar_reco
                SV_level_var_array['abspZstar_reco'][nSV] = v_pZstar_reco.Mag()
                SV_level_var_array['absnormalVector_reco'][nSV] = normalvector_reco.Mag()
                SV_level_var_array['beta_reco'][nSV] = v_pZstar_reco.Angle(normalvector_reco)
                                
                muonmass2 = 0.1056 *0.1056
                electronmass2 = 0.0005* 0.0005
                pionmass2 = 0.1396*0.1396
                
                ZstarP4_muoncase = TLorentzVector()
                ZstarP4_electroncase = TLorentzVector()

                track1E_muoncase =TMath.Sqrt(pow(TLV1.Px(),2)+pow(TLV1.Py(),2)+pow(TLV1.Pz(),2) + muonmass2)
                track2E_muoncase =TMath.Sqrt(pow(TLV2.Px(),2)+pow(TLV2.Py(),2)+pow(TLV2.Pz(),2) + muonmass2)

                ZstarETot_muoncase = track1E_muoncase + track2E_muoncase
                ZstarP4_muoncase.SetPxPyPzE(TLV1.Px()+TLV2.Px(), TLV1.Py()+TLV2.Py(), TLV1.Pz()+TLV2.Pz(), ZstarETot_muoncase)			
                mZstar_reco_muoncase = sqrt(ZstarP4_muoncase*ZstarP4_muoncase)
                mtransverse2_reco_muoncase = mZstar_reco_muoncase*mZstar_reco_muoncase + ((v_pZstar_reco.Cross(normalvector_reco))*(v_pZstar_reco.Cross(normalvector_reco)))

                SV_level_var_array['mZstar_reco_muon'][nSV] = mZstar_reco_muoncase
                SV_level_var_array['mtransverse2_reco_muon'][nSV] = mtransverse2_reco_muoncase
                theta = asin((v_pZstar_reco.Cross(normalvector_reco)).Mag()/(v_pZstar_reco.Mag()))
                SV_level_var_array['theta_reco'][nSV] = degrees(theta)
                SV_level_var_array['error_mtransverse2_reco_muon'][nSV] = pow(mZstar_reco_muoncase*mZstar_reco_muoncase+pow(v_pZstar_reco.Mag(),2)*pow(sin(theta),2),-1/2)*sin(theta)*cos(theta)*pow(v_pZstar_reco.Mag(),2)

                track1E_electroncase =TMath.Sqrt(pow(TLV1.Px(),2)+pow(TLV1.Py(),2)+pow(TLV1.Pz(),2) + muonmass2)
                track2E_electroncase =TMath.Sqrt(pow(TLV2.Px(),2)+pow(TLV2.Py(),2)+pow(TLV2.Pz(),2) + muonmass2)

                ZstarETot_electroncase = track1E_electroncase + track2E_electroncase
                ZstarP4_electroncase.SetPxPyPzE(TLV1.Px()+TLV2.Px(), TLV1.Py()+TLV2.Py(), TLV1.Pz()+TLV2.Pz(), ZstarETot_electroncase)			
                mZstar_reco_electroncase = sqrt(ZstarP4_electroncase*ZstarP4_electroncase)
                mtransverse2_reco_electroncase = mZstar_reco_electroncase*mZstar_reco_electroncase + ((v_pZstar_reco.Cross(normalvector_reco))*(v_pZstar_reco.Cross(normalvector_reco)))
                
                SV_level_var_array['mZstar_reco_electron'][nSV] = mZstar_reco_electroncase
                SV_level_var_array['mtransverse2_reco_electron'][nSV] = mtransverse2_reco_electroncase
                
                ip_Low = IPcalculator(trackLow, primaryvertices[0])
                ip_High = IPcalculator(trackHigh, primaryvertices[0])
                SV_level_var_array['IPsignificance_Low'][nSV] = ip_Low.getIPsignificance()
                SV_level_var_array['IPxyz_Low'][nSV] = ip_Low.getIP()
                SV_level_var_array['IPxy_Low'][nSV] = ip_Low.getDxy()
                SV_level_var_array['IPz_Low'][nSV] = ip_Low.getDz()
                SV_level_var_array['log10IPsignificance_Low'][nSV] = ROOT.TMath.Log10(ip_Low.getIPsignificance())
                SV_level_var_array['log10IPxyz_Low'][nSV] = ROOT.TMath.Log10(ip_Low.getIP())
                SV_level_var_array['log10IPxy_Low'][nSV] = ROOT.TMath.Log10( ip_Low.getDxy())
                SV_level_var_array['log10IPz_Low'][nSV] = ROOT.TMath.Log10(ip_Low.getDz())

                SV_level_var_array['IPsignificance_High'][nSV] = ip_High.getIPsignificance()
                SV_level_var_array['IPxyz_High'][nSV] = ip_High.getIP()
                SV_level_var_array['IPxy_High'][nSV] = ip_High.getDxy()
                SV_level_var_array['IPz_High'][nSV] = ip_High.getDz()
                SV_level_var_array['log10IPsignificance_High'][nSV] = ROOT.TMath.Log10(ip_High.getIPsignificance())
                SV_level_var_array['log10IPxyz_High'][nSV] = ROOT.TMath.Log10(ip_High.getIP())
                SV_level_var_array['log10IPxy_High'][nSV] = ROOT.TMath.Log10(ip_High.getDxy())
                SV_level_var_array['log10IPz_High'][nSV] = ROOT.TMath.Log10(ip_High.getDz())

                
                minipPU_Low = None
                minivPU_Low= -1
                minIPsignificancePU_Low = 999
                
                minipPU_High = None
                minivPU_High= -1
                minIPsignificancePU_High = 999
                
                for iv, v in enumerate(primaryvertices[1:]):
                    thisipPU_Low = IPcalculator(trackLow, v)
                    thisIPsignificancePU_Low = thisipPU_Low.getIPsignificance()
                    if thisIPsignificancePU_Low < minIPsignificancePU_Low:
                        minipPU_Low = thisipPU_Low
                        minivPU_Low = iv
                        minIPsignificancePU_Low = thisIPsignificancePU_Low
                        
                    thisipPU_High = IPcalculator(trackHigh, v)
                    thisIPsignificancePU_High = thisipPU_High.getIPsignificance()
                    if thisIPsignificancePU_High < minIPsignificancePU_High:
                        minipPU_High = thisipPU_High
                        minivPU_High = iv
                        minIPsignificancePU_High = thisIPsignificancePU_High
                if not minivPU_Low == -1:
                    #SV_level_var_array['idxpvPU'][nSV] = minivPU+1
                    SV_level_var_array['IPsignificancePU_Low'][nSV] = minipPU_Low.getIPsignificance()
                    SV_level_var_array['IPxyzPU_Low'][nSV] = minipPU_Low.getIP()
                    SV_level_var_array['IPxyPU_Low'][nSV] = minipPU_Low.getDxy()
                    SV_level_var_array['IPzPU_Low'][nSV] = minipPU_Low.getDz()
                    SV_level_var_array['log10IPsignificancePU_Low'][nSV] = ROOT.TMath.Log10(minipPU_Low.getIPsignificance())
                    SV_level_var_array['log10IPxyzPU_Low'][nSV] = ROOT.TMath.Log10(minipPU_Low.getIP())
                    SV_level_var_array['log10IPxyPU_Low'][nSV] = ROOT.TMath.Log10(minipPU_Low.getDxy())
                    SV_level_var_array['log10IPzPU_Low'][nSV] = ROOT.TMath.Log10(minipPU_Low.getDz())

                
                if not minivPU_High == -1:
                    #SV_level_var_array['idxpvPU'][nSV] = minivPU+1
                    SV_level_var_array['IPsignificancePU_High'][nSV] = minipPU_High.getIPsignificance()
                    SV_level_var_array['IPxyzPU_High'][nSV] = minipPU_High.getIP()
                    SV_level_var_array['IPxyPU_High'][nSV] = minipPU_High.getDxy()
                    SV_level_var_array['IPzPU_High'][nSV] = minipPU_High.getDz()
                    SV_level_var_array['log10IPsignificancePU_High'][nSV] = ROOT.TMath.Log10(minipPU_High.getIPsignificance())
                    SV_level_var_array['log10IPxyzPU_High'][nSV] = ROOT.TMath.Log10(minipPU_High.getIP())
                    SV_level_var_array['log10IPxyPU_High'][nSV] = ROOT.TMath.Log10(minipPU_High.getDxy())
                    SV_level_var_array['log10IPzPU_High'][nSV] = ROOT.TMath.Log10(minipPU_High.getDz())

                
                mounMatch_Low, mounDR_Low = matchToMuon(trackLow, muons)
                if mounDR_Low < 0.01 : 
                    SV_level_var_array['muonMatched_Low'][nSV] = 1
                    SV_level_var_array['numberOfChambers_Low'][nSV] = muons[mounMatch_Low].numberOfChambers()
                    SV_level_var_array['numberOfMatchedStations_Low'][nSV] = muons[mounMatch_Low].numberOfMatchedStations()
                    #SV_level_var_array['numberOfSegments_Low'][nSV] = muons[mounMatch_Low].numberOfSegments()
                    SV_level_var_array['isGlobalMuon_Low'][nSV] = muons[mounMatch_Low].isGlobalMuon()

                    
                    if muons[mounMatch_Low].isTrackerMuon() and not muons[mounMatch_Low].innerTrack()==None:
                        
                        SV_level_var_array['isTrackerMuon_Low'][nSV] = 1
                        SV_level_var_array['normalizedChi2Muon_Low'][nSV] = muons[mounMatch_Low].innerTrack().normalizedChi2()
                        SV_level_var_array['trackerLayersWithMeasurementMuon_Low'][nSV] = muons[mounMatch_Low].innerTrack().hitPattern().trackerLayersWithMeasurement()
                        SV_level_var_array['pixelLayersWithMeasurementMuon_Low'][nSV] = muons[mounMatch_Low].innerTrack().hitPattern().pixelLayersWithMeasurement()
                        SV_level_var_array['dxyPVMuon_Low'][nSV] = abs(muons[mounMatch_Low].innerTrack().dxy(pv_pos))
                        SV_level_var_array['dzPVMuon_Low'][nSV] = abs(muons[mounMatch_Low].innerTrack().dz(pv_pos))	
                        if ((muons[mounMatch_Low].innerTrack().hitPattern().trackerLayersWithMeasurement() > 10) and (muons[mounMatch_Low].innerTrack().hitPattern().pixelLayersWithMeasurement() > 2) and (muons[mounMatch_Low].innerTrack().normalizedChi2() < 1.8) and  (abs(muons[mounMatch_Low].innerTrack().dxy(pv_pos)) < 3) and (abs(muons[mounMatch_Low].innerTrack().dz(pv_pos)) < 20)):							
                            SV_level_var_array['isSoftMuon_Low'][nSV] = 1
                            hasSoftMuon = 1	
                        else: SV_level_var_array['isSoftMuon_Low'][nSV] = 0
                    else: SV_level_var_array['isSoftMuon_Low'][nSV] = 0
                                            
                else: 
                    SV_level_var_array['muonMatched_Low'][nSV] = 0
                    SV_level_var_array['isSoftMuon_Low'][nSV] = 0

                mounMatch_High, mounDR_High = matchToMuon(trackHigh, muons)
                if mounDR_High < 0.01 : 
                    SV_level_var_array['muonMatched_High'][nSV] = 1
                    SV_level_var_array['numberOfChambers_High'][nSV] = muons[mounMatch_High].numberOfChambers()
                    SV_level_var_array['numberOfMatchedStations_High'][nSV] = muons[mounMatch_High].numberOfMatchedStations()
                    #SV_level_var_array['numberOfSegments_High'][nSV] = muons[mounMatch_High].numberOfSegments()
                    SV_level_var_array['isGlobalMuon_High'][nSV] = muons[mounMatch_High].isGlobalMuon()
                    #SV_level_var_array['isGoodMuon_High'][nSV] = muons[mounMatch_High].isGoodMuon()
                    
                    if muons[mounMatch_High].isTrackerMuon() and not muons[mounMatch_High].innerTrack()==None:

                        SV_level_var_array['isTrackerMuon_High'][nSV] = 1
                        SV_level_var_array['normalizedChi2Muon_High'][nSV] = muons[mounMatch_High].innerTrack().normalizedChi2()					
                        SV_level_var_array['trackerLayersWithMeasurementMuon_High'][nSV] = muons[mounMatch_High].innerTrack().hitPattern().trackerLayersWithMeasurement()
                        SV_level_var_array['pixelLayersWithMeasurementMuon_High'][nSV] = muons[mounMatch_High].innerTrack().hitPattern().pixelLayersWithMeasurement()	
                        SV_level_var_array['dxyPVMuon_High'][nSV] = abs(muons[mounMatch_High].innerTrack().dxy(pv_pos))
                        SV_level_var_array['dzPVMuon_High'][nSV] = abs(muons[mounMatch_High].innerTrack().dz(pv_pos))	
                        if ((muons[mounMatch_High].innerTrack().hitPattern().trackerLayersWithMeasurement() > 10) and (muons[mounMatch_High].innerTrack().hitPattern().pixelLayersWithMeasurement() > 2) and (muons[mounMatch_High].innerTrack().normalizedChi2() < 1.8) and  (abs(muons[mounMatch_High].innerTrack().dxy(pv_pos)) < 3) and (abs(muons[mounMatch_High].innerTrack().dz(pv_pos)) < 20)):	
                            SV_level_var_array['isSoftMuon_High'][nSV] = 1	
                            hasSoftMuon = 1
                        else: SV_level_var_array['isSoftMuon_High'][nSV] = 0
                    else: SV_level_var_array['isSoftMuon_High'][nSV] = 0
                                            
                else: 
                    SV_level_var_array['muonMatched_High'][nSV] = 0
                    SV_level_var_array['isSoftMuon_High'][nSV] = 0
                    
                SV_level_var_array['eta_Low'][nSV] = trackLow.eta()
                SV_level_var_array['pt_Low'][nSV] = trackLow.pt()
                SV_level_var_array['log10PttrackerrorPttrack_Low'][nSV] =ROOT.TMath.Log10(fabs((trackLow.ptError())/(trackLow.pt())))
                SV_level_var_array['log10dxy_Low'][nSV] = ROOT.TMath.Log10(fabs(trackLow.dxy(pv_pos)))
                SV_level_var_array['log10dz_Low'][nSV] = ROOT.TMath.Log10(fabs(trackLow.dz(pv_pos)))
                SV_level_var_array['log10dxyerrorDxy_Low'][nSV] = ROOT.TMath.Log10(fabs(trackLow.dxyError()/trackLow.dxy()))
                SV_level_var_array['log10dzerrorDz_Low'][nSV] = ROOT.TMath.Log10(fabs(trackLow.dzError()/trackLow.dz()))
                SV_level_var_array['nvalidhits_Low'][nSV] = trackLow.numberOfValidHits()
                SV_level_var_array['absChi2_Low'][nSV] = abs(trackLow.normalizedChi2())
                SV_level_var_array['mvaSingle_Low'][nSV] =filesWithMVA[ifile][ievent][idxLow]
                SV_level_var_array['quality_Low'][nSV] = 10
                for i in range(8):
                    if trackLow.quality(i): SV_level_var_array['quality_Low'][nSV] = i
                SV_level_var_array['trackiso_Low'][nSV] , SV_level_var_array['trackdrmin_Low'][nSV], SV_level_var_array['tracknumneighbours_Low'][nSV] = calcIso_track(trackLow,tracks, pv_pos, False)
                SV_level_var_array['trackisoLoose_Low'][nSV] , SV_level_var_array['trackdrminLoose_Low'][nSV], SV_level_var_array['tracknumneighboursLoose_Low'][nSV] = calcIso_track(trackLow,tracks, pv_pos, True)

                
                SV_level_var_array['eta_High'][nSV] = trackHigh.eta()
                SV_level_var_array['pt_High'][nSV] = trackHigh.pt()
                SV_level_var_array['log10PttrackerrorPttrack_High'][nSV] =ROOT.TMath.Log10(fabs((trackHigh.ptError())/(trackHigh.pt())))
                SV_level_var_array['log10dxy_High'][nSV] = ROOT.TMath.Log10(fabs(trackHigh.dxy(pv_pos)))
                SV_level_var_array['log10dxyerrorDxy_High'][nSV] = ROOT.TMath.Log10(fabs(trackHigh.dxyError()/trackHigh.dxy()))
                SV_level_var_array['log10dzerrorDz_High'][nSV] = ROOT.TMath.Log10(fabs(trackHigh.dzError()/trackHigh.dz()))
                SV_level_var_array['log10dz_High'][nSV] = ROOT.TMath.Log10(fabs(trackHigh.dz(pv_pos)))
                SV_level_var_array['nvalidhits_High'][nSV] = trackHigh.numberOfValidHits()
                SV_level_var_array['absChi2_High'][nSV] = abs(trackHigh.normalizedChi2())
                SV_level_var_array['mvaSingle_High'][nSV] = filesWithMVA[ifile][ievent][idxHigh]
                SV_level_var_array['quality_High'][nSV] = 10

                for i in range(8):
                    if trackHigh.quality(i): SV_level_var_array['quality_High'][nSV] = i
                SV_level_var_array['trackiso_High'][nSV] , SV_level_var_array['trackdrmin_High'][nSV], SV_level_var_array['tracknumneighbours_High'][nSV] = calcIso_track(trackHigh,tracks, pv_pos, False)
                SV_level_var_array['trackisoLoose_High'][nSV] , SV_level_var_array['trackdrminLoose_High'][nSV], SV_level_var_array['tracknumneighboursLoose_High'][nSV] = calcIso_track(trackHigh,tracks, pv_pos, True)
                
                SV_level_var_array['jetrelpt_High'][nSV], SV_level_var_array['jetdrmin_High'][nSV], SV_level_var_array['jetnum_High'][nSV] = calcIso_jet(trackHigh,jets, pv_pos, False)
                SV_level_var_array['jetrelpt_Low'][nSV], SV_level_var_array['jetdrmin_Low'][nSV], SV_level_var_array['jetnum_Low'][nSV]  = calcIso_jet(trackLow,jets, pv_pos, False)

                SV_level_var_array['deltaPhi'][nSV] = abs(TLV1.DeltaPhi(TLV2))
                SV_level_var_array['deltaR'][nSV] = TLV1.DeltaR(TLV2)
                SV_level_var_array['deltaEta'][nSV] = abs(matchingTrk[0].eta()-matchingTrk[1].eta())
                SV_level_var_array['invMass'][nSV] = (TLV1+TLV2).M()
                SV_level_var_array['sumCharge'][nSV] = (trackLow.charge()+trackHigh.charge())
                SV_level_var_array['vectorSumPt'][nSV] = (TLV1+TLV2).Pt()
                SV_level_var_array['vectorSumPxy'][nSV] = sqrt(pow((TLV1+TLV2).Px(),2)+pow((TLV1+TLV2).Py(),2))            

                ######################################
                #### "gen match of SV constituent for background SVs"
                ######################################
                if 'signal' in options.tag:

                    idxGP = [-1,-1]
                    pdgIds = [-1, -1]
                    pdgIdsMother = [-1, -1]
                    drminGP = [-1, -1]
                    tlvMother = [None, None]
                    hasEWancestors = [-1, -1]
                    hasEWancestors_new = [-1, -1]
                    #relatives = [[], []]
                    numDaughtersOfMother = [-1, -1]
                    
                    idxGP_new = [-1,-1]
                    pdgIds_new = [-1, -1]
                    pdgIdsMother_new = [-1, -1]
                    drminGP_new = [-1, -1]
                    dxyzmin = [-1, -1]
                    tlvMother_new = [None, None]
                    #relatives_new = [[], []]
                    numDaughtersOfMother_new = [-1, -1]
                    ignoreIndices = []
                    
                    if not isSignal and  nSV != signalIdx:
                        
                        if secondary.numberOfDaughters()> 2: continue
                        
                        ### first element is always high PT track	
                        ### if 1st is higher in pt put first at first position
                        if secondary.daughter(0).pt() > secondary.daughter(1).pt():
                                    idxGP[0],  pdgIds[0], pdgIdsMother[0], drminGP[0], tlvMother[0], hasEWancestors[0], numDaughtersOfMother[0] = findMinDr_ancestors(secondary.daughter(0),genparticles, ignoreIndices)
                                    if drminGP[0] < 0.01 and idxGP[0] not in ignoreIndices: ignoreIndices.append(idxGP[0])
                                    idxGP_new[0], dxyzmin[0], pdgIds_new[0], pdgIdsMother_new[0], drminGP_new[0], tlvMother_new[0], hasEWancestors_new[0], numDaughtersOfMother_new[0]= findMatch_ancestor_new(secondary.daughter(0),genparticles, ignoreIndices)
                                    if dxyzmin[0] < 0.03 and drminGP_new[0] < 0.01 and idxGP_new[0] not in ignoreIndices: ignoreIndices.append(idxGP_new[0])
                                    
                                    idxGP[1],  pdgIds[1], pdgIdsMother[1], drminGP[1], tlvMother[1], hasEWancestors[1], numDaughtersOfMother[1]  = findMinDr_ancestors(secondary.daughter(1),genparticles,  ignoreIndices)
                                    idxGP_new[1], dxyzmin[1], pdgIds_new[1], pdgIdsMother_new[1], drminGP_new[1], tlvMother_new[1], hasEWancestors_new[1], numDaughtersOfMother_new[1]= findMatch_ancestor_new(secondary.daughter(1),genparticles, ignoreIndices)

                        ### else if 2nd is higher inpt, put 2nd at 1st position
                        else:
                            idxGP[1],  pdgIds[1], pdgIdsMother[1], drminGP[1] , tlvMother[1], hasEWancestors[1], numDaughtersOfMother[1] = findMinDr_ancestors(secondary.daughter(0),genparticles, ignoreIndices)
                            if drminGP[1] < 0.01 and idxGP[1] not in ignoreIndices: ignoreIndices.append(idxGP[1])
                            idxGP_new[1], dxyzmin[1], pdgIds_new[1], pdgIdsMother_new[1], drminGP_new[1] , tlvMother_new[1], hasEWancestors_new[1], numDaughtersOfMother_new[1] = findMatch_ancestor_new(secondary.daughter(0),genparticles, ignoreIndices)
                            if dxyzmin[1] < 0.03 and drminGP_new[1] < 0.01 and idxGP_new[1] not in ignoreIndices: ignoreIndices.append(idxGP_new[1])	
                                                        
                            idxGP[0],  pdgIds[0], pdgIdsMother[0], drminGP[0] , tlvMother[0], hasEWancestors[0], numDaughtersOfMother[0] = findMinDr_ancestors(secondary.daughter(1),genparticles, ignoreIndices)		
                            idxGP_new[0], dxyzmin[0], pdgIds_new[0], pdgIdsMother_new[0], drminGP_new[0] , tlvMother_new[0], hasEWancestors_new[0], numDaughtersOfMother_new[0] = findMatch_ancestor_new(secondary.daughter(1),genparticles, ignoreIndices)		

                            

                    if (drminGP[0] < -1 or (dxyzmin[0] < -1 and drminGP_new[0]<-1)): SV_level_var_array['hasGenMatch_High'][nSV] = -1 #higher Pt is signal
                    elif (drminGP[0] < 0.01 or (dxyzmin[0] < 0.03 and drminGP_new[0]<0.01)): SV_level_var_array['hasGenMatch_High'][nSV] = 1 #higherPt is gen matched (alwazs first value in dR [x,x]
                    else: SV_level_var_array['hasGenMatch_High'][nSV] = 0 #higher pt is not gen matched
                       
                    if (drminGP[1] < -1 or (dxyzmin[1] < -1 and drminGP_new[1]<-1)): SV_level_var_array['hasGenMatch_Low'][nSV] = -1
                    elif (drminGP[1] < 0.01 or (dxyzmin[1] < 0.03 and drminGP_new[1]<0.01)): SV_level_var_array['hasGenMatch_Low'][nSV] = 1
                    else: SV_level_var_array['hasGenMatch_Low'][nSV] = 0
                    
                    SV_level_var_array['numDaughtersOfMother_Low'][nSV] = abs(numDaughtersOfMother[1])
                    SV_level_var_array['numDaughtersOfMother_High'][nSV] = abs(numDaughtersOfMother[0])
                
                    if SV_level_var_array['hasGenMatch_High'][nSV] == 1 and SV_level_var_array['hasGenMatch_Low'][nSV] == 1:#both have matching gen
                        if (pdgIdsMother[0] != -1) and (pdgIdsMother[0] == pdgIdsMother[1]) and (drminGP[0] < 0.01 and drminGP[1] < 0.01): 			
                            if (abs(tlvMother[0].Eta()-tlvMother[1].Eta()) < 0.001):
                                if (abs(tlvMother[0].Pt()-tlvMother[1].Pt()) < 0.001):
                                    SV_level_var_array['events_hasGenMatchWithSameMother1'][nSV] = ievent
                                    SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 1	
                                    SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother[1])
                                    SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother[0])
                                    SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds[1])
                                    SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds[0])

                                else: 
                                    SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 0	
                                    SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother[1])
                                    SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother[0])
                                    SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds[1])
                                    SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds[0])

                            else: 
                                SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 0	
                                SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother[1])
                                SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother[0])						
                                SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds[1])
                                SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds[0])

                        elif (pdgIdsMother_new[0] != -1) and (pdgIdsMother_new[0]== pdgIdsMother_new[1]) and (dxyzmin[0] < 0.03 and drminGP_new[0]<0.01) and (dxyzmin[1] < 0.03 and drminGP_new[1]<0.01): 
                            #if (tlvMother[0].Charge()== tlvMother[1].Charge()) or (tlvMother_new[0].Charge()== tlvMother_new[1].Charge()): 					
                            if (abs(tlvMother_new[0].Eta()-tlvMother_new[1].Eta()) < 0.001):
                                if (abs(tlvMother_new[0].Pt()-tlvMother_new[1].Pt()) < 0.001):
                                    SV_level_var_array['events_hasGenMatchWithSameMother1'][nSV] = ievent
                                    SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 1	
                                    SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother_new[0])
                                    SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother_new[0])
                                    SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds_new[1])
                                    SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds_new[0])
                                    SV_level_var_array['numDaughtersOfMother_Low'][nSV] = abs(numDaughtersOfMother_new[1])

                                else: 
                                    SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 0	
                                    SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother_new[1])
                                    SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother_new[0])
                                    SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds_new[1])
                                    SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds_new[0])

                            else: 
                                SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 0	
                                SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother_new[1])
                                SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother_new[0])
                                SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds_new[1])
                                SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds_new[0])

                        elif (((pdgIdsMother[0] == -1) or (pdgIdsMother[1] == -1)) and (drminGP[0] < 0.01 and drminGP[1] < 0.01)): 
                            SV_level_var_array['events_hasGenMatchWithSameMotherm1'][nSV] = ievent
                            SV_level_var_array['hasGenMatchWithSameMother'][nSV] = -1 #signal SV
                            SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds[1])
                            SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds[0])
                            
                        elif (((pdgIdsMother_new[0] == -1) or (pdgIdsMother_new[1] == -1)) and  ((dxyzmin[0] < 0.03 and drminGP_new[0]<0.01) and (dxyzmin[1] < 0.03 and drminGP_new[1]<0.01))): 
                            SV_level_var_array['events_hasGenMatchWithSameMotherm1'][nSV] = ievent
                            SV_level_var_array['hasGenMatchWithSameMother'][nSV] = -1 #signal SV
                            SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds_new[1])
                            SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds_new[0])
                            #exit

                        else: 
                            SV_level_var_array['hasGenMatchWithSameMother'][nSV] = 0
                            SV_level_var_array['pdgIDMother_Low'][nSV] = abs(pdgIdsMother[1])
                            SV_level_var_array['pdgIDMother_High'][nSV] = abs(pdgIdsMother[0])	
                            SV_level_var_array['pdgID_Low'][nSV] = abs(pdgIds[1])
                            SV_level_var_array['pdgID_High'][nSV] = abs(pdgIds[0])

                    SV_level_var_array['hasEWancestor_Low'][nSV] = hasEWancestors[1]+hasEWancestors_new[1]
                    SV_level_var_array['hasEWancestor_High'][nSV] = hasEWancestors[0]+hasEWancestors_new[0]

                    if not leptons[0] == None and not leptons[1]==None and not theChi01 == None and not theChi02 == None:
                        SV_level_var_array['mtransverse2_hybrid'][nSV] = mZstar_reco_muoncase*mZstar_reco_muoncase + ((v_pZstar_reco.Cross(normalvector))*(v_pZstar_reco.Cross(normalvector)))
                        if isSignal and  nSV == signalIdx:
                            event_level_var_array['res_vx'][0] = secondary.vx() - theChi01.vx()
                            event_level_var_array['res_vy'][0] = secondary.vy() - theChi01.vy()
                            event_level_var_array['res_vz'][0] = secondary.vz() - theChi01.vz()
                            event_level_var_array['res_PVx'][0] = (pv_pos.x() - theChi02.vx())
                            event_level_var_array['res_PVy'][0] = (pv_pos.y() - theChi02.vy())
                            event_level_var_array['res_PVz'][0] = (pv_pos.z() - theChi02.vz())
                            event_level_var_array['res_dx'][0] = (secondary.vx() - pv_pos.x()) - (theChi01.vx() - theChi02.vx())
                            event_level_var_array['res_dy'][0] = (secondary.vy() - pv_pos.x()) - (theChi01.vy() - theChi02.vx())
                            event_level_var_array['res_dz'][0] = (secondary.vz() - pv_pos.x()) - (theChi01.vz() - theChi02.vx())
                            event_level_var_array['res_theta'][0] = (degrees(asin((v_pZstar_reco.Cross(normalvector_reco)).Mag()/(v_pZstar_reco.Mag()))))-(degrees(asin((v_pZstar.Cross(normalvector)).Mag()/(v_pZstar.Mag()))))
                            event_level_var_array['gen_theta'][0] = (degrees(asin((v_pZstar.Cross(normalvector)).Mag()/(v_pZstar.Mag()))))
                            event_level_var_array['reco_theta'][0] = (degrees(asin((v_pZstar_reco.Cross(normalvector_reco)).Mag()/(v_pZstar_reco.Mag()))))
                            event_level_var_array['res_deltaEta'][0] = (summedLeptons.Eta()-PVSV.Eta()) - (summedTracks.Eta()-PVVtx.Eta())
                            event_level_var_array['res_deltaPhi'][0] = deltaPhi(summedLeptons.Phi(), PVSV.Phi()) - deltaPhi(summedTracks.Phi(), PVVtx.Phi())
                            event_level_var_array['res_eta'][0] = - PVVtx.Eta()
                            event_level_var_array['res_phi'][0] = - PVVtx.Phi()
                            event_level_var_array['res_ncrossn'][0] = (normalvector.Cross(normalvector_reco)).Mag()
                            event_level_var_array['res_alphan'][0] = asin((normalvector.Cross(normalvector_reco).Mag()))
                            event_level_var_array['res_mtransverse2'][0] = mtransverse2 - mtransverse2_reco_muoncase
                            event_level_var_array['res_mtransverse'][0] = sqrt(mtransverse2)- sqrt(mtransverse2_reco_muoncase)    



                numsvsfinalpreselection += 1

        event_level_var_array['n_sv'][0] = numsvsfinalpreselection
        event_level_var_array['numSVs'][0] = numSVs
        event_level_var_array['n_sv_total'][0] = n_sv_total

        '''
        ###############################################################################################
        # get track-level info
        ###############################################################################################
        '''

        # TODO: add (ch)PF iso with pT cut

        tracksforiso0 = np.array([(t.pt(), t.eta(), t.phi()) for t in tracks
                                  if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=0.)])
        tracksforiso1 = np.array([(t.pt(), t.eta(), t.phi()) for t in tracks
                                  if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=1.)])
        tracksforiso5 = np.array([(t.pt(), t.eta(), t.phi()) for t in tracks
                                  if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=5.)])
        tracksforiso10 = np.array([(t.pt(), t.eta(), t.phi()) for t in tracks
                                   if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=10.)])

        jetsforiso0 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                if passesPreselection_iso_jet(j, pt_threshold=0.)])
        jetsforiso10 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                if passesPreselection_iso_jet(j, pt_threshold=10.)])
        jetsforiso20 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                if passesPreselection_iso_jet(j, pt_threshold=20.)])
        jetsforiso30 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                if passesPreselection_iso_jet(j, pt_threshold=30.)])

        jetsforisoNoLepton15 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                         if passesPreselection_iso_jet(j, pt_threshold=15.) and ij not in isleptonjetidxlist])

        bjetsforisoLoose15 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                       if passesPreselection_iso_jet(j, pt_threshold=15.) and btaglist[ij] > loosewp])
        bjetsforisoLoose30 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                       if passesPreselection_iso_jet(j, pt_threshold=30.) and btaglist[ij] > loosewp])
        bjetsforisoMedium15 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                       if passesPreselection_iso_jet(j, pt_threshold=15.) and btaglist[ij] > mediumwp])
        bjetsforisoMedium30 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                       if passesPreselection_iso_jet(j, pt_threshold=30.) and btaglist[ij] > mediumwp])
        bjetsforisoTight15 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                       if passesPreselection_iso_jet(j, pt_threshold=15.) and btaglist[ij] > tightwp])
        bjetsforisoTight30 = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                       if passesPreselection_iso_jet(j, pt_threshold=30.) and btaglist[ij] > tightwp])

        # see https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideParticleFlow#Output
        neutralhadrons0 = [p for p in pfcands if p.particleId() == 5]
        neutralhadrons1 = [p for p in pfcands if p.particleId() == 5 and p.pt() > 1.]
        neutralhadrons5 = [p for p in pfcands if p.particleId() == 5 and p.pt() > 5.]
        neutralhadrons10 = [p for p in pfcands if p.particleId() == 5 and p.pt() > 10.]

        if 'genmatchtracks' in options.tag or 'genmatchalltracks' in options.tag:
            genparticlesformatching = [gp for gp in genparticles if gp.status() == 1]

        numtracksbasicpreselection = 0
        numtracksfinalpreselection = 0

        i = 0
        for itrack, track in enumerate(tracks):

            if not passesPreselection_basic_track(track): continue

            numtracksbasicpreselection += 1

            # TODO: adapt preselection
            if not abs(track.dz(pv_pos)) < 1: continue
            jetiso30, jetisomulti30, jetdrmin30, jetisobtag30, jetminv30 = calcIso_jet_new(track, jetsforiso30, isTrack=True, btagvalues=btagvalues)
            if not jetdrmin30 > 0.4: continue

            numtracksfinalpreselection += 1

            track_level_var_array['track_random'][i] = random.randrange(10)

            track_level_var_array['track_charge'][i] = track.charge()

            track_level_var_array['track_px'][i] = track.px()
            track_level_var_array['track_py'][i] = track.py()
            track_level_var_array['track_pz'][i] = track.pz()

            track_level_var_array['track_pt'][i] = track.pt()

            track_level_var_array['track_ptError'][i] = track.ptError()
            track_level_var_array['track_log10(ptError)'][i] = ROOT.TMath.Log10(track.ptError())

            track_level_var_array['track_ptError/pt'][i] = track.ptError()/track.pt()
            track_level_var_array['track_log10(ptError/pt)'][i] = ROOT.TMath.Log10(track.ptError()/track.pt())

            track_level_var_array['track_eta'][i] = track.eta()
            track_level_var_array['track_phi'][i] = track.phi()

            track_level_var_array['track_etaError'][i] = track.etaError()
            track_level_var_array['track_phiError'][i] = track.phiError()

            ispfcand = False
            thepfcand = None
            for pfc in pfcands:
                if pfc.charge() == 0: continue
                if pfc.trackRef().isNull(): continue
                if pfc.trackRef().get() == track:
                    ispfcand = True
                    thepfcand = pfc
                    break

            if ispfcand:
                track_level_var_array['track_isPfCand'][i] = 1
                track_level_var_array['track_pfCandPt'][i] = thepfcand.pt()
                track_level_var_array['track_pfCandEta'][i] = thepfcand.eta()
                track_level_var_array['track_pfCandPhi'][i] = thepfcand.phi()
                track_level_var_array['track_pfCandPdgId'][i] = thepfcand.pdgId()
                track_level_var_array['track_pfCandParticleId'][i] = thepfcand.particleId()
                track_level_var_array['track_pfCandEnergy'][i] = thepfcand.energy()
                track_level_var_array['track_pfCandEcalEnergy'][i] = thepfcand.ecalEnergy()
                track_level_var_array['track_pfCandHcalEnergy'][i] = thepfcand.hcalEnergy()
            else:
                track_level_var_array['track_isPfCand'][i] = 0
                track_level_var_array['track_pfCandPt'][i] = -1
                track_level_var_array['track_pfCandEta'][i] = -1
                track_level_var_array['track_pfCandPhi'][i] = -1
                track_level_var_array['track_pfCandPdgId'][i] = -1
                track_level_var_array['track_pfCandParticleId'][i] = -1
                track_level_var_array['track_pfCandEnergy'][i] = -1
                track_level_var_array['track_pfCandEcalEnergy'][i] = -1
                track_level_var_array['track_pfCandHcalEnergy'][i] = -1


            assPV = -1
            for ipv in range(n_pv):
                if track in tracksByPV[ipv]:
                    assPV = ipv
                    break
            track_level_var_array['track_associatedPV'][i] = assPV
            if assPV > -1:
                assPV_pos = primaryvertices[assPV].position()
            else:  # not consistent but ok
                assPV = 0
                assPV_pos = pv_pos

            track_level_var_array['track_distPVAssPVxy'][i] = ROOT.TMath.Sqrt(pow(assPV_pos.x() - pv_pos.x(), 2)
                                                                              + pow(assPV_pos.y() - pv_pos.y(), 2))
            track_level_var_array['track_distPVAssPVz'][i] = abs(assPV_pos.z() - pv_pos.z())

            ip = IPcalculator(track, primaryvertices[0])
            track_level_var_array['track_IPsig'][i] = ip.getIPsignificance()
            track_level_var_array['track_IPxyz'][i] = ip.getIP()
            track_level_var_array['track_IPxy'][i] = ip.getDxy()
            track_level_var_array['track_IPz'][i] = ip.getDz()
            track_level_var_array['track_log10(IPsig)'][i] = ROOT.TMath.Log10(ip.getIPsignificance())
            track_level_var_array['track_log10(IPxyz)'][i] = ROOT.TMath.Log10(ip.getIP())
            track_level_var_array['track_log10(IPxy)'][i] = ROOT.TMath.Log10(ip.getDxy())
            track_level_var_array['track_log10(IPz)'][i] = ROOT.TMath.Log10(ip.getDz())

            minipPU = None
            minivPU = -1
            minIPsignificancePU = float('inf')
            for iv, v in enumerate(primaryvertices[1:]):
                thisipPU = IPcalculator(track, v)
                thisIPsignificancePU = thisipPU.getIPsignificance()
                if thisIPsignificancePU < minIPsignificancePU:
                    minipPU = thisipPU
                    minivPU = iv
                    minIPsignificancePU = thisIPsignificancePU

            if not minivPU == -1:
                track_level_var_array['track_associatedPU'][i] = minivPU+1
                track_level_var_array['track_IPsigPU'][i] = minipPU.getIPsignificance()
                track_level_var_array['track_IPxyzPU'][i] = minipPU.getIP()
                track_level_var_array['track_IPxyPU'][i] = minipPU.getDxy()
                track_level_var_array['track_IPzPU'][i] = minipPU.getDz()
                track_level_var_array['track_log10(IPsigPU)'][i] = ROOT.TMath.Log10(minipPU.getIPsignificance())
                track_level_var_array['track_log10(IPxyzPU)'][i] = ROOT.TMath.Log10(minipPU.getIP())
                track_level_var_array['track_log10(IPxyPU)'][i] = ROOT.TMath.Log10(minipPU.getDxy())
                track_level_var_array['track_log10(IPzPU)'][i] = ROOT.TMath.Log10(minipPU.getDz())
            else:
                track_level_var_array['track_associatedPU'][i] = -1
                track_level_var_array['track_IPsigPU'][i] = -1
                track_level_var_array['track_IPxyzPU'][i] = -1
                track_level_var_array['track_IPxyPU'][i] = -1
                track_level_var_array['track_IPzPU'][i] = -1
                track_level_var_array['track_log10(IPsigPU)'][i] = -10
                track_level_var_array['track_log10(IPxyzPU)'][i] = -10
                track_level_var_array['track_log10(IPxyPU)'][i] = -10
                track_level_var_array['track_log10(IPzPU)'][i] = -10


            ipAssPV = IPcalculator(track, primaryvertices[assPV])
            track_level_var_array['track_IPsigAssPV'][i] = ipAssPV.getIPsignificance()
            track_level_var_array['track_IPxyzAssPV'][i] = ipAssPV.getIP()
            track_level_var_array['track_IPxyAssPV'][i] = ipAssPV.getDxy()
            track_level_var_array['track_IPzAssPV'][i] = ipAssPV.getDz()
            track_level_var_array['track_log10(IPsigAssPV)'][i] = ROOT.TMath.Log10(ipAssPV.getIPsignificance())
            track_level_var_array['track_log10(IPxyzAssPV)'][i] = ROOT.TMath.Log10(ipAssPV.getIP())
            track_level_var_array['track_log10(IPxyAssPV)'][i] = ROOT.TMath.Log10(ipAssPV.getDxy())
            track_level_var_array['track_log10(IPzAssPV)'][i] = ROOT.TMath.Log10(ipAssPV.getDz())

            minipPUAssPV = None
            minivPUAssPV = -1
            minIPsignificancePUAssPV = float('inf')
            for iv, v in enumerate(primaryvertices[:assPV] + primaryvertices[assPV+1:]):
                thisipPUAssPV = IPcalculator(track, v)
                thisIPsignificancePUAssPV = thisipPUAssPV.getIPsignificance()
                if thisIPsignificancePUAssPV < minIPsignificancePUAssPV:
                    minipPUAssPV = thisipPUAssPV
                    minivPUAssPV = iv
                    minIPsignificancePUAssPV = thisIPsignificancePUAssPV

            if not minivPUAssPV == -1:
                track_level_var_array['track_associatedPUAssPV'][i] = minivPUAssPV if minivPUAssPV < assPV else minivPUAssPV+1
                track_level_var_array['track_IPsigPUAssPV'][i] = minipPUAssPV.getIPsignificance()
                track_level_var_array['track_IPxyzPUAssPV'][i] = minipPUAssPV.getIP()
                track_level_var_array['track_IPxyPUAssPV'][i] = minipPUAssPV.getDxy()
                track_level_var_array['track_IPzPUAssPV'][i] = minipPUAssPV.getDz()
                track_level_var_array['track_log10(IPsigPUAssPV)'][i] = ROOT.TMath.Log10(minipPUAssPV.getIPsignificance())
                track_level_var_array['track_log10(IPxyzPUAssPV)'][i] = ROOT.TMath.Log10(minipPUAssPV.getIP())
                track_level_var_array['track_log10(IPxyPUAssPV)'][i] = ROOT.TMath.Log10(minipPUAssPV.getDxy())
                track_level_var_array['track_log10(IPzPUAssPV)'][i] = ROOT.TMath.Log10(minipPUAssPV.getDz())
            else:
                track_level_var_array['track_associatedPUAssPV'][i] = -1
                track_level_var_array['track_IPsigPUAssPV'][i] = -1
                track_level_var_array['track_IPxyzPUAssPV'][i] = -1
                track_level_var_array['track_IPxyPUAssPV'][i] = -1
                track_level_var_array['track_IPzPUAssPV'][i] = -1
                track_level_var_array['track_log10(IPsigPUAssPV)'][i] = -10
                track_level_var_array['track_log10(IPxyzPUAssPV)'][i] = -10
                track_level_var_array['track_log10(IPxyPUAssPV)'][i] = -10
                track_level_var_array['track_log10(IPzPUAssPV)'][i] = -10

            track_level_var_array['track_dxy0'][i] = abs(track.dxy())
            track_level_var_array['track_dz0'][i] = abs(track.dz())
            track_level_var_array['track_log10(dxy0)'][i] = ROOT.TMath.Log10(abs(track.dxy()))
            track_level_var_array['track_log10(dz0)'][i] = ROOT.TMath.Log10(abs(track.dz()))

            track_level_var_array['track_dxyNoAbs'][i] = track.dxy(pv_pos)
            track_level_var_array['track_dzNoAbs'][i] = track.dz(pv_pos)
            track_level_var_array['track_dxySign'][i] = track.dxy(pv_pos) / abs(track.dxy(pv_pos))
            track_level_var_array['track_dzSign'][i] = track.dz(pv_pos) / abs(track.dz(pv_pos))

            track_level_var_array['track_dxy'][i] = abs(track.dxy(pv_pos))
            track_level_var_array['track_dz'][i] = abs(track.dz(pv_pos))
            track_level_var_array['track_log10(dxy)'][i] = ROOT.TMath.Log10(abs(track.dxy(pv_pos)))
            track_level_var_array['track_log10(dz)'][i] = ROOT.TMath.Log10(abs(track.dz(pv_pos)))

            dxyhandmade, dzhandmade = handmadeDxyDz(track, pv_pos)
            track_level_var_array['track_dxyHandmade'][i] = dxyhandmade
            track_level_var_array['track_dzHandmade'][i] = dzhandmade
            track_level_var_array['track_log10(dxyHandmade)'][i] = ROOT.TMath.Log10(dxyhandmade)
            track_level_var_array['track_log10(dzHandmade)'][i] = ROOT.TMath.Log10(dzhandmade)

            mindxyPU = float('inf')
            mindzPU = float('inf')
            for v in primaryvertices[1:]:
                thisdxyPU = abs(track.dxy(v.position()))
                if thisdxyPU < mindxyPU:
                    mindxyPU = thisdxyPU
                thisdzPU = abs(track.dz(v.position()))
                if thisdzPU < mindzPU:
                    mindzPU = thisdzPU
            track_level_var_array['track_dxyPU'][i] = mindxyPU
            track_level_var_array['track_dzPU'][i] = mindzPU
            track_level_var_array['track_log10(dxyPU)'][i] = ROOT.TMath.Log10(mindxyPU)
            track_level_var_array['track_log10(dzPU)'][i] = ROOT.TMath.Log10(mindzPU)
            
            
            track_level_var_array['track_dxyNoAbsAssPV'][i] = track.dxy(assPV_pos)
            track_level_var_array['track_dzNoAbsAssPV'][i] = track.dz(assPV_pos)
            track_level_var_array['track_dxySignAssPV'][i] = track.dxy(assPV_pos) / abs(track.dxy(assPV_pos))
            track_level_var_array['track_dzSignAssPV'][i] = track.dz(assPV_pos) / abs(track.dz(assPV_pos))

            track_level_var_array['track_dxyAssPV'][i] = abs(track.dxy(assPV_pos))
            track_level_var_array['track_dzAssPV'][i] = abs(track.dz(assPV_pos))
            track_level_var_array['track_log10(dxyAssPV)'][i] = ROOT.TMath.Log10(abs(track.dxy(assPV_pos)))
            track_level_var_array['track_log10(dzAssPV)'][i] = ROOT.TMath.Log10(abs(track.dz(assPV_pos)))

            dxyhandmade, dzhandmade = handmadeDxyDz(track, assPV_pos)
            track_level_var_array['track_dxyHandmadeAssPV'][i] = dxyhandmade
            track_level_var_array['track_dzHandmadeAssPV'][i] = dzhandmade
            track_level_var_array['track_log10(dxyHandmadeAssPV)'][i] = ROOT.TMath.Log10(dxyhandmade)
            track_level_var_array['track_log10(dzHandmadeAssPV)'][i] = ROOT.TMath.Log10(dzhandmade)

            mindxyPUAssPV = float('inf')
            mindzPUAssPV = float('inf')
            for v in primaryvertices[:assPV] + primaryvertices[assPV+1:]:
                thisdxyPUAssPV = abs(track.dxy(v.position()))
                if thisdxyPUAssPV < mindxyPUAssPV:
                    mindxyPUAssPV = thisdxyPUAssPV
                thisdzPUAssPV = abs(track.dz(v.position()))
                if thisdzPUAssPV < mindzPUAssPV:
                    mindzPUAssPV = thisdzPUAssPV
            track_level_var_array['track_dxyPUAssPV'][i] = mindxyPUAssPV
            track_level_var_array['track_dzPUAssPV'][i] = mindzPUAssPV
            track_level_var_array['track_log10(dxyPUAssPV)'][i] = ROOT.TMath.Log10(mindxyPUAssPV)
            track_level_var_array['track_log10(dzPUAssPV)'][i] = ROOT.TMath.Log10(mindzPUAssPV)

            track_level_var_array['track_dxyError'][i] = abs(track.dxyError())
            track_level_var_array['track_dzError'][i] = abs(track.dzError())
            track_level_var_array['track_log10(dxyError)'][i] = ROOT.TMath.Log10(abs(track.dxyError()))
            track_level_var_array['track_log10(dzError)'][i] = ROOT.TMath.Log10(abs(track.dzError()))

            track_level_var_array['track_pfAbsIso'][i], track_level_var_array['track_pfRelIso'][i], track_level_var_array['track_drminPf'][i], track_level_var_array['track_numneighboursPf'][i] = calcIso_pf_or_track_new(track, pfcandsforiso0, subtractObject=ispfcand)
            track_level_var_array['track_chPfAbsIso'][i], track_level_var_array['track_chPfRelIso'][i], track_level_var_array['track_drminChPf'][i], track_level_var_array['track_numneighboursChPf'][i] = calcIso_pf_or_track_new(track, chpfcandsforiso0, subtractObject=(ispfcand and abs(track.dxy(pv_pos)) < 0.1 and abs(track.dz(pv_pos)) < 0.1))

            track_level_var_array['track_tkAbsIso0'][i], track_level_var_array['track_tkRelIso0'][i], track_level_var_array['track_drminTrack0'][i], track_level_var_array['track_numneighboursTrack0'][i] = calcIso_pf_or_track_new(track, tracksforiso0, subtractObject=passesPreselection_iso_track(track, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=0.))
            track_level_var_array['track_tkAbsIso1'][i], track_level_var_array['track_tkRelIso1'][i], track_level_var_array['track_drminTrack1'][i], track_level_var_array['track_numneighboursTrack1'][i] = calcIso_pf_or_track_new(track, tracksforiso1, subtractObject=passesPreselection_iso_track(track, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=1.))
            track_level_var_array['track_tkAbsIso5'][i], track_level_var_array['track_tkRelIso5'][i], track_level_var_array['track_drminTrack5'][i], track_level_var_array['track_numneighboursTrack5'][i] = calcIso_pf_or_track_new(track, tracksforiso5, subtractObject=passesPreselection_iso_track(track, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=5.))
            track_level_var_array['track_tkAbsIso10'][i], track_level_var_array['track_tkRelIso10'][i], track_level_var_array['track_drminTrack10'][i], track_level_var_array['track_numneighboursTrack10'][i] = calcIso_pf_or_track_new(track, tracksforiso10, subtractObject=passesPreselection_iso_track(track, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=10.))

            track_level_var_array['track_jetIso0'][i], track_level_var_array['track_jetIsoMulti0'][i], track_level_var_array['track_drminJet0'][i], track_level_var_array['track_btagJet0'][i], track_level_var_array['track_minvJet0'][i] = calcIso_jet_new(track, jetsforiso0, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_jetIso10'][i], track_level_var_array['track_jetIsoMulti10'][i], track_level_var_array['track_drminJet10'][i], track_level_var_array['track_btagJet10'][i], track_level_var_array['track_minvJet10'][i] = calcIso_jet_new(track, jetsforiso10, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_jetIso15'][i], track_level_var_array['track_jetIsoMulti15'][i], track_level_var_array['track_drminJet15'][i], track_level_var_array['track_btagJet15'][i], track_level_var_array['track_minvJet15'][i] = calcIso_jet_new(track, jetsforiso15, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_jetIso20'][i], track_level_var_array['track_jetIsoMulti20'][i], track_level_var_array['track_drminJet20'][i], track_level_var_array['track_btagJet20'][i], track_level_var_array['track_minvJet20'][i] = calcIso_jet_new(track, jetsforiso20, isTrack=True, btagvalues=btagvalues)

            track_level_var_array['track_jetIso30'][i] = jetiso30
            track_level_var_array['track_jetIsoMulti30'][i] = jetisomulti30
            track_level_var_array['track_drminJet30'][i] = jetdrmin30
            track_level_var_array['track_btagJet30'][i] = jetisobtag30
            track_level_var_array['track_minvJet30'][i] = jetminv30

            track_level_var_array['track_jetIsoNoLepton15'][i], track_level_var_array['track_jetIsoMultiNoLepton15'][i], track_level_var_array['track_drminJetNoLepton15'][i], track_level_var_array['track_btagJetNoLepton15'][i], track_level_var_array['track_minvJetNoLepton15'][i] = calcIso_jet_new(track, jetsforisoNoLepton15, isTrack=True, btagvalues=btagvalues)
            
            track_level_var_array['track_bjetLooseIso15'][i], track_level_var_array['track_bjetLooseIsoMulti15'][i], track_level_var_array['track_drminBjetLoose15'][i], track_level_var_array['track_btagBjetLoose15'][i], track_level_var_array['track_minvBjetLoose15'][i] = calcIso_jet_new(track, bjetsforisoLoose15, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_bjetLooseIso30'][i], track_level_var_array['track_bjetLooseIsoMulti30'][i], track_level_var_array['track_drminBjetLoose30'][i], track_level_var_array['track_btagBjetLoose30'][i], track_level_var_array['track_minvBjetLoose30'][i] = calcIso_jet_new(track, bjetsforisoLoose30, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_bjetMediumIso15'][i], track_level_var_array['track_bjetMediumIsoMulti15'][i], track_level_var_array['track_drminBjetMedium15'][i], track_level_var_array['track_btagBjetMedium15'][i], track_level_var_array['track_minvBjetMedium15'][i] = calcIso_jet_new(track, bjetsforisoMedium15, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_bjetMediumIso30'][i], track_level_var_array['track_bjetMediumIsoMulti30'][i], track_level_var_array['track_drminBjetMedium30'][i], track_level_var_array['track_btagBjetMedium30'][i], track_level_var_array['track_minvBjetMedium30'][i] = calcIso_jet_new(track, bjetsforisoMedium30, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_bjetTightIso15'][i], track_level_var_array['track_bjetTightIsoMulti15'][i], track_level_var_array['track_drminBjetTight15'][i], track_level_var_array['track_btagBjetTight15'][i], track_level_var_array['track_minvBjetTight15'][i] = calcIso_jet_new(track, bjetsforisoTight15, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['track_bjetTightIso30'][i], track_level_var_array['track_bjetTightIsoMulti30'][i], track_level_var_array['track_drminBjetTight30'][i], track_level_var_array['track_btagBjetTight30'][i], track_level_var_array['track_minvBjetTight30'][i] = calcIso_jet_new(track, bjetsforisoTight30, isTrack=True, btagvalues=btagvalues)


            trackTlv = ROOT.TLorentzVector()
            trackTlv.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), 0.1396)
            for suffix, neutralhadroncollection in [('0', neutralhadrons0), ('1', neutralhadrons1), ('5', neutralhadrons5), ('10', neutralhadrons10)]:
                absisoneutralhadron = 0
                drminneutralhadron = 9
                closestneutralhadronTlv = ROOT.TLorentzVector()
                for p in neutralhadroncollection:
                    dr = deltaR(track.eta(), p.eta(), track.phi(), p.phi())
                    if dr < 0.3:
                        absisoneutralhadron += p.energy()
                    if dr < drminneutralhadron:
                        drminneutralhadron = dr
                        closestneutralhadronTlv.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), 0.1350)
                track_level_var_array['track_neHadAbsIso' + suffix][i] = absisoneutralhadron
                track_level_var_array['track_drminNeHad' + suffix][i] = drminneutralhadron
                if drminneutralhadron < 9:
                    track_level_var_array['track_invmNeHad' + suffix][i] = (closestneutralhadronTlv + trackTlv).M()
                else:
                    track_level_var_array['track_invmNeHad' + suffix][i] = 0

            drminphoton = 9
            for p in photons:
                dr = deltaR(track.eta(), p.eta(), track.phi(), p.phi())
                if dr < drminphoton:
                    drminphoton = dr
            track_level_var_array['track_drminPhoton'][i] = drminphoton

            drminelectron = 9
            for e in electrons:
                dr = deltaR(track.eta(), e.eta(), track.phi(), e.phi())
                if dr < drminelectron:
                    drminelectron = dr
            track_level_var_array['track_drminElectron'][i] = drminelectron

            drminmuon = 9
            for m in muons:
                dr = deltaR(track.eta(), m.eta(), track.phi(), m.phi())
                if dr < drminmuon:
                    drminmuon = dr
            track_level_var_array['track_drminMuon'][i] = drminmuon

            for itaucolleciton, taucollection in enumerate([tauswithdiscriminators, [t for t in tauswithdiscriminators if t[0].pt() > 20]]):

                istauleadpfchhadcand = 0
                drmintau = 9
                closesttaumvadiscr = -1
                closesttaudecaymode = -1
                taudr3wp = 0
                taudr4wp = 0
                taudr5wp = 0
                for t in taucollection:

                    try:
                        leadpfchhadcand = t[0].leadPFChargedHadrCand().trackRef().get()
                        leadpfchhadcand.pt()
                    except:
                        leadpfchhadcand = None

                    if leadpfchhadcand is not None:
                        if deltaR(leadpfchhadcand.eta(), track.eta(), leadpfchhadcand.phi(), track.phi()) < 0.001:
                            istauleadpfchhadcand = 1

                    dr = deltaR(track.eta(), t[0].eta(), track.phi(), t[0].phi())
                    if dr < drmintau:
                        drmintau = dr
                        closesttaumvadiscr = t[2]
                        closesttaudecaymode = t[0].decayMode()
                    if dr < 0.5:
                        itauwp = 1
                        while itauwp < 7 and t[2+itauwp] > 0.5:
                            if taudr5wp < itauwp: taudr5wp = itauwp
                            itauwp += 1
                        if dr < 0.4:
                            itauwp = 1
                            while itauwp < 7 and t[2+itauwp] > 0.5:
                                if taudr4wp < itauwp: taudr4wp = itauwp
                                itauwp += 1
                            if dr < 0.3:
                                itauwp = 1
                                while itauwp < 7 and t[2+itauwp] > 0.5:
                                    if taudr3wp < itauwp: taudr3wp = itauwp
                                    itauwp += 1

                suffix = '0'
                if itaucolleciton == 1: suffix = '20'

                track_level_var_array['track_isTauLeadPfChHadCand' + suffix][i] = istauleadpfchhadcand
                track_level_var_array['track_drminTau' + suffix][i] = drmintau
                track_level_var_array['track_mvaDiscrTau' + suffix][i] = closesttaumvadiscr
                track_level_var_array['track_decayModeTau' + suffix][i] = closesttaudecaymode
                track_level_var_array['track_dr3highestWpTau' + suffix][i] = taudr3wp
                track_level_var_array['track_dr4highestWpTau' + suffix][i] = taudr4wp
                track_level_var_array['track_dr5highestWpTau' + suffix][i] = taudr5wp

            track_level_var_array['track_detaLeadingJet'][i] = abs(track.eta() - jets[idxhighestptjet].eta())
            track_level_var_array['track_dphiLeadingJet'][i] = deltaPhi(track.phi(), jets[idxhighestptjet].phi())

            track_level_var_array['track_dphiMet'][i] = deltaPhi(track.phi(), met.phi())
            track_level_var_array['track_dphiMetPca'][i], _, _ = handmadeDphiMetPCA(track, pv_pos, met)

            track_level_var_array['track_chi2'][i] = track.normalizedChi2()

            quality = 0
            if track.quality(track.qualityByName('loose')): quality = 1
            if track.quality(track.qualityByName('tight')): quality = 2
            if track.quality(track.qualityByName('highPurity')): quality = 3
            track_level_var_array['track_quality'][i] = quality

            track_level_var_array['track_numValidHits'][i] = track.numberOfValidHits()
            track_level_var_array['track_numLostHits'][i] = track.numberOfLostHits()


            hasGenMatch = -1
            tmingen = -1
            genmatchpdgid = -1
            genmatchmotherpdgid = -1
            genmatchpt = -1
            genmatchmotherpt = -1
            genmatchstatus = -1
            genmatchmotherstatus = -1
            genmatchishardprocess = -1
            genmatchmotherishardprocess = -1
            genmatchisfromhardprocess = -1
            genmatchisprompt = -1
            genmatchisdirecthadrondecayproduct = -1
            genmatchisdirecttaudecayproduct = -1
            genmatchmotheristhetau = 0

            gentaujetmatchdrmin = 9
            gentaujetmatchpt = -1

            dothegenmatch = False

            if 'genmatchtracks' in options.tag:
                if ievent % 10 == 0 and met.pt() > metthresholdtrackgenmatch: dothegenmatch = True

            if 'genmatchalltracks' in options.tag:
                if met.pt() > metthresholdtrackgenmatch: dothegenmatch = True

            if dothegenmatch:

                idxgentaujet, gentaujetmatchdrmin = findMatch_gen_old_easy(track, gentaujets)
                if gentaujetmatchdrmin < 9:
                    gentaujetmatchpt = gentaujets[idxgentaujet].pt()

                idxgen, dxyzmingen, tmingen, drmingen = findMatch_gen_new(track, genparticlesformatching)

                hasGenMatch = 0

                if drmingen < matchingDrThreshold and dxyzmingen < matchingDxyzThreshold:

                    hasGenMatch = 1

                    genmatch = genparticlesformatching[idxgen]

                else:

                    idxgenold, drmingenold = findMatch_gen_old_easy(track, genparticlesformatching)

                    if drmingenold < 0.02:

                        hasGenMatch = 2

                        genmatch = genparticlesformatching[idxgenold]

                if hasGenMatch > 0:

                    mm = 0
                    while not genmatch.statusFlags().isFirstCopy():
                        for idx in range(genmatch.numberOfMothers()):
                            if genmatch.pdgId() == genmatch.mother(idx).pdgId():
                                genmatch = genmatch.mother(idx)
                                break
                        mm += 1
                        if mm > 100: break

                    genmatchmother = genmatch.mother(0)

                    genmatchpdgid = genmatch.pdgId()
                    genmatchmotherpdgid = genmatchmother.pdgId()

                    genmatchpt = genmatch.pt()
                    genmatchmotherpt = genmatchmother.pt()

                    genmatchstatus = genmatch.status()
                    genmatchmotherstatus = genmatchmother.status()

                    genmatchishardprocess = genmatch.statusFlags().fromHardProcess()
                    genmatchmotherishardprocess = genmatchmother.statusFlags().fromHardProcess()

                    genmatchisfromhardprocess = 0
                    while genmatchmother.numberOfMothers() > 0:
                        if genmatchmother.statusFlags().fromHardProcess():
                            genmatchisfromhardprocess = 1
                            break
                        genmatchmother = genmatchmother.mother(0)

                    genmatchisprompt = genmatch.statusFlags().isPrompt()
                    genmatchisdirecthadrondecayproduct = genmatch.statusFlags().isDirectHadronDecayProduct()
                    genmatchisdirecttaudecayproduct = genmatch.statusFlags().isDirectTauDecayProduct()

                    if genmatchmother == thetau: genmatchmotheristhetau = 1

            track_level_var_array['track_hasGenMatch'][i] = hasGenMatch
            track_level_var_array['track_genMatchTmin'][i] = tmingen
            track_level_var_array['track_genMatchPdgId'][i] = genmatchpdgid
            track_level_var_array['track_genMatchMotherPdgId'][i] = genmatchmotherpdgid
            track_level_var_array['track_genMatchPt'][i] = genmatchpt
            track_level_var_array['track_genMatchMotherPt'][i] = genmatchmotherpt
            track_level_var_array['track_genMatchStatus'][i] = genmatchstatus
            track_level_var_array['track_genMatchMotherStatus'][i] = genmatchmotherstatus
            track_level_var_array['track_genMatchIsHardProcess'][i] = genmatchishardprocess
            track_level_var_array['track_genMatchMotherIsHardProcess'][i] = genmatchmotherishardprocess
            track_level_var_array['track_genMatchIsFromHardProcess'][i] = genmatchisfromhardprocess
            track_level_var_array['track_genMatchIsPrompt'][i] = genmatchisprompt
            track_level_var_array['track_genMatchIsDirectHadronDecayProduct'][i] = genmatchisdirecthadrondecayproduct
            track_level_var_array['track_genMatchIsDirectTauDecayProduct'][i] = genmatchisdirecttaudecayproduct
            track_level_var_array['track_genMatchMotherIsTheTau'][i] = genmatchmotheristhetau
            track_level_var_array['track_genMatchMotherTauDecay'][i] = decayWtau

            track_level_var_array['track_drminGenTauJet'][i] = gentaujetmatchdrmin
            track_level_var_array['track_genTauJetPt'][i] = gentaujetmatchpt

            issignaltrack = 0
            if itrack == matchedTrackIdxCharginoPion1 or itrack == matchedTrackIdxCharginoPion2: issignaltrack = 1
            track_level_var_array['track_isSignalTrack'][i] = issignaltrack

            issusytrack = 0
            if itrack in susytracks: issusytrack = 1
            track_level_var_array['track_isSusyTrack'][i] = issusytrack

            susytrackmother = 0
            susytrackpdgid = 0
            if issusytrack:
                susytrackmother = susytracks[itrack][0]
                susytrackpdgid = susytracks[itrack][1]
            track_level_var_array['track_susyTrackPdgIdMother'][i] = susytrackmother
            track_level_var_array['track_susyTrackPdgId'][i] = susytrackpdgid

            i += 1

        event_level_var_array['n_track_total'][0] = len(tracks)
        event_level_var_array['n_track_basic'][0] = numtracksbasicpreselection
        event_level_var_array['n_track'][0] = numtracksfinalpreselection

        event_level_var_array['cutflow'][0] = cutflow

        tEvent.Fill()

        nEventsPerFile += 1

'''
###############################################################################################
# write histos and tree to file
###############################################################################################
'''

if saveOutputFile:

    fout.cd()

    # write histos

    hCutflow.Write()

    hMetptRaw.Write()
    hMetptBeforeLeptonCleaning.Write()
    hMetpt.Write()
    hMindphimetjets.Write()
    hMtmetleadingjet.Write()
    hNjetsbtagmedium.Write()
    hNumphotons.Write()
    hNumleptons.Write()
    hNumtaus.Write()

    hNumjets.Write()
    hNumjets30.Write()
    hNumjets50.Write()
    hNumjets100.Write()
    hNumjets200.Write()

    hBtagjets.Write()

    hZllLeptonPt.Write()
    hZllDrTrack.Write()
    hZllDrPfc.Write()
    hZllDrJet.Write()

    hNPVsPerEvent.Write()
    hPV0x.Write()
    hPV0y.Write()
    hPV0z.Write()
    hPVsx.Write()
    hPVsy.Write()
    hPVsz.Write()

    for v in jetIDvars:
        jetIDhistos[v[0] + 'all'].Write()
        jetIDhistos[v[0] + 'pass'].Write()

    fout.Write('', ROOT.TObject.kWriteDelete)

    print 'just created ' + fout.GetName()

    '''
    ###############################################################################################
    # write json file
    ###############################################################################################
    '''

    if 'data' in options.tag:

        print 'runs'
        print runs

        if len(runs) > 0:
            runs_compacted = {}
            for run in runs:
                if run not in runs_compacted:
                    runs_compacted[run] = []
                for lumisec in runs[run]:
                    if len(runs_compacted[run]) > 0 and lumisec == runs_compacted[run][-1][-1]+1:
                        runs_compacted[run][-1][-1] = lumisec
                    else:
                        runs_compacted[run].append([lumisec, lumisec])

            print 'runs_compacted'
            print runs_compacted

            json_content = json.dumps(runs_compacted)
            with open(fout.GetName().replace('.root', '.json'), 'w') as fo:
                fo.write(json_content)

            print 'just created ' + fout.GetName().replace('.root', '.json')

        else:

            with open(fout.GetName().replace('.root', '.json'), 'w') as fo:
                fo.write(' ')

            print 'just created empty ' + fout.GetName().replace('.root', '.json')

    fout.Close()
