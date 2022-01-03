#! /usr/bin/env python2


"""

Runs over AOD files and writes file with histos and tree.

----------------------------------------------------------------------
python plantTrees.py inputFiles="file1, file2,..." tag="tag1 tag2 ..."
----------------------------------------------------------------------

tags:
-----

local -> don't open file via xrootd
redirinfn -> use infn redirector for xrootd
redirfnal -> use fnal redirector for xrootd
data -> check trigger flags and runnum/lumisec, save json file
geninfoZ -> save GEN info for DY process
geninfoW -> save GEN info for W boson production process
cleanleptons -> perform DY cleaning
signal -> signal GEN variables and FastSim correction for MET/JEC
pmssm -> save pMSSM IDs
noleptonveto -> don't veto events with leptons
nodphimetjetsveto -> don't veto events with dphi(MET, jets) < 0.5
genmatchtracks -> try to find GEN match to every 10. track
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

from DataFormats.FWLite import Events, Lumis, Handle
from FWCore.ParameterSet.VarParsing import VarParsing

from commons import *


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


nMaxEventsPerFile = 100000
nMaxTracksPerEvent = 10000

# TODO: check if saveoutputfile and if test
saveOutputFile = True
isTest = False
neventsTest = 10  # number of events to run over in case of test
printevery = 100  # print event number for every Xth event

# TODO: check thresholds for "new" matching
matchingDrThreshold = 0.05
matchingDxyzThreshold = 0.2

# TODO: check met threshold for event selection
metthreshold = 200
metthresholdtrackgenmatch = 200  # only relevant if "genmatch(all)tracks" in options.tag


'''
###############################################################################################
# get input arguments
###############################################################################################
'''

options = VarParsing('python')
options.parseArguments()


'''
###############################################################################################
# define output file with histos and variables for tree
###############################################################################################
'''

if True:

    if saveOutputFile:
        nameout = 'output'
        if isTest: nameout += '_test'
        fout = ROOT.TFile(nameout + '.root', 'recreate')

        fout.cd()

    tEvent = ROOT.TTree('tEvent', 'tEvent')


    event_level_var_names = []

    if 'pmssm' in options.tag:
        event_level_var_names += [('pMSSMid1', 'F'), ('pMSSMid2', 'F')]

    var_names_gen_signal = [
        ('deltamFILE', 'F'), ('chipmmFILE', 'F')
        ,('numC1', 'I'), ('numN2', 'I'), ('numN1', 'I')
        ,('hasChargino', 'I'), ('hasPion', 'I'), ('hasMatchedTrack', 'I')

        ,('deltamGEN', 'F'), ('chipmmGEN', 'F')
        ,('chipmptGEN', 'F'), ('chipmetaGEN', 'F'), ('chipmphiGEN', 'F')
        ,('chidecaylengthXY', 'F'), ('chidecaylengthZ', 'F'),('chidecaylength3D', 'F')
        ,('log10(chidecaylengthXY)', 'F'), ('log10(chidecaylengthZ)', 'F'), ('log10(chidecaylength3D)', 'F')
        ,('chipmnumdaughters', 'I')

        ,('deltamN2GEN', 'F'), ('chiN2mGEN', 'F')
        ,('chiN2ptGEN', 'F'), ('chiN2etaGEN', 'F'), ('chiN2phiGEN', 'F')
        ,('chidecaylengthXYN2', 'F'), ('chidecaylengthZN2', 'F'), ('chidecaylength3DN2', 'F')
        ,('log10(chidecaylengthXYN2)', 'F'), ('log10(chidecaylengthZN2)', 'F'), ('log10(chidecaylength3DN2)', 'F')
        ,('chiN2numdaughters', 'I')

        ,('pionptGEN', 'F'), ('pionetaGEN', 'F'), ('pionphiGEN', 'F'), ('pionchargeGEN', 'F')
        ,('tminmatching', 'F'), ('dxyzmin', 'F'), ('drmin', 'F')
        ,('dxyzminrandom', 'F'), ('drminrandom', 'F')
        ,('drminold', 'F'), ('drminoldrandom', 'F')
        ,('matchedTrackIdxCharginoPion1', 'F'), ('matchedTrackIdxCharginoPion2', 'F')

        ,('numchidaughters', 'I')
        ]
    event_level_var_names += var_names_gen_signal

    var_names_gen_background = [
        ('genmetpt', 'F'), ('genmetphi', 'F')
        ,('genht', 'F'), ('genhtmiss', 'F')
        ,('pTneutrinosum', 'F')

        ,('numW', 'I')
        ,('ptW', 'F'), ('etaW', 'F'), ('phiW', 'F')
        ,('numWDaughters', 'I'), ('ptWneutrino', 'F'), ('decayWtau', 'F')

        ,('numZgamma', 'I')
        ,('ptZgamma', 'F'), ('etaZgamma', 'F'), ('phiZgamma', 'F')
        ,('numZgammaDaughters', 'I'), ('ptsumZgammaNeutrinos', 'F')
        ]
    event_level_var_names += var_names_gen_background

    var_names_cleaning = [
        ('electronsCleaned', 'I'), ('muonsCleaned', 'I')
        ,('invmCleaning', 'F'), ('zptCleaning', 'F')
        ,('l1ptCleaning', 'F'), ('l2ptCleaning', 'F')
        ,('l1etaCleaning', 'F'), ('l2etaCleaning', 'F')
        ,('l1phiCleaning', 'F'), ('l2phiCleaning', 'F')
        ,('l1absisodbetaCleaning', 'F'), ('l2absisodbetaCleaning', 'F')
        ,('l1relisodbetaCleaning', 'F'), ('l2relisodbetaCleaning', 'F')
        ,('metptBeforeCleaning', 'F'), ('metphiBeforeCleaning', 'F')
        ]
    event_level_var_names += var_names_cleaning

    var_names_event = [
        ('cutflow', 'I'), ('randomevent', 'I')
        ,('numpvs', 'I'), ('rho', 'F')
        ,('metpt', 'F'), ('metphi', 'F')
        ,('nofastsimcorrmetpt', 'F'), ('nofastsimcorrmetphi', 'F')
        ,('ht', 'F'), ('ht5', 'F'), ('htmiss', 'F')
        ,('numgenjets', 'I')
        ,('numbadjets', 'I'), ('minetaabsbadjets', 'F'), ('numbadjetsEventVeto', 'I')
        ,('numbadjetsLepVeto', 'I'), ('minetaabsbadjetsLepVeto', 'F'), ('numbadjetsLepVetoEventVeto', 'I')
        ,('ptleadingjet', 'F'), ('etaleadingjet', 'F'), ('phileadingjet', 'F')
        ,('numjets', 'I'), ('numjets30', 'I'), ('numjets50','I'), ('numjets100', 'I'), ('numjets200', 'I')
        ,('njetsbtagloose', 'I'), ('njetsbtaglooseTIGHT', 'I')
        ,('njetsbtagmedium', 'I'), ('njetsbtagmediumTIGHT', 'I')
        ,('njetsbtagtight', 'I'), ('njetsbtagtightTIGHT', 'I')
        ,('mtmetleadingjet', 'F')
        ,('mindphimetjets', 'F')
        ,('numphotons', 'I'), ('numphotonsiso', 'I')
        ,('numpfleptons', 'I'), ('numpfleptonsiso', 'I')
        ,('numelectrons', 'I'), ('numelectronsiso', 'I')
        ,('nummuons', 'I'), ('nummuonsiso', 'I')
        ,('numleptons', 'I'), ('numleptonsiso', 'I')
        ,('numtaus', 'I'), ('numtausvloose', 'I'), ('numtausloose', 'I'), ('numtausmedium', 'I'), ('numtaustight', 'I'), ('numtausvtight', 'I'), ('numtausvvtight', 'I')
        ,('numtrackstotal', 'I'), ('numtracksbasicpreselection', 'I'), ('numtracksfinalpreselection', 'I')
        ]
    event_level_var_names += var_names_event

    var_names_data = [
        ('runnum', 'F'), ('lumisec', 'F'), ('eventnum', 'F')
        ]
    event_level_var_names += var_names_data

    event_level_var_array = {}
    for n in event_level_var_names:
        event_level_var_array[n[0]] = array(n[1].lower(), [0])
        tEvent.Branch(nice_string(n[0]), event_level_var_array[n[0]], nice_string(n[0]) + '/' + n[1])


    # TODO: implement the correct MET filters for data:
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2
    if 'era16' in options.tag:
        pass
    elif 'era17' in options.tag:
        pass
    elif 'era18' in options.tag:
        pass
    else:
        raise NotImplementedError('MET filters: era unknown or not specified')

    trigger_flags = [
        'globalSuperTightHalo2016Filter'
        ,'HBHENoiseFilter'
        ,'HBHEIsoNoiseFilter'
        ,'eeBadScFilter'
        ,'EcalDeadCellTriggerPrimitiveFilter'
        ,'ecalBadCalibFilter'
        ,'BadChargedHadronFilter'
        ,'BadPFMuonFilter'
        ,'globalTightHalo2016Filter'
        ,'PrimaryVertexFilter'
        ,'CSCTightHaloFilter'
        ]

    for tf in trigger_flags:
        event_level_var_array[tf] = array('i', [0])
        tEvent.Branch(tf, event_level_var_array[tf], tf + '/I')


    trigger_hlt = [
        'HLT_PFMET90_PFMHT90_IDTight_v'
        , 'HLT_PFMET100_PFMHT100_IDTight_BeamHaloCleaned_v'
        , 'HLT_PFMET100_PFMHT100_IDTight_v'
        , 'HLT_PFMET110_PFMHT110_IDTight_v'
        , 'HLT_PFMET120_PFMHT120_IDTight_v'
        , 'triggerfired'
        ]

    for t_hlt in trigger_hlt:
        event_level_var_array[t_hlt] = array('i', [0])
        tEvent.Branch(t_hlt, event_level_var_array[t_hlt], t_hlt + '/I')


    var_names_chidaughter = [
        ('motherchidaughter', 'I'), ('pdgIdchidaughter', 'F')
        ,('ptchidaughter', 'F'), ('etachidaughter', 'F'), ('phichidaughter', 'F')
        ,('hasmatchedtrackchidaughter', 'I')
        ]

    chidaughter_var_array = {}
    for n in var_names_chidaughter:
        chidaughter_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), chidaughter_var_array[n[0]], nice_string(n[0]) + '[numchidaughters]/F')


    var_names_zdaughter = [
        ('pdgIdZdaughter', 'F')
        ,('ptZdaughter', 'F'), ('etaZdaughter', 'F'), ('phiZdaughter', 'F')
        ]

    zdaughter_var_array = {}
    for n in var_names_zdaughter:
        zdaughter_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), zdaughter_var_array[n[0]], nice_string(n[0]) + '[numZgammaDaughters]/F')


    var_names_wdaughter = [
        ('pdgIdWdaughter', 'F')
        ,('ptWdaughter', 'F'), ('etaWdaughter', 'F'), ('phiWdaughter', 'F')
        ]

    wdaughter_var_array = {}
    for n in var_names_wdaughter:
        wdaughter_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), wdaughter_var_array[n[0]], nice_string(n[0]) + '[numWDaughters]/F')

    var_names_pv = [
        ('idxpv', 'I'), ('numtrackspv', 'F')
        ,('xpv', 'F'), ('ypv', 'F'), ('zpv', 'F')
        ]

    pv_var_array = {}
    for n in var_names_pv:
        pv_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), pv_var_array[n[0]], nice_string(n[0]) + '[numpvs]/F')


    var_names_genjet = [
        ('ptgenjet', 'F'), ('etagenjet', 'F'), ('phigenjet', 'F'), ('massgenjet', 'F')
        ]

    genjet_var_array = {}
    for n in var_names_genjet:
        genjet_var_array[n[0]] = array('f', 1000*[0.])
        tEvent.Branch(nice_string(n[0]), genjet_var_array[n[0]], nice_string(n[0]) + '[numgenjets]/F')


    var_names_jet = [
        ('pxjet', 'F'), ('pyjet', 'F'), ('pzjet', 'F')
        ,('ptjet', 'F'), ('energyjet', 'F')
        ,('etajet', 'F'), ('phijet', 'F')
        ,('nconstituentsjet', 'I')
        ,('btagjet', 'F')
        ,('drminleptonjet', 'F'), ('ptclosestleptonjet', 'F'), ('isleptonjet', 'F')
        ,('drmingenjetjet', 'F'), ('ptclosestgenjetjet', 'F'), ('isgenjetjet', 'F')
        ]

    jet_var_array = {}
    for n in var_names_jet:
        jet_var_array[n[0]] = array('f', 1000*[0.])
        tEvent.Branch(nice_string(n[0]), jet_var_array[n[0]], nice_string(n[0]) + '[numjets]/F')


    var_names_photon = [
        ('pxphoton', 'F'), ('pyphoton', 'F'), ('pzphoton', 'F')
        ,('ptphoton', 'F'), ('energyphoton', 'F')
        ,('etaphoton', 'F'), ('phiphoton', 'F')
        ,('pfabsisophoton', 'F'), ('pfabsisominiphoton', 'F')
        ,('chpfabsisophoton', 'F'), ('chpfabsisominiphoton', 'F')
        ,('jetisophoton', 'F'), ('jetisomultiphoton', 'F'), ('jetdrminphoton', 'F')
        ,('chhadisophoton', 'F'), ('neuhadisophoton', 'F'), ('photisophoton', 'F')
        ,('absisophoton', 'F'), ('relisophoton', 'F')
        ]

    photon_var_array = {}
    for n in var_names_photon:
        photon_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), photon_var_array[n[0]], nice_string(n[0]) + '[numphotons]/F')


    var_names_electron = [
        ('chargeelectron', 'I')
        ,('pxelectron', 'F'), ('pyelectron', 'F'), ('pzelectron', 'F')
        ,('ptelectron', 'F'), ('energyelectron', 'F')
        ,('etaelectron', 'F'), ('phielectron', 'F')
        ,('dzelectron', 'F'), ('dxyelectron', 'F')
        ,('pfabsisoelectron', 'F'), ('pfabsisominielectron', 'F')
        ,('chpfabsisoelectron', 'F'), ('chpfabsisominielectron', 'F')
        ,('jetisoelectron', 'F'), ('jetisomultielectron', 'F'), ('jetdrminelectron', 'F')
        ,('chhadisoelectron', 'F'), ('challisoelectron', 'F')
        ,('neuhadisoelectron', 'F'), ('photisoelectron', 'F')
        ,('puchhadisoelectron', 'F')
        ,('absisodbetaelectron', 'F'), ('relisodbetaelectron', 'F')
        ]

    electron_var_array = {}
    for n in var_names_electron:
        electron_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), electron_var_array[n[0]], nice_string(n[0]) + '[numelectrons]/F')


    var_names_muon = [
        ('chargemuon', 'I')
        ,('pxmuon', 'F'), ('pymuon', 'F'), ('pzmuon', 'F')
        ,('ptmuon', 'F'), ('energymuon', 'F')
        ,('etamuon', 'F'), ('phimuon', 'F')
        ,('dzmuon', 'F'), ('dxymuon', 'F')
        ,('pfabsisomuon', 'F'), ('pfabsisominimuon', 'F')
        ,('chpfabsisomuon', 'F'), ('chpfabsisominimuon', 'F')
        ,('jetisomuon', 'F'), ('jetisomultimuon', 'F'), ('jetdrminmuon', 'F')
        ,('chhadisomuon', 'F'), ('challisomuon', 'F')
        ,('neuhadisomuon', 'F'), ('photisomuon', 'F')
        ,('puchhadisomuon', 'F')
        ,('absisodbetamuon', 'F'), ('relisodbetamuon', 'F')
        ]

    muon_var_array = {}
    for n in var_names_muon:
        muon_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), muon_var_array[n[0]], nice_string(n[0]) + '[nummuons]/F')


    var_names_tau = [
        ('chargetau', 'I')
        ,('pxtau', 'F'), ('pytau', 'F'), ('pztau', 'F')
        ,('pttau', 'F'), ('energytau', 'F')
        ,('etatau', 'F'), ('phitau', 'F')
        ,('dztau', 'F'), ('dxytau', 'F')
        ,('chhadisotau', 'F'), ('photisotau', 'F')
        ,('decaymodetau', 'F'), ('decaymodefindingtau', 'F'), ('mvadiscrtau', 'F')
        ,('isvloosetau', 'I'), ('isloosetau', 'I'), ('ismediumtau', 'I'), ('istighttau', 'I'), ('isvtighttau', 'I'), ('isvvtighttau', 'I')
        ,('tauleadpfchhadcandpt', 'F'), ('tauleadpfchhadcandeta', 'F'), ('tauleadpfchhadcandphi', 'F')
        ,('genmatchpdgidtau', 'F'), ('genmatchdrtau', 'F')
        ]

    tau_var_array = {}
    for n in var_names_tau:
        tau_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), tau_var_array[n[0]], nice_string(n[0]) + '[numtaus]/F')


    var_names_pflepton = [
        ('chargepflepton', 'I')
        ,('pxpflepton', 'F'), ('pypflepton', 'F'), ('pzpflepton', 'F')
        ,('ptpflepton', 'F'), ('energypflepton', 'F')
        ,('etapflepton', 'F'), ('phipflepton', 'F')
        ,('dzpflepton', 'F'), ('dxypflepton', 'F')
        ,('pfrelisopflepton', 'F'), ('pfrelisominipflepton', 'F')
        ,('chpfrelisopflepton', 'F'), ('chpfrelisominipflepton', 'F')
        ,('jetisopflepton', 'F'), ('jetisomultipflepton', 'F'), ('jetdrminpflepton', 'F')
        ,('pdgidpflepton', 'F')
        ]

    pflepton_var_array = {}
    for n in var_names_pflepton:
        pflepton_var_array[n[0]] = array('f', 100*[0.])
        tEvent.Branch(nice_string(n[0]), pflepton_var_array[n[0]], nice_string(n[0]) + '[numpfleptons]/F')


    track_level_var_names = [
        ('randomtrack', 'I')

        ,('charge', 'F')
        ,('pxtrack', 'F'), ('pytrack', 'F'), ('pztrack', 'F')
        ,('pttrack', 'F')
        ,('pttrackerror', 'F'), ('log10(pttrackerror)', 'F')
        ,('pttrackerror/pttrack', 'F'), ('log10(pttrackerror/pttrack)', 'F')
        ,('eta', 'F'), ('phi', 'F')
        ,('etaerror', 'F'), ('phierror', 'F')

        ,('associatedpv', 'I')

        ,('idxpvPU', 'I')

        ,('IPsignificance', 'F'), ('IPxyz', 'F'), ('IPxy', 'F'), ('IPz', 'F')
        ,('log10(IPsignificance)', 'F'), ('log10(IPxyz)', 'F'), ('log10(IPxy)', 'F'), ('log10(IPz)', 'F')

        ,('IPsignificancePU', 'F'), ('IPxyzPU', 'F'), ('IPxyPU', 'F'), ('IPzPU', 'F')
        ,('log10(IPsignificancePU)', 'F'), ('log10(IPxyzPU)', 'F'), ('log10(IPxyPU)', 'F'), ('log10(IPzPU)', 'F')

        ,('dxy0', 'F'), ('dz0', 'F')
        ,('log10(dxy0)', 'F'), ('log10(dz0)', 'F')

        ,('dxynoabs', 'F')
        ,('dznoabs', 'F')

        ,('dxy', 'F'), ('dxyhandmade', 'F'), ('dxyclosestpv', 'F'), ('dxyclosestpvPU', 'F')
        ,('dz', 'F'), ('dzhandmade', 'F'), ('dzclosestpv', 'F'), ('dzclosestpvPU', 'F')
        ,('log10(dxy)', 'F'), ('log10(dxyhandmade)', 'F'), ('log10(dxyclosestpv)', 'F'), ('log10(dxyclosestpvPU)', 'F')
        ,('log10(dz)', 'F'), ('log10(dzhandmade)', 'F'), ('log10(dzclosestpv)', 'F'), ('log10(dzclosestpvPU)', 'F')

        ,('dxyerror', 'F'), ('dzerror', 'F')
        ,('log10(dxyerror)', 'F'), ('log10(dzerror)', 'F')

        ,('trackabsiso', 'F'), ('trackreliso', 'F'), ('trackdrmin', 'F'), ('tracknumneighbours', 'I')
        ,('trackabsisotight', 'F'), ('trackrelisotight', 'F'), ('trackdrmintight', 'F'), ('tracknumneighbourstight', 'I')
        ,('pfabsiso', 'F'), ('pfreliso', 'F'), ('pfdrmin', 'F'), ('pfnumneighbours', 'I')
        ,('chpfabsiso', 'F'), ('chpfreliso', 'F'), ('chpfdrmin', 'F'), ('chpfnumneighbours', 'I')
        ,('jetiso', 'F'), ('jetisomulti', 'F'), ('jetdrmin', 'F'), ('jetisobtag', 'F')
        ,('jetisotightNoLepton', 'F'), ('jetisomultitightNoLepton', 'F'), ('jetdrmintightNoLepton', 'F'), ('jetisobtagtightNoLepton', 'F')
        ,('jetisotight', 'F'), ('jetisomultitight', 'F'), ('jetdrmintight', 'F'), ('jetisobtagtight', 'F')
        ,('jetisomeditight', 'F'), ('jetisomultimeditight', 'F'), ('jetdrminmeditight', 'F'), ('jetisobtagmeditight', 'F')
        ,('jetisomedium', 'F'), ('jetisomultimedium', 'F'), ('jetdrminmedium', 'F'), ('jetisobtagmedium', 'F')
        ,('jetisoloose', 'F'), ('jetisomultiloose', 'F'), ('jetdrminloose', 'F'), ('jetisobtagloose', 'F')

        ,('drminneutralhadron', 'F'), ('invmclosestneutralhadrontrack', 'F')
        ,('drminphoton', 'F'), ('drminelectron', 'F'), ('drminmuon', 'F')
        ,('istauleadpfchhadcand', 'I'), ('drmintau', 'F'), ('closesttaumvadiscr', 'F'), ('closesttaudecaymode', 'F')
        ,('taudr3wp', 'I'), ('taudr4wp', 'I'), ('taudr5wp', 'I')

        ,('detahighestptjet', 'F'), ('dphihighestptjet', 'F')
        ,('dphimet', 'F'), ('dphimetpca', 'F')

        ,('chi2', 'F')
        ,('quality', 'I')
        ,('nvalidhits', 'I'), ('nlosthits', 'I')

        ,('issignaltrack', 'I'), ('issusytrack', 'I'), ('susytrackmother', 'I'), ('susytrackpdgid', 'I')

        ,('hasGenMatch', 'I'), ('genmatchtmin', 'F')
        ,('genmatchpdgid', 'F'), ('genmatchpt', 'F'), ('genmatchstatus', 'F'), ('genmatchishardprocess', 'F'), ('genmatchisfromhardprocess', 'F')
        ,('genmatchisprompt', 'F'), ('genmatchisdirecthadrondecayproduct', 'F'), ('genmatchisdirecttaudecayproduct', 'F')
        ,('genmatchmotherpdgid', 'F'), ('genmatchmotherpt', 'F'), ('genmatchmotherstatus', 'F'), ('genmatchmotherishardprocess', 'F')
        ,('genmatchmotheristhetau', 'I'), ('genmatchmothertaudecay', 'F')
        ]

    track_level_var_array = {}
    for n in track_level_var_names:
        track_level_var_array[n[0]] = array('f', nMaxTracksPerEvent*[0.])
        tEvent.Branch(nice_string(n[0]), track_level_var_array[n[0]], nice_string(n[0]) + '[numtracksfinalpreselection]/F')


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
# define... stuff
###############################################################################################
'''

if True:

    jettype = 'AK4PFchs'

    if 'era16_07Aug17' in options.tag:

        # https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2016Analysis#Re_reco_datasets_07Aug17
        goldenjson = {'273158': [[1, 1279]], '273302': [[1, 459]], '273402': [[100, 292]], '273403': [[1, 53]], '273404': [[1, 18]], '273405': [[2, 25]], '273406': [[1, 112]], '273408': [[1, 6]], '273409': [[1, 309]], '273410': [[1, 90]], '273411': [[1, 29]], '273425': [[62, 352], [354, 733]], '273446': [[1, 33]], '273447': [[1, 113], [115, 412]], '273448': [[1, 391]], '273449': [[1, 214]], '273450': [[1, 214], [219, 647]], '273492': [[71, 71], [73, 282], [284, 325], [327, 338]], '273493': [[1, 233]], '273494': [[1, 192]], '273502': [[73, 256], [258, 318], [320, 813], [815, 1064]], '273503': [[1, 598]], '273554': [[77, 437]], '273555': [[1, 173]], '273725': [[83, 252], [254, 2545]], '273728': [[1, 100]], '273730': [[1, 1814], [1820, 2126]], '274094': [[108, 332]], '274146': [[1, 67]], '274157': [[105, 534]], '274159': [[1, 43]], '274160': [[1, 207]], '274161': [[1, 516]], '274172': [[31, 95]], '274198': [[81, 191]], '274199': [[1, 623]], '274200': [[1, 678]], '274240': [[1, 40], [42, 82]], '274241': [[1, 1152], [1161, 1176]], '274244': [[1, 607]], '274250': [[1, 701]], '274251': [[1, 546]], '274283': [[2, 19]], '274284': [[1, 210]], '274286': [[1, 154]], '274314': [[97, 97], [99, 158]], '274315': [[1, 424]], '274316': [[1, 959]], '274317': [[1, 3]], '274319': [[1, 225]], '274335': [[60, 1003]], '274336': [[1, 14]], '274337': [[3, 17]], '274338': [[1, 698]], '274339': [[1, 29], [31, 31], [33, 33], [35, 93]], '274344': [[1, 632]], '274345': [[1, 170]], '274382': [[94, 144]], '274387': [[88, 439]], '274388': [[1, 1820]], '274420': [[94, 268]], '274421': [[1, 342]], '274422': [[1, 2207]], '274440': [[92, 493]], '274441': [[1, 431]], '274442': [[1, 752]], '274954': [[37, 37], [39, 57]], '274955': [[1, 91]], '274968': [[1, 1192]], '274969': [[1, 1003]], '274970': [[1, 47]], '274971': [[1, 905]], '274998': [[64, 782]], '274999': [[1, 1241]], '275000': [[1, 136]], '275001': [[1, 1781], [1786, 2061]], '275059': [[78, 81], [105, 137]], '275066': [[1, 96]], '275067': [[1, 392]], '275068': [[1, 915]], '275073': [[1, 517]], '275074': [[1, 442], [444, 647]], '275124': [[106, 106], [108, 431]], '275125': [[1, 989]], '275282': [[91, 180]], '275283': [[1, 132]], '275284': [[1, 74]], '275290': [[96, 143]], '275291': [[1, 347]], '275292': [[1, 121]], '275293': [[1, 142], [144, 201]], '275309': [[55, 617]], '275310': [[1, 1929]], '275311': [[1, 1253]], '275319': [[141, 282]], '275337': [[1, 427]], '275338': [[1, 520]], '275344': [[76, 356]], '275345': [[1, 353]], '275370': [[81, 365]], '275371': [[1, 22], [28, 569]], '275375': [[127, 1449]], '275376': [[1, 2667], [2669, 3096]], '275657': [[1, 105]], '275658': [[1, 337]], '275659': [[1, 17]], '275761': [[1, 9]], '275767': [[1, 4]], '275772': [[1, 56]], '275773': [[1, 7]], '275774': [[1, 311], [315, 315]], '275776': [[1, 140]], '275777': [[1, 300]], '275778': [[1, 305]], '275782': [[1, 131], [133, 762]], '275832': [[1, 367]], '275833': [[1, 53], [56, 115], [117, 251]], '275834': [[1, 297]], '275835': [[1, 13]], '275836': [[1, 429], [431, 1163], [1166, 1170], [1184, 1293]], '275837': [[1, 186], [198, 726]], '275847': [[1, 2263]], '275886': [[73, 109]], '275890': [[1, 1393]], '275911': [[62, 298], [300, 354], [356, 440]], '275912': [[1, 289]], '275913': [[1, 475]], '275918': [[1, 318], [348, 361]], '275920': [[5, 463]], '275921': [[1, 2], [4, 5], [17, 20]], '275923': [[3, 53], [63, 64], [66, 126]], '275931': [[1, 14], [19, 89]], '275963': [[82, 139], [141, 172]], '276092': [[74, 149]], '276097': [[1, 507]], '276242': [[1, 7], [18, 61], [72, 1664]], '276243': [[1, 15], [18, 480], [482, 611]], '276244': [[3, 1202]], '276282': [[75, 534], [537, 1142]], '276283': [[3, 1087]], '276315': [[40, 175], [178, 217]], '276317': [[3, 138]], '276318': [[3, 103], [106, 570]], '276355': [[1, 33]], '276361': [[1, 161], [169, 208], [210, 800], [802, 833]], '276363': [[1, 140], [142, 238], [242, 1482]], '276384': [[2, 1117]], '276437': [[63, 224], [227, 1074], [1076, 2190]], '276454': [[1, 527]], '276458': [[1, 341]], '276495': [[87, 268]], '276501': [[4, 221], [223, 2547]], '276502': [[2, 741]], '276525': [[88, 469], [471, 1606], [1626, 2893]], '276527': [[1, 214]], '276528': [[4, 394]], '276542': [[74, 857]], '276543': [[1, 638], [643, 952]], '276544': [[2, 161]], '276545': [[2, 110], [117, 213]], '276581': [[79, 444]], '276582': [[1, 871]], '276583': [[1, 52]], '276584': [[1, 2]], '276585': [[1, 238], [241, 242], [245, 246]], '276586': [[2, 658], [680, 773]], '276587': [[1, 1006]], '276653': [[72, 550]], '276655': [[1, 593], [595, 1106]], '276659': [[1, 127], [129, 252]], '276775': [[96, 1260]], '276776': [[1, 1823]], '276794': [[1, 885]], '276807': [[66, 220]], '276808': [[1, 875]], '276810': [[1, 287]], '276811': [[1, 1270], [1272, 2563]], '276831': [[64, 755], [761, 2702]], '276834': [[1, 720]], '276870': [[78, 1354], [1356, 3108], [3111, 3258], [3260, 3484]], '276935': [[79, 184], [186, 838], [842, 906]], '276940': [[70, 213]], '276946': [[1, 27]], '276947': [[1, 89], [91, 126], [135, 141]], '276948': [[1, 474]], '276950': [[1, 2353]], '277069': [[81, 265], [267, 390]], '277070': [[1, 309], [311, 1059]], '277071': [[1, 82], [90, 178]], '277072': [[1, 253], [256, 466]], '277073': [[1, 90]], '277076': [[1, 3], [5, 7], [9, 35], [38, 1037]], '277087': [[204, 1191]], '277094': [[1, 161], [164, 584]], '277096': [[1, 1309], [1311, 2086]], '277112': [[1, 155]], '277126': [[42, 59]], '277127': [[1, 438], [440, 902]], '277148': [[83, 190], [193, 700]], '277166': [[77, 186], [188, 431]], '277168': [[1, 1708], [1711, 1822], [1824, 2223]], '277180': [[88, 228]], '277194': [[113, 139], [144, 497], [500, 1115], [1117, 1312], [1320, 1749], [1754, 2067], [2070, 2070]], '277305': [[62, 744]], '277420': [[84, 84], [86, 291], [293, 346]], '277981': [[82, 83], [85, 163]], '277991': [[1, 98]], '277992': [[1, 260], [262, 312]], '278017': [[77, 97], [99, 213], [215, 512], [514, 589]], '278018': [[1, 263], [265, 422], [424, 615], [617, 627], [642, 1011], [1020, 1181]], '278167': [[87, 394], [397, 1153], [1155, 1660], [1662, 1707], [1709, 2258]], '278175': [[1, 88]], '278193': [[77, 231]], '278239': [[76, 339], [341, 558], [560, 740]], '278240': [[1, 64], [70, 113], [115, 1121], [1123, 1296], [1299, 1309]], '278273': [[75, 110]], '278274': [[1, 18], [20, 85]], '278288': [[67, 81]], '278289': [[1, 42], [44, 52]], '278290': [[1, 11]], '278308': [[87, 216], [219, 587], [589, 680], [683, 1200], [1217, 1410], [1413, 1848], [1880, 1880]], '278310': [[1, 32], [34, 709]], '278315': [[73, 254], [256, 661], [663, 767]], '278345': [[84, 500], [503, 831]], '278346': [[1, 117]], '278349': [[1, 401], [403, 612], [632, 633]], '278366': [[1, 453]], '278406': [[85, 360], [362, 1682]], '278509': [[91, 1557]], '278769': [[75, 104]], '278770': [[1, 767]], '278801': [[48, 85]], '278802': [[1, 17]], '278803': [[1, 87], [91, 133], [135, 297], [299, 323]], '278804': [[1, 4]], '278805': [[3, 26], [30, 167], [170, 193], [196, 280], [283, 284], [288, 288]], '278808': [[1, 445], [447, 462], [464, 1793]], '278820': [[17, 1533]], '278822': [[1, 1627]], '278873': [[70, 129]], '278874': [[1, 273], [275, 478]], '278875': [[1, 210], [212, 834]], '278923': [[55, 467]], '278957': [[79, 227]], '278962': [[68, 408]], '278963': [[1, 23], [25, 175]], '278969': [[70, 511], [514, 1051], [1053, 1291], [1293, 1397], [1399, 1460]], '278975': [[1, 475], [477, 745], [747, 850]], '278976': [[1, 20]], '278986': [[71, 199]], '279024': [[82, 382]], '279029': [[1, 260]], '279071': [[71, 244]], '279080': [[68, 224]], '279115': [[118, 524]], '279116': [[38, 485]], '279479': [[86, 190]], '279588': [[100, 1259]], '279653': [[77, 77], [82, 261]], '279654': [[1, 108], [110, 1231], [1285, 1299]], '279656': [[1, 43]], '279658': [[1, 689], [691, 713]], '279667': [[68, 1033]], '279681': [[77, 104]], '279682': [[1, 29], [33, 34], [37, 38]], '279683': [[1, 26]], '279684': [[1, 22]], '279685': [[1, 93], [95, 209]], '279691': [[71, 113]], '279694': [[1, 2235]], '279715': [[71, 474], [476, 477], [480, 480], [511, 511], [523, 691]], '279716': [[1, 860], [875, 1528], [1530, 1653]], '279760': [[68, 578], [585, 728]], '279766': [[1, 1689]], '279767': [[1, 776]], '279794': [[77, 1100]], '279823': [[61, 395]], '279841': [[75, 398], [407, 2122]], '279844': [[72, 295]], '279887': [[79, 221], [225, 397]], '279931': [[84, 628], [630, 743], [746, 801], [803, 1043], [1045, 3022]], '279966': [[79, 441]], '279975': [[70, 190], [192, 253], [256, 281], [283, 709], [734, 1121]], '279993': [[85, 156]], '279994': [[1, 47]], '280013': [[1, 25]], '280015': [[1, 39], [41, 56], [59, 554], [560, 580]], '280016': [[1, 149]], '280017': [[1, 608]], '280018': [[1, 1281]], '280020': [[1, 45]], '280024': [[1, 427]], '280187': [[4, 60]], '280188': [[1, 245]], '280191': [[1, 781], [783, 866], [869, 900]], '280194': [[1, 238]], '280242': [[1, 411], [414, 627]], '280249': [[1, 486], [488, 1433]], '280251': [[1, 165], [167, 372]], '280327': [[49, 85]], '280330': [[1, 857]], '280349': [[1, 247], [252, 623], [626, 626]], '280363': [[1, 359]], '280364': [[1, 370], [372, 617], [619, 619], [621, 1090], [1102, 1363]], '280383': [[64, 65]], '280384': [[2, 34]], '280385': [[1, 519], [523, 569], [574, 1187], [1189, 1533], [1536, 2022]], '281613': [[101, 128], [130, 130], [133, 133], [135, 139], [143, 256], [258, 903]], '281639': [[1, 132]], '281641': [[1, 319]], '281693': [[1, 2191]], '281707': [[99, 982], [1000, 1065]], '281726': [[1, 288]], '281727': [[1, 1605]], '281797': [[125, 2176]], '281975': [[1, 215]], '281976': [[1, 2166]], '282033': [[82, 117]], '282034': [[1, 33]], '282035': [[1, 40]], '282037': [[1, 457], [459, 1862]], '282092': [[92, 222], [624, 2276]], '282708': [[1, 8]], '282710': [[1, 2], [8, 8]], '282712': [[1, 1], [10, 68]], '282730': [[89, 164]], '282731': [[1, 172]], '282732': [[1, 69]], '282733': [[1, 177]], '282734': [[1, 327]], '282735': [[1, 642], [645, 1232], [1235, 1823]], '282800': [[1, 377]], '282807': [[1, 326]], '282814': [[1, 1843]], '282842': [[1, 80]], '282917': [[117, 157], [159, 191]], '282918': [[1, 51]], '282919': [[1, 243]], '282922': [[1, 131]], '282923': [[1, 17], [19, 30], [32, 36], [38, 39], [41, 86], [88, 224]], '283042': [[1, 6]], '283043': [[1, 105], [108, 519]], '283049': [[82, 93]], '283050': [[1, 212]], '283052': [[1, 111]], '283059': [[1, 125], [127, 451]], '283270': [[76, 573], [576, 1502], [1504, 1888], [1890, 1912]], '283283': [[4, 1668], [1670, 1748]], '283305': [[79, 85]], '283306': [[1, 289]], '283307': [[1, 153], [156, 456]], '283308': [[1, 547], [549, 571], [573, 895], [897, 948]], '283353': [[80, 822]], '283358': [[1, 243], [245, 981]], '283359': [[1, 428]], '283407': [[82, 114]], '283408': [[1, 27], [29, 2088], [2098, 2125], [2203, 2416], [2528, 2542]], '283416': [[49, 151], [154, 245]], '283453': [[83, 537]], '283469': [[74, 74]], '283478': [[76, 303], [324, 969]], '283548': [[145, 288]], '283680': [[1, 81]], '283681': [[1, 17]], '283682': [[1, 384]], '283685': [[1, 314]], '283820': [[67, 1548]], '283830': [[1, 722]], '283834': [[1, 67], [69, 82]], '283835': [[1, 14], [16, 112]], '283865': [[1, 1177]], '283876': [[65, 211], [215, 724]], '283877': [[1, 1496]], '283884': [[349, 504], [509, 756]], '283885': [[1, 1723]], '283933': [[88, 232]], '283934': [[1, 784], [793, 870], [875, 1245], [1267, 1291]], '283946': [[85, 1448], [1450, 1462]], '283964': [[1, 388]], '284006': [[73, 390]], '284014': [[1, 266]], '284025': [[110, 157]], '284029': [[1, 112]], '284035': [[1, 360]], '284036': [[1, 140], [143, 348]], '284037': [[1, 340]], '284038': [[1, 55]], '284039': [[1, 30]], '284040': [[1, 33]], '284041': [[1, 44]], '284042': [[1, 129]], '284043': [[1, 205], [210, 224]], '284044': [[1, 30]]}

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        if 'data' in options.tag:  # data
            jet_energy_corrections = [
                [1, 276811, 'Summer16_07Aug2017BCD_V11_DATA'],
                [276831, 278801, 'Summer16_07Aug2017EF_V11_DATA'],
                [278802, float('inf'), 'Summer16_07Aug2017GH_V11_DATA']]
            DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'signal' in options.tag:  # FastSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer16_FastSimV1_MC/Summer16_FastSimV1_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
        else:  # FullSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

    elif 'era16_UL' in options.tag:

        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions16/13TeV/Legacy_2016/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage
        goldenjson = {'273158': [[1, 1283]],  '273302': [[1, 459]],  '273402': [[100, 292]],  '273403': [[1, 68]],  '273404': [[1, 22]],  '273405': [[2, 34]],  '273406': [[1, 125]],  '273408': [[1, 9]],  '273409': [[1, 317]],  '273410': [[1, 99]],  '273411': [[1, 29]],  '273425': [[62, 352], [354, 742]],  '273446': [[1, 48]],  '273447': [[1, 113], [115, 420]],  '273448': [[1, 396]],  '273449': [[1, 216]],  '273450': [[1, 214], [219, 647]],  '273492': [    [71, 282],    [284, 325],    [327, 343]  ],  '273493': [[1, 241]],  '273494': [[1, 192]],  '273502': [    [73, 256],    [258, 318],    [320, 813],    [815, 1077]  ],  '273503': [[1, 598]],  '273554': [[77, 444]],  '273555': [[1, 173]],  '273725': [[83, 252], [254, 2556]],  '273728': [[1, 112]],  '273730': [[1, 2126]],  '274094': [[108, 332]],  '274146': [[1, 73]],  '274157': [[105, 537]],  '274159': [[1, 47]],  '274160': [[1, 214]],  '274161': [[1, 516]],  '274172': [[31, 95]],  '274198': [[81, 192]],  '274199': [[1, 630]],  '274200': [[1, 678]],  '274240': [[1, 40], [42, 86]],  '274241': [[1, 1180]],  '274244': [[1, 607]],  '274250': [[1, 704]],  '274251': [[1, 546]],  '274283': [[1, 22]],  '274284': [[1, 215]],  '274286': [[1, 154]],  '274314': [[97, 165]],  '274315': [[1, 432]],  '274316': [[1, 974]],  '274317': [[1, 22]],  '274319': [[1, 225]],  '274335': [[60, 1013]],  '274336': [[1, 20]],  '274337': [[1, 22]],  '274338': [[1, 703]],  '274339': [[1, 99]],  '274344': [[1, 639]],  '274345': [[1, 170]],  '274382': [[94, 144]],  '274387': [[88, 447]],  '274388': [[1, 1820]],  '274420': [[93, 268]],  '274421': [[1, 356]],  '274422': [[1, 2207]],  '274440': [[92, 498]],  '274441': [[1, 443]],  '274442': [[1, 752]],  '274954': [[37, 60]],  '274955': [[1, 98]],  '274968': [[1, 1206]],  '274969': [[1, 1007]],  '274970': [[1, 47]],  '274971': [[1, 905]],  '274998': [[63, 794]],  '274999': [[1, 1257]],  '275000': [[1, 140]],  '275001': [[1, 1781], [1786, 2061]],  '275059': [[78, 81], [105, 150]],  '275066': [[1, 111]],  '275067': [[1, 396]],  '275068': [[1, 922]],  '275073': [[1, 523]],  '275074': [[1, 647]],  '275124': [[106, 433]],  '275125': [[1, 989]],  '275282': [[90, 185]],  '275283': [[1, 137]],  '275284': [[1, 74]],  '275290': [[96, 150]],  '275291': [[1, 356]],  '275292': [[1, 125]],  '275293': [[1, 142], [144, 201]],  '275309': [[55, 627]],  '275310': [[1, 1939]],  '275311': [[1, 1253]],  '275319': [[141, 292]],  '275337': [[1, 433]],  '275338': [[1, 520]],  '275344': [[76, 368]],  '275345': [[1, 353]],  '275370': [[81, 371]],  '275371': [[1, 569]],  '275375': [[127, 1453]],  '275376': [[1, 3096]],  '275657': [[1, 111]],  '275658': [[1, 344]],  '275659': [[1, 17]],  '275761': [[1, 11]],  '275767': [[1, 8]],  '275772': [[1, 61]],  '275773': [[1, 21]],  '275774': [[1, 317]],  '275776': [[1, 149]],  '275777': [[1, 304]],  '275778': [[1, 319]],  '275782': [[1, 131], [133, 762]],  '275832': [[1, 371]],  '275833': [    [1, 53],    [56, 115],    [117, 254]  ],  '275834': [[1, 303]],  '275835': [[1, 20]],  '275836': [    [1, 429],    [431, 1163],    [1166, 1170],    [1184, 1306]  ],  '275837': [[1, 186], [198, 726]],  '275847': [[1, 2263]],  '275886': [[70, 109]],  '275890': [[1, 1393]],  '275911': [    [62, 298],    [300, 354],    [356, 445]  ],  '275912': [[1, 303]],  '275913': [[1, 484]],  '275918': [[1, 318], [348, 361]],  '275920': [[1, 472]],  '275921': [[1, 32]],  '275923': [[1, 127]],  '275931': [[1, 89]],  '275963': [[82, 139], [141, 172]],  '276092': [[74, 153]],  '276097': [[1, 507]],  '276242': [[1, 7], [18, 61], [72, 1669]],  '276243': [[1, 15], [18, 627]],  '276244': [[1, 1202]],  '276282': [[75, 534], [537, 1151]],  '276283': [[1, 1087]],  '276315': [[40, 175], [178, 227]],  '276317': [[1, 147]],  '276318': [[1, 103], [106, 576]],  '276355': [[1, 34]],  '276361': [    [1, 161],    [169, 208],    [210, 800],    [802, 844]  ],  '276363': [[1, 238], [242, 1489]],  '276384': [[1, 1117]],  '276437': [    [63, 224],    [227, 1074],    [1076, 2190]  ],  '276454': [[1, 527]],  '276458': [[1, 341]],  '276495': [[87, 279]],  '276501': [[1, 221], [223, 2556]],  '276502': [[1, 741]],  '276525': [[87, 1606], [1626, 2904]],  '276527': [[1, 214]],  '276528': [[1, 394]],  '276542': [[74, 858]],  '276543': [[1, 961]],  '276544': [[1, 163]],  '276545': [[1, 110], [117, 213]],  '276581': [[79, 447]],  '276582': [[1, 873]],  '276583': [[1, 60]],  '276584': [[1, 2]],  '276585': [[1, 253]],  '276586': [[1, 658], [680, 781]],  '276587': [[1, 1006]],  '276653': [[72, 562]],  '276655': [[1, 593], [595, 1114]],  '276659': [[1, 127], [129, 252]],  '276775': [[96, 1269]],  '276776': [[1, 1823]],  '276794': [[1, 885]],  '276807': [[66, 227]],  '276808': [[1, 883]],  '276810': [[1, 291]],  '276811': [[1, 2563]],  '276831': [[64, 2711]],  '276834': [[1, 729]],  '276870': [    [78, 1354],    [1356, 3108],    [3111, 3258],    [3260, 3484]  ],  '276935': [    [79, 184],    [186, 838],    [842, 906]  ],  '276940': [[70, 214]],  '276946': [[1, 34]],  '276947': [[1, 150]],  '276948': [[1, 481]],  '276950': [[1, 2353]],  '277069': [[81, 394]],  '277070': [[1, 1063]],  '277071': [[1, 82], [90, 178]],  '277072': [[1, 253], [256, 484]],  '277073': [[1, 98]],  '277076': [    [1, 3],    [5, 7],    [9, 35],    [38, 1037]  ],  '277087': [[204, 1191]],  '277094': [[1, 161], [164, 598]],  '277096': [[1, 2093]],  '277112': [[1, 155]],  '277126': [[42, 60]],  '277127': [[1, 438], [440, 902]],  '277148': [[83, 715]],  '277166': [[77, 433]],  '277168': [[1, 2223]],  '277180': [[88, 228]],  '277194': [[113, 139], [144, 2070]],  '277305': [[62, 744]],  '277420': [[84, 346]],  '277981': [[82, 83], [85, 163]],  '277991': [[1, 98]],  '277992': [[1, 260], [262, 312]],  '278017': [    [77, 97],    [99, 213],    [215, 512],    [514, 600]  ],  '278018': [    [1, 263],    [265, 627],    [642, 1011],    [1020, 1181]  ],  '278167': [[87, 1660], [1662, 2260]],  '278175': [[1, 88]],  '278193': [[77, 231]],  '278239': [[76, 754]],  '278240': [    [1, 64],    [70, 113],    [115, 1309]  ],  '278273': [[75, 114]],  '278274': [[1, 85]],  '278288': [[67, 84]],  '278289': [[1, 42], [44, 60]],  '278290': [[1, 11]],  '278308': [    [87, 216],    [219, 1200],    [1217, 1848],    [1876, 1885]  ],  '278310': [[1, 709]],  '278315': [    [73, 254],    [256, 661],    [663, 767]  ],  '278345': [[84, 500], [503, 833]],  '278346': [[1, 150]],  '278349': [    [1, 401],    [403, 612],    [630, 639]  ],  '278366': [[1, 453]],  '278406': [[85, 360], [362, 1682]],  '278509': [[91, 1557]],  '278769': [[75, 111]],  '278770': [[1, 767]],  '278801': [[48, 92]],  '278802': [[1, 21]],  '278803': [[1, 330]],  '278804': [[1, 8]],  '278805': [[1, 26], [30, 291]],  '278808': [    [1, 445],    [447, 462],    [464, 1793]  ],  '278820': [[17, 1540]],  '278822': [[1, 1627]],  '278873': [[70, 136]],  '278874': [[1, 484]],  '278875': [[1, 834]],  '278923': [[55, 467]],  '278957': [[79, 227]],  '278962': [[68, 418]],  '278963': [[1, 23], [25, 175]],  '278969': [    [70, 1051],    [1053, 1291],    [1293, 1465]  ],  '278975': [[1, 857]],  '278976': [[1, 20]],  '278986': [[71, 199]],  '279024': [[82, 382]],  '279029': [[1, 260]],  '279071': [[71, 244]],  '279080': [[68, 224]],  '279115': [[118, 524]],  '279116': [[38, 485]],  '279479': [[86, 190]],  '279588': [[100, 1259]],  '279653': [[77, 77], [82, 268]],  '279654': [    [1, 108],    [110, 1231],    [1285, 1307]  ],  '279656': [[1, 43], [82, 87]],  '279658': [[1, 713]],  '279667': [[68, 1033]],  '279681': [[77, 111]],  '279682': [[1, 47]],  '279683': [[1, 34]],  '279684': [[1, 34]],  '279685': [[1, 93], [95, 209]],  '279691': [[71, 124]],  '279694': [[1, 2235]],  '279715': [    [71, 474],    [476, 477],    [480, 480],    [511, 511],    [523, 691]  ],  '279716': [    [1, 860],    [875, 1528],    [1530, 1653]  ],  '279760': [    [68, 578],    [585, 728],    [798, 806]  ],  '279766': [[1, 1694]],  '279767': [[1, 776]],  '279794': [[77, 1100]],  '279823': [[61, 395]],  '279841': [[75, 398], [407, 2122]],  '279844': [[72, 304]],  '279887': [[79, 397]],  '279931': [    [84, 628],    [630, 801],    [803, 1043],    [1045, 3022]  ],  '279966': [[79, 441]],  '279975': [    [70, 190],    [192, 253],    [256, 281],    [283, 709],    [734, 1121]  ],  '279993': [[85, 163]],  '279994': [[1, 59]],  '280013': [[1, 34]],  '280015': [    [1, 39],    [41, 56],    [59, 554],    [560, 584]  ],  '280016': [[1, 163]],  '280017': [[1, 613]],  '280018': [[1, 1282]],  '280020': [[1, 47]],  '280024': [[1, 427]],  '280187': [[4, 70]],  '280188': [[1, 253]],  '280191': [[1, 781], [783, 909]],  '280194': [[1, 238]],  '280242': [[1, 411], [414, 639]],  '280249': [[1, 1437]],  '280251': [[1, 165], [167, 372]],  '280327': [[49, 98]],  '280330': [[1, 870]],  '280349': [[1, 247], [252, 639]],  '280363': [[1, 367]],  '280364': [    [1, 619],    [621, 1090],    [1102, 1363]  ],  '280383': [[64, 73]],  '280384': [[1, 47]],  '280385': [[1, 519], [523, 2022]],  '281613': [[101, 903]],  '281639': [[1, 136]],  '281641': [[1, 319]],  '281693': [[1, 2191]],  '281707': [    [99, 982],    [1000, 1065],    [1087, 1089]  ],  '281726': [[1, 291]],  '281727': [[1, 1605]],  '281797': [[125, 2176]],  '281975': [[1, 215]],  '281976': [[1, 2166]],  '282033': [[82, 124]],  '282034': [[1, 35]],  '282035': [[1, 47]],  '282037': [[1, 457], [459, 1862]],  '282092': [[92, 222], [624, 2276]],  '282708': [[1, 8]],  '282710': [[1, 9]],  '282712': [[1, 1], [10, 68]],  '282730': [[89, 171]],  '282731': [[1, 176]],  '282732': [[1, 73]],  '282733': [[1, 178]],  '282734': [[1, 330]],  '282735': [[1, 1823]],  '282800': [[1, 382]],  '282807': [[1, 330]],  '282814': [[1, 1843]],  '282842': [[1, 80]],  '282917': [[117, 201]],  '282918': [[1, 59]],  '282919': [[1, 243]],  '282922': [[1, 137]],  '282923': [    [1, 17],    [19, 30],    [32, 86],    [88, 229]  ],  '283042': [[1, 10]],  '283043': [[1, 519]],  '283049': [[82, 98]],  '283050': [[1, 227]],  '283052': [[1, 124]],  '283059': [[1, 458]],  '283270': [[76, 1913]],  '283283': [[1, 1748]],  '283305': [[79, 93]],  '283306': [[1, 291]],  '283307': [[1, 461]],  '283308': [    [1, 547],    [549, 571],    [573, 948]  ],  '283353': [[80, 832]],  '283358': [[1, 243], [245, 986]],  '283359': [[1, 428]],  '283407': [[82, 124]],  '283408': [    [1, 2125],    [2203, 2416],    [2528, 2543]  ],  '283416': [[49, 245]],  '283453': [[83, 537]],  '283469': [[74, 74]],  '283478': [[76, 303], [324, 973]],  '283548': [[144, 291]],  '283680': [[1, 87]],  '283681': [[1, 23]],  '283682': [[1, 389]],  '283685': [[1, 314]],  '283820': [[67, 1552]],  '283830': [[1, 729]],  '283834': [[1, 85]],  '283835': [[1, 112]],  '283865': [[1, 1177]],  '283876': [[65, 736]],  '283877': [[1, 1496]],  '283884': [[349, 756]],  '283885': [[1, 1723]],  '283933': [[88, 240]],  '283934': [    [1, 784],    [793, 870],    [875, 1245],    [1267, 1291]  ],  '283946': [[85, 1462]],  '283964': [[1, 388]],  '284006': [[73, 394]],  '284014': [[1, 266]],  '284025': [[109, 162]],  '284029': [[1, 112]],  '284035': [[1, 369]],  '284036': [[1, 356]],  '284037': [[1, 343]],  '284038': [[1, 60]],  '284039': [[1, 34]],  '284040': [[1, 35]],  '284041': [[1, 47]],  '284042': [[1, 137]],  '284043': [[1, 227]],  '284044': [[1, 30]]}

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        # and https://twiki.cern.ch/twiki/bin/view/CMS/PdmVDataReprocessingUL2016
        if 'data' in options.tag:  # data
            jet_energy_corrections = [
                [1, 276811, 'Summer19UL16APV_RunBCD_V7_DATA'],
                [276831, 278807, 'Summer19UL16APV_RunEF_V7_DATA'],
                [278769, float('inf'), 'Summer19UL16_RunFGH_V7_DATA']]
            DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'signal' in options.tag:  # FastSim
            raise NotImplementedError('no JECs for UL FastSim')
        else:  # FullSim
            if 'era16_UL_APV' in options.tag:
                jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer19UL16_V7_MC/Summer19UL16_V7_MC',
                                   ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
            else:
                jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Summer19UL16APV_V7_MC/Summer19UL16APV_V7_MC',
                                   ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

    elif 'era17_17Nov2017' in options.tag:

        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions17/13TeV/ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage
        goldenjson = {"297050": [[12, 137], [193, 776]], "297056": [[12, 203]], "297057": [[1, 4], [14, 105], [112, 377], [385, 418], [424, 509], [516, 906]], "297099": [[24, 62]], "297100": [[1, 15], [21, 369], [375, 381]], "297101": [[1, 668], [673, 697], [700, 856], [862, 937], [943, 1101]], "297113": [[1, 204], [211, 252]], "297114": [[1, 99], [106, 161]], "297175": [[1, 85]], "297176": [[11, 120], [125, 214]], "297177": [[1, 162]], "297178": [[1, 54], [59, 334], [342, 749], [754, 967], [972, 1037], [1043, 1264], [1272, 1282], [1290, 1385]], "297215": [[1, 47]], "297218": [[1, 27]], "297219": [[1, 80], [85, 281], [288, 579], [585, 916], [921, 1429], [1436, 2004], [2010, 2638]], "297224": [[10, 19], [24, 138]], "297225": [[1, 32]], "297227": [[9, 192]], "297292": [[1, 125], [130, 131], [136, 667], [675, 753]], "297293": [[1, 121], [127, 150]], "297296": [[1, 236], [240, 401], [406, 418], [425, 497]], "297308": [[1, 44]], "297359": [[39, 70], [164, 180]], "297411": [[32, 737], [740, 800], [807, 950]], "297424": [[32, 149]], "297425": [[1, 107], [112, 157]], "297426": [[1, 28], [34, 84], [90, 111]], "297429": [[1, 72]], "297430": [[1, 199]], "297431": [[1, 49], [55, 64], [71, 188]], "297432": [[1, 112]], "297433": [[1, 159]], "297434": [[1, 161]], "297435": [[1, 94]], "297467": [[50, 138]], "297468": [[1, 74]], "297469": [[1, 4], [9, 70]], "297483": [[37, 68], [71, 201], [206, 214]], "297484": [[1, 47], [53, 208], [214, 214]], "297485": [[1, 16], [23, 253], [258, 299], [302, 314], [321, 420]], "297486": [[1, 74], [79, 598], [603, 625]], "297487": [[1, 433], [439, 491], [495, 603], [609, 613]], "297488": [[1, 73], [80, 424]], "297503": [[5, 275], [282, 559], [566, 606], [612, 635], [642, 772], [777, 779]], "297504": [[1, 41], [125, 136]], "297505": [[1, 394]], "297557": [[8, 28], [67, 113], [119, 167], [173, 174], [180, 394]], "297558": [[9, 266]], "297562": [[1, 69], [120, 369]], "297563": [[1, 254], [260, 264]], "297598": [[17, 17], [22, 33]], "297599": [[1, 169], [211, 225], [230, 312], [319, 385], [395, 407]], "297603": [[1, 420]], "297604": [[1, 126], [131, 272], [279, 375], [381, 407]], "297605": [[1, 6], [13, 20], [24, 89], [95, 223], [257, 407]], "297606": [[1, 94], [99, 231]], "297620": [[32, 318]], "297656": [[64, 116], [123, 135], [140, 230], [269, 307], [313, 330], [341, 388], [393, 433]], "297665": [[1, 153], [159, 209], [214, 279]], "297666": [[1, 11], [17, 81], [86, 121]], "297670": [[21, 34]], "297674": [[3, 102], [108, 188]], "297675": [[1, 123], [129, 239], [244, 328], [334, 467], [470, 471]], "297722": [[55, 160], [165, 353]], "297723": [[1, 13], [51, 222]], "298996": [[33, 216]], "298997": [[1, 37], [47, 47]], "299000": [[4, 77]], "299042": [[33, 55]], "299061": [[38, 355]], "299062": [[1, 163], [166, 303]], "299064": [[7, 85]], "299065": [[13, 248], [251, 342]], "299067": [[1, 459]], "299096": [[2, 97]], "299149": [[29, 470]], "299178": [[37, 56], [58, 111]], "299180": [[5, 98]], "299184": [[1, 561]], "299185": [[1, 120]], "299327": [[1, 72]], "299329": [[1, 172]], "299368": [[37, 175]], "299369": [[1, 303]], "299370": [[1, 7], [47, 705]], "299380": [[34, 227]], "299381": [[1, 45]], "299394": [[5, 33]], "299395": [[1, 187]], "299396": [[1, 81]], "299420": [[2, 50]], "299443": [[145, 164]], "299450": [[39, 88]], "299477": [[39, 42], [82, 87]], "299478": [[1, 175]], "299479": [[1, 123]], "299480": [[1, 6], [8, 715]], "299481": [[1, 196], [199, 236], [260, 479], [487, 940], [943, 1037], [1061, 1257]], "299593": [[95, 177], [179, 896]], "299594": [[1, 317]], "299595": [[1, 134], [138, 138]], "299597": [[3, 91], [93, 540]], "299649": [[151, 332]], "300087": [[36, 59], [61, 126], [128, 216], [218, 239]], "300105": [[1, 21]], "300106": [[1, 74]], "300107": [[1, 28], [30, 47]], "300117": [[35, 67]], "300122": [[46, 730], [735, 924], [927, 1295]], "300123": [[1, 384], [387, 612]], "300155": [[35, 1229]], "300156": [[1, 72]], "300157": [[9, 1107]], "300226": [[43, 448]], "300233": [[43, 162]], "300234": [[1, 59]], "300235": [[1, 187]], "300236": [[11, 187]], "300237": [[1, 713], [716, 717]], "300238": [[30, 58], [62, 329]], "300239": [[1, 145], [148, 167], [171, 213]], "300240": [[1, 7], [11, 46], [51, 362]], "300280": [[52, 56], [61, 69], [73, 150], [155, 165], [178, 198], [207, 222], [226, 251], [255, 268], [275, 345], [349, 370], [381, 548], [553, 607], [617, 639], [663, 691]], "300281": [[3, 8]], "300282": [[1, 9], [13, 59], [73, 92], [97, 114], [142, 151], [156, 186]], "300283": [[1, 34]], "300284": [[1, 22], [38, 47], [50, 82], [90, 98], [108, 130], [133, 152], [156, 250], [260, 414], [420, 561], [568, 585], [590, 680], [691, 751]], "300364": [[27, 46]], "300365": [[1, 20]], "300366": [[1, 21]], "300367": [[1, 20]], "300368": [[1, 20]], "300369": [[1, 20]], "300370": [[1, 20]], "300371": [[1, 20]], "300372": [[1, 8]], "300373": [[1, 21]], "300374": [[1, 21]], "300375": [[1, 93]], "300389": [[1, 1], [4, 5], [8, 8], [11, 20], [23, 39], [60, 149]], "300390": [[2, 21]], "300391": [[1, 21]], "300392": [[1, 21]], "300393": [[1, 20]], "300394": [[1, 21]], "300395": [[1, 20]], "300396": [[1, 20]], "300397": [[1, 20]], "300398": [[1, 20]], "300399": [[1, 20]], "300400": [[1, 677]], "300401": [[19, 673]], "300459": [[40, 332]], "300461": [[1, 98]], "300462": [[1, 97]], "300463": [[1, 124]], "300464": [[1, 103], [126, 265]], "300466": [[1, 650]], "300467": [[1, 563]], "300497": [[26, 175]], "300514": [[38, 150]], "300515": [[1, 838], [957, 1013]], "300516": [[1, 111]], "300517": [[1, 8], [103, 623]], "300558": [[8, 548]], "300560": [[1, 640], [645, 844]], "300574": [[15, 111]], "300575": [[1, 82]], "300576": [[7, 123], [125, 1206]], "300631": [[41, 49], [63, 66], [75, 226]], "300632": [[1, 21]], "300633": [[1, 447]], "300635": [[1, 23], [26, 176]], "300636": [[1, 335], [338, 1572]], "300673": [[41, 47], [49, 49], [52, 56], [59, 66]], "300674": [[1, 33]], "300675": [[1, 33]], "300676": [[1, 26]], "300742": [[56, 343]], "300777": [[21, 509]], "300780": [[3, 341]], "300785": [[1, 549], [552, 750], [752, 1201], [1219, 1272]], "300806": [[36, 214]], "300811": [[6, 508]], "300812": [[1, 59]], "300816": [[6, 161]], "300817": [[1, 33], [36, 74], [80, 383], [410, 493]], "301046": [[162, 223]], "301141": [[25, 31]], "301142": [[1, 897]], "301161": [[36, 805]], "301165": [[1, 145]], "301179": [[35, 59]], "301180": [[1, 97]], "301183": [[3, 10], [13, 303]], "301281": [[38, 157]], "301283": [[3, 886]], "301298": [[45, 949]], "301323": [[35, 474], [477, 990]], "301330": [[22, 353]], "301359": [[33, 319]], "301384": [[1, 476]], "301391": [[38, 214]], "301392": [[1, 627]], "301393": [[2, 18]], "301396": [[1, 33]], "301397": [[1, 228], [231, 517], [519, 728]], "301398": [[1, 9]], "301399": [[1, 108]], "301417": [[50, 367]], "301447": [[86, 96], [99, 400], [404, 512]], "301448": [[1, 329]], "301449": [[1, 404]], "301450": [[1, 173]], "301461": [[28, 581]], "301472": [[35, 830]], "301475": [[1, 18]], "301476": [[1, 844]], "301519": [[42, 250]], "301524": [[1, 110], [117, 263]], "301529": [[1, 49]], "301530": [[1, 110]], "301531": [[1, 394]], "301532": [[1, 611]], "301567": [[14, 372]], "301627": [[57, 943]], "301664": [[28, 445]], "301665": [[1, 294], [319, 487]], "301694": [[36, 102]], "301912": [[43, 52], [101, 422]], "301913": [[1, 58]], "301914": [[1, 350]], "301941": [[31, 568]], "301959": [[30, 1938]], "301960": [[1, 147]], "301970": [[6, 123]], "301984": [[17, 317]], "301985": [[1, 367]], "301986": [[1, 381]], "301987": [[1, 1128]], "301997": [[37, 407]], "301998": [[1, 1704]], "302019": [[34, 86]], "302026": [[24, 53], [66, 72]], "302029": [[1, 98]], "302031": [[1, 401], [403, 446], [448, 675], [678, 818]], "302033": [[1, 40], [44, 46]], "302034": [[1, 20]], "302037": [[18, 20]], "302038": [[10, 10]], "302040": [[1, 174]], "302041": [[1, 72]], "302042": [[1, 523]], "302043": [[1, 228]], "302131": [[71, 943]], "302159": [[33, 140]], "302163": [[32, 671], [674, 1230]], "302165": [[1, 85]], "302166": [[1, 16]], "302225": [[54, 133], [136, 923]], "302228": [[58, 78], [81, 293]], "302229": [[1, 457]], "302240": [[1, 960]], "302262": [[37, 471]], "302263": [[1, 1250]], "302277": [[15, 17], [22, 192], [194, 391]], "302279": [[1, 71]], "302280": [[1, 152]], "302322": [[33, 870]], "302328": [[42, 722]], "302337": [[27, 162]], "302342": [[19, 72]], "302343": [[1, 98]], "302344": [[3, 482]], "302350": [[1, 136]], "302388": [[27, 157], [164, 717]], "302392": [[45, 407]], "302393": [[1, 887]], "302448": [[21, 312], [317, 442], [445, 483], [486, 1926]], "302472": [[28, 808]], "302473": [[1, 368], [398, 406]], "302474": [[1, 305]], "302475": [[1, 7]], "302476": [[1, 259]], "302479": [[30, 222], [225, 340]], "302484": [[8, 176]], "302485": [[1, 922]], "302492": [[10, 21], [23, 59]], "302493": [[1, 7]], "302494": [[1, 618]], "302509": [[73, 92]], "302513": [[37, 89]], "302522": [[29, 46]], "302523": [[1, 59]], "302525": [[1, 677], [747, 778]], "302526": [[1, 582]], "302548": [[40, 124]], "302551": [[1, 7]], "302553": [[1, 188]], "302554": [[1, 7]], "302555": [[1, 11]], "302563": [[40, 46]], "302565": [[1, 7]], "302572": [[6, 291]], "302573": [[1, 693], [730, 1285]], "302596": [[47, 534], [545, 705], [710, 986]], "302597": [[1, 1054]], "302634": [[37, 73], [75, 123], [125, 129], [133, 165], [168, 175], [177, 216], [218, 358], [361, 375], [378, 404], [407, 423], [425, 503], [505, 578], [581, 594], [596, 638]], "302635": [[1, 22], [24, 28], [30, 39], [41, 53], [55, 132], [134, 144], [146, 265], [267, 271], [274, 344], [347, 357], [359, 375], [378, 384], [386, 414], [416, 494], [497, 608], [611, 634], [637, 684], [687, 706], [708, 724], [726, 901], [904, 954], [957, 982], [984, 1072], [1075, 1124], [1126, 1129], [1132, 1206], [1209, 1234], [1236, 1291]], "302651": [[1, 149]], "302654": [[1, 317]], "302661": [[1, 72]], "302663": [[1, 706]], "303825": [[1, 180]], "303832": [[54, 1334], [1338, 1913]], "303838": [[54, 54], [83, 2044]], "303885": [[60, 2052]], "303948": [[55, 1678]], "303998": [[58, 319]], "303999": [[1, 751]], "304000": [[1, 56]], "304062": [[54, 2014]], "304119": [[71, 138], [143, 150]], "304120": [[1, 253]], "304125": [[1, 1769]], "304144": [[76, 2596], [2598, 2656]], "304158": [[165, 1750], [1752, 2087]], "304169": [[50, 1714], [1731, 1733]], "304170": [[1, 620]], "304199": [[10, 18]], "304200": [[1, 321]], "304204": [[55, 607]], "304209": [[52, 98], [100, 133], [135, 157], [176, 253], [255, 477]], "304291": [[56, 85]], "304292": [[1, 1125], [1183, 1779], [1781, 1811]], "304333": [[74, 1653]], "304354": [[82, 295]], "304366": [[44, 1387], [1390, 1396], [1399, 1402], [1404, 1407], [1409, 1412], [1414, 1416], [1419, 1421], [1424, 1873]], "304446": [[40, 92], [110, 111]], "304447": [[1, 534], [540, 1644]], "304451": [[1, 60]], "304505": [[60, 86]], "304506": [[1, 370]], "304507": [[1, 239]], "304508": [[1, 1324]], "304562": [[52, 56], [60, 848]], "304616": [[52, 223], [227, 740], [747, 1002]], "304625": [[73, 536]], "304626": [[1, 8]], "304654": [[53, 704]], "304655": [[1, 1194]], "304661": [[53, 67], [69, 143], [147, 173], [175, 198], [237, 240]], "304662": [[1, 150]], "304663": [[1, 689]], "304671": [[51, 1193]], "304672": [[1, 60]], "304737": [[69, 149]], "304738": [[1, 1681]], "304739": [[3, 16]], "304740": [[1, 278]], "304776": [[49, 98]], "304777": [[1, 431], [438, 510]], "304778": [[4, 1300]], "304797": [[28, 87], [91, 306], [308, 377], [385, 1202], [1205, 2950]], "305044": [[3, 203], [302, 306], [309, 310], [313, 313], [318, 330]], "305045": [[1, 873]], "305046": [[1, 667], [671, 686]], "305059": [[63, 518], [520, 575]], "305062": [[1, 8]], "305063": [[1, 35]], "305064": [[1, 2045]], "305081": [[52, 1107]], "305112": [[68, 1527]], "305113": [[9, 72]], "305114": [[1, 526]], "305178": [[69, 124]], "305179": [[1, 21]], "305180": [[1, 9]], "305181": [[1, 8]], "305182": [[1, 8]], "305183": [[1, 231], [262, 266]], "305184": [[1, 8]], "305186": [[1, 112], [120, 422]], "305188": [[1, 1002]], "305202": [[74, 132], [136, 729]], "305204": [[1, 1229]], "305207": [[52, 1077]], "305208": [[1, 372]], "305234": [[52, 99]], "305236": [[1, 23]], "305237": [[1, 16], [18, 1147]], "305247": [[57, 433]], "305248": [[1, 957]], "305252": [[1, 548]], "305282": [[75, 207]], "305310": [[60, 157], [163, 458]], "305311": [[1, 153]], "305312": [[1, 227]], "305313": [[1, 741]], "305314": [[1, 404]], "305336": [[36, 241]], "305338": [[1, 107]], "305341": [[1, 503]], "305349": [[1, 34]], "305350": [[1, 21]], "305351": [[1, 868]], "305358": [[91, 231], [233, 253]], "305364": [[50, 147]], "305365": [[1, 668], [676, 832]], "305366": [[1, 721], [724, 756], [769, 934], [936, 1254]], "305376": [[71, 168]], "305377": [[9, 1292], [1294, 1383], [1386, 1525]], "305405": [[44, 536], [573, 575]], "305406": [[1, 394], [401, 520], [528, 535], [540, 1475]], "305440": [[20, 291]], "305441": [[1, 121]], "305516": [[46, 518], [558, 639]], "305517": [[1, 163]], "305518": [[1, 1134]], "305586": [[53, 583]], "305589": [[1, 691]], "305590": [[1, 500], [517, 1020]], "305636": [[60, 339], [342, 667], [671, 2390]], "305766": [[55, 902]], "305809": [[56, 197]], "305814": [[85, 689], [692, 978], [980, 1074], [1077, 1912]], "305821": [[59, 830]], "305832": [[87, 266]], "305840": [[1, 1144]], "305842": [[1, 862]], "305862": [[81, 705]], "305898": [[70, 780]], "305902": [[53, 521]], "305967": [[1, 32]], "306029": [[63, 96]], "306030": [[1, 110]], "306036": [[60, 63]], "306037": [[1, 49]], "306038": [[1, 139]], "306041": [[1, 320]], "306042": [[1, 371]], "306048": [[1, 140]], "306049": [[1, 358]], "306051": [[1, 415]], "306091": [[422, 629]], "306092": [[1, 588], [593, 976]], "306095": [[1, 300]], "306121": [[57, 152]], "306122": [[1, 127]], "306125": [[1, 756], [770, 2642], [2667, 3007]], "306126": [[1, 497]], "306134": [[53, 84]], "306135": [[1, 1095]], "306138": [[1, 1298]], "306139": [[1, 1112]], "306153": [[78, 165]], "306154": [[1, 251], [253, 691], [709, 1233]], "306155": [[1, 1440]], "306169": [[1, 745]], "306170": [[1, 22]], "306171": [[1, 503]], "306418": [[1, 33], [35, 75]], "306419": [[1, 62]], "306420": [[1, 108]], "306422": [[9, 126]], "306423": [[1, 333]], "306432": [[1, 339]], "306454": [[13, 101]], "306455": [[1, 11]], "306456": [[1, 237], [239, 787]], "306457": [[1, 31]], "306458": [[1, 17], [20, 35], [37, 41], [43, 47], [49, 53], [56, 60], [62, 66], [68, 72], [74, 77], [79, 83], [85, 89], [93, 102], [104, 108], [110, 114], [116, 120], [122, 126], [129, 139], [141, 145], [147, 151], [153, 166], [169, 173], [175, 179], [181, 185], [187, 191], [193, 197], [200, 210], [212, 216], [218, 222], [225, 229], [231, 235], [237, 241], [243, 247], [249, 249], [252, 256], [258, 268]], "306459": [[1, 512], [514, 2275]], "306460": [[1, 73]]}

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        if 'data' in options.tag:  # data
            raise NotImplementedError('no JECs yet for 2017 data')
            # TODO: implement JECs for 2017 data
            # jet_energy_corrections = []
            # DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'signal' in options.tag:  # FastSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Fall17_FastSimV1_MC/Fall17_FastSimV1_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
        else:  # FullSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Fall17_17Nov2017_V32_MC/Fall17_17Nov2017_V32_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

    elif 'era18_17Sep2018' in options.tag:

        # https://cms-service-dqmdc.web.cern.ch/CAF/certification/Collisions18/13TeV/ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
        # from https://twiki.cern.ch/twiki/bin/viewauth/CMS/DCUserPage
        goldenjson = {"315257": [[1, 88], [91, 92]], "315259": [[1, 172]], "315264": [[32, 261]], "315265": [[4, 58]], "315267": [[1, 244]], "315270": [[1, 633]], "315322": [[23, 118], [122, 1354]], "315339": [[37, 654]], "315357": [[44, 732], [736, 770], [780, 831]], "315361": [[40, 619]], "315363": [[1, 35], [37, 47], [49, 67], [69, 80], [82, 90]], "315366": [[10, 61], [67, 750]], "315420": [[28, 920], [924, 942], [954, 1748]], "315488": [[42, 843]], "315489": [[1, 653], [672, 709]], "315490": [[1, 24]], "315506": [[13, 100]], "315510": [[1, 345]], "315512": [[1, 1122]], "315543": [[55, 171]], "315555": [[22, 97]], "315556": [[1, 26]], "315557": [[1, 279]], "315640": [[46, 87]], "315641": [[1, 4]], "315642": [[1, 92]], "315644": [[1, 184]], "315645": [[1, 40], [47, 390], [395, 565], [567, 594]], "315646": [[1, 1033]], "315647": [[1, 58]], "315648": [[1, 110]], "315689": [[24, 1127], [1180, 1186]], "315690": [[10, 654]], "315702": [[38, 113]], "315703": [[1, 545]], "315704": [[1, 61]], "315705": [[1, 700]], "315713": [[35, 359], [374, 385], [400, 1123]], "315721": [[33, 50], [56, 626]], "315741": [[34, 92]], "315764": [[37, 309]], "315770": [[39, 332]], "315784": [[29, 33], [40, 156], [158, 161]], "315785": [[1, 198], [201, 305]], "315786": [[1, 72]], "315790": [[1, 716], [718, 922]], "315800": [[41, 621]], "315801": [[1, 344]], "315840": [[33, 1154]], "315973": [[39, 240], [262, 914]], "315974": [[1, 71]], "316058": [[42, 405]], "316059": [[1, 321], [323, 567]], "316060": [[1, 935]], "316061": [[1, 23], [194, 206]], "316062": [[1, 4]], "316082": [[37, 407]], "316110": [[1, 210]], "316111": [[1, 48]], "316113": [[1, 64]], "316114": [[1, 777], [779, 1562]], "316153": [[1, 770]], "316186": [[38, 81]], "316187": [[1, 1091], [1093, 1100], [1207, 2077]], "316199": [[33, 1197]], "316200": [[1, 10]], "316201": [[1, 498]], "316202": [[1, 403]], "316216": [[25, 466]], "316217": [[1, 264]], "316218": [[1, 1008]], "316219": [[1, 283]], "316239": [[38, 626]], "316240": [[1, 1224]], "316241": [[1, 325]], "316271": [[36, 121]], "316361": [[22, 124], [126, 131], [133, 135], [137, 137], [139, 142], [144, 145], [147, 147], [149, 159], [161, 174], [176, 178], [180, 189], [191, 197], [199, 208], [210, 223]], "316362": [[1, 208], [210, 212], [214, 225], [227, 242], [244, 269], [271, 319], [332, 392], [394, 395], [397, 402], [404, 404], [406, 410], [412, 412], [414, 418], [420, 428], [430, 450]], "316363": [[1, 39], [41, 49]], "316377": [[19, 19], [21, 40]], "316378": [[1, 29]], "316379": [[1, 70]], "316380": [[1, 708], [714, 1213]], "316455": [[36, 71]], "316457": [[1, 1454]], "316469": [[17, 444]], "316470": [[1, 476]], "316472": [[1, 70], [76, 333]], "316505": [[44, 205], [207, 921], [923, 1364]], "316569": [[20, 703], [742, 1945]], "316590": [[17, 526]], "316613": [[49, 241]], "316615": [[1, 338]], "316666": [[1, 981]], "316667": [[1, 197]], "316700": [[46, 346], [388, 397]], "316701": [[1, 479]], "316702": [[1, 388]], "316715": [[33, 45]], "316716": [[1, 181]], "316717": [[1, 192]], "316718": [[1, 311]], "316719": [[1, 91], [100, 144]], "316720": [[1, 182]], "316721": [[1, 15]], "316722": [[1, 751]], "316723": [[1, 64]], "316758": [[11, 1609]], "316766": [[51, 1920], [1922, 2199]], "316876": [[34, 38], [40, 644]], "316877": [[1, 164], [171, 401]], "316879": [[1, 156]], "316928": [[40, 188]], "316985": [[33, 503]], "316993": [[44, 254]], "316994": [[1, 14]], "316995": [[1, 623]], "317080": [[41, 66]], "317087": [[43, 177], [213, 222], [257, 852]], "317089": [[1, 1003]], "317182": [[47, 63], [65, 1424]], "317212": [[36, 175]], "317213": [[1, 375]], "317279": [[43, 508]], "317291": [[34, 824]], "317292": [[1, 330]], "317297": [[1, 283], [347, 760]], "317319": [[44, 182]], "317320": [[1, 326], [333, 411], [413, 1827]], "317338": [[66, 107]], "317339": [[1, 163]], "317340": [[1, 418]], "317382": [[58, 128]], "317383": [[1, 58]], "317391": [[39, 46]], "317392": [[1, 1116], [1119, 1900]], "317435": [[1, 1397]], "317438": [[1, 68], [71, 309]], "317475": [[33, 89], [105, 115]], "317478": [[1, 23]], "317484": [[1, 448], [467, 514], [519, 545]], "317488": [[1, 844]], "317527": [[41, 1487]], "317591": [[43, 334]], "317626": [[40, 2045]], "317640": [[29, 829]], "317641": [[1, 1390]], "317648": [[45, 139]], "317649": [[1, 621]], "317650": [[1, 1304]], "317661": [[35, 1256]], "317663": [[1, 858]], "317683": [[83, 402]], "317696": [[38, 682]], "318733": [[1, 33]], "318828": [[54, 123]], "318872": [[16, 287]], "318874": [[1, 320]], "318876": [[1, 161]], "318877": [[1, 615]], "319077": [[52, 92]], "319337": [[48, 2240]], "319347": [[40, 690]], "319348": [[1, 37]], "319349": [[1, 148]], "319449": [[35, 559], [562, 734]], "319450": [[1, 287], [290, 683]], "319456": [[138, 346]], "319459": [[1, 78]], "319486": [[38, 103]], "319503": [[1, 317]], "319524": [[36, 1459]], "319526": [[1, 282]], "319528": [[1, 259]], "319579": [[41, 3168]], "319625": [[17, 206]], "319639": [[31, 1509]], "319656": [[51, 310]], "319657": [[1, 167]], "319658": [[1, 225]], "319659": [[1, 87]], "319678": [[36, 294]], "319687": [[46, 90]], "319697": [[47, 482], [490, 490]], "319698": [[1, 312]], "319756": [[44, 1966]], "319840": [[41, 388]], "319841": [[1, 167]], "319847": [[49, 51]], "319848": [[1, 53]], "319849": [[1, 492]], "319851": [[1, 4]], "319853": [[1, 40], [47, 262]], "319854": [[1, 225]], "319908": [[1, 40], [43, 53]], "319909": [[1, 7]], "319910": [[1, 983]], "319912": [[1, 59]], "319913": [[1, 56]], "319914": [[1, 32]], "319915": [[1, 416]], "319941": [[43, 298]], "319942": [[1, 50]], "319950": [[38, 205]], "319991": [[46, 882]], "319992": [[1, 264]], "319993": [[1, 955]], "320002": [[52, 192]], "320006": [[1, 34], [36, 341]], "320010": [[1, 330]], "320011": [[1, 302]], "320012": [[1, 99]], "320023": [[17, 292]], "320024": [[1, 410]], "320025": [[1, 113]], "320026": [[1, 204]], "320038": [[43, 663]], "320039": [[1, 30]], "320040": [[1, 737]], "320059": [[1, 105]], "320060": [[1, 42]], "320061": [[1, 49]], "320062": [[1, 21]], "320063": [[1, 64]], "320064": [[1, 200]], "320065": [[1, 920]], "320673": [[35, 901]], "320674": [[1, 599]], "320688": [[49, 531]], "320712": [[39, 242]], "320757": [[51, 382]], "320804": [[46, 1274]], "320807": [[1, 7]], "320809": [[1, 716]], "320821": [[41, 221]], "320822": [[1, 523]], "320823": [[1, 360]], "320824": [[1, 1051]], "320838": [[93, 357]], "320840": [[1, 471]], "320841": [[1, 205]], "320853": [[41, 369]], "320854": [[1, 125]], "320855": [[1, 565]], "320856": [[1, 159]], "320857": [[1, 272]], "320858": [[1, 230]], "320859": [[1, 40]], "320887": [[49, 321]], "320888": [[1, 26]], "320916": [[2, 25]], "320917": [[1, 1926]], "320920": [[1, 178]], "320933": [[40, 214]], "320934": [[1, 831]], "320936": [[1, 407]], "320941": [[1, 93]], "320980": [[44, 142]], "320995": [[26, 214]], "320996": [[1, 380]], "321004": [[39, 188]], "321005": [[1, 61]], "321006": [[1, 162]], "321007": [[1, 831]], "321009": [[1, 85]], "321010": [[1, 342]], "321011": [[1, 213]], "321012": [[1, 35], [190, 201]], "321051": [[58, 1179]], "321055": [[1, 302], [304, 326], [328, 340], [368, 759]], "321067": [[39, 225], [232, 639]], "321068": [[1, 715]], "321069": [[1, 313]], "321119": [[45, 214]], "321121": [[1, 47]], "321122": [[1, 395]], "321124": [[1, 819]], "321126": [[1, 493]], "321134": [[33, 70]], "321138": [[1, 741]], "321140": [[1, 798]], "321149": [[35, 1424], [1426, 1476], [1478, 1553], [1558, 1576], [1578, 1588], [1591, 1743]], "321165": [[1, 8]], "321166": [[1, 10]], "321167": [[1, 141], [143, 143], [145, 510], [512, 552], [554, 691], [693, 923]], "321177": [[38, 74], [77, 214], [216, 232], [234, 247], [249, 321], [323, 365], [367, 455]], "321178": [[5, 78]], "321218": [[49, 962]], "321219": [[1, 934]], "321221": [[1, 40]], "321230": [[41, 124]], "321231": [[1, 59]], "321232": [[1, 30]], "321233": [[1, 727]], "321262": [[1, 4]], "321283": [[48, 357]], "321294": [[1, 62]], "321295": [[1, 307], [309, 316], [318, 384], [390, 394], [396, 604], [606, 616], [619, 646], [649, 690], [693, 754]], "321296": [[1, 24], [34, 41], [44, 67]], "321305": [[20, 2600], [2605, 2651]], "321311": [[1, 10]], "321312": [[1, 768]], "321313": [[1, 408]], "321393": [[1, 127], [134, 148]], "321396": [[1, 1475]], "321397": [[1, 365]], "321414": [[31, 1283]], "321415": [[1, 804]], "321431": [[30, 189]], "321432": [[1, 47]], "321433": [[1, 125]], "321434": [[1, 642]], "321436": [[1, 710]], "321457": [[43, 451], [453, 1888]], "321461": [[1, 149]], "321475": [[50, 518], [526, 2084]], "321710": [[1, 57]], "321712": [[1, 2], [16, 54], [57, 115], [117, 263]], "321730": [[2, 257], [259, 291]], "321732": [[1, 127], [129, 181], [185, 189], [192, 245], [248, 252], [254, 373], [375, 381], [386, 386], [389, 392], [395, 424], [426, 432], [434, 448], [450, 452], [454, 459], [467, 586], [589, 680], [682, 686], [689, 903], [905, 973], [975, 1448]], "321735": [[1, 146]], "321755": [[33, 361], [363, 470], [472, 473], [475, 487], [489, 729]], "321758": [[1, 47], [49, 75], [77, 121], [128, 130], [146, 148], [151, 155], [161, 165], [168, 189]], "321760": [[1, 171], [175, 205], [207, 238], [240, 258], [260, 420], [422, 520], [526, 586], [588, 593], [598, 602], [604, 607], [613, 716], [719, 721], [727, 788], [794, 818], [822, 824], [828, 830], [834, 836], [840, 841], [845, 855]], "321773": [[11, 14], [25, 35], [39, 52], [54, 79]], "321774": [[1, 12], [14, 52], [54, 119]], "321775": [[1, 12], [14, 14]], "321776": [[1, 12], [15, 19], [30, 45]], "321777": [[1, 81], [83, 169], [174, 176], [192, 207]], "321778": [[8, 150]], "321780": [[1, 332], [336, 338], [342, 346], [351, 357], [359, 360], [362, 371], [374, 383], [392, 412], [414, 420], [422, 493], [496, 499], [502, 503], [505, 508], [517, 518]], "321781": [[6, 37], [53, 56], [58, 66], [69, 69], [77, 180], [186, 209], [212, 265], [269, 274], [276, 290], [293, 312], [316, 410], [412, 427]], "321813": [[32, 352]], "321815": [[1, 23]], "321817": [[1, 536]], "321818": [[1, 690]], "321820": [[1, 214]], "321831": [[25, 781]], "321832": [[1, 389], [403, 510]], "321833": [[1, 407]], "321834": [[1, 333]], "321879": [[39, 47], [50, 52], [55, 68], [71, 73], [77, 89], [93, 95], [99, 111], [114, 116], [120, 132], [136, 138], [141, 154], [157, 159], [163, 175], [178, 181], [185, 197], [200, 202], [207, 218], [222, 356]], "321880": [[1, 41], [44, 132]], "321887": [[54, 948]], "321908": [[43, 472]], "321909": [[1, 208], [210, 1654]], "321917": [[4, 156], [164, 808]], "321919": [[1, 6]], "321933": [[43, 232], [235, 326]], "321960": [[18, 47]], "321961": [[1, 354]], "321973": [[37, 746], [748, 968], [972, 1253]], "321975": [[1, 866]], "321988": [[45, 996], [1106, 1486]], "321990": [[1, 471]], "322013": [[14, 22]], "322014": [[1, 17]], "322022": [[42, 185], [201, 1805]], "322040": [[32, 70]], "322057": [[38, 58]], "322068": [[51, 724]], "322079": [[39, 200], [216, 393], [409, 428]], "322106": [[48, 871]], "322113": [[48, 159]], "322118": [[1, 516], [530, 874]], "322179": [[43, 820], [823, 1783]], "322201": [[39, 266]], "322204": [[1, 280], [282, 301], [303, 331], [337, 1143]], "322222": [[1, 526]], "322252": [[42, 1586]], "322317": [[48, 101]], "322319": [[1, 163]], "322322": [[1, 170], [267, 1205]], "322324": [[1, 416]], "322332": [[37, 1055]], "322348": [[40, 1505]], "322355": [[36, 137]], "322356": [[1, 779]], "322381": [[45, 577]], "322407": [[46, 582]], "322430": [[46, 501]], "322431": [[59, 1166]], "322480": [[60, 408]], "322492": [[1, 1386]], "322510": [[37, 45]], "322599": [[43, 294]], "322602": [[1, 69], [72, 72]], "322603": [[1, 10]], "322605": [[1, 280]], "322617": [[1, 601]], "322625": [[41, 484], [492, 1167]], "322633": [[1, 249]], "323414": [[1, 46]], "323423": [[1, 136]], "323470": [[38, 172], [176, 218], [223, 266]], "323471": [[1, 238]], "323472": [[1, 64]], "323473": [[1, 227]], "323474": [[1, 355]], "323475": [[1, 77]], "323487": [[42, 177], [184, 498]], "323488": [[1, 514], [555, 734], [738, 793]], "323492": [[1, 33]], "323493": [[1, 144]], "323495": [[1, 187]], "323524": [[25, 561]], "323525": [[1, 91], [97, 1126]], "323526": [[1, 248], [253, 466]], "323693": [[38, 151]], "323696": [[1, 257]], "323702": [[1, 808]], "323725": [[18, 346]], "323726": [[1, 60]], "323727": [[1, 83], [88, 677], [682, 813], [819, 822], [826, 987]], "323755": [[27, 815], [818, 823], [826, 826], [828, 830], [833, 861], [864, 964]], "323775": [[38, 81], [84, 171]], "323778": [[1, 934]], "323790": [[45, 948]], "323794": [[1, 68]], "323841": [[46, 510]], "323857": [[1, 357]], "323940": [[49, 1567]], "323954": [[1, 77]], "323976": [[31, 85]], "323978": [[1, 73]], "323980": [[1, 202]], "323983": [[1, 188]], "323997": [[1, 498]], "324021": [[44, 819]], "324022": [[1, 554]], "324077": [[54, 710], [712, 753]], "324201": [[20, 834], [837, 1385]], "324202": [[1, 240]], "324205": [[1, 163]], "324206": [[1, 149]], "324207": [[1, 34]], "324209": [[1, 142]], "324237": [[33, 236]], "324245": [[23, 1681]], "324293": [[39, 1440], [1442, 2176], [2178, 2342]], "324315": [[1, 200], [203, 204]], "324318": [[1, 332]], "324420": [[1, 625]], "324729": [[1, 193]], "324747": [[63, 1139]], "324764": [[1, 150]], "324765": [[1, 481]], "324769": [[1, 328]], "324772": [[1, 165]], "324785": [[77, 664]], "324791": [[1, 1217]], "324835": [[40, 230], [302, 369]], "324840": [[1, 96]], "324841": [[1, 1347]], "324846": [[1, 151], [154, 517]], "324878": [[62, 111], [113, 175], [180, 1800]], "324897": [[30, 170]], "324970": [[1, 425], [428, 598], [606, 632], [634, 1529], [1532, 2195]], "324980": [[39, 917], [919, 954], [956, 968], [1005, 1042], [1044, 2340]], "324997": [[29, 150]], "324998": [[1, 368]], "324999": [[1, 14]], "325000": [[1, 371]], "325001": [[1, 105], [108, 171], [173, 595]], "325022": [[45, 1594]], "325057": [[42, 383]], "325097": [[40, 96]], "325098": [[1, 8]], "325099": [[1, 394]], "325100": [[1, 254]], "325101": [[1, 462], [464, 485]], "325110": [[1, 21]], "325117": [[1, 533]], "325159": [[48, 266]], "325168": [[1, 21]], "325169": [[1, 23]], "325170": [[1, 692], [694, 1205]], "325172": [[1, 267], [269, 485]]}

        # from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
        if 'data' in options.tag:  # data
            raise NotImplementedError('no JECs yet for 2018 data')
            # TODO: implement JECs for 2018 data
            # jet_energy_corrections = []
            # DataJECs = DataJEC(jet_energy_corrections, jettype)
        elif 'signal' in options.tag:  # FastSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Autumn18_FastSimV1_MC/Autumn18_FastSimV1_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
        else:  # FullSim
            jecAK4 = createJEC('/nfs/dust/cms/user/wolfmor/JECs/Autumn18_V19_MC/Autumn18_V19_MC',
                               ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)

    else:

        raise NotImplementedError('JECs: era unknown or not specified')


'''
###############################################################################################
# define handles and labels
###############################################################################################
'''

if 'pmssm' in options.tag:
    handle_lumis = Handle('GenLumiInfoHeader')

if 'signal' not in options.tag and 'pmssm' not in options.tag:

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
label_taudiscriminatorDM = ('hpsPFTauDiscriminationByDecayModeFindingNewDMs')

# this discriminator (MVArun2v1DBnewDM) seems to be the one referred to as "MVA2017v2"???
# see https://cms-nanoaod-integration.web.cern.ch/integration/master/mc102X_doc.html#Tau

handle_taudiscriminatorMVA_VLoose = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_VLoose = ('hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVA_Loose = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_Loose = ('hpsPFTauDiscriminationByLooseIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVA_Medium = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_Medium = ('hpsPFTauDiscriminationByMediumIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVA_Tight = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_Tight = ('hpsPFTauDiscriminationByTightIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVA_VTight = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_VTight = ('hpsPFTauDiscriminationByVTightIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVA_VVTight = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA_VVTight = ('hpsPFTauDiscriminationByVVTightIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVAraw = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVAraw = ('hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw')

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


runs = {}
lastlumi = -1
lastrun = -1

print ''
print 'n input files: ' + str(len(options.inputFiles))

for f in options.inputFiles:

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

    try:
        nevents = events.size()
    except:
        print 'skipping file ' + f
        continue

    print '### with ' + str(nevents) + ' events'
    print '### printing every ' + str(printevery) + '. event'

    if saveOutputFile: fout.cd()

    if isTest: nevents = neventsTest

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

        if ievent >= nevents:
            print 'nevents boundary'
            break
        if nEventsPerFile >= nMaxEventsPerFile:
            print 'nMaxEventsPerFile boundary'
            break

        if 'local' not in options.tag:
            if fin.IsZombie() or not fin.IsOpen():
                print 'file not usable'
                sys.exit(1)

        if saveOutputFile and ievent % 100 == 0: fout.Write('', ROOT.TObject.kWriteDelete)

        if ievent % printevery == 0: print 'analyzing event %d of %d' % (ievent, min(nevents, nMaxEventsPerFile))


        hCutflow.Fill(0)
        cutflow = 0

        random.seed()
        event_level_var_array['randomevent'][0] = random.randrange(10)

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
        ###########################################################################################veto

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

        hCutflow.Fill(1)
        cutflow = 1

        event_level_var_array['runnum'][0] = runnum
        event_level_var_array['lumisec'][0] = lumisec
        event_level_var_array['eventnum'][0] = eventnum

        numC1 = -1
        numN2 = -1
        numN1 = -1

        chipmmFILE = -1
        deltamFILE = -1

        chipmmGEN = -1
        deltamGEN = -1

        chipmptGEN = -1
        chipmetaGEN = -999
        chipmphiGEN = -999

        chiN2mGEN = -1
        deltamN2GEN = -1

        chiN2ptGEN = -1
        chiN2etaGEN = -999
        chiN2phiGEN = -999

        pionptGEN = -1
        pionetaGEN = -999
        pionphiGEN = -999
        pionchargeGEN = -999

        hasChargino = 0
        hasPion = 0
        hasMatchedTrack = 0

        decaylength3D = -1
        decaylengthXY = -1
        decaylengthZ = -1

        decaylength3DN2 = -1
        decaylengthXYN2 = -1
        decaylengthZN2 = -1

        chipmnumdaughters = 0
        chiN2numdaughters = 0
        numchidaughters = 0

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
        if not len(primaryvertices) > 0: continue
        pv_pos = primaryvertices[0].position()

        event.getByLabel(label_tracks, handle_tracks)
        tracks = handle_tracks.product()
        if not len(tracks) > 0: continue

        ###########################################################################################veto

        hCutflow.Fill(2)
        cutflow = 2

        '''
        ###############################################################################################
        # trigger
        ###############################################################################################
        '''

        trigger_flags_accept = {}
        for tf in trigger_flags:
            trigger_flags_accept[tf] = -1

        allfine = True

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

        ###########################################################################################veto

        hCutflow.Fill(3)
        cutflow = 3

        trigger_hlt_accept = {}
        for t_hlt in trigger_hlt:
            if t_hlt == 'triggerfired': continue
            trigger_hlt_accept[t_hlt] = -1

        triggerfired = 0

        if 'signal' not in options.tag and 'pmssm' not in options.tag:

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

        tauswithdiscriminators = [
            (tau, taudiscriminatorDM.value(itau)
                , taudiscriminatorMVAraw.value(itau)
                , taudiscriminatorMVA_VLoose.value(itau)
                , taudiscriminatorMVA_Loose.value(itau)
                , taudiscriminatorMVA_Medium.value(itau)
                , taudiscriminatorMVA_Tight.value(itau)
                , taudiscriminatorMVA_VTight.value(itau)
                , taudiscriminatorMVA_VVTight.value(itau)
                , taudiscriminatorElectronRej.value(itau)
                , taudiscriminatorMuonRej.value(itau))
            for itau, tau in enumerate(taus)
        ]


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
                                  if t[3] > 0.5  # vloose working point
                                  and t[9] > 0.5  # electron rejection
                                  and t[10] > 0.5  # muon rejection
                                  and t[0].pt() > 20
                                  and abs(t[0].eta()) < 2.3]

        '''
        ###############################################################################################
        # GEN info
        ###############################################################################################
        '''

        numZgamma = 0
        ptZgamma = -1
        etaZgamma = -999
        phiZgamma = -999
        numZgammaDaughters = 0
        ptsumZgammaNeutrinos = -1
        if 'geninfoZ' in options.tag:

            Zgammas = [gp for gp in genparticles if gp.isLastCopy() and gp.statusFlags().fromHardProcess()
                       and (abs(gp.pdgId()) == 22 or abs(gp.pdgId()) == 23)]

            numZgamma = len(Zgammas)

            if numZgamma == 1:  # else there's something funny going on...

                Zgamma = Zgammas[0]

                ptZgamma = Zgamma.pt()
                etaZgamma = Zgamma.eta()
                phiZgamma = Zgamma.phi()

                numZgammaDaughters = Zgamma.numberOfDaughters()

                ZgammaNeutrinos = ROOT.TLorentzVector()

                for i in range(numZgammaDaughters):

                    daughter = Zgamma.daughter(i)
                    if not abs(daughter.pdgId()) == 15:
                        daughter = getLastCopyStatusOne(daughter)
                    else:
                        daughter = getLastCopy(daughter)

                    if daughter is None: continue

                    zdaughter_var_array['pdgIdZdaughter'][i] = daughter.pdgId()
                    zdaughter_var_array['ptZdaughter'][i] = daughter.pt()
                    zdaughter_var_array['etaZdaughter'][i] = daughter.eta()
                    zdaughter_var_array['phiZdaughter'][i] = daughter.phi()

                    if abs(daughter.pdgId()) == 12 or abs(daughter.pdgId()) == 14 or abs(daughter.pdgId()) == 16:
                        nTlv = ROOT.TLorentzVector(daughter.px(), daughter.py(), daughter.pz(), daughter.energy())
                        ZgammaNeutrinos += nTlv

                ptsumZgammaNeutrinos = ZgammaNeutrinos.Pt()

        event_level_var_array['numZgamma'][0] = numZgamma
        event_level_var_array['ptZgamma'][0] = ptZgamma
        event_level_var_array['etaZgamma'][0] = etaZgamma
        event_level_var_array['phiZgamma'][0] = phiZgamma
        event_level_var_array['numZgammaDaughters'][0] = numZgammaDaughters
        event_level_var_array['ptsumZgammaNeutrinos'][0] = ptsumZgammaNeutrinos


        numW = 0
        ptW = -1
        etaW = -999
        phiW = -999
        numWDaughters = 0
        ptWneutrino = -1
        decayWtau = -1
        thetau = None
        if 'geninfoW' in options.tag:

            Ws = [gp for gp in genparticles if gp.isLastCopy() and gp.statusFlags().fromHardProcess()
                  and abs(gp.pdgId()) == 24]

            numW = len(Ws)

            if numW == 1:  # else there's something funny going on...

                W = Ws[0]

                ptW = W.pt()
                etaW = W.eta()
                phiW = W.phi()

                numWDaughters = W.numberOfDaughters()

                ptWneutrino = -1
                decayWtau = -1

                for i in range(numWDaughters):

                    daughter = W.daughter(i)
                    if not abs(daughter.pdgId()) == 15:
                        daughter = getLastCopyStatusOne(daughter)
                    else:
                        daughter = getLastCopy(daughter)

                    if daughter is None: continue

                    wdaughter_var_array['pdgIdWdaughter'][i] = daughter.pdgId()
                    wdaughter_var_array['ptWdaughter'][i] = daughter.pt()
                    wdaughter_var_array['etaWdaughter'][i] = daughter.eta()
                    wdaughter_var_array['phiWdaughter'][i] = daughter.phi()

                    if abs(daughter.pdgId()) == 15:

                        thetau = daughter

                        nprongs = 0
                        nneutral = 0

                        for k in range(daughter.numberOfDaughters()):

                            taudaughterpdgid = abs(daughter.daughter(k).pdgId())

                            if taudaughterpdgid in [12, 14, 16]:  # skip the neutrinos
                                continue
                            elif taudaughterpdgid in [11, 13]:  # charged lepton
                                decayWtau = 10 + taudaughterpdgid
                            elif taudaughterpdgid == 111:  # neutral pion
                                nneutral += 1
                            else:  # has to be a charged hadron
                                nprongs += 1

                        if decayWtau not in [21, 23]:
                            decayWtau = 5 * (nprongs - 1) + nneutral

                    if abs(daughter.pdgId()) == 12 or abs(daughter.pdgId()) == 14 or abs(daughter.pdgId()) == 16:
                        ptWneutrino = daughter.pt()

        event_level_var_array['numW'][0] = numW
        event_level_var_array['ptW'][0] = ptW
        event_level_var_array['etaW'][0] = etaW
        event_level_var_array['phiW'][0] = phiW
        event_level_var_array['numWDaughters'][0] = numWDaughters
        event_level_var_array['ptWneutrino'][0] = ptWneutrino
        event_level_var_array['decayWtau'][0] = decayWtau


        '''
        ###############################################################################################
        # jetID and JECs
        ###############################################################################################
        '''

        jetsP4Raw = []
        jetsP4Corr = []
        jetsIdxGood = []
        numBadJets = 0
        minetaabsbadjets = 10
        numBadJetsEventVeto = 0
        numBadJetsLepVeto = 0
        minetaabsbadjetsLepVeto = 10
        numBadJetsLepVetoEventVeto = 0
        for ijet, jet in enumerate(jets):

            jetP4Raw = ROOT.TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy())
            jetsP4Raw.append(jetP4Raw)

            if 'data' in options.tag:
                correction = getJEC(DataJECs.jecAK4(runnum), jetP4Raw, jet.jetArea(), rho, len(primaryvertices))
            else:
                correction = getJEC(jecAK4, jetP4Raw, jet.jetArea(), rho, len(primaryvertices))

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
                    jetIDhistos[v[0] + 'all'].Fill(globals()[v[0].replace('jet','')])
                    if goodJet: jetIDhistos[v[0] + 'pass'].Fill(globals()[v[0].replace('jet','')])

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

        ###########################################################################################veto

        hCutflow.Fill(4)
        cutflow = 4

        event_level_var_array['numbadjets'][0] = numBadJets
        event_level_var_array['minetaabsbadjets'][0] = minetaabsbadjets
        event_level_var_array['numbadjetsEventVeto'][0] = numBadJetsEventVeto

        event_level_var_array['numbadjetsLepVeto'][0] = numBadJetsLepVeto
        event_level_var_array['minetaabsbadjetsLepVeto'][0] = minetaabsbadjetsLepVeto
        event_level_var_array['numbadjetsLepVetoEventVeto'][0] = numBadJetsLepVetoEventVeto

        jets = [j for ij, j in enumerate(jets) if ij in jetsIdxGood and abs(j.eta()) < 5.]

        if not len(jets) > 0: continue

        ###########################################################################################veto

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
        l1eta = -999
        l2eta = -999
        l1phi = -999
        l2phi = -999
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
                sys.exit(0)

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

        ###########################################################################################veto

        hCutflow.Fill(5)
        cutflow = 5

        event_level_var_array['electronsCleaned'][0] = electronsCleaned
        event_level_var_array['muonsCleaned'][0] = muonsCleaned
        event_level_var_array['invmCleaning'][0] = invm
        event_level_var_array['zptCleaning'][0] = zpt
        event_level_var_array['l1ptCleaning'][0] = l1pt
        event_level_var_array['l2ptCleaning'][0] = l2pt
        event_level_var_array['l1etaCleaning'][0] = l1eta
        event_level_var_array['l2etaCleaning'][0] = l2eta
        event_level_var_array['l1phiCleaning'][0] = l1phi
        event_level_var_array['l2phiCleaning'][0] = l2phi
        event_level_var_array['l1absisodbetaCleaning'][0] = l1absisodbeta
        event_level_var_array['l1relisodbetaCleaning'][0] = l1relisodbeta
        event_level_var_array['l2absisodbetaCleaning'][0] = l2absisodbeta
        event_level_var_array['l2relisodbetaCleaning'][0] = l2relisodbeta
        event_level_var_array['metptBeforeCleaning'][0] = metptBeforeCleaning
        event_level_var_array['metphiBeforeCleaning'][0] = metphiBeforeCleaning

        '''
        ###############################################################################################
        # GEN MET and HT(miss) and FastSim MET correction
        ###############################################################################################
        '''

        pTneutrinosum = 0
        genmetpt = -1
        genmetphi = -10
        genht = 0
        genhtmiss = -1
        if 'data' not in options.tag:

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
                if abs(genjet.eta()) < 2.4 and genjet.pt() > 30: genht += genjet.pt()
                if abs(genjet.eta()) < 5 and genjet.pt() > 30:
                    genjetTlv = ROOT.TLorentzVector(genjet.px(), genjet.py(), genjet.pz(), genjet.energy())
                    genhtmissTlv -= genjetTlv
            genhtmiss = genhtmissTlv.Pt()

        nofastsimcorrmetpt = met.pt()
        nofastsimcorrmetphi = met.phi()

        if 'signal' in options.tag:

            met.setP4(ROOT.Math.LorentzVector('ROOT::Math::PxPyPzE4D<double>')(0.5 * (genmet.px() + met.px())
                                                                               , 0.5 * (genmet.py() + met.py())
                                                                               , 0
                                                                               , 0.5 * (genmet.energy() + met.energy())))

        event_level_var_array['pTneutrinosum'][0] = pTneutrinosum
        event_level_var_array['genmetpt'][0] = genmetpt
        event_level_var_array['genmetphi'][0] = genmetphi
        event_level_var_array['genht'][0] = genht
        event_level_var_array['genhtmiss'][0] = genhtmiss
        event_level_var_array['nofastsimcorrmetpt'][0] = nofastsimcorrmetpt
        event_level_var_array['nofastsimcorrmetphi'][0] = nofastsimcorrmetphi
        event_level_var_array['metphi'][0] = met.phi()
        event_level_var_array['metpt'][0] = met.pt()

        hMetpt.Fill(met.pt())

        '''
        ###############################################################################################
        # get GEN info for signal and match tracks
        ###############################################################################################
        '''

        matchedTrackIdxCharginoPion1 = -1
        matchedTrackIdxCharginoPion2 = -1

        tminmatching = -1
        dxyzmin = -1
        drmin = -1
        dxyzminrandom = -1
        drminrandom = -1

        drminold = -1
        drminoldrandom = -1

        susytracks = {}

        if 'signal' in options.tag:

            chipmmFILE = float(re.search(r'mChipm(.*?)GeV', f).group(1))
            deltamFILE = float(re.search(r'dm(.*?)GeV', f).group(1).replace('p','.'))

            C1s = [gp for gp in genparticles if gp.isLastCopy() and abs(gp.pdgId()) == 1000024]
            N2s = [gp for gp in genparticles if gp.isLastCopy() and abs(gp.pdgId()) == 1000023]
            N1s = [gp for gp in genparticles if gp.isLastCopy() and abs(gp.pdgId()) == 1000022]

            numC1 = len(C1s)
            numN2 = len(N2s)
            numN1 = len(N1s)

            c1daughters = []
            n2daughters = []

            for gp in C1s:

                decaylength3D = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                + pow(gp.vy() - gp.daughter(0).vy(), 2)
                                                + pow(gp.vz() - gp.daughter(0).vz(), 2))

                if not decaylength3D > 0: continue

                hasChargino += 1

                chipmmGEN = round(gp.mass(), 2)

                deltamGEN = round((gp.mass() - gp.daughter(0).mass()), 2)

                chipmptGEN = gp.pt()
                chipmetaGEN = gp.eta()
                chipmphiGEN = gp.phi()

                decaylengthXY = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                + pow(gp.vy() - gp.daughter(0).vy(), 2))

                decaylengthZ = abs(gp.vz() - gp.daughter(0).vz())

                c1daughters += findDaughters(gp)

            for c1d in c1daughters:

                if abs(c1d.pdgId()) == 211:

                    pion = c1d

                    hasPion += 1

                    chidaughterpdgid = pion.pdgId()

                    pionptGEN = pion.pt()
                    pionetaGEN = pion.eta()
                    pionphiGEN = pion.phi()
                    pionchargeGEN = pion.charge()

                    idxold, drminold = findMatch_track_old(pion, tracks)
                    _, drminoldrandom = findMatch_track_old_random(pion, tracks)

                    idx, dxyzmin, tminmatching, drmin = findMatch_track_new(pion, tracks)
                    _, dxyzminrandom, _, drminrandom = findMatch_track_new_random(pion, tracks)

                    if not idx == -1:
                        if drmin < matchingDrThreshold and dxyzmin < matchingDxyzThreshold:
                            if not hasMatchedTrack: matchedTrackIdxCharginoPion1 = idx
                            else: matchedTrackIdxCharginoPion2 = idx
                            hasMatchedTrack += 1

            for gp in N2s:

                decaylength3DN2 = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                  + pow(gp.vy() - gp.daughter(0).vy(), 2)
                                                  + pow(gp.vz() - gp.daughter(0).vz(), 2))

                if not decaylength3DN2 > 0: continue

                chiN2mGEN = round(gp.mass(), 2)

                deltamN2GEN = round((gp.mass() - gp.daughter(0).mass()), 2)

                chiN2ptGEN = gp.pt()
                chiN2etaGEN = gp.eta()
                chiN2phiGEN = gp.phi()

                decaylengthXYN2 = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(), 2)
                                                  + pow(gp.vy() - gp.daughter(0).vy(), 2))

                decaylengthZN2 = abs(gp.vz() - gp.daughter(0).vz())

                n2daughters += findDaughters(gp)

            chipmnumdaughters = len(c1daughters)
            chiN2numdaughters = len(n2daughters)
            numchidaughters = chipmnumdaughters + chiN2numdaughters

            i = 0

            for c1d in c1daughters:

                chidaughter_var_array['motherchidaughter'][i] = 1
                chidaughter_var_array['pdgIdchidaughter'][i] = c1d.pdgId()
                chidaughter_var_array['ptchidaughter'][i] = c1d.pt()
                chidaughter_var_array['etachidaughter'][i] = c1d.eta()
                chidaughter_var_array['phichidaughter'][i] = c1d.phi()

                chidaughter_var_array['hasmatchedtrackchidaughter'][i] = -1

                if c1d.charge() != 0:

                    idxC1, dxyzminC1, _, drminC1 = findMatch_track_new(c1d, tracks)

                    if not idxC1 == -1:
                        if drminC1 < matchingDrThreshold and dxyzminC1 < matchingDxyzThreshold:
                            susytracks[idxC1] = (1, c1d.pdgId())
                            chidaughter_var_array['hasmatchedtrackchidaughter'][i] = idxC1

                i += 1

            for n2d in n2daughters:

                chidaughter_var_array['motherchidaughter'][i] = 2
                chidaughter_var_array['pdgIdchidaughter'][i] = n2d.pdgId()
                chidaughter_var_array['ptchidaughter'][i] = n2d.pt()
                chidaughter_var_array['etachidaughter'][i] = n2d.eta()
                chidaughter_var_array['phichidaughter'][i] = n2d.phi()

                chidaughter_var_array['hasmatchedtrackchidaughter'][i] = -1

                if n2d.charge() != 0:

                    idxN2, dxyzminN2, _, drminN2 = findMatch_track_new(n2d, tracks)

                    if not idxN2 == -1:
                        if drminN2 < matchingDrThreshold and dxyzminN2 < matchingDxyzThreshold:
                            susytracks[idxN2] = (2, n2d.pdgId())
                            chidaughter_var_array['hasmatchedtrackchidaughter'][i] = idxN2

                i += 1

        event_level_var_array['chipmmFILE'][0] = chipmmFILE
        event_level_var_array['deltamFILE'][0] = deltamFILE

        event_level_var_array['numC1'][0] = numC1
        event_level_var_array['numN2'][0] = numN2
        event_level_var_array['numN1'][0] = numN1

        event_level_var_array['chipmmGEN'][0] = chipmmGEN
        event_level_var_array['deltamGEN'][0] = deltamGEN

        event_level_var_array['chiN2mGEN'][0] = chiN2mGEN
        event_level_var_array['deltamN2GEN'][0] = deltamN2GEN

        event_level_var_array['chipmptGEN'][0] = chipmptGEN
        event_level_var_array['chipmetaGEN'][0] = chipmetaGEN
        event_level_var_array['chipmphiGEN'][0] = chipmphiGEN

        event_level_var_array['chiN2ptGEN'][0] = chiN2ptGEN
        event_level_var_array['chiN2etaGEN'][0] = chiN2etaGEN
        event_level_var_array['chiN2phiGEN'][0] = chiN2phiGEN

        event_level_var_array['pionptGEN'][0] = pionptGEN
        event_level_var_array['pionetaGEN'][0] = pionetaGEN
        event_level_var_array['pionphiGEN'][0] = pionphiGEN
        event_level_var_array['pionchargeGEN'][0] = pionchargeGEN

        event_level_var_array['hasChargino'][0] = hasChargino
        event_level_var_array['hasPion'][0] = hasPion
        event_level_var_array['hasMatchedTrack'][0] = hasMatchedTrack

        event_level_var_array['tminmatching'][0] = tminmatching
        event_level_var_array['dxyzmin'][0] = dxyzmin
        event_level_var_array['drmin'][0] = drmin
        event_level_var_array['dxyzminrandom'][0] = dxyzminrandom
        event_level_var_array['drminrandom'][0] = drminrandom
        event_level_var_array['drminold'][0] = drminold
        event_level_var_array['drminoldrandom'][0] = drminoldrandom

        event_level_var_array['matchedTrackIdxCharginoPion1'][0] = matchedTrackIdxCharginoPion1
        event_level_var_array['matchedTrackIdxCharginoPion2'][0] = matchedTrackIdxCharginoPion2

        event_level_var_array['chidecaylength3D'][0] = decaylength3D
        event_level_var_array['chidecaylengthXY'][0] = decaylengthXY
        event_level_var_array['chidecaylengthZ'][0] = decaylengthZ
        if decaylength3D > 0:
            event_level_var_array['log10(chidecaylength3D)'][0] = ROOT.TMath.Log10(decaylength3D)
            event_level_var_array['log10(chidecaylengthXY)'][0] = ROOT.TMath.Log10(decaylengthXY)
            event_level_var_array['log10(chidecaylengthZ)'][0] = ROOT.TMath.Log10(decaylengthZ)
        else:
            event_level_var_array['log10(chidecaylength3D)'][0] = -1
            event_level_var_array['log10(chidecaylengthXY)'][0] = -1
            event_level_var_array['log10(chidecaylengthZ)'][0] = -1

        event_level_var_array['chidecaylength3DN2'][0] = decaylength3DN2
        event_level_var_array['chidecaylengthXYN2'][0] = decaylengthXYN2
        event_level_var_array['chidecaylengthZN2'][0] = decaylengthZN2
        if decaylength3DN2 > 0:
            event_level_var_array['log10(chidecaylength3DN2)'][0] = ROOT.TMath.Log10(decaylength3DN2)
            event_level_var_array['log10(chidecaylengthXYN2)'][0] = ROOT.TMath.Log10(decaylengthXYN2)
            event_level_var_array['log10(chidecaylengthZN2)'][0] = ROOT.TMath.Log10(decaylengthZN2)
        else:
            event_level_var_array['log10(chidecaylength3DN2)'][0] = -1
            event_level_var_array['log10(chidecaylengthXYN2)'][0] = -1
            event_level_var_array['log10(chidecaylengthZN2)'][0] = -1

        event_level_var_array['chipmnumdaughters'][0] = chipmnumdaughters
        event_level_var_array['chiN2numdaughters'][0] = chiN2numdaughters
        event_level_var_array['numchidaughters'][0] = numchidaughters


        '''
        ###############################################################################################
        # get event-level info
        ###############################################################################################
        '''

        event_level_var_array['numpvs'][0] = len(primaryvertices)
        event_level_var_array['rho'][0] = rho


        numjets = len(jets)
        numjets30 = 0
        numjets50 = 0
        numjets100 = 0
        numjets200 = 0
        ht = 0
        ht5 = 0
        htmissTlv = ROOT.TLorentzVector()
        idxhighestptjet = 0
        for ijet, jet in enumerate(jets):

            if abs(jet.eta()) < 2.4 and jet.pt() > 30: ht += jet.pt()
            if abs(jet.eta()) < 5. and jet.pt() > 30:
                ht5 += jet.pt()
                jetTlv = ROOT.TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy())
                htmissTlv -= jetTlv

            jetpt = jet.pt()
            if jetpt > 30: numjets30 += 1
            if jetpt > 50: numjets50 += 1
            if jetpt > 100: numjets100 += 1
            if jetpt > 200: numjets200 += 1

            if jetpt > jets[idxhighestptjet].pt(): idxhighestptjet = ijet

        if len(jets) > 0:
            event_level_var_array['ptleadingjet'][0] = jets[idxhighestptjet].pt()
            event_level_var_array['etaleadingjet'][0] = jets[idxhighestptjet].eta()
            event_level_var_array['phileadingjet'][0] = jets[idxhighestptjet].phi()
        else:
            event_level_var_array['ptleadingjet'][0] = -1
            event_level_var_array['etaleadingjet'][0] = -1
            event_level_var_array['phileadingjet'][0] = -1
        event_level_var_array['numjets'][0] = numjets
        event_level_var_array['numjets30'][0] = numjets30
        event_level_var_array['numjets50'][0] = numjets50
        event_level_var_array['numjets100'][0] = numjets100
        event_level_var_array['numjets200'][0] = numjets200
        event_level_var_array['ht'][0] = ht
        event_level_var_array['ht5'][0] = ht5
        event_level_var_array['htmiss'][0] = htmissTlv.Pt()

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
            event_level_var_array['mindphimetjets'][0] = min(dphimetjets)
            hMindphimetjets.Fill(min(dphimetjets))
        else:
            event_level_var_array['mindphimetjets'][0] = -1
            hMindphimetjets.Fill(-1)


        hNPVsPerEvent.Fill(len(primaryvertices))

        hPV0x.Fill(primaryvertices[0].x())
        hPV0y.Fill(primaryvertices[0].y())
        hPV0z.Fill(primaryvertices[0].z())

        for i in range(1, len(primaryvertices)):

            hPVsx.Fill(primaryvertices[i].x())
            hPVsy.Fill(primaryvertices[i].y())
            hPVsz.Fill(primaryvertices[i].z())


        '''
        ###############################################################################################
        # event selection
        ###############################################################################################
        '''

        if not met.pt() > metthreshold: continue

        ###########################################################################################veto

        hCutflow.Fill(6)
        cutflow = 6

        if not numjets100 > 0: continue

        ###########################################################################################veto

        hCutflow.Fill(7)
        cutflow = 7

        if 'nodphimetjetsveto' not in options.tag:
            if len(dphimetjets) > 0:
                if not min(dphimetjets) > 0.5: continue
            else:
                continue

        ###########################################################################################veto

        hCutflow.Fill(8)
        cutflow = 8


        # chpfcandsforiso = [p for p in pfcands if passesPreselection_iso_pfc(p, pv_pos, dz_threshold=0.1)]
        chpfcandsforiso = np.array([(p.pt(), p.eta(), p.phi()) for p in pfcands if passesPreselection_iso_pfc(p, pv_pos, dz_threshold=0.1)])
        pfcandsforiso = np.array([(p.pt(), p.eta(), p.phi()) for p in pfcands])
        jetsforiso = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets])


        event_level_var_array['numelectrons'][0] = len(electrons)

        numelectronsiso = 0
        for ie, e in enumerate(electrons):

            electron_var_array['chargeelectron'][ie] = e.charge()
            electron_var_array['pxelectron'][ie] = e.px()
            electron_var_array['pyelectron'][ie] = e.py()
            electron_var_array['pzelectron'][ie] = e.pz()
            electron_var_array['ptelectron'][ie] = e.pt()
            electron_var_array['energyelectron'][ie] = e.energy()
            electron_var_array['etaelectron'][ie] = e.eta()
            electron_var_array['phielectron'][ie] = e.phi()
            electron_var_array['dzelectron'][ie] = abs(e.gsfTrack().dz(pv_pos))
            electron_var_array['dxyelectron'][ie] = abs(e.gsfTrack().dxy(pv_pos))
            electron_var_array['pfabsisoelectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, pfcandsforiso)
            electron_var_array['pfabsisominielectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, pfcandsforiso, isMini=True)
            electron_var_array['chpfabsisoelectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, chpfcandsforiso)
            electron_var_array['chpfabsisominielectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, chpfcandsforiso, isMini=True)
            electron_var_array['jetisoelectron'][ie], electron_var_array['jetisomultielectron'][ie], electron_var_array['jetdrminelectron'][ie], _ = calcIso_jet_new(e, jetsforiso, isTrack=False, btagvalues=btagvalues)
            electron_var_array['chhadisoelectron'][ie] = e.pfIsolationVariables().sumChargedHadronPt
            electron_var_array['challisoelectron'][ie] = e.pfIsolationVariables().sumChargedParticlePt
            electron_var_array['neuhadisoelectron'][ie] = e.pfIsolationVariables().sumNeutralHadronEt
            electron_var_array['photisoelectron'][ie] = e.pfIsolationVariables().sumPhotonEt
            electron_var_array['puchhadisoelectron'][ie] = e.pfIsolationVariables().sumPUPt
            electron_var_array['absisodbetaelectron'][ie] = calcIso_dBeta(e.pfIsolationVariables())
            electron_var_array['relisodbetaelectron'][ie] = calcIso_dBeta(e.pfIsolationVariables()) / e.pt()

            if calcIso_dBeta(e.pfIsolationVariables()) / e.pt() < 0.2: numelectronsiso += 1

        event_level_var_array['numelectronsiso'][0] = numelectronsiso


        event_level_var_array['nummuons'][0] = len(muons)

        nummuonsiso = 0
        for im, m in enumerate(muons):

            muon_var_array['chargemuon'][im] = m.charge()
            muon_var_array['pxmuon'][im] = m.px()
            muon_var_array['pymuon'][im] = m.py()
            muon_var_array['pzmuon'][im] = m.pz()
            muon_var_array['ptmuon'][im] = m.pt()
            muon_var_array['energymuon'][im] = m.energy()
            muon_var_array['etamuon'][im] = m.eta()
            muon_var_array['phimuon'][im] = m.phi()
            muon_var_array['dzmuon'][im] = abs(m.muonBestTrack().dz(pv_pos))
            muon_var_array['dxymuon'][im] = abs(m.muonBestTrack().dxy(pv_pos))
            muon_var_array['pfabsisomuon'][im], _, _, _ = calcIso_pf_or_track_new(m, pfcandsforiso)
            muon_var_array['pfabsisominimuon'][im], _, _, _ = calcIso_pf_or_track_new(m, pfcandsforiso, isMini=True)
            muon_var_array['chpfabsisomuon'][im], _, _, _ = calcIso_pf_or_track_new(m, chpfcandsforiso)
            muon_var_array['chpfabsisominimuon'][im], _, _, _ = calcIso_pf_or_track_new(m, chpfcandsforiso, isMini=True)
            muon_var_array['jetisomuon'][im], muon_var_array['jetisomultimuon'][im], muon_var_array['jetdrminmuon'][im], _ = calcIso_jet_new(m, jetsforiso, isTrack=False, btagvalues=btagvalues)
            muon_var_array['chhadisomuon'][im] = m.pfIsolationR03().sumChargedHadronPt
            muon_var_array['challisomuon'][im] = m.pfIsolationR03().sumChargedParticlePt
            muon_var_array['neuhadisomuon'][im] = m.pfIsolationR03().sumNeutralHadronEt
            muon_var_array['photisomuon'][im] = m.pfIsolationR03().sumPhotonEt
            muon_var_array['puchhadisomuon'][im] = m.pfIsolationR03().sumPUPt
            muon_var_array['absisodbetamuon'][im] = calcIso_dBeta(m.pfIsolationR03())
            muon_var_array['relisodbetamuon'][im] = calcIso_dBeta(m.pfIsolationR03()) / m.pt()

            if calcIso_dBeta(m.pfIsolationR03()) / m.pt() < 0.2: nummuonsiso += 1

        event_level_var_array['nummuonsiso'][0] = nummuonsiso

        event_level_var_array['numleptons'][0] = len(electrons) + len(muons)
        event_level_var_array['numleptonsiso'][0] = numelectronsiso + nummuonsiso

        hNumleptons.Fill(numelectronsiso+nummuonsiso)

        if 'noleptonveto' not in options.tag:
            if not (numelectronsiso+nummuonsiso) == 0: continue

        ###########################################################################################veto

        hCutflow.Fill(9)
        cutflow = 9


        event_level_var_array['numgenjets'][0] = 0

        if 'data' not in options.tag:

            event_level_var_array['numgenjets'][0] = len(genjets)

            for igj, gj in enumerate(genjets):

                genjet_var_array['ptgenjet'][igj] = gj.pt()
                genjet_var_array['etagenjet'][igj] = gj.eta()
                genjet_var_array['phigenjet'][igj] = gj.phi()
                genjet_var_array['massgenjet'][igj] = gj.mass()


        event_level_var_array['numphotons'][0] = len(photons)

        numphotonsiso = 0
        for ip, p in enumerate(photons):

            photon_var_array['pxphoton'][ip] = p.px()
            photon_var_array['pyphoton'][ip] = p.py()
            photon_var_array['pzphoton'][ip] = p.pz()
            photon_var_array['ptphoton'][ip] = p.pt()
            photon_var_array['energyphoton'][ip] = p.energy()
            photon_var_array['etaphoton'][ip] = p.eta()
            photon_var_array['phiphoton'][ip] = p.phi()
            photon_var_array['pfabsisophoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, pfcandsforiso, dontSubtractObject=True)
            photon_var_array['pfabsisominiphoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, pfcandsforiso, isMini=True, dontSubtractObject=True)
            photon_var_array['chpfabsisophoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, chpfcandsforiso, dontSubtractObject=True)
            photon_var_array['chpfabsisominiphoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, chpfcandsforiso, isMini=True, dontSubtractObject=True)
            photon_var_array['jetisophoton'][ip], photon_var_array['jetisomultiphoton'][ip], photon_var_array['jetdrminphoton'][ip], _ = calcIso_jet_new(p, jetsforiso, isTrack=False, btagvalues=btagvalues)
            photon_var_array['chhadisophoton'][ip] = p.chargedHadronIso()
            photon_var_array['neuhadisophoton'][ip] = p.neutralHadronIso()
            photon_var_array['photisophoton'][ip] = p.photonIso()
            absisophoton = p.chargedHadronIso() + p.neutralHadronIso() + p.photonIso()
            photon_var_array['absisophoton'][ip] = absisophoton
            relisophoton = absisophoton / p.pt()
            photon_var_array['relisophoton'][ip] = relisophoton

            if relisophoton < 0.2: numphotonsiso += 1

        event_level_var_array['numphotonsiso'][0] = numphotonsiso
        hNumphotons.Fill(numphotonsiso)


        pfleptons = [l for l in pfcands
                     if ((abs(l.pdgId()) == 11 or abs(l.pdgId()) == 13)
                         and l.pt() > 10
                         and abs(l.eta()) < 2.4
                         and not (1.4442 < abs(l.eta()) < 1.566))]

        event_level_var_array['numpfleptons'][0] = len(pfleptons)

        numpfleptonsiso = 0
        for il, l in enumerate(pfleptons):

            pflepton_var_array['chargepflepton'][il] = l.charge()
            pflepton_var_array['pdgidpflepton'][il] = l.pdgId()
            pflepton_var_array['pxpflepton'][il] = l.px()
            pflepton_var_array['pypflepton'][il] = l.py()
            pflepton_var_array['pzpflepton'][il] = l.pz()
            pflepton_var_array['ptpflepton'][il] = l.pt()
            pflepton_var_array['energypflepton'][il] = l.energy()
            pflepton_var_array['etapflepton'][il] = l.eta()
            pflepton_var_array['phipflepton'][il] = l.phi()
            if l.trackRef().isNull():
                pflepton_var_array['dzpflepton'][il] = -1
                pflepton_var_array['dxypflepton'][il] = -1
            else:
                pflepton_var_array['dzpflepton'][il] = abs(l.trackRef().get().dz(pv_pos))
                pflepton_var_array['dxypflepton'][il] = abs(l.trackRef().get().dxy(pv_pos))
            _, pfrelisopflepton, _, _ = calcIso_pf_or_track_new(l, pfcandsforiso)
            pflepton_var_array['pfrelisopflepton'][il] = pfrelisopflepton
            _, pflepton_var_array['pfrelisominipflepton'][il], _, _ = calcIso_pf_or_track_new(l, pfcandsforiso, isMini=True)
            _, pflepton_var_array['chpfrelisopflepton'][il], _, _ = calcIso_pf_or_track_new(l, chpfcandsforiso)
            _, pflepton_var_array['chpfrelisominipflepton'][il], _, _ = calcIso_pf_or_track_new(l, chpfcandsforiso, isMini=True)
            pflepton_var_array['jetisopflepton'][il], pflepton_var_array['jetisomultipflepton'][il], pflepton_var_array['jetdrminpflepton'][il], _ = calcIso_jet_new(l, jetsforiso, isTrack=False, btagvalues=btagvalues)

            if pfrelisopflepton < 0.2: numpfleptonsiso += 1

        event_level_var_array['numpfleptonsiso'][0] = numpfleptonsiso


        event_level_var_array['numtaus'][0] = len(tauswithdiscriminators)

        if len(tauswithdiscriminators) > 0:
            genparticlesfortaumatching = [gp for gp in genparticles if gp.status() == 1 and abs(gp.pdgId()) in [11, 13, 15]]

        numtausvloose = 0
        numtausloose = 0
        numtausmedium = 0
        numtaustight = 0
        numtausvtight = 0
        numtausvvtight = 0
        for it, t in enumerate(tauswithdiscriminators):

            tau_var_array['chargetau'][it] = t[0].charge()
            tau_var_array['pxtau'][it] = t[0].px()
            tau_var_array['pytau'][it] = t[0].py()
            tau_var_array['pztau'][it] = t[0].pz()
            tau_var_array['pttau'][it] = t[0].pt()
            tau_var_array['energytau'][it] = t[0].energy()
            tau_var_array['etatau'][it] = t[0].eta()
            tau_var_array['phitau'][it] = t[0].phi()
            if t[0].leadPFChargedHadrCand().trackRef().get() and t[0].leadPFChargedHadrCand().trackRef().isNonnull():
                leadpfchhadcand = t[0].leadPFChargedHadrCand().trackRef().get()
                tau_var_array['dztau'][it] = leadpfchhadcand.dz(pv_pos)
                tau_var_array['dxytau'][it] = leadpfchhadcand.dxy(pv_pos)
                tau_var_array['tauleadpfchhadcandpt'][it] = leadpfchhadcand.pt()
                tau_var_array['tauleadpfchhadcandeta'][it] = leadpfchhadcand.eta()
                tau_var_array['tauleadpfchhadcandphi'][it] = leadpfchhadcand.phi()
            else:
                tau_var_array['dztau'][it] = -1
                tau_var_array['dxytau'][it] = -1
                tau_var_array['tauleadpfchhadcandpt'][it] = 0.
                tau_var_array['tauleadpfchhadcandeta'][it] = 0.
                tau_var_array['tauleadpfchhadcandphi'][it] = 0.
            tau_var_array['chhadisotau'][it] = t[0].isolationPFChargedHadrCandsPtSum()
            tau_var_array['photisotau'][it] = t[0].isolationPFGammaCandsEtSum()
            tau_var_array['decaymodetau'][it] = t[0].decayMode()
            tau_var_array['decaymodefindingtau'][it] = t[1]
            tau_var_array['mvadiscrtau'][it] = t[2]

            if t[3] > 0.5:
                tau_var_array['isvloosetau'][it] = t[3]
                numtausvloose += 1
            if t[4] > 0.5:
                tau_var_array['isloosetau'][it] = t[4]
                numtausloose += 1
            if t[5] > 0.5:
                tau_var_array['ismediumtau'][it] = t[5]
                numtausmedium += 1
            if t[6] > 0.5:
                tau_var_array['istighttau'][it] = t[6]
                numtaustight += 1
            if t[7] > 0.5:
                tau_var_array['isvtighttau'][it] = t[7]
                numtausvtight += 1
            if t[8] > 0.5:
                tau_var_array['isvvtighttau'][it] = t[8]
                numtausvvtight += 1

            idxtaumatch, drmintaumatch = findMatch_gen_old_easy(t[0], genparticlesfortaumatching)
            tau_var_array['genmatchdrtau'][it] = drmintaumatch
            if drmintaumatch < 0.5:
                tau_var_array['genmatchpdgidtau'][it] = genparticlesfortaumatching[idxtaumatch].pdgId()
            else:
                tau_var_array['genmatchpdgidtau'][it] = -1.

        hNumtaus.Fill(numtausvloose)
        event_level_var_array['numtausvloose'][0] = numtausvloose
        event_level_var_array['numtausloose'][0] = numtausloose
        event_level_var_array['numtausmedium'][0] = numtausmedium
        event_level_var_array['numtaustight'][0] = numtaustight
        event_level_var_array['numtausvtight'][0] = numtausvtight
        event_level_var_array['numtausvvtight'][0] = numtausvvtight


        tracksByPV = {}
        for ipv, pv in enumerate(primaryvertices):

            pv_var_array['idxpv'][ipv] = ipv
            pv_var_array['numtrackspv'][ipv] = pv.tracksSize()
            pv_var_array['xpv'][ipv] = pv.position().x()
            pv_var_array['ypv'][ipv] = pv.position().y()
            pv_var_array['zpv'][ipv] = pv.position().z()
            tracksByPV[ipv] = [pv.trackRefAt(i).get() for i in range(pv.tracksSize())]


        btagvalues = []
        isleptonjetidxlist = []
        for ijet, jet in enumerate(jets):

            jet_var_array['ptjet'][ijet] = jet.pt()
            jet_var_array['etajet'][ijet] = jet.eta()
            jet_var_array['phijet'][ijet] = jet.phi()
            jet_var_array['pxjet'][ijet] = jet.px()
            jet_var_array['pyjet'][ijet] = jet.py()
            jet_var_array['pzjet'][ijet] = jet.pz()
            jet_var_array['energyjet'][ijet] = jet.energy()
            jet_var_array['nconstituentsjet'][ijet] = jet.numberOfDaughters()

            dRmin = 0.1
            thisbtag = -2.
            for ib in range(nbtags):
                dR = deltaR(handle_btag.product().key(ib).get().eta(), jet.eta(), handle_btag.product().key(ib).get().phi(), jet.phi())
                if dR < dRmin:
                    dRmin = dR
                    thisbtag = max(-1.0, handle_btag.product().value(ib))

            jet_var_array['btagjet'][ijet] = thisbtag
            btagvalues.append((thisbtag, jet.pt(), jet.eta()))
            hBtagjets.Fill(thisbtag)

            drminleptonjet = 10.
            ptclosestleptonjet = 0.
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

            jet_var_array['drminleptonjet'][ijet] = drminleptonjet
            jet_var_array['ptclosestleptonjet'][ijet] = ptclosestleptonjet
            jet_var_array['isleptonjet'][ijet] = isleptonjet

            drmingenjetjet = 10.
            ptclosestgenjetjet = 0.
            isgenjetjet = 0.

            if 'data' not in options.tag:

                idx_genjet, drmingenjetjet = findMatch_gen_old_easy(jet, genjets)

                if drmingenjetjet < 0.2:
                    ptclosestgenjetjet = genjets[idx_genjet].pt()
                    isgenjetjet = 1.

            jet_var_array['drmingenjetjet'][ijet] = drmingenjetjet
            jet_var_array['ptclosestgenjetjet'][ijet] = ptclosestgenjetjet
            jet_var_array['isgenjetjet'][ijet] = isgenjetjet

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

        event_level_var_array['njetsbtagloose'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > loosewp and jetpt > 30 and abs(jeteta) < 2.4)])
        event_level_var_array['njetsbtaglooseTIGHT'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > loosewp and jetpt > 15 and abs(jeteta) < 2.4)])

        njetsbtagmedium = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > mediumwp and jetpt > 30 and abs(jeteta) < 2.4)])
        hNjetsbtagmedium.Fill(njetsbtagmedium)
        event_level_var_array['njetsbtagmedium'][0] = njetsbtagmedium
        event_level_var_array['njetsbtagmediumTIGHT'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > mediumwp and jetpt > 15 and abs(jeteta) < 2.4)])
        
        event_level_var_array['njetsbtagtight'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > tightwp and jetpt > 30 and abs(jeteta) < 2.4)])
        event_level_var_array['njetsbtagtightTIGHT'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > tightwp and jetpt > 15 and abs(jeteta) < 2.4)])

        mtmetleadingjet = ROOT.TMath.Sqrt(2 * met.pt() * jets[idxhighestptjet].pt()
                                          * (1 - ROOT.TMath.Cos(deltaPhi(met.phi(), jets[idxhighestptjet].phi()))))
        hMtmetleadingjet.Fill(mtmetleadingjet)
        event_level_var_array['mtmetleadingjet'][0] = mtmetleadingjet


        if phifirsttrack == tracks[0].phi() and etafirsttrack == tracks[0].eta():
            print 'suspicious... better get out of here!'
            sys.exit(1)
        phifirsttrack = tracks[0].phi()
        etafirsttrack = tracks[0].eta()


        '''
        ###############################################################################################
        # get track-level info
        ###############################################################################################
        '''

        tracksforiso = np.array([(t.pt(), t.eta(), t.phi()) for t in tracks
                                 if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1, dxy_threshold=1000., pt_threshold=0.)])
        tracksforisotight = np.array([(t.pt(), t.eta(), t.phi()) for t in tracks
                                      if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1, dxy_threshold=0.1, pt_threshold=1.)])

        # jetsforisotightNoLepton = [j for ij, j in enumerate(jets) if passesPreselection_iso_jet(j, pt_threshold=30) and ij not in isleptonjetidxlist]
        # jetsforisotight = [j for j in jets if passesPreselection_iso_jet(j, pt_threshold=30)]
        # jetsforisomeditight = [j for j in jets if passesPreselection_iso_jet(j, pt_threshold=20)]
        # jetsforisomedium = [j for j in jets if passesPreselection_iso_jet(j, pt_threshold=15)]
        # jetsforisoloose = [j for j in jets if passesPreselection_iso_jet(j, pt_threshold=10)]

        jetsforisotightNoLepton = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for ij, j in enumerate(jets)
                                            if passesPreselection_iso_jet(j, pt_threshold=30) and ij not in isleptonjetidxlist])
        jetsforisotight = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                    if passesPreselection_iso_jet(j, pt_threshold=30)])
        jetsforisomeditight = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                        if passesPreselection_iso_jet(j, pt_threshold=20)])
        jetsforisomedium = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                     if passesPreselection_iso_jet(j, pt_threshold=15)])
        jetsforisoloose = np.array([(j.pt(), j.eta(), j.phi(), j.energy(), j.numberOfDaughters()) for j in jets
                                    if passesPreselection_iso_jet(j, pt_threshold=10)])

        neutralhadrons = [p for p in pfcands if p.particleId() == 5]

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
            jetisotight, jetisomultitight, jetdrmintight, jetisobtagtight = calcIso_jet_new(track, jetsforisotight, isTrack=True, btagvalues=btagvalues)
            if not jetdrmintight > 0.4: continue

            numtracksfinalpreselection += 1

            track_level_var_array['randomtrack'][i] = random.randrange(10)

            track_level_var_array['charge'][i] = track.charge()

            track_level_var_array['pxtrack'][i] = track.px()
            track_level_var_array['pytrack'][i] = track.py()
            track_level_var_array['pztrack'][i] = track.pz()

            track_level_var_array['pttrack'][i] = track.pt()
            track_level_var_array['pttrackerror/pttrack'][i] = track.ptError()/track.pt()
            track_level_var_array['log10(pttrackerror/pttrack)'][i] = ROOT.TMath.Log10(track.ptError()/track.pt())

            track_level_var_array['pttrackerror'][i] = track.ptError()
            track_level_var_array['log10(pttrackerror)'][i] = ROOT.TMath.Log10(track.ptError())

            track_level_var_array['etaerror'][i] = track.etaError()
            track_level_var_array['phierror'][i] = track.phiError()

            track_level_var_array['eta'][i] = track.eta()
            track_level_var_array['phi'][i] = track.phi()

            track_level_var_array['associatedpv'][i] = -1
            for ipv in range(len(primaryvertices)):
                if track in tracksByPV[ipv]: track_level_var_array['associatedpv'][i] = ipv

            ip = IPcalculator(track, primaryvertices[0])
            track_level_var_array['IPsignificance'][i] = ip.getIPsignificance()
            track_level_var_array['IPxyz'][i] = ip.getIP()
            track_level_var_array['IPxy'][i] = ip.getDxy()
            track_level_var_array['IPz'][i] = ip.getDz()
            track_level_var_array['log10(IPsignificance)'][i] = ROOT.TMath.Log10(ip.getIPsignificance())
            track_level_var_array['log10(IPxyz)'][i] = ROOT.TMath.Log10(ip.getIP())
            track_level_var_array['log10(IPxy)'][i] = ROOT.TMath.Log10(ip.getDxy())
            track_level_var_array['log10(IPz)'][i] = ROOT.TMath.Log10(ip.getDz())

            minipPU = None
            minivPU = -1
            minIPsignificancePU = 999
            for iv, v in enumerate(primaryvertices[1:]):
                thisipPU = IPcalculator(track, v)
                thisIPsignificancePU = thisipPU.getIPsignificance()
                if thisIPsignificancePU < minIPsignificancePU:
                    minipPU = thisipPU
                    minivPU = iv
                    minIPsignificancePU = thisIPsignificancePU

            if not minivPU == -1:
                track_level_var_array['idxpvPU'][i] = minivPU+1
                track_level_var_array['IPsignificancePU'][i] = minipPU.getIPsignificance()
                track_level_var_array['IPxyzPU'][i] = minipPU.getIP()
                track_level_var_array['IPxyPU'][i] = minipPU.getDxy()
                track_level_var_array['IPzPU'][i] = minipPU.getDz()
                track_level_var_array['log10(IPsignificancePU)'][i] = ROOT.TMath.Log10(minipPU.getIPsignificance())
                track_level_var_array['log10(IPxyzPU)'][i] = ROOT.TMath.Log10(minipPU.getIP())
                track_level_var_array['log10(IPxyPU)'][i] = ROOT.TMath.Log10(minipPU.getDxy())
                track_level_var_array['log10(IPzPU)'][i] = ROOT.TMath.Log10(minipPU.getDz())

            track_level_var_array['dxy0'][i] = abs(track.dxy())
            track_level_var_array['dz0'][i] = abs(track.dz())
            track_level_var_array['log10(dxy0)'][i] = ROOT.TMath.Log10(abs(track.dxy()))
            track_level_var_array['log10(dz0)'][i] = ROOT.TMath.Log10(abs(track.dz()))

            track_level_var_array['dxynoabs'][i] = track.dxy(pv_pos)
            track_level_var_array['dznoabs'][i] = track.dz(pv_pos)

            track_level_var_array['dxy'][i] = abs(track.dxy(pv_pos))
            track_level_var_array['dz'][i] = abs(track.dz(pv_pos))
            track_level_var_array['log10(dxy)'][i] = ROOT.TMath.Log10(abs(track.dxy(pv_pos)))
            track_level_var_array['log10(dz)'][i] = ROOT.TMath.Log10(abs(track.dz(pv_pos)))

            dxyhandmade, dzhandmade = handmadeDxyDz(track, pv_pos)
            track_level_var_array['dxyhandmade'][i] = dxyhandmade
            track_level_var_array['dzhandmade'][i] = dzhandmade
            track_level_var_array['log10(dxyhandmade)'][i] = ROOT.TMath.Log10(dxyhandmade)
            track_level_var_array['log10(dzhandmade)'][i] = ROOT.TMath.Log10(dzhandmade)

            mindxy = 999
            mindz = 999
            for v in primaryvertices:
                thisdxy = abs(track.dxy(v.position()))
                if thisdxy < mindxy:
                    mindxy = thisdxy
                thisdz = abs(track.dz(v.position()))
                if thisdz < mindz:
                    mindz = thisdz
            track_level_var_array['dxyclosestpv'][i] = mindxy
            track_level_var_array['dzclosestpv'][i] = mindz
            track_level_var_array['log10(dxyclosestpv)'][i] = ROOT.TMath.Log10(mindxy)
            track_level_var_array['log10(dzclosestpv)'][i] = ROOT.TMath.Log10(mindz)

            mindxyPU = 999
            mindzPU = 999
            for v in primaryvertices[1:]:
                thisdxyPU = abs(track.dxy(v.position()))
                if thisdxyPU < mindxyPU:
                    mindxyPU = thisdxyPU
                thisdzPU = abs(track.dz(v.position()))
                if thisdzPU < mindzPU:
                    mindzPU = thisdzPU
            track_level_var_array['dxyclosestpvPU'][i] = mindxyPU
            track_level_var_array['dzclosestpvPU'][i] = mindzPU
            track_level_var_array['log10(dxyclosestpvPU)'][i] = ROOT.TMath.Log10(mindxyPU)
            track_level_var_array['log10(dzclosestpvPU)'][i] = ROOT.TMath.Log10(mindzPU)

            track_level_var_array['dxyerror'][i] = abs(track.dxyError())
            track_level_var_array['dzerror'][i] = abs(track.dzError())
            track_level_var_array['log10(dxyerror)'][i] = ROOT.TMath.Log10(abs(track.dxyError()))
            track_level_var_array['log10(dzerror)'][i] = ROOT.TMath.Log10(abs(track.dzError()))

            dontSubtractTrackPt = False
            if abs(track.dz(pv_pos)) >= 0.1 or abs(track.dxy(pv_pos)) >= 0.1 or track.pt() <= 1.: dontSubtractTrackPt = True
            track_level_var_array['trackabsisotight'][i], track_level_var_array['trackrelisotight'][i], track_level_var_array['trackdrmintight'][i], track_level_var_array['tracknumneighbourstight'][i] = calcIso_pf_or_track_new(track, tracksforisotight, dontSubtractObject=dontSubtractTrackPt)

            dontSubtractTrackPt = False
            if abs(track.dz(pv_pos)) >= 0.1: dontSubtractTrackPt = True
            track_level_var_array['trackabsiso'][i], track_level_var_array['trackreliso'][i], track_level_var_array['trackdrmin'][i], track_level_var_array['tracknumneighbours'][i] = calcIso_pf_or_track_new(track, tracksforiso, dontSubtractObject=dontSubtractTrackPt)

            track_level_var_array['pfabsiso'][i], track_level_var_array['pfreliso'][i], track_level_var_array['pfdrmin'][i], track_level_var_array['pfnumneighbours'][i] = calcIso_pf_or_track_new(track, pfcandsforiso)
            track_level_var_array['chpfabsiso'][i], track_level_var_array['chpfreliso'][i], track_level_var_array['chpfdrmin'][i], track_level_var_array['chpfnumneighbours'][i] = calcIso_pf_or_track_new(track, chpfcandsforiso, dontSubtractObject=dontSubtractTrackPt)

            track_level_var_array['jetiso'][i], track_level_var_array['jetisomulti'][i], track_level_var_array['jetdrmin'][i], track_level_var_array['jetisobtag'][i] = calcIso_jet_new(track, jetsforiso, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['jetisoloose'][i], track_level_var_array['jetisomultiloose'][i], track_level_var_array['jetdrminloose'][i], track_level_var_array['jetisobtagloose'][i] = calcIso_jet_new(track, jetsforisoloose, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['jetisomeditight'][i], track_level_var_array['jetisomultimeditight'][i], track_level_var_array['jetdrminmeditight'][i], track_level_var_array['jetisobtagmeditight'][i] = calcIso_jet_new(track, jetsforisomeditight, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['jetisomedium'][i], track_level_var_array['jetisomultimedium'][i], track_level_var_array['jetdrminmedium'][i], track_level_var_array['jetisobtagmedium'][i] = calcIso_jet_new(track, jetsforisomedium, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['jetisotightNoLepton'][i], track_level_var_array['jetisomultitightNoLepton'][i], track_level_var_array['jetdrmintightNoLepton'][i], track_level_var_array['jetisobtagtightNoLepton'][i] = calcIso_jet_new(track, jetsforisotightNoLepton, isTrack=True, btagvalues=btagvalues)
            track_level_var_array['jetisotight'][i] = jetisotight
            track_level_var_array['jetisomultitight'][i] = jetisomultitight
            track_level_var_array['jetdrmintight'][i] = jetdrmintight
            track_level_var_array['jetisobtagtight'][i] = jetisobtagtight


            drminneutralhadron = 10
            closestneutralhadronTlv = ROOT.TLorentzVector()
            trackTlv = ROOT.TLorentzVector()
            trackTlv.SetPtEtaPhiM(track.pt(), track.eta(), track.phi(), 0.1396)
            for p in neutralhadrons:
                dr = deltaR(track.eta(), p.eta(), track.phi(), p.phi())
                if dr < drminneutralhadron:
                    drminneutralhadron = dr
                    closestneutralhadronTlv.SetPtEtaPhiM(p.pt(), p.eta(), p.phi(), 0.1350)
            track_level_var_array['drminneutralhadron'][i] = drminneutralhadron
            if drminneutralhadron < 10:
                track_level_var_array['invmclosestneutralhadrontrack'][i] = (closestneutralhadronTlv + trackTlv).M()
            else:
                track_level_var_array['invmclosestneutralhadrontrack'][i] = 0

            drminphoton = 10
            for p in photons:
                dr = deltaR(track.eta(), p.eta(), track.phi(), p.phi())
                if dr < drminphoton:
                    drminphoton = dr
            track_level_var_array['drminphoton'][i] = drminphoton

            drminelectron = 10
            for e in electrons:
                dr = deltaR(track.eta(), e.eta(), track.phi(), e.phi())
                if dr < drminelectron:
                    drminelectron = dr
            track_level_var_array['drminelectron'][i] = drminelectron

            drminmuon = 10
            for m in muons:
                dr = deltaR(track.eta(), m.eta(), track.phi(), m.phi())
                if dr < drminmuon:
                    drminmuon = dr
            track_level_var_array['drminmuon'][i] = drminmuon

            istauleadpfchhadcand = 0
            drmintau = 10
            closesttaumvadiscr = -1
            closesttaudecaymode = -1
            taudr3wp = 0
            taudr4wp = 0
            taudr5wp = 0
            for t in tauswithdiscriminators:

                leadpfchhadcand = t[0].leadPFChargedHadrCand().trackRef().get()
                if leadpfchhadcand:
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

            track_level_var_array['istauleadpfchhadcand'][i] = istauleadpfchhadcand
            track_level_var_array['drmintau'][i] = drmintau
            track_level_var_array['closesttaumvadiscr'][i] = closesttaumvadiscr
            track_level_var_array['closesttaudecaymode'][i] = closesttaudecaymode
            track_level_var_array['taudr3wp'][i] = taudr3wp
            track_level_var_array['taudr4wp'][i] = taudr4wp
            track_level_var_array['taudr5wp'][i] = taudr5wp

            track_level_var_array['detahighestptjet'][i] = abs(track.eta() - jets[idxhighestptjet].eta())
            track_level_var_array['dphihighestptjet'][i] = deltaPhi(track.phi(), jets[idxhighestptjet].phi())

            track_level_var_array['dphimet'][i] = deltaPhi(track.phi(), met.phi())
            track_level_var_array['dphimetpca'][i], _, _ = handmadeDphiMetPCA(track, pv_pos, met)

            track_level_var_array['chi2'][i] = track.normalizedChi2()

            quality = 0
            if track.quality(track.qualityByName('loose')): quality = 1
            if track.quality(track.qualityByName('tight')): quality = 2
            if track.quality(track.qualityByName('highPurity')): quality = 3
            track_level_var_array['quality'][i] = quality

            track_level_var_array['nvalidhits'][i] = track.numberOfValidHits()
            track_level_var_array['nlosthits'][i] = track.numberOfLostHits()


            hasGenMatch = -1
            tmingen = -999
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

            dothegenmatch = False

            if 'genmatchtracks' in options.tag:
                if ievent % 10 == 0 and met.pt() > 250: dothegenmatch = True

            if 'genmatchalltracks' in options.tag:
                if met.pt() > metthresholdtrackgenmatch: dothegenmatch = True

            if dothegenmatch:

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

            track_level_var_array['hasGenMatch'][i] = hasGenMatch
            track_level_var_array['genmatchtmin'][i] = tmingen
            track_level_var_array['genmatchpdgid'][i] = genmatchpdgid
            track_level_var_array['genmatchmotherpdgid'][i] = genmatchmotherpdgid
            track_level_var_array['genmatchpt'][i] = genmatchpt
            track_level_var_array['genmatchmotherpt'][i] = genmatchmotherpt
            track_level_var_array['genmatchstatus'][i] = genmatchstatus
            track_level_var_array['genmatchmotherstatus'][i] = genmatchmotherstatus
            track_level_var_array['genmatchishardprocess'][i] = genmatchishardprocess
            track_level_var_array['genmatchmotherishardprocess'][i] = genmatchmotherishardprocess
            track_level_var_array['genmatchisfromhardprocess'][i] = genmatchisfromhardprocess
            track_level_var_array['genmatchisprompt'][i] = genmatchisprompt
            track_level_var_array['genmatchisdirecthadrondecayproduct'][i] = genmatchisdirecthadrondecayproduct
            track_level_var_array['genmatchisdirecttaudecayproduct'][i] = genmatchisdirecttaudecayproduct
            track_level_var_array['genmatchmotheristhetau'][i] = genmatchmotheristhetau
            track_level_var_array['genmatchmothertaudecay'][i] = decayWtau

            issignaltrack = 0
            if itrack == matchedTrackIdxCharginoPion1 or itrack == matchedTrackIdxCharginoPion2: issignaltrack = 1
            track_level_var_array['issignaltrack'][i] = issignaltrack

            issusytrack = 0
            if itrack in susytracks: issusytrack = 1
            track_level_var_array['issusytrack'][i] = issusytrack

            susytrackmother = 0
            susytrackpdgid = 0
            if issusytrack:
                susytrackmother = susytracks[itrack][0]
                susytrackpdgid = susytracks[itrack][1]
            track_level_var_array['susytrackmother'][i] = susytrackmother
            track_level_var_array['susytrackpdgid'][i] = susytrackpdgid

            i += 1

        event_level_var_array['numtrackstotal'][0] = len(tracks)
        event_level_var_array['numtracksbasicpreselection'][0] = numtracksbasicpreselection
        event_level_var_array['numtracksfinalpreselection'][0] = numtracksfinalpreselection

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
