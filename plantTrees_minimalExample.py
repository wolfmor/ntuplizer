#! /usr/bin/env python2


"""

Runs over AOD files and writes file with histos and tree.

----------------------------------------------------------------------
python plantTrees.py inputFiles="file1, file2,..." tag="tag1 tag2 ..."

minimal example: 
python plantTrees.py inputFiles="/pnfs/desy.de/cms/tier2/store/user/sbein/CommonSamples/RadiativeMu_2016Fast/v2/higgsino94x_susyall_mChipm115GeV_dm0p77GeV_pu35_part22of25.root" tag="test, local, era16_07Aug17, Signal, skipSVs"
python plantTrees.py inputFiles="/nfs/dust/cms/user/wolfmor/testsamples/ZJetsToNuNu_Zpt-200toInf/B044CEA0-F8C9-E611-8F67-0CC47AD990C4.root" tag="test, local, era16_07Aug17"
----------------------------------------------------------------------

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
from commons_Alex import findMatch_tracktrack_new, findMinDr_track, calcIso_vtx, matchToMuon, calcIso_track, calcIso_jet, findMinDr_ancestors, findMatch_ancestor_new, findMinDr



'''
###############################################################################################
# get input arguments
###############################################################################################
'''

options = VarParsing('python')
options.parseArguments()

if 'test' in options.tag: isTest = True
else: isTest = False
nEventsTest = 3 # number of events that are analyzed in case of test
printevery = 1


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
    
    for ifile, f in enumerate(options.inputFiles):
        
        vertexfile = localpath+f.split("/")[-2]+"_"+f.split("/")[-1]

        '''
        ###############################################################################################
        # get sv-level info
        ###############################################################################################
        '''

            
        ### get sv files
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
            for nSV, secondary in enumerate(secondaryVertices):

                filesWithSV[ifile][v_ievent][nSV] = secondary.vx()
                #filesWithSV[ifile][v_ievent][nSV] = secondary
                print 'file', ifile, 'event', v_ievent, 'SV loop', nSV, secondary.vx()

                
                
        print "----------Finished loop over SV files----------------"


for ifile , afile in enumerate(filesWithSV):
    #print 'File loop ', ifile, afile
    for ievent, event in enumerate(filesWithSV[ifile]):
        
        if isTest and ievent >= nEventsTest:
            print 'nEventsTest boundary'
            break
        #print 'Event loop', ievent, event
        #continue
        
        for isv, sv, in enumerate(filesWithSV[ifile][ievent]):
            print 'file', ifile, 'event', ievent, 'SV loop', isv, sv
            #print 'file', ifile, 'event', ievent, 'SV loop', isv, sv.vx()
    
sys.exit()

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


    nEventsPerFile = 0

    phifirsttrack = 0
    etafirsttrack = 0
        

    '''
    ###############################################################################################
    # event loop
    ###############################################################################################
    '''

    #print 'File loop ', ifile, filesWithSV[ifile]
    


    for ievent, event in enumerate(events):

    
        if isTest and ievent >= nEventsTest:
            print 'nEventsTest boundary'
            break
        
        #print 'Event loop', ievent, filesWithSV[ifile][ievent]
        
        
        if 'local' not in options.tag:
            if fin.IsZombie() or not fin.IsOpen():
                print 'file not usable'
                sys.exit(1)

 

        if ievent % printevery == 0: print 'analyzing event %d of %d' % (ievent, nevents)

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
        # get sv-level info
        ###############################################################################################
        '''
        


        if not 'skipSVs' in options.tag: 

            
            numSVs = -1
            n_sv_total = -1
            numsvsfinalpreselection = 0  
            
            n_sv_total = len(filesWithSV[ifile][ievent]) 

            for isv, sv, in enumerate(filesWithSV[ifile][ievent]):
                print 'file', ifile, 'event', ievent, 'SV loop', isv, sv.vx()
            
            continue



        tEvent.Fill()

        nEventsPerFile += 1
