#! /usr/bin/env python


"""Methods used by other scripts.
"""


import sys
from glob import glob
from array import array
from math import ceil, cosh
import numpy as np
import scipy.optimize
from scipy.linalg import block_diag

from ROOT import gROOT, gSystem, FWLiteEnabler, TFile, TH1F, TMath, TLorentzVector, TH2F, TTree

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') # white background

gSystem.Load("libFWCoreFWLite.so")
gSystem.Load("libDataFormatsFWLite.so")
FWLiteEnabler.enable()

from DataFormats.FWLite import Events, Handle

'''Check if electron passes ID (without isolation requirement)
https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/Heppy/python/physicsobjects/Electron.py
https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_2016_data_for
'''
def passesEleID(e, pv_pos, wp='veto'):
	
	if wp == 'veto':
		
		if abs(e.eta()) < 1.4442:  # barrel
			
			if e.full5x5_sigmaIetaIeta() >= 0.0115: return False
			if abs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00749: return False
			if abs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.228: return False
			if e.hadronicOverEm() >= 0.356: return False
			if e.ecalEnergy() > 0:
				if abs(1.0/e.ecalEnergy() - e.eSuperClusterOverP() / e.ecalEnergy()) >= 0.299: return False
			else: 
				return False
			#if e.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS) > 2 : return False
			###if not e.passConversionVeto() >= : return False  # not implemented
			#if abs(e.gsfTrack().dxy(pv_pos)) >= 0.05: return False
			#if abs(e.gsfTrack().dz(pv_pos)) >= 0.1: return False
						
			return True
		
		elif abs(e.eta()) > 1.566:  # endcap
			
			if e.full5x5_sigmaIetaIeta() >= 0.037: return False
			if abs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00895: return False
			if abs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.213: return False
			if e.hadronicOverEm() >= 0.211: return False
			if e.ecalEnergy() > 0:
				if abs(1.0/e.ecalEnergy() - e.eSuperClusterOverP() / e.ecalEnergy()) >= 0.15: return False
			else: 
				return False
			#if e.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS) > 3 : return False
# 			##if not e.passConversionVeto() >= : return False  # not implemented
			#if abs(e.gsfTrack().dxy(pv_pos)) >= 0.1: return False
			#if abs(e.gsfTrack().dz(pv_pos)) >= 0.2: return False
						
			return True
		
		else:
			
			return False
		
	else:
		
		print('working point not implemented')
		sys.exit(0)


'''Check if muon passes ID (without isolation requirement)
https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/Heppy/python/physicsobjects/Muon.py
'''
def passesMuonID(m, wp='loose'):
	
	if wp == 'loose':
		
		if not m.isPFMuon(): return False
		if not (m.isGlobalMuon() or m.isTrackerMuon()): return False
		
		return True
		
	elif wp == 'medium':
		
		print('working point not implemented')  # don't know how to access segment compatibility...
		sys.exit(0)
		
		if not passesMuonID(m, 'loose'): return False
		if m.innerTrack().validFraction() <= 0.8 : return False
		
		goodGlobalMuon = m.isGlobalMuon() \
			and m.globalTrack().normalizedChi2() < 3 \
			and m.combinedQuality().chi2LocalPosition < 12 \
			and m.combinedQuality().trkKink < 20 \
# 			and ROOT.reco.Muon.segmentCompatibility(m) > 0.303
		
# 		if not (goodGlobalMuon or ROOT.reco.Muon.segmentCompatibility(m) > 0.451): return False
		
		return True
		
	else:
		
		print('working point not implemented')
		sys.exit(0)



'''Check if photon passes ID (without isolation requirement)
https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_2016_data_for
'''
def passesPhotID(p, wp='loose'):
	
	if wp == 'loose':
		
		if abs(p.eta()) < 1.4442:  # barrel
			
			if p.hadTowOverEm() >= 0.0597: return False
			if p.sigmaIetaIeta() >= 0.01031: return False
						
			return True
		
		elif abs(p.eta()) > 1.566:  # endcap
			
			if p.hadTowOverEm() >= 0.0481: return False
			if p.sigmaIetaIeta() >= 0.03013: return False
						
			return True
		
		else:
			
			return False
		
	else:
		
		print('working point not implemented')
		sys.exit(0)

'''Calculates delta-beta corrected isolation.
https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
'''
def calcIso_dBeta(pfIsoVars):
	
	neutralIso = pfIsoVars.sumNeutralHadronEt + pfIsoVars.sumPhotonEt
	neutralIsoCorrected = neutralIso - 0.5 * pfIsoVars.sumPUPt
	
	return pfIsoVars.sumChargedHadronPt + max(0., neutralIsoCorrected)



###############################################################################################

'''Converts strings to be used as branch names.
'''
def nice_string(ugly_string):
	return ugly_string.replace('(','_').replace(')','_').replace('/','_D_').replace('*','_T_')


"""Angle between two angles, returns value between -pi and +pi.
"""
def angleBetween(phi1, phi2):

	phi = TMath.ATan2(TMath.Sin(phi1) + TMath.Sin(phi2), TMath.Cos(phi1) + TMath.Cos(phi2))
	
	while phi >= TMath.Pi(): phi -= 2 * TMath.Pi()
	while phi < -TMath.Pi(): phi += 2 * TMath.Pi()
	
	return phi

'''Calculates delta-beta corrected isolation.
https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
'''
def calcIso_dBeta(pfIsoVars):
	
	neutralIso = pfIsoVars.sumNeutralHadronEt + pfIsoVars.sumPhotonEt
	neutralIsoCorrected = neutralIso - 0.5 * pfIsoVars.sumPUPt
	
	return pfIsoVars.sumChargedHadronPt + max(0., neutralIsoCorrected)


"""Adds two angles, returns value between -pi and +pi.
"""
def addPhi(phi1, phi2):
	
	sumphi = phi1 + phi2
	
	while sumphi >= TMath.Pi(): sumphi -= 2 * TMath.Pi()
	while sumphi < -TMath.Pi(): sumphi += 2 * TMath.Pi()
	
	return sumphi


"""Subtracts two angles, returns value between -pi and +pi.
"""
def deltaPhi(phi1, phi2):
	
	dphi = phi1 - phi2
	
	while dphi >= TMath.Pi(): dphi -= 2 * TMath.Pi()
	while dphi < -TMath.Pi(): dphi += 2 * TMath.Pi()
	
	return dphi


"""Returns deltaR.
"""
def deltaR(eta1, eta2, phi1, phi2):

	deta = eta1 - eta2
	dphi = deltaPhi(phi1, phi2)
	
	return TMath.Sqrt(deta * deta + dphi * dphi)
	

"""Calculates isolation-variables for a track inside jets.
"""
def calcIso_jet(one, many, pv_pos, isTight):

	pt_threshold = 15
	conesize = 0.4
	
	if isTight:
		pt_threshold = 30
		conesize = 0.4
	
	ptsum = 0
	num = 0
	
	dRmin = 10
	for m in many:
		
		if not passesPreselection_iso_jet(m, pt_threshold): continue
		
		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


"""Calculates isolation-variables for a PF inside PFs.
"""
def calcIso_pf(one, many, pv_pos, isTight):

	pt_threshold = 0
	conesize = 0.3
	
	if isTight:
		pt_threshold = 0.3
		conesize = 0.1
	
	ptsum = -one.pt()
	num = -1
	
	if not passesPreselection_iso_pfc(one, pv_pos, pt_threshold):
		ptsum = 0
		num = 0
	
	dRmin = 999
	for m in many:
		
		if not passesPreselection_iso_pfc(m, pv_pos, pt_threshold): continue

		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


"""Calculates isolation-variables for a track inside tracks.
"""
def calcIso_track(one, many, pv_pos, isTight):
	
	dz_threshold = 1
	dxy_threshold = 1
	pt_threshold = 0
	conesize = 0.3
	
	if isTight:
		dz_threshold = 10.
		dxy_threshold = 10.
		pt_threshold = 0.
		conesize = 0.3
	
	ptsum = -one.pt()
	num = -1
	
	if not passesPreselection_iso_track(one, pv_pos, dz_threshold, dxy_threshold, pt_threshold):			
		ptsum = 0
		num = 0

	dRmin = 10	
	for m in many:
		
		if not passesPreselection_iso_track(m, pv_pos, dz_threshold, dxy_threshold, pt_threshold): continue
		
		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num
	
"""Calculates isolation-variables for a vtx inside vtxs.
"""
def calcIso_vtx(one, many):

	conesize = 0.3
	
	ptsum = -one.pt()
	num = -1

	dRmin = 10	
	for m in many:

		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


"""Calculates isolation-variables for a gen-particle inside tracks.
"""
def calcIso_gen(one, many, pv_pos, isTight):
	
	dz_threshold = 1
	dxy_threshold = 1
	pt_threshold = 0
	conesize = 0.3
	
	if isTight:
		dz_threshold = 0.15
		dxy_threshold = 0.15
		pt_threshold = 1
		conesize = 0.3
	
	oneTlv = TLorentzVector()
	oneTlv.SetPxPyPzE(one.px(), one.py(), one.pz(), one.energy())

	ptsum = 0
	num = 0
	
	dRmin = 999
	for m in many:
		
		if not passesPreselection_iso_track(m, pv_pos, dz_threshold, dxy_threshold, pt_threshold): continue
		
		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


"""Returns last copy of a gen. particle.
"""
def getLastCopy(gp):
		
	i = 0
	while not gp.status() == 1:
		for idx in range(gp.numberOfDaughters()):
			if gp.pdgId() == gp.daughter(idx).pdgId():
				gp = gp.daughter(idx)
				break
		i += 1
		if i > 100: return None
		
	return gp


"""Returns lepton and neutrino from chi1pm decay to chi10.
"""
def findLeptonNeutrino(gp):
		
	leptonDaughterFound = False
	neutrinoDaughterFound = False
	
	lepton = None
	neutrino = None
	
	for i in range(gp.numberOfDaughters()):
		
		pdgIdDaughter = abs(gp.daughter(i).pdgId())
			
		if pdgIdDaughter == 11 or pdgIdDaughter == 13:
			
			leptonDaughterFound = True
			lepton = gp.daughter(i)
			
		if pdgIdDaughter == 12 or pdgIdDaughter == 14:
			
			neutrinoDaughterFound = True
			neutrino = gp.daughter(i)
			
		if leptonDaughterFound and neutrinoDaughterFound: break
		
	if not (leptonDaughterFound and neutrinoDaughterFound): return None, None
	if not abs(lepton.pdgId()) + 1 == abs(neutrino.pdgId()): return None, None
	if not lepton.charge() * gp.charge() > 0: return None, None
	
	lepton = getLastCopy(lepton)
	neutrino = getLastCopy(neutrino)
	
	return lepton, neutrino
	
	
"""Returns pion from chi1pm decay to chi10.
"""
def findPion(gp):
	
	pionDaughterFound = False
	
	pion = None
	
	for i in range(gp.numberOfDaughters()):
		
		pdgIdDaughter = abs(gp.daughter(i).pdgId())
		
		if pdgIdDaughter == 211:
			
			pionDaughterFound = True
			pion = gp.daughter(i)
			
		if pionDaughterFound: break
		
	if not pionDaughterFound: return None
	if not pion.charge() * gp.charge() > 0: return None
	
	pion = getLastCopy(pion)
				
	return pion	
	

"""Returns lepton pair from chi20 decay to chi10.
"""
def findLeptons(gp):
	
	oneleptonDaughterFound = False
	
	twoleptonDaughtersFound = False
	
	lepton1 = None
	lepton2 = None
	
	for i in range(gp.numberOfDaughters()):
		
		pdgIdDaughter = abs(gp.daughter(i).pdgId())
			
		if pdgIdDaughter == 11 or pdgIdDaughter == 13:
			
			if oneleptonDaughterFound:
				twoleptonDaughtersFound = True
				lepton2 = gp.daughter(i)
			else:
				oneleptonDaughterFound = True
				lepton1 = gp.daughter(i)
				
		if twoleptonDaughtersFound: break

	if not twoleptonDaughtersFound: return None, None
	
	if not abs(lepton1.pdgId()) == abs(lepton2.pdgId()): return None, None
	
	if not lepton1.charge() * lepton2.charge() < 0: return None, None
	
	lepton1 = getLastCopy(lepton1)
	
	lepton2 = getLastCopy(lepton2)
	
	return lepton1, lepton2
	
"""Returns invariant mass from the TLorentzVectors of a pair of electrons
"""
def invariantMass(electron1, electron2):

	electronTlv1 = TLorentzVector()
	
	electronTlv2 = TLorentzVector()
	
	electronTlv1.SetPxPyPzE(electron1.px(), electron1.py(), electron1.pz(), electron1.energy())	
	
	electronTlv2.SetPxPyPzE(electron2.px(), electron2.py(), electron2.pz(), electron2.energy())	

	invMass = (electronTlv1+electronTlv2).M()
	
	return invMass
"""Returns invariant mass from the TLorentzVectors of a pair of tracks
"""
def invariantMassTracks(track1, track2):

	trkTlv1 = TLorentzVector()
	
	trkTlv2 = TLorentzVector()
	
	#trkTlv1.SetPxPyPzE(track1.px(), track1.py(), track1.pz(), track1.pt()*np.cosh(track1.eta()))	
	
	#trkTlv2.SetPxPyPzE(track2.px(), track2.py(), track2.pz(), track2.pt()*np.cosh(track2.eta()))
	trkTlv1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), 0.13957018)
		
	trkTlv2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), 0.13957018)		
	
	invMass = (trkTlv1 +trkTlv2).M()
	
	return invMass
	
"""Returns invariant mass from the TLorentzVectors of a pair of tracks assuming pion mass
"""
def invariantMassTracksAsPions(track1, track2):

	trkTlv1 = TLorentzVector()
	
	trkTlv2 = TLorentzVector()
	
	trkTlv1.SetPtEtaPhiM(track1.pt(), track1.eta(), track1.phi(), 0.13957018)
		
	trkTlv2.SetPtEtaPhiM(track2.pt(), track2.eta(), track2.phi(), 0.13957018)	
	
	invMass = (trkTlv1 +trkTlv2).M()
	
	return invMass

	
"""Defines basic preselection for PFs.
"""
def passesPreselection_basic_pfc(pfc):
	
	if pfc.trackRef().isNull(): return False
	
	if not passesPreselection_basic_track(pfc.trackRef().get()): return False
	
	return True


"""Defines preselection used for calculation of isolation-variables for PFs.
"""
def passesPreselection_iso_pfc(pfc, pv_pos, pt_threshold):
	
	if pfc.trackRef().isNull(): return False
	
	if pfc.trackRef().get().numberOfValidHits() == 0: return False
	
	if pfc.trackRef().get().ndof() == 0: return False
	
	if pfc.trackRef().get().charge() == 0: return False
	
	if abs(pfc.trackRef().get().dz(pv_pos)) > 1: return False
	
	if pfc.pt() <= pt_threshold: return False
	
	return True


"""Defines preselection used for BDT for PFs.
"""
def passesPreselection_final_pfc(pfc, pv_pos):

	if not passesPreselection_basic_pfc(pfc): return False
	
	if pfc.pt() <= 0.7: return False
	
	return True


"""Defines basic preselection for tracks.
"""
def passesPreselection_basic_track(track):
	
	if track.numberOfValidHits() == 0: return False
	
	if track.ndof() == 0: return False
	
	if track.charge() == 0: return False
	
	if track.pt() > 5: return False
	
	return True


"""Defines preselection used for calculation of isolation-variables for tracks.
"""
def passesPreselection_iso_track(track, pv_pos, dz_threshold, dxy_threshold, pt_threshold):
		
	if track.numberOfValidHits() == 0: return False
	
	if track.ndof() == 0: return False
	
	if track.charge() == 0: return False
	
	if abs(track.dz(pv_pos)) >= dz_threshold: return False
	
	if abs(track.dxy(pv_pos)) >= dxy_threshold: return False
	
	if track.pt() <= pt_threshold: return False
			
	return True


"""Defines preselection used for BDT for tracks.
"""
def passesPreselection_final_track(track, pv_pos):
	
	if not passesPreselection_basic_track(track): return False
	
	return True


"""Defines preselection used for calculation of isolation-variables for jets.
"""
def passesPreselection_iso_jet(jet, pt_threshold):
	
	if jet.pt() <= pt_threshold: return False

	return True


"""Finds matching track for lepton via dxyz.
"""
def findMatch_track_new(lepton, tracks):
	
	dxyzmin = 10
	tmin = 10
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
		
	
		for itrack, track in enumerate(tracks):
			
			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * lepton.charge() > 0: continue
			
			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
					
			if not abs(track.eta() - lepton.eta()) < 0.1: continue
			
			if not abs(deltaPhi(track.phi(), lepton.phi())) < 1.57: continue
			
			res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(track, lepton.vertex()))  # other minimization method?
			
			dxyz = distance(res.x, track, lepton.vertex())
			
			if dxyz < dxyzmin:
					
				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(lepton.eta(), track.eta(), lepton.phi(), addPhi(track.phi(), res.x[0]))
				idx = itrack
				
	return idx, dxyzmin, tmin, drmin

	
"""Finds maching gp for tracks in SV; loops over all gen particles per SV daughter, uses classic delta R
"""
def findMinDr_gp_to_SVDaughter(svdaughter, gps):

	pdgId = -1
	pdgIdMother = -1
	drmin = 10
	idx = -1
	
	tlvMother = TLorentzVector()	
	
	if svdaughter.pt() > 0:
	
		for igp, gp in enumerate(gps):
			
			#if not gp.isLastCopy(): continue

			if not gp.charge()== svdaughter.charge(): continue
							
			if not abs(gp.eta() - svdaughter.eta()) < 0.1: continue
			
			if not abs(deltaPhi(gp.phi(), svdaughter.phi())) < 1.57: continue
			
			dr = deltaR(svdaughter.eta(), gp.eta(), gp.phi(), svdaughter.phi())
			
				
			if dr < drmin:
					
				drmin = dr
				idx = igp
				pdgId = gp.pdgId()
				pdgIdMother = gp.mother().pdgId()

				tlvMother.SetPxPyPzE(gp.mother().px(),gp.mother().py(),gp.mother().pz(),gp.mother().energy())
				
	return idx,  pdgId, pdgIdMother, drmin, tlvMother
	
"""Foo
"""
def findMinDr_ancestors(svdaughter, gps, ignoreIdx=[None]):

	pdgId = -1
	pdgIdMother = -1
	drmin = 999
	idx = -1
	
	tlvMother = TLorentzVector()	
	genmatchmother = None
	
	relatives = [None, None]
	charge = 0
	
	numDaughtersOfMother = 0
	
	if svdaughter.pt() > 0:
	
		for igp, gp in enumerate(gps):
				
			#if igp == 797: print " cant use ", ignoreIdx
			#if igp == 797: print " checked gp ", igp," shall ignore this gp", (igp in ignoreIdx), 
			
			if igp in ignoreIdx: continue
			
			#if igp == 797: print " but I dont"
			
			if not gp.isLastCopy(): continue

			if not gp.charge() == svdaughter.charge(): continue
							
			if not abs(gp.eta() - svdaughter.eta()) < 0.1: continue
			
			if not abs(deltaPhi(gp.phi(), svdaughter.phi())) < 1.57: continue
			
			dr = deltaR(svdaughter.eta(), gp.eta(), gp.phi(), svdaughter.phi())
			
				
			if dr < drmin:
					
				drmin = dr
				idx = igp
				pdgId = gp.pdgId()
				pdgIdMother = gp.mother().pdgId()
				genmatchmother = gp.mother()
				tlvMother.SetPxPyPzE(gp.mother().px(),gp.mother().py(),gp.mother().pz(),gp.mother().energy())
				charge = gp.charge()
				
				#relatives[0] = gp
				#relatives[1] = gp.mother()
				numDaughtersOfMother = gp.mother().numberOfDaughters()
	
				
	hasEWancestors = 0
	
	#numDaughtersOfMother = 0
	#if dxyzmin < 0.03: print "match has dxyzmin" , dxyzmin , " charge SV" , svdaughter.charge(), " gp ", charge, " same? ", charge == svdaughter.charge()
	if not (genmatchmother == None):
		iteration = 0
		
		while genmatchmother.numberOfMothers() > 0:			
			# if iteration == 0:
				# d = 0
				# while d < genmatchmother.numberOfDaughters():
					# if genmatchmother.daughter(d).isLastCopy():
						# numDaughtersOfMother +=1
					# d += 1	
					
			if (genmatchmother.pdgId())>1000000:
				hasEWancestors = 1
				#print " I think this is EW", genmatchmother.pdgId()
				break
			
			genmatchmother = genmatchmother.mother(0)
			#relatives.append(genmatchmother)
			iteration += 1
			
			
	#return idx,  pdgId, pdgIdMother, drmin, tlvMother, hasEWancestors, relatives
	return idx,  pdgId, pdgIdMother, drmin, tlvMother, hasEWancestors, numDaughtersOfMother
	
"""Finds matching gp for tracks in SV; loops over all gen particles per SV daughter
"""
def findMatch_gp_new(svdaughter, gps):
	
	dxyzmin = 10
	pdgId = -1
	pdgIdMother = -1
	drmin = 10
	idx = -1
	
	tlvMother = TLorentzVector()	
	if svdaughter.pt() > 0:
	
		for igp, gp in enumerate(gps):
			
			#if not gp.isLastCopy(): continue

			if not gp.charge()== svdaughter.charge(): continue
							
			if not abs(gp.eta() - svdaughter.eta()) < 0.1: continue
			
			if not abs(deltaPhi(gp.phi(), svdaughter.phi())) < 1.57: continue
			
			res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(svdaughter, gp.vertex()))  # other minimization method?
			
			dxyz = distance(res.x, svdaughter, gp.vertex())
			
			if dxyz < dxyzmin:
					
				dxyzmin = dxyz
				drmin = deltaR(svdaughter.eta(), gp.eta(), gp.phi(), addPhi(svdaughter.phi(), res.x[0]))
				idx = igp
				pdgId = gp.pdgId()
				pdgIdMother = gp.mother().pdgId()
				tlvMother.SetPxPyPzE(gp.mother().px(),gp.mother().py(),gp.mother().pz(),gp.mother().energy())
				
	#return idx, dxyzmin, pdgId, pdgIdMother, drmin
	return idx, dxyzmin, pdgId, pdgIdMother, drmin, tlvMother
	
"""Fooo
"""
def findMatch_ancestor_new(svdaughter, gps, ignoreIdx=[None]):
	
	dxyzmin = 10
	pdgId = -1
	pdgIdMother = -1
	drmin = 10
	idx = -1
	genmatchmother = None
	numDaughtersOfMother = 0
	relatives = [None, None]
	charge = 0
	tlvMother = TLorentzVector()	
	if svdaughter.pt() > 0:
	
		for igp, gp in enumerate(gps):
			
			if igp in ignoreIdx: continue
			
			if not gp.isLastCopy(): continue

			if not gp.charge()== svdaughter.charge(): continue
										
			if not abs(gp.eta() - svdaughter.eta()) < 0.1: continue
			
			if not abs(deltaPhi(gp.phi(), svdaughter.phi())) < 1.57: continue
			
			res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(svdaughter, gp.vertex()))  # other minimization method?
			
			dxyz = distance(res.x, svdaughter, gp.vertex())
			
			if dxyz < dxyzmin:
					
				dxyzmin = dxyz
				drmin = deltaR(svdaughter.eta(), gp.eta(), gp.phi(), addPhi(svdaughter.phi(), res.x[0]))
				idx = igp
				pdgId = gp.pdgId()
				pdgIdMother = gp.mother().pdgId()
				tlvMother.SetPxPyPzE(gp.mother().px(),gp.mother().py(),gp.mother().pz(),gp.mother().energy())
				genmatchmother = gp.mother()
				charge = gp.charge()
				
				#relatives[0] = gp				
				#relatives[1] = gp.mother()
				numDaughtersOfMother  = gp.mother().numberOfDaughters()
	
				
	hasEWancestors = 0
	
	#numDaughtersOfMother = 0
	#if dxyzmin < 0.03: print "match has dxyzmin" , dxyzmin , " charge SV" , svdaughter.charge(), " gp ", charge, " same? ", charge == svdaughter.charge()
	if not (genmatchmother == None):
		iteration = 0
		
		while genmatchmother.numberOfMothers() > 0:			
			# if iteration == 0:
				# d = 0
				# while d < genmatchmother.numberOfDaughters():
					# if genmatchmother.daughter(d).isLastCopy():
						# numDaughtersOfMother +=1
					# d += 1	
					
			if (genmatchmother.pdgId())>1000000:
				hasEWancestors = 1
				#print " I think this is EW", genmatchmother.pdgId()
				break
			
			genmatchmother = genmatchmother.mother(0)
			#relatives.append(genmatchmother)
			iteration += 1
		
	#print "idx ",idx,   "dxyzmin",dxyzmin, "pdgId ", pdgId, "pdgIdMother ",pdgIdMother, "drmin ",drmin, "hasEWancestors ", hasEWancestors
	#return idx,  dxyzmin,pdgId, pdgIdMother, drmin, tlvMother, hasEWancestors, relatives
	return idx,  dxyzmin,pdgId, pdgIdMother, drmin, tlvMother, hasEWancestors, numDaughtersOfMother


"""Finds matching tracks in SV for GP; loops over all SVs and all tracks per Sv for each gen  particle (opposite of fn=infMatch_GP_new)
"""
def findMatch_sv_new(gp, secondaries):
	
	dxyzmin = [10, 10]
	dxyz = [999, 999]
	drmin = [10, 10]
	drminNew = [10, 10]
	idx = [-1,-1]
	
	if gp.pt() > 0:
			
		gpDaughter1 = gp.daughter(0)
		gpDaughter2 = gp.daughter(1)
		
		for isv, sv in enumerate(secondaries):
			
			for isvDaughter in range(sv.numberOfDaughters()):
				
				svDaughter = sv.daughter(isvDaughter)
											
				if not (gpDaughter1.pt()- svDaughter.pt())/(gpDaughter1.pt()) < 0.2: continue
				
				#print  ' Sv no.', isv, 'daughter no.', isvDaughter, 'stage 1'
							
				#if not abs(gpDaughter1.eta() - svDaughter.eta()) < 0.1: continue
			
				#print  ' Sv no.', isv, 'daughter no.', isvDaughter, 'stage 2'
			
				if not abs(deltaPhi(gpDaughter1.phi(), svDaughter.phi())) < 1.57: continue
								
				#print  ' Sv no.', isv, 'daughter no.', isvDaughter, 'stage 3'
				
				#print ' Sv no.', isv, 'daughter no.', isvDaughter, 'sv daughter eta', svDaughter.eta(), 'gpDaughter1.phi()', gpDaughter1.phi()
			
				res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(svDaughter, gpDaughter1.vertex()))  # other minimization method?
			
				dxyz[0] = distance(res.x, svDaughter, gpDaughter1.vertex())	
				
				#print  ' Sv no.', isv, 'daughter no.', isvDaughter, 'dxyz[0]', dxyz[0], 'dxyzmin[0]', dxyzmin[0]
			
				if (dxyz[0] < dxyzmin[0]):
					
					drmin[0] = deltaR(svDaughter.eta(), gpDaughter1.eta(), gpDaughter1.phi(), svDaughter.phi())
					drminNew[0] = deltaR(svDaughter.eta(), gpDaughter1.eta(), gpDaughter1.phi(), addPhi(svDaughter.phi(), res.x[0]))

					idx[0] = isv
					
					dxyzmin[0] = dxyz[0]
					
					if isvDaughter == 1: other = 0
					elif isvDaughter == 0: other = 1
											
					if not (gpDaughter2.pt()- sv.daughter(other).pt())/(gpDaughter2.pt()) < 0.2: continue
							
					#if not abs(gpDaughter2.eta() - sv.daughter(other).eta()) < 0.1: continue
			
					if not abs(deltaPhi(gpDaughter2.phi(), sv.daughter(other).phi())) < 1.57: continue
			
					res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(sv.daughter(other), gpDaughter2.vertex()))  # other minimization method?
			
					dxyzmin[1] = distance(res.x, sv.daughter(other), gpDaughter2.vertex())				
					
					drmin[1] = deltaR(sv.daughter(other).eta(), gpDaughter2.eta(), gpDaughter2.phi(), sv.daughter(other).phi())
					drminNew[1] = deltaR(sv.daughter(other).eta(), gpDaughter2.eta(), gpDaughter2.phi(), addPhi(sv.daughter(other).phi(), res.x[0]))

					idx[1] = isv					
				
	return idx, dxyzmin, drmin, drminNew
	
"""Finds matching gp for SV to gp.vtx
"""
def findMatch_gpToSV_new(secondary, gps):
	
	dxyzmin = 10
	pdgId = -1
	drmin = 10
	idx = -1
	
	
	for igp, gp in enumerate(gps):
			
		#if not gp.isLastCopy(): continue

		if not gp.charge()== secondary.charge(): continue
		
		if gp.numberOfDaughters() < 1 : continue
		
		genVtx = gp.daughter(0).vertex() 
		
		sV = secondary.vertex()							
			
		dxyz = TMath.Sqrt(pow(genVtx.x()-sV.x(), 2) + pow(genVtx.y()-sV.y(), 2) + pow(genVtx.z()-sV.z(), 2))
		
		#print "found neutral gp with daughter", genVtx.x(), sV.x(), dxyz
		
		dr = deltaR(secondary.eta(), gp.eta(), secondary.phi(), gp.phi())
		
		if dr < drmin:
			drmin = dr
			idx = igp
			pdgId = gp.pdgId()
				
		if dxyz < dxyzmin:
					
			dxyzmin = dxyz
				
	return idx, dxyzmin, pdgId, drmin

"""Finds matching Sv for gen via delta R.
"""
def findMinDr_GP(lepton, tracks, threshold):
	
	drmin = 10
	idx = -1
	matchingTrack = None
	match = False
	
	if lepton.pt() > 0:
	
		for itrack, track in enumerate(tracks):
			
			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * lepton.charge() > 0: continue
			
			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
					
			if not abs(track.eta() - lepton.eta()) < 0.1: continue
			
			if not abs(deltaPhi(track.phi(), lepton.phi())) < 1.57: continue
			
			leptonTlv = TLorentzVector()
			leptonTlv.SetPxPyPzE(lepton.px(),lepton.py(),lepton.pz(),lepton.energy())	
			trkTlv = TLorentzVector()
			trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*np.cosh(track.eta()))
			
			dr = trkTlv.DeltaR(leptonTlv)
			
			if dr < drmin:
					
				drmin = dr
				idx = itrack
				matchingTrack = track
				

	if drmin < threshold: match = True
	return match, idx, drmin, matchingTrack


"""Finds matching track for track via dxyz.
"""
def findMatch_tracktrack_new(aTrack, tracks):
	
	dxyzmin = 999
	tmin = 10
	drmin = 999
	idx = -1
	
	if aTrack.pt() > 0:
	
		for itrack, track in enumerate(tracks):
			
			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * aTrack.charge() > 0: continue
			##
			if not abs(track.pt() - aTrack.pt()) / aTrack.pt() < 0.2: continue
					
			if not abs(track.eta() - aTrack.eta()) < 0.1: continue
			
			if not abs(deltaPhi(track.phi(), aTrack.phi())) < 1.57: continue
			
			res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(track, aTrack.vertex()))  # other minimization method?
			
			dxyz = distance(res.x, track, aTrack.vertex())
			
			if dxyz < dxyzmin:
					
				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(aTrack.eta(), track.eta(), aTrack.phi(), addPhi(track.phi(), res.x[0]))
				idx = itrack
				
	return idx, dxyzmin, tmin, drmin
	
"""Finds matching track for lepton via delta R.
"""
def findMinDr(lepton, tracks, threshold):
	
	drmin = 10
	idx = -1
	matchingTrack = None
	match = False
	
	if lepton.pt() > 0:
	
		for itrack, track in enumerate(tracks):
			
			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * lepton.charge() > 0: continue
			
			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
					
			if not abs(track.eta() - lepton.eta()) < 0.1: continue
			
			if not abs(deltaPhi(track.phi(), lepton.phi())) < 1.57: continue
			
			leptonTlv = TLorentzVector()
			leptonTlv.SetPxPyPzE(lepton.px(),lepton.py(),lepton.pz(),lepton.energy())	
			trkTlv = TLorentzVector()
			trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*np.cosh(track.eta()))
			
			dr = trkTlv.DeltaR(leptonTlv)
			
			if dr < drmin:
					
				drmin = dr
				idx = itrack
				matchingTrack = track
				

	if drmin < threshold: match = True
	return match, idx, drmin, matchingTrack
	
"""Finds matching track for lepton via delta R.
"""
def findMinDr_track(aTrack, tracks, threshold):
	
	drmin = 10
	idx = -1
	matchingTrack = None
	match = False
	
	if aTrack.pt() > 0:
	
		for itrack, track in enumerate(tracks):
			
			if not passesPreselection_basic_track(track): continue
			#print "1"
			
			if not track.charge() * aTrack.charge() > 0: continue
			#print "2"
			
			#if not abs(track.pt() - aTrack.pt()) / aTrack.pt() < 0.2: continue
			#print "3"
					
			#if not abs(track.eta() - aTrack.eta()) < 0.1: continue
			#print "4"
			
			#if not abs(deltaPhi(track.phi(), aTrack.phi())) < 1.57: continue
			#print "5fin"
			
			aTrackTlv = TLorentzVector()
			aTrackTlv.SetPxPyPzE(aTrack.px(),aTrack.py(),aTrack.pz(),aTrack.pt()*np.cosh(aTrack.eta()))	
			trkTlv = TLorentzVector()
			trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*np.cosh(track.eta()))
			
			dr = trkTlv.DeltaR(aTrackTlv)
			
			if dr < drmin:
					
				drmin = dr
				idx = itrack
				matchingTrack = track
				

	if drmin < threshold: match = True
	return match, idx, drmin, matchingTrack


"""Finds matching track for lepton via dxyz, random with opposite charge.
"""
def findMatch_track_new_random(lepton, tracks):
	
	dxyzmin = 10
	tmin = 10
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
	
		for itrack, track in enumerate(tracks):
			
			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * lepton.charge() < 0: continue
			
			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
			
			if not abs(track.eta() - lepton.eta()) < 0.1: continue
			
			if not abs(deltaPhi(track.phi(), lepton.phi())) < 1.57: continue
			
			res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(track, lepton.vertex()))  # other minimization method?
			
			dxyz = distance(res.x, track, lepton.vertex())
			
			if dxyz < dxyzmin:
				
				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(lepton.eta(), track.eta(), lepton.phi(), addPhi(track.phi(), res.x[0]))
				idx = itrack
			
	return idx, dxyzmin, tmin, drmin


"""Finds matching pfc for lepton via deltaEta and taking distance in 3D into account.
"""
def findMatch_pfc_new(lepton, pfcands):
	
	dxyzmin = 10
	tmin = 10
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
		
		for ipfc, pfc in enumerate(pfcands):
			
			if not passesPreselection_basic_pfc(pfc): continue
			
			if not pfc.charge() * lepton.charge() > 0: continue
			
			if not abs(pfc.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
			
			if not abs(pfc.eta() - lepton.eta()) < 0.1: continue
			
			if not abs(deltaPhi(pfc.phi(), lepton.phi())) < 1.57: continue
			
			res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(pfc.trackRef().get(), lepton.vertex()))  # other minimization method?
			
			dxyz = distance(res.x, pfc.trackRef().get(), lepton.vertex())
			
			if dxyz < dxyzmin:
				
				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(lepton.eta(), pfc.eta(), lepton.phi(), addPhi(pfc.phi(), res.x[0]))
				idx = ipfc
				
	return idx, dxyzmin, tmin, drmin

"""Calculates the pointing angle betweenn the line from PV to SV and the vectoriel combined momentum of a track pair coming from the SV.
"""
def pointingAngle(track1, track2, pv_pos, sv_pos):
	
	trk1Tlv = TLorentzVector()
	trk1Tlv.SetPxPyPzE(track1.px(),track1.py(),track1.pz(),track1.pt()*np.cosh(track1.eta()))
	trk2Tlv = TLorentzVector()
	trk2Tlv.SetPxPyPzE(track2.px(),track2.py(),track2.pz(),track2.pt()*np.cosh(track2.eta()))	
	combinedMomentum =trk2Tlv-trk1Tlv
	
	PVtoSV = [sv_pos[0]-pv_pos.x(),sv_pos[1]-pv_pos.y(),sv_pos[2]-pv_pos.z()]
	
	pointingAngle=np.arctan2(PVtoSV[0] -combinedMomentum.X(),  PVtoSV[1]-combinedMomentum.Y())
	
	return pointingAngle
	

"""Finds a common vertex for two tracks if the minimin distance is smaller than 0.1 cm and returns the vertex' distance to the PV.
"""
def vertexFinder(track1, track2, pv_pos):
	
	mindist, tmin1, tmin2, p1, p2 = minDistanceTrackTrack(track1, track2)
	
	vertex = np.array([0, 0, 0])
	angletoorigin = -10
	angletopv = -10
	disttopv3D = -1
	disttopvXY = -1
	disttopvZ = -1
	if mindist < 0.2:
		
		vertex = (p1 + p2) / 2
		
		disttopv3D = np.linalg.norm(vertex - np.array([pv_pos.x(), pv_pos.y(), pv_pos.z()]))
		disttopvXY = np.linalg.norm(vertex[:2] - np.array([pv_pos.x(), pv_pos.y()]))
		disttopvZ = np.linalg.norm(vertex[2:] - np.array([pv_pos.z()]))
		
		angletoorigin = np.arctan2(vertex[1], vertex[0])
		angletopv = np.arctan2(vertex[1] - pv_pos.y(), vertex[0] - pv_pos.x())
	
	return mindist, tmin1, tmin2, vertex[0], vertex[1], vertex[2], angletoorigin, angletopv, disttopv3D, disttopvXY, disttopvZ
	
"""Finds the minimin distance in R phi of two tracks
"""
def clostestApproachInRPhi(track1, track2, pv_pos):
	
	mindist, tmin1, tmin2, p1, p2 = minDistanceTrackTrack(track1, track2)
	
	vertex = np.array([0, 0, 0])
	angletoorigin = -10
	angletopv = -10
	disttopv3D = -1
	disttopvXY = -1
	disttopvZ = -1
	if mindist < 0.2:
		
		vertex = (p1 + p2) / 2
		
		disttopv3D = np.linalg.norm(vertex - np.array([pv_pos.x(), pv_pos.y(), pv_pos.z()]))
		disttopvXY = np.linalg.norm(vertex[:2] - np.array([pv_pos.x(), pv_pos.y()]))
		disttopvZ = np.linalg.norm(vertex[2:] - np.array([pv_pos.z()]))
		
		angletoorigin = np.arctan2(vertex[1], vertex[0])
		angletopv = np.arctan2(vertex[1] - pv_pos.y(), vertex[0] - pv_pos.x())
	
	return mindist, tmin1, tmin2, vertex[0], vertex[1], vertex[2], angletoorigin, angletopv, disttopv3D, disttopvXY, disttopvZ
	
"""Calculates the new track phi at the PCA
"""
def addPhi(phiALT, tmin):
        
        phiNEU = phiALT + tmin
        
        while phiNEU >= TMath.Pi(): phiNEU -= 2 * TMath.Pi()
        while phiNEU < -TMath.Pi(): phiNEU += 2 * TMath.Pi()
        
        return phiNEU
        
'''Adds two angles, returns value between -pi and +pi.
'''
def addPhiNew(phi1, phi2):
	
	sumphi = phi1 + phi2
	
	while sumphi >= TMath.Pi(): sumphi -= 2 * TMath.Pi()
	while sumphi < -TMath.Pi(): sumphi += 2 * TMath.Pi()
	
	return sumphi
	
""" Corrects the track px, py to be evaluated at the PCA (with new phi)
"""
def trackPxPyAtPCA(track, tmin):
	
	phiNeu = addPhi(track.phi(), tmin)
	pX = track.pt() * np.cos(phiNeu)
	pY = track.pt() * np.sin(phiNeu)
	return pX, pY


"""Gives minimum distance between two tracks and the corresponding points of closest approach on each track.
"""
def minDistanceTrackTrack(track1, track2):
	
	res = scipy.optimize.minimize(distanceTrackTrack, x0=[0.0, 0.0], bounds=((-1.57, 1.57), (-1.57, 1.57)), args=(track1, track2))
	
	tmin1 = res.x[0]
	tmin2 = res.x[1]
	
	d, p1, p2 = getDistanceAndPoints(res.x, track1, track2)
	
	return d, tmin1, tmin2, p1, p2


"""Distance between a specific point along a track (given by t1) and a specific point along another track (given by t2). 
Additionally returns the points.
"""
def getDistanceAndPoints(tarray, track1, track2):
	
	t1 = np.array([tarray[0]])
	t2 = np.array([tarray[1]])
	
	p1 = helix(t1, track1.phi(), track1.eta(), track1.charge(), track1.pt(), track1.vx(), track1.vy(), track1.vz())
	p2 = helix(t2, track2.phi(), track2.eta(), track2.charge(), track2.pt(), track2.vx(), track2.vy(), track2.vz())
	
	d = np.linalg.norm(p1 - p2)
	
	return d, p1, p2


"""Distance between a specific point along a track (given by t1) and a specific point along another track (given by t2).
"""
def distanceTrackTrack(tarray, track1, track2):
	
	t1 = np.array([tarray[0]])
	t2 = np.array([tarray[1]])
	
	d = np.linalg.norm(helix(t1, track1.phi(), track1.eta(), track1.charge(), track1.pt(), track1.vx(), track1.vy(), track1.vz())
					 - helix(t2, track2.phi(), track2.eta(), track2.charge(), track2.pt(), track2.vx(), track2.vy(), track2.vz()))
	
	return d


"""Helix parametrization in 3D.
"""
def helix(t, phi, eta, q, pt, vx, vy, vz):
	
	r = 87.78  # radius [cm] for particle with pT=1GeV in B=3.8T
	
	x = vx + r * q * pt * (np.sin(phi) - np.sin(phi + t))
	
	y = vy + r * q * pt * (-np.cos(phi) + np.cos(phi + t))
	
	z = vz - t * r * q * pt / np.tan(2 * np.arctan(np.exp(-eta)))
	
	return np.array([x[0], y[0], z[0]])



def distanceNew(t, track, (vtx_x, vtx_y, vtx_z)):

	return np.linalg.norm(helixNew(t, track) - np.array([vtx_x, vtx_y, vtx_z]))
 

"""Distance between a specific point along a track (given by t) and a vertex in 3D.
"""
def distance(t, track, vertex):
	
	d = np.linalg.norm(helix(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())
                    - np.array([vertex.x(), vertex.y(), vertex.z()]))
	
	return d


"""Distance between a specific point along a track (given by t) and a vertex in XY.
"""
def distanceXY(t, track, vertex):
	
	d = np.linalg.norm(helix(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())[:2]
                    - np.array([vertex.x(), vertex.y()]))
	
	return d


"""Distance between a specific point along a track (given by t) and a vertex in Z.
"""
def distanceZ(t, track, vertex):
	
	d = np.linalg.norm(helix(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())[2:]
                    - np.array([vertex.z()]))
	
	return d


"""Impact parameters in XY and Z at PCA to PV in 3D (calculated with helix extrapolation).
"""
def handmadeDxyDz(track, pv_pos):
	
	res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(track, pv_pos))
	
	dxy = distanceXY(res.x, track, pv_pos)
	dz = distanceZ(res.x, track, pv_pos)
	
	return dxy, dz


"""Impact parameters in XY and Z at PCA to PV in XY (calculated with helix extrapolation).
"""
def handmadeDxyDzTransversePCA(track, pv_pos):
	
	res = scipy.optimize.minimize(distanceXY, x0=0.0, bounds=((-1.57, 1.57),), args=(track, pv_pos))
	
	dxy = distanceXY(res.x, track, pv_pos)
	dz = distanceZ(res.x, track, pv_pos)
	
	return dxy, dz

	
"""dPhi between MET and vector from PV to PCA to PV of a track
"""
def handmadeDphiMetPCA(track, pv_pos, met):
	
	#gives the helix parameter value for which the distance is minimal
	res = scipy.optimize.minimize(distance, x0=0.0, bounds=((-1.57, 1.57),), args=(track, pv_pos))
	
	#insert found helix parameter and track parameters to get PCA
	pca = helix(res.x, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())
	
	#compute phi angle of vector from PV to PCA
	phipca = np.arctan2(pca[1] - pv_pos.y(), pca[0] - pv_pos.x())
	
	#compute dPhi
	dphi = deltaPhi(phipca, met.phi())
	
	return dphi, phipca, res.x[0]
	

"""Finds matching track for lepton via dR.
"""
def findMatch_track_old(lepton, tracks):
	
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
		
		leptonTlv = TLorentzVector()
		leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
	
		for itrack, track in enumerate(tracks):
			
			if track.numberOfValidHits() == 0: continue
			if track.ndof() == 0: continue
			if track.charge() == 0: continue
		
#			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * lepton.charge() > 0: continue
			
#			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
					
			trackTlv = TLorentzVector()
			trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*TMath.CosH(track.eta()))
			
			dr = trackTlv.DeltaR(leptonTlv)
			
			if dr < drmin:
				
				drmin = dr
				idx = itrack
				
	return idx, drmin
	
"""Finds matching muon for track via dR.
"""
def matchToMuon(track, muons):
	
	drmin = 10
	idx = -1
	
	if track.pt() > 0:
		
		trackTlv = TLorentzVector()
		trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*TMath.CosH(track.eta()))
	
		for imuon, muon in enumerate(muons):
			
			#if muon.numberOfValidHits() == 0: continue
			#if muon.ndof() == 0: continue
			if muon.charge() == 0: continue
		
#			if not passesPreselection_basic_muon(muon): continue
			
			if not muon.charge() * track.charge() > 0: continue
			
#			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
					
			muonTlv = TLorentzVector()
			muonTlv.SetPxPyPzE(muon.px(), muon.py(), muon.pz(), muon.pt()*TMath.CosH(muon.eta()))
			
			dr = muonTlv.DeltaR(trackTlv)
			
			if dr < drmin:
				
				drmin = dr
				idx = imuon
				
	return idx, drmin


"""Finds matching track for lepton via dR, random with opposite charge.
"""
def findMatch_track_old_random(lepton, tracks):
	
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
		
		leptonTlv = TLorentzVector()
		leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
	
		for itrack, track in enumerate(tracks):
						
			if track.numberOfValidHits() == 0: continue
			if track.ndof() == 0: continue
			if track.charge() == 0: continue
			
#			if not passesPreselection_basic_track(track): continue
			
			if not track.charge() * lepton.charge() < 0: continue
			
#			if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
					
			trackTlv = TLorentzVector()
			trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*TMath.CosH(track.eta()))
			
			dr = trackTlv.DeltaR(leptonTlv)
			
			if dr < drmin:
				
				drmin = dr
				idx = itrack
				
	return idx, drmin


"""Finds matching pfc for lepton via deltaR.
"""
def findMatch_pfc_old(lepton, pfcands):
	
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
		
		leptonTlv = TLorentzVector()
		leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
		
		for ipfc, pfc in enumerate(pfcands):
			
			if pfc.trackRef().isNull(): continue
			if pfc.trackRef().get().numberOfValidHits() == 0: continue
			if pfc.trackRef().get().ndof() == 0: continue
			if pfc.trackRef().get().charge() == 0: continue
			
#			if not passesPreselection_basic_pfc(pfc): continue
			
			if not pfc.charge() * lepton.charge() > 0: continue
			
#			if not abs(pfc.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
			
			pfcTlv = TLorentzVector()
			pfcTlv.SetPxPyPzE(pfc.px(), pfc.py(), pfc.pz(), pfc.energy())
			
			dr = pfcTlv.DeltaR(leptonTlv)
			
			if dr < drmin:
				
				drmin = dr
				idx = ipfc
				
	return idx, drmin


"""Finds matching jet for lepton via deltaR.
"""
def findMatch_jet_old(lepton, jets):
	
	drmin = 10
	idx = -1
		
	for ijet, jet in enumerate(jets):
		
		dr = deltaR(lepton.eta(), jet.eta(), lepton.phi(), jet.phi())
		
		if dr < drmin:
			
			drmin = dr
			idx = ijet
				
	return idx, drmin


'''Helix parametrization using CMSSW track parameters.
Particles with positive (negative) charge travel along their helix in direction of positive (negative) t parameter.
To get the azimuthal angle phi for a given t parameter, -t has to be added to phi0 (same for positive and negative charges).
Apart from phi, the track parameters don't change.
Mostly taken from Lucas Wiens (Higgs -> tau tau CP analysis: CMS AN-19-192),
but simplified trigonometric phi term and respecting track charge (differnt direction of rotation).
'''
def helixNew(t, track):

		# conversion factor to convert from T to GeV/(e * cm) = c [m/s] * 10^-9 * 10^-2
		# 1 T = 1 V s / m^2 = 1 eV s / (e m^2) = 1 s/m * eV / (e m) = c * eV / (e m)
		B = 3.8 * 0.0029979

		qoverp = track.parameter(0)
		lmbd = track.parameter(1)
		phi = track.parameter(2)

		x = track.referencePoint().x() + np.cos(lmbd) / B / qoverp * (np.sin(phi) - np.sin(phi - t))

		y = track.referencePoint().y() + np.cos(lmbd) / B / qoverp * (-np.cos(phi) + np.cos(phi - t))

		z =  track.referencePoint().z() + t * np.sin(lmbd) / B / qoverp

		return np.array([x[0], y[0], z[0]])


###############################################################################################

##########################################################

class IPcalculator:

	def __init__(self, track, vertex, verbose=False):

		self.verbose = verbose

		self.track = track

		self.vtx_x = vertex.x()
		self.vtx_y = vertex.y()
		self.vtx_z = vertex.z()
		self.vtx_covariance = np.array([[vertex.covariance(i,j) for j in range(3)] for i in range(3)])

		# 0:qoverp, 1:lambda, 2:phi, 3:dxy, 4:dsz
		self.track_parameters = np.array([track.parameter(i) for i in range(5)])
		self.track_covariance = np.array([[track.covariance(i,j) for j in range(5)] for i in range(5)])

		self.IPvector, self.tmin = self.calculateIPvector()
		self.IPcovariance = self.calculateIPcovariance()

		self.PCA = helixNew(self.tmin, self.track)
		self.newPhi = addPhiNew(self.track.phi(), -1.*self.tmin)

		self.dxyz = np.linalg.norm(self.IPvector)
		self.dxy = np.linalg.norm(self.IPvector[:2])
		self.dz = np.linalg.norm(self.IPvector[2:])

		if self.dxyz == 0: self.IPvectorNorm = self.IPvector
		else: self.IPvectorNorm = self.IPvector / self.dxyz

		self.IPerror = np.sqrt(np.dot(self.IPvectorNorm, np.dot(self.IPcovariance, self.IPvectorNorm)))

		if self.verbose:
			print 'vertex', (self.vtx_x, self.vtx_y, self.vtx_z)

			print 'track parameters'
			print self.track_parameters
			print 'track covariance matrix'
			print self.track_covariance

			print 'IPvector', self.IPvector
			print 'IPvectorNorm', self.IPvectorNorm

			print 'IPcovariance', self.IPcovariance
			print 'IPerror', self.IPerror

	def calculateIPvector(self):

		res = scipy.optimize.minimize(distanceNew, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(self.track, (self.vtx_x, self.vtx_y, self.vtx_z)))

		return helixNew(res.x, self.track) - np.array([self.vtx_x, self.vtx_y, self.vtx_z]), res.x

	def calculateIPcovariance(self):

		# conversion factor to convert from T to GeV/(e * cm) = c [m/s] * 10^-9 * 10^-2
		# 1 T = 1 V s / m^2 = 1 eV s / (e m^2) = 1 s/m * eV / (e m) = c * eV / (e m)
		B = 3.8 * 0.0029979

		qoverp = self.track_parameters[0]
		lmbd = self.track_parameters[1]
		phi = self.track_parameters[2]
		dxy = self.track_parameters[3]
		dsz = self.track_parameters[4]

		HelixJacobian = np.zeros((3,5))
		# derivatives of the helix coordinates with respect to
		# qoverp
		HelixJacobian[0][0] = - (np.cos(lmbd) * (np.sin(phi) - np.sin(phi - self.tmin))) / (B * qoverp**2.)
		HelixJacobian[1][0] = - (np.cos(lmbd) * (-np.cos(phi) + np.cos(phi - self.tmin))) / (B * qoverp**2.)
		HelixJacobian[2][0] = - (np.sin(lmbd) * self.tmin) / (B * qoverp**2.)
		# lambda
		HelixJacobian[0][1] = - (np.sin(lmbd) * (np.sin(phi) - np.sin(phi - self.tmin))) / (B * qoverp)
		HelixJacobian[1][1] = - (np.sin(lmbd) * (-np.cos(phi) + np.cos(phi - self.tmin))) / (B * qoverp)
		HelixJacobian[2][1] = (dsz * np.tan(lmbd)) / (np.cos(lmbd)) + (np.cos(lmbd) * self.tmin) / (B * qoverp)
		# phi
		HelixJacobian[0][2] = - dxy * np.cos(phi) + (np.cos(lmbd) * (np.cos(phi) - np.cos(phi - self.tmin))) / (B * qoverp)
		HelixJacobian[1][2] = - dxy * np.sin(phi) + (np.cos(lmbd) * (np.sin(phi) - np.sin(phi - self.tmin))) / (B * qoverp)
		HelixJacobian[2][2] = 0.
		# dxy
		HelixJacobian[0][3] = - np.sin(phi)
		HelixJacobian[1][3] = np.cos(phi)
		HelixJacobian[2][3] = 0.
		# dsz
		HelixJacobian[0][4] = 0.
		HelixJacobian[1][4] = 0.
		HelixJacobian[2][4] = 1. / np.cos(lmbd)

		HelixCovariance = np.dot(HelixJacobian, np.dot(self.track_covariance, HelixJacobian.T))

		HelixAndPVcovariance = block_diag(HelixCovariance, self.vtx_covariance)

		IPjacobian = np.zeros((3,6))
		# derivatives of the IP with respect to
		# helix
		IPjacobian[0][0] = 1.
		IPjacobian[1][1] = 1.
		IPjacobian[2][2] = 1.
		# PV
		IPjacobian[0][3] = -1.
		IPjacobian[1][4] = -1.
		IPjacobian[2][5] = -1.

		if self.verbose:
			print 'HelixJacobian'
			print HelixJacobian
			print 'HelixCovariance'
			print HelixCovariance
			print 'HelixAndPVcovariance'
			print HelixAndPVcovariance
			print 'IPjacobian'
			print IPjacobian

		return np.dot(IPjacobian, np.dot(HelixAndPVcovariance, IPjacobian.T))

	def getPCAonTrack(self):
		return self.PCA[0], self.PCA[1], self.PCA[2]

	def getPhiAtPCA(self):
		return self.newPhi[0]

	def getIPvector(self):
		return self.IPvector

	def getDxy(self):
		return self.dxy

	def getDz(self):
		return self.dz

	def getDxyDz(self):
		return self.dxy, self.dz

	def getIP(self):
		return self.dxyz

	def getIPsignificance(self):
		if self.IPerror == 0: return 0.
		return self.getIP() / self.IPerror
