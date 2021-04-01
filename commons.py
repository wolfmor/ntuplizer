#! /usr/bin/env python


'''Methods used by other scripts.
'''


import sys
from glob import glob
from array import array
from math import ceil
import numpy as np
import scipy.optimize
from scipy.linalg import block_diag

import ROOT

from ROOT import gROOT, gSystem, FWLiteEnabler, TMath, TLorentzVector

gROOT.SetBatch()        # don't pop up canvases
gROOT.SetStyle('Plain') # white background

gSystem.Load('libFWCoreFWLite.so')
gSystem.Load('libDataFormatsFWLite.so')
FWLiteEnabler.enable()

from DataFormats.FWLite import Events, Handle


'''Converts strings to be used as branch names.
'''
def nice_string(ugly_string):
	return ugly_string.replace('(','_').replace(')','_').replace('/','_D_').replace('*','_T_')


'''Implements dummy physics object class from TLorentzVector.
'''
class Dummy:
	
	def __init__(self, P4):
		
		self._pt = P4.Pt()
		self._eta = P4.PseudoRapidity()
		self._phi = P4.Phi()
		self._px = P4.Px()
		self._py = P4.Py()
		self._pz = P4.Pz()
		self._energy = P4.Energy()

	def pt(self):
		return self._pt
	
	def eta(self):
		return self._eta
	
	def phi(self):
		return self._phi
	
	def px(self):
		return self._px
	
	def py(self):
		return self._py
	
	def pz(self):
		return self._pz
	
	def energy(self):
		return self._energy

###############################################################################################

'''Angle between two angles, returns value between -pi and +pi.
'''
def angleBetween(phi1, phi2):

	phi = TMath.ATan2(TMath.Sin(phi1) + TMath.Sin(phi2), TMath.Cos(phi1) + TMath.Cos(phi2))
	
	while phi >= TMath.Pi(): phi -= 2 * TMath.Pi()
	while phi < -TMath.Pi(): phi += 2 * TMath.Pi()
	
	return phi


'''Adds two angles, returns value between -pi and +pi.
'''
def addPhi(phi1, phi2):
	
	sumphi = phi1 + phi2
	
	while sumphi >= TMath.Pi(): sumphi -= 2 * TMath.Pi()
	while sumphi < -TMath.Pi(): sumphi += 2 * TMath.Pi()
	
	return sumphi


'''Subtracts two angles, returns value between -pi and +pi.
'''
def deltaPhi(phi1, phi2):
	
	dphi = phi1 - phi2
	
	while dphi >= TMath.Pi(): dphi -= 2 * TMath.Pi()
	while dphi < -TMath.Pi(): dphi += 2 * TMath.Pi()
	
	return dphi


'''Returns deltaR.
'''
def deltaR(eta1, eta2, phi1, phi2):

	deta = eta1 - eta2
	dphi = deltaPhi(phi1, phi2)
	
	return TMath.Sqrt(deta * deta + dphi * dphi)

###############################################################################################

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
			if e.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS) > 2 : return False
# 			if not e.passConversionVeto() >= : return False  # not implemented
			if abs(e.gsfTrack().dxy(pv_pos)) >= 0.05: return False
			if abs(e.gsfTrack().dz(pv_pos)) >= 0.1: return False
						
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
			if e.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS) > 3 : return False
# 			if not e.passConversionVeto() >= : return False  # not implemented
			if abs(e.gsfTrack().dxy(pv_pos)) >= 0.1: return False
			if abs(e.gsfTrack().dz(pv_pos)) >= 0.2: return False
						
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

###############################################################################################

'''Calculates delta-beta corrected isolation.
https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
'''
def calcIso_dBeta(pfIsoVars):
	
	neutralIso = pfIsoVars.sumNeutralHadronEt + pfIsoVars.sumPhotonEt
	neutralIsoCorrected = neutralIso - 0.5 * pfIsoVars.sumPUPt
	
	return pfIsoVars.sumChargedHadronPt + max(0., neutralIsoCorrected)


'''Calculates isolation-variables for an object inside jets.
'''
def calcIso_jet_new(one, many, isTrack=False):
	
	jetiso = 0
	jetisomulti = 0
	
	dRmin = 999
	closestjet = None
	for m in many:
		
		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin:
			dRmin = dR
			closestjet = m
			
	if dRmin < 0.4:
		
		if isTrack: oneTlv = TLorentzVector(one.px(), one.py(), one.pz(), one.pt()*TMath.CosH(one.eta()))
		else: oneTlv = TLorentzVector(one.px(), one.py(), one.pz(), one.energy())
		jetTlv = TLorentzVector(closestjet.px(), closestjet.py(), closestjet.pz(), closestjet.energy())
		
		jetiso = (jetTlv - oneTlv).Pt()
		jetisomulti = closestjet.numberOfDaughters()
		
	return jetiso, jetisomulti, dRmin


'''Calculates isolation-variables for an object inside PFs or tracks.
'''
def calcIso_pf_or_track_new(one, many, isMini=False, dontSubtractObject=False):

	conesize = 0.3
	
	if isMini:
		if one.pt() <= 50: conesize = 0.2
		elif one.pt() <= 200: conesize = 10.0/one.pt()
		else: conesize = 0.05
	
	ptsum = -one.pt()
	num = -1
	if dontSubtractObject:
		ptsum = 0
		num = 0
	
	dRmin = 999
	for m in many:

		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum, ptsum/one.pt(), dRmin, num


'''Calculates isolation-variables for a track inside jets.
'''
def calcIso_jet(one, many, pv_pos, isTight):

#	pt_threshold = 10
	conesize = 0.5
	
	if isTight:
#		pt_threshold = 30
		conesize = 0.5
	
	ptsum = 0
	num = 0
	
	dRmin = 10
	for m in many:
		
		# expects only jets that pass preselection !!!
# 		if not passesPreselection_iso_jet(m, pt_threshold): continue
		
		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


'''Calculates isolation-variables for a PF inside PFs.
'''
def calcIso_pf(one, many, pv_pos, isTight):

#	pt_threshold = 0
	conesize = 0.3
	
	if isTight:
#		pt_threshold = 1
		conesize = 0.2
	
	ptsum = -one.pt()
	num = -1
	
#	if not passesPreselection_iso_pfc(one, pv_pos, pt_threshold):
#		ptsum = 0
#		num = 0
	
	dRmin = 999
	for m in many:
		
		# no preselection
#		if not passesPreselection_iso_pfc(m, pv_pos, pt_threshold): continue

		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


'''Calculates isolation-variables for a track inside tracks.
'''
def calcIso_track(one, many, pv_pos, isTight):
	
	dz_threshold = 1
	dxy_threshold = 1
	pt_threshold = 0
	conesize = 0.3
	
	if isTight:
		dz_threshold = 0.15
		dxy_threshold = 0.15
		pt_threshold = 1
		conesize = 0.3
	
	ptsum = -one.pt()
	num = -1
	
	if not passesPreselection_iso_track(one, pv_pos, dz_threshold, dxy_threshold, pt_threshold):			
		ptsum = 0
		num = 0

	dRmin = 10	
	for m in many:
		
		# expects only tracks that pass preselection !!!
#		if not passesPreselection_iso_track(m, pv_pos, dz_threshold, dxy_threshold, pt_threshold): continue
		
		dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
		
		if dR < dRmin and dR > 0.001:
			dRmin = dR
			
		if dR < conesize:
			ptsum += m.pt()
			num += 1
			
	return ptsum/one.pt(), dRmin, num


'''Calculates isolation-variables for a gen-particle inside tracks.
'''
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

###############################################################################################

'''Returns last copy of a gen. particle with status 1 = final state.
'''
def getLastCopyStatusOne(gp):
		
	i = 0
	while not gp.status() == 1:
		for idx in range(gp.numberOfDaughters()):
			if gp.pdgId() == gp.daughter(idx).pdgId():
				gp = gp.daughter(idx)
				break
		i += 1
		if i > 100: return None
		
	return gp


'''Returns last copy of a gen. particle.
'''
def getLastCopy(gp):
		
	i = 0
	while not gp.isLastCopy():
		for idx in range(gp.numberOfDaughters()):
			if gp.pdgId() == gp.daughter(idx).pdgId():
				gp = gp.daughter(idx)
				break
		i += 1
		if i > 100: return None
		
	return gp


'''Returns lepton and neutrino from chi1pm decay to chi10.
'''
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
	
	lepton = getLastCopyStatusOne(lepton)
	neutrino = getLastCopyStatusOne(neutrino)
	
	return lepton, neutrino
	
	
'''Returns pion from chi1pm decay to chi10.
'''
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
	
	pion = getLastCopyStatusOne(pion)
				
	return pion	
	

'''Returns lepton pair from chi20 decay to chi10.
'''
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
	
	lepton1 = getLastCopyStatusOne(lepton1)
	lepton2 = getLastCopyStatusOne(lepton2)
	
	return lepton1, lepton2


'''Returns list of status 1 (grand)daughters.
'''
def findDaughters(gp):
	
	daughters = []
	
	for i in range(gp.numberOfDaughters()):
		
		daughter = getLastCopy(gp.daughter(i))
		
		if daughter.status() == 1:
			
			daughters.append(daughter)
					
		else:
			
			for granddaughter in findDaughters(daughter):
				
				daughters.append(granddaughter)
	
	return daughters

###############################################################################################

'''Defines basic preselection for PFs.
'''
def passesPreselection_basic_pfc(pfc):
	
	if pfc.trackRef().isNull(): return False
	
	if not passesPreselection_basic_track(pfc.trackRef().get()): return False
	
	return True


'''Defines preselection used for calculation of isolation-variables for PFs.
'''
def passesPreselection_iso_pfc(pfc, pv_pos, dz_threshold):
	
	if pfc.charge() == 0: return False
	
	if pfc.trackRef().isNull(): return False
	
	if pfc.trackRef().get().numberOfValidHits() == 0: return False
	
	if pfc.trackRef().get().ndof() == 0: return False
	
	if pfc.trackRef().get().charge() == 0: return False
	
	if abs(pfc.trackRef().get().dz(pv_pos)) >= dz_threshold: return False
	
	return True


'''Defines preselection used for BDT for PFs.
'''
def passesPreselection_final_pfc(pfc, pv_pos):

	if not passesPreselection_basic_pfc(pfc): return False
	
	if pfc.pt() <= 0.7: return False
	
	return True


'''Defines basic preselection for tracks.
'''
def passesPreselection_basic_track(track):
	
	if track.numberOfValidHits() == 0: return False
	
	if track.ndof() == 0: return False
	
	if track.charge() == 0: return False
	
	if track.pt() > 5: return False
	
	return True


'''Defines preselection used for calculation of isolation-variables for tracks.
'''
def passesPreselection_iso_track(track, pv_pos, dz_threshold, dxy_threshold, pt_threshold):
		
	if track.numberOfValidHits() == 0: return False
	
	if track.ndof() == 0: return False
	
	if track.charge() == 0: return False
	
	if abs(track.dz(pv_pos)) >= dz_threshold: return False
	if abs(track.dxy(pv_pos)) >= dxy_threshold: return False

	if track.pt() <= pt_threshold: return False
			
	return True


'''Defines preselection used for BDT for tracks.
'''
def passesPreselection_final_track(track, pv_pos):
	
	if not passesPreselection_basic_track(track): return False
	
	return True


'''Defines preselection used for calculation of isolation-variables for jets.
'''
def passesPreselection_iso_jet(jet, pt_threshold):
	
	if jet.pt() <= pt_threshold: return False

	return True

###############################################################################################

'''Finds a common vertex for two tracks if the minimin distance is smaller than 0.1 cm and returns the vertex' distance to the PV.
'''
def vertexFinder(track1, track2, pv_pos):
	
	mindist, tmin1, tmin2, p1, p2 = minDistanceTrackTrack(track1, track2)
	
	vertex = np.array([0, 0, 0])
	angletoorigin = -10
	angletopv = -10
	disttopv3D = -1
	disttopvXY = -1
	disttopvZ = -1
	if mindist < 0.1:
		
		vertex = (p1 + p2) / 2
		
		disttopv3D = np.linalg.norm(vertex - np.array([pv_pos.x(), pv_pos.y(), pv_pos.z()]))
		disttopvXY = np.linalg.norm(vertex[:2] - np.array([pv_pos.x(), pv_pos.y()]))
		disttopvZ = np.linalg.norm(vertex[2:] - np.array([pv_pos.z()]))
		
		angletoorigin = np.arctan2(vertex[1], vertex[0])
		angletopv = np.arctan2(vertex[1] - pv_pos.y(), vertex[0] - pv_pos.x())
	
	return mindist, tmin1, tmin2, vertex[0], vertex[1], vertex[2], angletoorigin, angletopv, disttopv3D, disttopvXY, disttopvZ


'''Gives minimum distance between two tracks and the corresponding points of closest approach on each track.
'''
def minDistanceTrackTrack(track1, track2):
	
	res = scipy.optimize.minimize(distanceTrackTrack, x0=np.array([0.0, 0.0]), bounds=((-1.57, 1.57), (-1.57, 1.57)), args=(track1, track2))
	
	tmin1 = res.x[0]
	tmin2 = res.x[1]
	
	d, p1, p2 = getDistanceAndPoints(res.x, track1, track2)
	
	return d, tmin1, tmin2, p1, p2


'''Distance between a specific point along a track (given by t1) and a specific point along another track (given by t2). 
Additionally returns the points.
'''
def getDistanceAndPoints(tarray, track1, track2):
	
	t1 = np.array([tarray[0]])
	t2 = np.array([tarray[1]])
	
	p1 = helixOld(t1, track1.phi(), track1.eta(), track1.charge(), track1.pt(), track1.vx(), track1.vy(), track1.vz())
	p2 = helixOld(t2, track2.phi(), track2.eta(), track2.charge(), track2.pt(), track2.vx(), track2.vy(), track2.vz())
	
	d = np.linalg.norm(p1 - p2)
	
	return d, p1, p2


'''Distance between a specific point along a track (given by t1) and a specific point along another track (given by t2).
'''
def distanceTrackTrack(tarray, track1, track2):
	
	t1 = np.array([tarray[0]])
	t2 = np.array([tarray[1]])
	
	d = np.linalg.norm(helixOld(t1, track1.phi(), track1.eta(), track1.charge(), track1.pt(), track1.vx(), track1.vy(), track1.vz())
					 - helixOld(t2, track2.phi(), track2.eta(), track2.charge(), track2.pt(), track2.vx(), track2.vy(), track2.vz()))
	
	return d


'''Helix parametrization in 3D. Same result as helix(...) but not using CMSSW track parameters. 
'''
def helixOld(t, phi, eta, q, pt, vx, vy, vz):
	
	r = 87.78  # radius [cm] for particle with pT=1GeV in B=3.8T
	
	x = vx + r * q * pt * (np.sin(phi) - np.sin(phi - t))
	
	y = vy + r * q * pt * (-np.cos(phi) + np.cos(phi - t))
	
	z = vz + t * r * q * pt / np.tan(2 * np.arctan(np.exp(-eta)))
	
	return np.array([x[0], y[0], z[0]])


'''Distance between a specific point along a track (given by t) and a vertex in 3D.
'''
def distanceOld(t, track, vertex):
	
	d = np.linalg.norm(helixOld(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())
					- np.array([vertex.x(), vertex.y(), vertex.z()]))
	
	return d


'''Distance in XY between a specific point along a track (given by t) and a vertex.
'''
def distanceXYOld(t, track, vertex):
	
	d = np.linalg.norm(helixOld(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())[:2]
					- np.array([vertex.x(), vertex.y()]))
	
	return d


'''Distance in Z between a specific point along a track (given by t) and a vertex.
'''
def distanceZOld(t, track, vertex):
	
	d = np.linalg.norm(helixOld(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())[2:]
					- np.array([vertex.z()]))
	
	return d


'''Impact parameters in XY and Z at PCA to PV in 3D (calculated with helix extrapolation).
'''
def handmadeDxyDz(track, pv_pos):
	
	res = scipy.optimize.minimize(distanceOld, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, pv_pos))
	
	dxy = distanceXYOld(res.x, track, pv_pos)
	dz = distanceZOld(res.x, track, pv_pos)
	
	return dxy, dz


'''Impact parameters in XY and Z at PCA to PV in XY (calculated with helix extrapolation).
'''
def handmadeDxyDzTransversePCA(track, pv_pos):
	
	res = scipy.optimize.minimize(distanceXYOld, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, pv_pos))
	
	dxy = distanceXYOld(res.x, track, pv_pos)
	dz = distanceZOld(res.x, track, pv_pos)
	
	return dxy, dz

	
'''dPhi between MET and vector from PV to PCA to PV of a track
'''
def handmadeDphiMetPCA(track, pv_pos, met):
	
	#gives the helix parameter value for which the distance is minimal
	res = scipy.optimize.minimize(distanceOld, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, pv_pos))
	
	#insert found helix parameter and track parameters to get PCA
	pca = helixOld(res.x, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())
	
	#compute phi angle of vector from PV to PCA
	phipca = np.arctan2(pca[1] - pv_pos.y(), pca[0] - pv_pos.x())
	
	#compute dPhi
	dphi = deltaPhi(phipca, met.phi())
	
	return dphi, phipca, res.x[0]


def distance(t, track, (vtx_x, vtx_y, vtx_z)):

	return np.linalg.norm(helix(t, track) - np.array([vtx_x, vtx_y, vtx_z]))


'''Helix parametrization using CMSSW track parameters.

Particles with positive (negative) charge travel along their helix in direction of positive (negative) t parameter.
To get the azimuthal angle phi for a given t parameter, -t has to be added to phi0 (same for positive and negative charges).
Apart from phi, the track parameters don't change.

Mostly taken from Lucas Wiens (Higgs -> tau tau CP analysis: CMS AN-19-192),
but simplified trigonometric phi term and respecting track charge (differnt direction of rotation).
'''
def helix(t, track):

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

'''Finds matching track for lepton/pion/... via dxyz.
'''
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

			res = scipy.optimize.minimize(distance, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z())))

			dxyz = distance(res.x, track, (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z()))

			if dxyz < dxyzmin:

				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(lepton.eta(), track.eta(), lepton.phi(), addPhi(track.phi(), -1. * res.x[0]))
				idx = itrack

	return idx, dxyzmin, tmin, drmin


'''Finds matching track for lepton via dxyz, random with opposite charge.
'''
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

			res = scipy.optimize.minimize(distance, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z())))

			dxyz = distance(res.x, track, (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z()))

			if dxyz < dxyzmin:

				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(lepton.eta(), track.eta(), lepton.phi(), addPhi(track.phi(), -1. * res.x[0]))
				idx = itrack

	return idx, dxyzmin, tmin, drmin


'''Finds matching pfc for lepton via dxyz.
'''
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

			res = scipy.optimize.minimize(distance, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(pfc.trackRef().get(), (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z())))

			dxyz = distance(res.x, pfc.trackRef().get(), (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z()))

			if dxyz < dxyzmin:

				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(lepton.eta(), pfc.eta(), lepton.phi(), addPhi(pfc.phi(), -1. * res.x[0]))
				idx = ipfc

	return idx, dxyzmin, tmin, drmin


'''Finds matching track for lepton via dR.
'''
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


'''Finds matching track for lepton via dR, random with opposite charge.
'''
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


'''Finds matching pfc for lepton via deltaR.
'''
def findMatch_pfc_old(lepton, pfcands):
	
	drmin = 10
	idx = -1
	
	if lepton.pt() > 0:
		
		leptonTlv = TLorentzVector()
		leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
		
		for ipfc, pfc in enumerate(pfcands):
			
# 			if pfc.trackRef().isNull(): continue
# 			if pfc.trackRef().get().numberOfValidHits() == 0: continue
# 			if pfc.trackRef().get().ndof() == 0: continue
# 			if pfc.trackRef().get().charge() == 0: continue
			
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


'''Finds matching jet for lepton via deltaR.
'''
def findMatch_jet_old(lepton, jets):
	
	drmin = 10
	idx = -1
		
	for ijet, jet in enumerate(jets):
		
		dr = deltaR(lepton.eta(), jet.eta(), lepton.phi(), jet.phi())
		
		if dr < drmin:
			
			drmin = dr
			idx = ijet
				
	return idx, drmin


'''Finds matching gp for track via dxyz.
'''
def findMatch_gen_new(track, genparticles):

	dxyzmin = 10
	tmin = 10
	drmin = 10
	idx = -1

	if track.pt() > 0:

		for igp, gp in enumerate(genparticles):

			if not track.charge() * gp.charge() > 0: continue

			if not abs(track.pt() - gp.pt()) / track.pt() < 0.2: continue

			if not abs(track.eta() - gp.eta()) < 0.1: continue

			if not abs(deltaPhi(track.phi(), gp.phi())) < 1.57: continue

			res = scipy.optimize.minimize(distance, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, (gp.vertex().x(), gp.vertex().y(), gp.vertex().z())))

			dxyz = distance(res.x, track, (gp.vertex().x(), gp.vertex().y(), gp.vertex().z()))

			if dxyz < dxyzmin:

				dxyzmin = dxyz
				tmin = res.x[0]
				drmin = deltaR(gp.eta(), track.eta(), gp.phi(), addPhi(track.phi(), -1. * res.x[0]))
				idx = igp

	return idx, dxyzmin, tmin, drmin


'''Finds matching gp for track via dR.
'''
def findMatch_gen_old(track, genparticles):
	
	drmin = 10
	idx = -1
	
	if track.pt() > 0:
		
		trackTlv = TLorentzVector()
		trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*TMath.CosH(track.eta()))
	
		for igp, gp in enumerate(genparticles):
			
			if not track.charge() * gp.charge() > 0: continue
					
			gpTlv = TLorentzVector()
			gpTlv.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy())
			
			dr = trackTlv.DeltaR(gpTlv)
			
			if dr < drmin:
				
				drmin = dr
				idx = igp
				
	return idx, drmin


'''Finds matching gp for track via dR.
'''
def findMatch_gen_old_easy(track, genparticles):
	
	drmin = 10
	idx = -1
	
	for igp, gp in enumerate(genparticles):
		
		dr = deltaR(track.eta(), gp.eta(), track.phi(), gp.phi())
		
		if dr < drmin:
				
			drmin = dr
			idx = igp
				
	return idx, drmin

###############################################################################################


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

		self.PCA = helix(self.tmin, self.track)
		self.newPhi = addPhi(self.track.phi(), -1.*self.tmin)

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

		res = scipy.optimize.minimize(distance, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(self.track, (self.vtx_x, self.vtx_y, self.vtx_z)))

		return helix(res.x, self.track) - np.array([self.vtx_x, self.vtx_y, self.vtx_z]), res.x

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
