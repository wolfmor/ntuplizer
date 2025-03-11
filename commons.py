#! /usr/bin/env python2


"""Helper functions
"""

import numpy as np
import scipy.optimize
from scipy.linalg import block_diag
from numba import jit

import ROOT

ROOT.gROOT.SetBatch()


def nice_string(ugly_string):
    """Converts strings to be used as branch names."""
    return ugly_string.replace('(', '_').replace(')', '_').replace('/', '_D_').replace('*', '_T_')


class Dummy:
    """Implements dummy physics object class from a ROOT.TLorentzVector."""

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


class DummyJet:
    """Implements dummy jet class."""

    def __init__(self, jet, systs):

        self._pt_raw = jet.pt()
        self._eta = jet.eta()
        self._phi = jet.phi()
        self._px_raw = jet.px()
        self._py_raw = jet.py()
        self._pz_raw = jet.pz()
        self._energy_raw = jet.energy()
        self._numberOfDaughters = jet.numberOfDaughters()

        self._JECdn = systs[0][0]
        self._JEC = systs[0][1]
        self._JECup = systs[0][2]

        self._JERdn = systs[1][0]
        self._JER = systs[1][1]
        self._JERup = systs[1][2]

        self._resolution = systs[2]
        self._stochasticSmearing = systs[3]

    def pt(self, jec=0, jer=0, raw=False):
        if raw: return self._pt_raw
        corr = {-1: self._JECdn, 0: self._JEC, 1: self._JECup}[jec]
        smear = {-1: self._JERdn, 0: self._JER, 1: self._JERup}[jer]
        return self._pt_raw * corr * smear

    def eta(self):
        return self._eta

    def phi(self):
        return self._phi

    def px(self, jec=0, jer=0, raw=False):
        if raw: return self._px_raw
        corr = {-1: self._JECdn, 0: self._JEC, 1: self._JECup}[jec]
        smear = {-1: self._JERdn, 0: self._JER, 1: self._JERup}[jer]
        return self._px_raw * corr * smear

    def py(self, jec=0, jer=0, raw=False):
        if raw: return self._py_raw
        corr = {-1: self._JECdn, 0: self._JEC, 1: self._JECup}[jec]
        smear = {-1: self._JERdn, 0: self._JER, 1: self._JERup}[jer]
        return self._py_raw * corr * smear

    def pz(self, jec=0, jer=0, raw=False):
        if raw: return self._pz_raw
        corr = {-1: self._JECdn, 0: self._JEC, 1: self._JECup}[jec]
        smear = {-1: self._JERdn, 0: self._JER, 1: self._JERup}[jer]
        return self._pz_raw * corr * smear

    def energy(self, jec=0, jer=0, raw=False):
        if raw: return self._energy_raw
        corr = {-1: self._JECdn, 0: self._JEC, 1: self._JECup}[jec]
        smear = {-1: self._JERdn, 0: self._JER, 1: self._JERup}[jer]
        return self._energy_raw * corr * smear

    def numberOfDaughters(self):
        return self._numberOfDaughters


class DummyMet:
    """Implements dummy MET class."""

    def __init__(self, pt_list, phi_list):

        self._pt_raw, self._pt_nom = pt_list[0]
        self._pt_JECdn, self._pt_JECup = pt_list[1]
        self._pt_JERdn, self._pt_JERup = pt_list[2]

        self._phi_raw, self._phi_nom = phi_list[0]
        self._phi_JECdn, self._phi_JECup = phi_list[1]
        self._phi_JERdn, self._phi_JERup = phi_list[2]

    def pt(self, jec=0, jer=0, raw=False):
        if raw: return self._pt_raw
        if abs(jec) + abs(jer) == 0: return self._pt_nom
        elif abs(jec) + abs(jer) > 1:
            raise NotImplementedError('MET variation not implemented for simultaneous JEC and JER variation')
        else:
            if abs(jec) > 0: return {-1: self._pt_JECdn, 1: self._pt_JECup}[jec]
            if abs(jer) > 0: return {-1: self._pt_JERdn, 1: self._pt_JERup}[jer]

    def phi(self, jec=0, jer=0, raw=False):
        if raw: return self._phi_raw
        if abs(jec) + abs(jer) == 0: return self._phi_nom
        elif abs(jec) + abs(jer) > 1:
            raise NotImplementedError('MET variation not implemented for simultaneous JEC and JER variation')
        else:
            if abs(jec) > 0: return {-1: self._phi_JECdn, 1: self._phi_JECup}[jec]
            if abs(jer) > 0: return {-1: self._phi_JERdn, 1: self._phi_JERup}[jer]

    def __str__(self):
        return '\n'.join([
            '\t'.join(['', 'raw', '\tnominal', '\tJECdn', '\tJECup', '\tJERdn', '\tJERup']),
            '\t'.join(['MET p_T'] + [str(p) for p in [self._pt_raw, self._pt_nom, self._pt_JECdn, self._pt_JECup, self._pt_JERdn, self._pt_JERup]]),
            '\t'.join(['MET phi'] + [str(p) for p in [self._phi_raw, self._phi_nom, self._phi_JECdn, self._phi_JECup, self._phi_JERdn, self._phi_JERup]])
        ])


###############################################################################################


def angleBetween(phi1, phi2):
    """Angle between two angles, returns value between -pi and +pi.
    """

    phi = ROOT.TMath.ATan2(ROOT.TMath.Sin(phi1) + ROOT.TMath.Sin(phi2), ROOT.TMath.Cos(phi1) + ROOT.TMath.Cos(phi2))
    
    while phi >= ROOT.TMath.Pi():
        phi -= 2 * ROOT.TMath.Pi()

    while phi < -ROOT.TMath.Pi():
        phi += 2 * ROOT.TMath.Pi()
    
    return phi


@jit(nopython=True)
def addPhi(phi1, phi2):
    """Adds two angles, returns value between -pi and +pi.
    """
    
    sumphi = phi1 + phi2
    
    while sumphi >= np.pi:
        sumphi -= 2 * np.pi

    while sumphi < -np.pi:
        sumphi += 2 * np.pi
    
    return sumphi


@jit(nopython=True)
def deltaPhi(phi1, phi2):
    """Subtracts two angles, returns value between -pi and +pi.
    """
    
    dphi = phi1 - phi2
    
    while dphi >= np.pi:
        dphi -= 2 * np.pi

    while dphi < -np.pi:
        dphi += 2 * np.pi
    
    return dphi


@jit(nopython=True)
def deltaR(eta1, eta2, phi1, phi2):
    """Returns deltaR.
    """

    deta = eta1 - eta2
    dphi = deltaPhi(phi1, phi2)
    
    return np.sqrt(deta * deta + dphi * dphi)


###############################################################################################


def jetID(tag, eta, nhf, nef, chf, cef, mef, nconstituents, cm, nm, lepveto=True):

    # old: https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID13TeVRun2016
    # new: https://twiki.cern.ch/twiki/bin/view/CMS/JetID13TeVUL

    if 'era16' in tag:
        if abs(eta) <= 2.4:

            if lepveto:
                return nhf < 0.9 and \
                       nef < 0.9 and \
                       nconstituents > 1 and \
                       mef < 0.8 and \
                       chf > 0. and \
                       cm > 0. and \
                       cef < 0.8
            else:
                return nhf < 0.9 and \
                       nef < 0.9 and \
                       nconstituents > 1 and \
                       chf > 0. and \
                       cm > 0.

        elif 2.4 < abs(eta) <= 2.7:

            return nhf < 0.9 and \
                   nef < 0.99

            # if lepveto:
            #     return nhf < 0.9 and \
            #            nef < 0.9 and \
            #            nconstituents > 1 and \
            #            mef < 0.8
            # else:
            #     return nhf < 0.9 and \
            #            nef < 0.9 and \
            #            nconstituents > 1

        elif 2.7 < abs(eta) <= 3.:

            return nhf < 0.9 and \
                   0. < nef < 0.99 and \
                   nm > 1

        elif abs(eta) > 3.:

            return nhf > 0.2 and \
                   nef < 0.9 and \
                   nm > 10

        else:

            raise NotImplementedError('what is with the eta?')

    elif 'era17' in tag or 'era18' in tag:
        if abs(eta) <= 2.6:

            if lepveto:
                return nhf < 0.9 and \
                       nef < 0.9 and \
                       nconstituents > 1 and \
                       mef < 0.8 and \
                       chf > 0. and \
                       cm > 0. and \
                       cef < 0.8
            else:
                return nhf < 0.9 and \
                       nef < 0.9 and \
                       nconstituents > 1 and \
                       chf > 0. and \
                       cm > 0.

        elif 2.6 < abs(eta) <= 2.7:

            if lepveto:
                return nhf < 0.9 and \
                       nef < 0.99 and \
                       mef < 0.8 and \
                       cm > 0. and \
                       cef < 0.8
            else:
                return nhf < 0.9 and \
                       nef < 0.99 and \
                       cm > 0.

        elif 2.7 < abs(eta) <= 3.:

            return 0.01 < nef < 0.99 and \
                   nm > 1

        elif abs(eta) > 3.:

            return nhf > 0.2 and \
                   nef < 0.9 and \
                   nm > 10

        else:

            raise NotImplementedError('what is with the eta?')

    # elif 'era18' in tag:
    #     if abs(eta) <= 2.6:
    #
    #         if lepveto:
    #             return nhf < 0.9 and \
    #                    nef < 0.9 and \
    #                    nconstituents > 1 and \
    #                    mef < 0.8 and \
    #                    chf > 0. and \
    #                    cm > 0. and \
    #                    cef < 0.8
    #         else:
    #             return nhf < 0.9 and \
    #                    nef < 0.9 and \
    #                    nconstituents > 1 and \
    #                    chf > 0. and \
    #                    cm > 0.
    #
    #     elif 2.6 < abs(eta) <= 2.7:
    #
    #         if lepveto:
    #             return nhf < 0.9 and \
    #                    nef < 0.99 and \
    #                    mef < 0.8 and \
    #                    cm > 0. and \
    #                    cef < 0.8
    #         else:
    #             return nhf < 0.9 and \
    #                    nef < 0.99 and \
    #                    cm > 0.
    #
    #     elif 2.7 < abs(eta) <= 3.:
    #
    #         return 0.02 < nef < 0.99 and \
    #                nm > 2
    #
    #     elif 3. < abs(eta) <= 5.:
    #
    #         return nhf > 0.2 and \
    #                nef < 0.9 and \
    #                nm > 10
    #
    #     elif abs(eta) > 5.:
    #
    #         return False
    #
    #     else:
    #
    #         raise NotImplementedError('what is with the eta?')

    else:

        raise NotImplementedError('JetID: era unknown or not specified')


def passesEleID(e, pv_pos, rho, wp='veto'):
    """Check if electron passes ID (without isolation requirement)
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/Heppy/python/physicsobjects/Electron.py
    https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedElectronIdentificationRun2#Working_points_for_94X_and_later
    "94X-V2 ID is the recommended ID for all 3 years of Run 2, ie, 2016 legacy rereco, 2017 rereco & UL and 2018 rereco."
    """
    
    if wp == 'veto':
        
        if abs(e.eta()) < 1.4442:  # barrel
            
            if e.full5x5_sigmaIetaIeta() >= 0.0126: return False
            if abs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00463: return False
            if abs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.148: return False
            if e.superCluster().isNonnull(): energy = e.superCluster().energy()
            else: energy = 0
            if energy > 0:
                if e.hadronicOverEm() >= 0.05 + 1.16/energy + 0.0324*rho/energy: return False
            else:
                return False
            # relIsoWithEA < 0.198+0.506/pT	--> don't know what this is
            if e.ecalEnergy() > 0:
                if abs(1.0 / e.ecalEnergy() - e.eSuperClusterOverP() / e.ecalEnergy()) >= 0.209: return False
            else: 
                return False
            if e.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS) > 2: return False
# 			if not e.passConversionVeto() >= : return False  # not implemented
            if abs(e.gsfTrack().dxy(pv_pos)) >= 0.05: return False
            if abs(e.gsfTrack().dz(pv_pos)) >= 0.1: return False
                        
            return True
        
        elif abs(e.eta()) > 1.566:  # endcap

            if e.full5x5_sigmaIetaIeta() >= 0.0457: return False
            if abs(e.deltaEtaSuperClusterTrackAtVtx()) >= 0.00814: return False
            if abs(e.deltaPhiSuperClusterTrackAtVtx()) >= 0.19: return False
            if e.superCluster().isNonnull(): energy = e.superCluster().energy()
            else: energy = 0
            if energy > 0:
                if e.hadronicOverEm() >= 0.05 + 2.54/energy + 0.183*rho/energy: return False
            else:
                return False
            # relIsoWithEA < 0.203+0.963/pT	--> don't know what this is
            if e.ecalEnergy() > 0:
                if abs(1.0 / e.ecalEnergy() - e.eSuperClusterOverP() / e.ecalEnergy()) >= 0.132: return False
            else:
                return False
            if e.gsfTrack().hitPattern().numberOfLostHits(ROOT.reco.HitPattern.MISSING_INNER_HITS) > 3: return False
# 			if not e.passConversionVeto() >= : return False  # not implemented
            if abs(e.gsfTrack().dxy(pv_pos)) >= 0.1: return False
            if abs(e.gsfTrack().dz(pv_pos)) >= 0.2: return False
                        
            return True
        
        else:
            
            return False
        
    else:

        raise NotImplementedError('working point not implemented')


def passesMuonID(m, wp='loose'):
    """Check if muon passes ID (without isolation requirement)
    https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2
    https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/Heppy/python/physicsobjects/Muon.py
    """
    
    if wp == 'loose':
        
        if not m.isPFMuon(): return False
        if not (m.isGlobalMuon() or m.isTrackerMuon()): return False
        
        return True
        
    elif wp == 'medium':
        
        raise NotImplementedError('muon working point not implemented')  # don't know how to access segment compatibility...
        
        if not passesMuonID(m, 'loose'): return False
        if m.innerTrack().validFraction() <= 0.8: return False
        
        goodGlobalMuon = m.isGlobalMuon() \
            and m.globalTrack().normalizedChi2() < 3 \
            and m.combinedQuality().chi2LocalPosition < 12 \
            and m.combinedQuality().trkKink < 20 \
            # and ROOT.reco.Muon.segmentCompatibility(m) > 0.303
        
# 		if not (goodGlobalMuon or ROOT.reco.Muon.segmentCompatibility(m) > 0.451): return False
        
        return True
        
    else:
        
        raise NotImplementedError('working point not implemented')


def passesPhotID(p, wp='loose'):
    """Check if photon passes ID (without isolation requirement)
    https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2#Working_points_for_94X_and_later
    "94X-V2 ID is the recommended ID for all 3 years of Run 2, ie, 2016 legacy rereco, 2017 rereco & UL and 2018 rereco."
    """
    
    if wp == 'loose':
        
        if abs(p.eta()) < 1.4442:  # barrel
            
            if p.hadTowOverEm() >= 0.04596: return False
            if p.sigmaIetaIeta() >= 0.0106: return False
                        
            return True
        
        elif abs(p.eta()) > 1.566:  # endcap
            
            if p.hadTowOverEm() >= 0.0590: return False
            if p.sigmaIetaIeta() >= 0.0272: return False
                        
            return True
        
        else:
            
            return False
        
    else:
        
        raise NotImplementedError('working point not implemented')


###############################################################################################

def calcIso_dBeta(pfIsoVars):
    """Calculates delta-beta corrected isolation.
    https://twiki.cern.ch/twiki/bin/view/CMS/EgammaPFBasedIsolationRun2
    """
    
    neutralIso = pfIsoVars.sumNeutralHadronEt + pfIsoVars.sumPhotonEt
    neutralIsoCorrected = neutralIso - 0.5 * pfIsoVars.sumPUPt
    
    return pfIsoVars.sumChargedHadronPt + max(0., neutralIsoCorrected)


# def calcIso_jet_new(one, many, isTrack=False, btagvalues=None):
#     """Calculates isolation variables for an object inside jets.
#     """
#
#     jetiso = 0
#     jetisomulti = 0
#     jetisobtag = -2.
#
#     dRmin = 999
#     closestjet = None
#     for m in many:
#
#         dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
#
#         if dR < dRmin:
#             dRmin = dR
#             closestjet = m
#
#     if dRmin < 0.4:
#
#         if isTrack: oneTlv = ROOT.TLorentzVector(one.px(), one.py(), one.pz(), one.pt()*ROOT.TMath.CosH(one.eta()))
#         else: oneTlv = ROOT.TLorentzVector(one.px(), one.py(), one.pz(), one.energy())
#         jetTlv = ROOT.TLorentzVector(closestjet.px(), closestjet.py(), closestjet.pz(), closestjet.energy())
#
#         jetiso = (jetTlv - oneTlv).Pt()
#         jetisomulti = closestjet.numberOfDaughters()
#         if btagvalues:
#             for (bt, jetpt, jeteta) in btagvalues:
#                 if np.allclose([jetpt, jeteta], [closestjet.pt(), closestjet.eta()]):
#                     jetisobtag = bt
#                     break
#
#     return jetiso, jetisomulti, dRmin, jetisobtag


# @jitclass({
#     '_px': float32,
#     '_py': float32,
#     '_pz': float32,
#     '_energy': float32,
# })
# class MyTLorentzVector:
#
#     def __init__(self, x, y, z, e, isPtEtaPhiE=False):
#         if isPtEtaPhiE:
#             self._px = x * np.cos(z)
#             self._py = x * np.sin(z)
#             self._pz = x * np.sinh(y)
#         else:
#             self._px = x
#             self._py = y
#             self._pz = z
#         self._energy = e
#
#     def px(self):
#         return self._px
#
#     def py(self):
#         return self._py
#
#     def pz(self):
#         return self._pz
#
#     def energy(self):
#         return self._energy
#
#     def __sub__(self, other):
#         return MyTLorentzVector(self.px() - other.px(),
#                                 self.py() - other.py(),
#                                 self.pz() - other.pz(),
#                                 self.energy() - other.energy())


def converter_jet(func):
    """Converts physics object "one" to tuple to be used in jit compiled function.
    Also accounts for case if many is empty (numba would fail).
    """

    def helper(one, many, isTrack=False, btagvalues=None):
        if len(many) > 0:
            if isTrack:
                return func(np.array((one.pt(), one.eta(), one.phi(), one.pt()*np.cosh(one.eta()))), many, isTrack, btagvalues)
            else:
                return func(np.array((one.pt(), one.eta(), one.phi(), one.energy())), many, isTrack, btagvalues)
        else:
            return 0, 0, 9., -2., -1.

    return helper


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
    
@converter_jet
@jit(nopython=True)
def calcIso_jet_new(one, many, isTrack=False, btagvalues=None):
    """Calculates isolation variables for an object inside jets.
    """

    jetiso = 0
    jetisomulti = 0
    jetisobtag = -2.
    minv = -1.

    dRmin = 9.
    closestjet = None
    for i in range(len(many)):

        dR = deltaR(one[1], many[i][1], one[2], many[i][2])

        if dR < dRmin:
            dRmin = dR
            closestjet = many[i]

    if dRmin < 9.:

        # oneTlv = MyTLorentzVector(one[0], one[1], one[2], one[3], isPtEtaPhiE=True)
        # jetTlv = MyTLorentzVector(closestjet[0], closestjet[1], closestjet[2], closestjet[3], isPtEtaPhiE=True)
        #
        # subTlv = jetTlv - oneTlv
        # jetiso = np.sqrt(subTlv.px()**2 + subTlv.py()**2)  # pT


        # pT of jet with track subtracted = sqrt(pX^2 + pY^2)
        jetiso = np.sqrt((closestjet[0] * np.cos(closestjet[2]) - one[0] * np.cos(one[2]))**2 +
                         (closestjet[0] * np.sin(closestjet[2]) - one[0] * np.sin(one[2]))**2)

        jetisomulti = closestjet[4]

        if btagvalues:
            for (bt, jetpt, jeteta) in btagvalues:
                # if np.allclose([jetpt, jeteta], [closestjet[0], closestjet[1]]):
                # handmade implementation because not in numba:
                if abs(jetpt - closestjet[0]) <= (1.e-8 + 1.e-5 * abs(closestjet[0])) and abs(jeteta - closestjet[1]) <= (1.e-8 + 1.e-5 * abs(closestjet[1])):
                    jetisobtag = bt
                    break

        minv = np.sqrt(2 * one[0] * closestjet[0] * (np.cosh(one[1] - closestjet[1]) - np.cos(one[2] - closestjet[2])))

    return jetiso, jetisomulti, dRmin, jetisobtag, minv


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
# def calcIso_pf_or_track_new(one, many, isMini=False, dontSubtractObject=False):
#     """Calculates isolation variables for an object inside PFs or tracks.
#     """
#
#     conesize = 0.3
#
#     if isMini:
#         if one.pt() <= 50: conesize = 0.2
#         elif one.pt() <= 200: conesize = 10.0/one.pt()
#         else: conesize = 0.05
#
#     ptsum = -one.pt()
#     num = -1
#     if dontSubtractObject:
#         ptsum = 0
#         num = 0
#
#     dRmin = 999
#     for m in many:
#
#         dR = deltaR(one.eta(), m.eta(), one.phi(), m.phi())
#
#         if dR < dRmin and dR > 0.001:
#             dRmin = dR
#
#         if dR < conesize:
#             ptsum += m.pt()
#             num += 1
#
#     return ptsum, ptsum/one.pt(), dRmin, num


def converter(func):
    """Converts physics object "one" to tuple to be used in jit compiled function.
    Also accounts for case if many is empty (numba would fail).
    """

    def helper(one, many, isMini=False, subtractObject=True):
        if len(many) > 0:
            return func(np.array((one.pt(), one.eta(), one.phi())), many, isMini, subtractObject)
        else:
            if subtractObject:
                return -one.pt(), -1, 9, -1
            else:
                return 0, 0, 9, 0


    return helper


def converter_withDrmin2nd(func):
    """Converts physics object "one" to tuple to be used in jit compiled function.
    Also accounts for case if many is empty (numba would fail).
    """

    def helper(one, many, isMini=False, subtractObject=True):
        if len(many) > 0:
            return func(np.array((one.pt(), one.eta(), one.phi())), many, isMini, subtractObject)
        else:
            if subtractObject:
                return -one.pt(), -1, 9, 9, -1
            else:
                return 0, 0, 9, 9, 0


    return helper


	
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


@converter
@jit(nopython=True)
def calcIso_pf_or_track_new(one, many, isMini=False, subtractObject=True):
    """Calculates isolation variables for an object inside PFs or tracks.
    Expects (list of) tuples of the form (pT, eta, phi).
    """

    conesize = 0.3

    if isMini:
        if one[0] <= 50: conesize = 0.2
        elif one[0] <= 200: conesize = 10.0/one[0]
        else: conesize = 0.05

    if subtractObject:
        ptsum = -one[0]
        num = -1
    else:
        ptsum = 0
        num = 0

    dRmin = 9

    for i in range(len(many)):

        dR = deltaR(one[1], many[i][1], one[2], many[i][2])

        if dRmin > dR > 0.001:
            dRmin = dR

        if dR < conesize:
            ptsum += many[i][0]
            num += 1

    return ptsum, ptsum/one[0], dRmin, num


@converter_withDrmin2nd
@jit(nopython=True)
def calcIso_pf_or_track_new_withDrmin2nd(one, many, isMini=False, subtractObject=True):
    """Calculates isolation variables for an object inside PFs or tracks.
    Expects (list of) tuples of the form (pT, eta, phi).
    """

    conesize = 0.3

    if isMini:
        if one[0] <= 50: conesize = 0.2
        elif one[0] <= 200: conesize = 10.0/one[0]
        else: conesize = 0.05

    if subtractObject:
        ptsum = -one[0]
        num = -1
    else:
        ptsum = 0
        num = 0

    dRmin = 9
    dRmin2nd = 9

    for i in range(len(many)):

        dR = deltaR(one[1], many[i][1], one[2], many[i][2])

        if dRmin > dR > 0.001:
            dRmin = dR
        elif dRmin2nd > dR > 0.001:
            dRmin2nd = dR

        if dR < conesize:
            ptsum += many[i][0]
            num += 1

    return ptsum, ptsum/one[0], dRmin, dRmin2nd, num


###############################################################################################

def getLastCopyStatusOne(gp):
    """Returns last copy of a gen. particle with status 1 = final state.
    """
        
    i = 0
    while not gp.status() == 1:
        for idx in range(gp.numberOfDaughters()):
            if gp.pdgId() == gp.daughter(idx).pdgId():
                gp = gp.daughter(idx)
                break
        i += 1
        if i > 100: return None
        
    return gp


def getLastCopy(gp):
    """Returns last copy of a gen. particle.
    """
        
    i = 0
    while not gp.isLastCopy():
        for idx in range(gp.numberOfDaughters()):
            if gp.pdgId() == gp.daughter(idx).pdgId():
                gp = gp.daughter(idx)
                break
        i += 1
        if i > 100: return None
        
    return gp


def findLeptonNeutrino(gp):
    """Returns lepton and neutrino from chi1pm decay to chi10.
    """
        
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
    

def findPion(gp):
    """Returns pion from chi1pm decay to chi10.
    """
    
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
    

def findLeptons(gp):
    """Returns lepton pair from chi20 decay to chi10.
    """
    
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


def findDaughters(gp):
    """Returns list of status 1 (grand)daughters.
    """
    
    daughters = []
    
    for i in range(gp.numberOfDaughters()):
        
        daughter = getLastCopy(gp.daughter(i))
        
        if daughter.status() == 1:
            
            daughters.append(daughter)
                    
        else:
            
            for granddaughter in findDaughters(daughter):
                
                daughters.append(granddaughter)
    
    return daughters



"""Finds matching track for lepton via delta R.
"""
def findMinDr_track(aTrack, tracks, threshold):

    drmin = 10
    idx = -1
    matchingTrack = None
    match = False

    if aTrack.pt() > 0:

        for itrack, track in enumerate(tracks):
            
            if track == None: continue
            
            if not passesPreselection_basic_track(track): continue

            
            if not track.charge() * aTrack.charge() > 0: continue
            
            aTrackTlv = ROOT.TLorentzVector()
            aTrackTlv.SetPxPyPzE(aTrack.px(),aTrack.py(),aTrack.pz(),aTrack.pt()*np.cosh(aTrack.eta()))	
            trkTlv = ROOT.TLorentzVector()
            trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*np.cosh(track.eta()))
            
            dr = trkTlv.DeltaR(aTrackTlv)
            
            if dr < drmin:
                    
                drmin = dr
                idx = itrack
                matchingTrack = track
                

    if drmin < threshold: match = True
    return match, idx, drmin, matchingTrack


def getTauDecayMode(tau, decaymode, decaymode_alt):

    nprongs = 0
    nneutral = 0

    nprongs_alt = 0
    nneutral_alt = 0

    for k in range(tau.numberOfDaughters()):

        taudaughterpdgid = abs(tau.daughter(k).pdgId())
        taudaughtercharge = abs(tau.daughter(k).charge())

        # see: https://twiki.cern.ch/twiki/bin/viewauth/CMS/TauIDRecommendationForRun2#Decay_Mode_Reconstruction
        if taudaughterpdgid in [12, 14, 16]:  # skip the neutrinos
            continue
        elif taudaughterpdgid in [11, 13]:  # charged lepton
            decaymode = 10 + taudaughterpdgid
            decaymode_alt = 10 + taudaughterpdgid
        else:

            if taudaughterpdgid == 111:  # neutral pion
                nneutral += 1
            else:  # has to be a charged hadron (but has it? somehow this gives 2-prong taus...)
                nprongs += 1

            if taudaughtercharge > 0:
                nprongs_alt += 1
            else:
                nneutral_alt += 1

    if decaymode not in [21, 23]:
        decaymode = 5 * (nprongs - 1) + nneutral
        decaymode_alt = 5 * (nprongs_alt - 1) + nneutral_alt

    return decaymode, decaymode_alt


###############################################################################################

def passesPreselection_basic_pfc(pfc):
    """Defines basic preselection for PFs.
    """
    
    if pfc.trackRef().isNull(): return False
    
    if not passesPreselection_basic_track(pfc.trackRef().get()): return False
    
    return True


def passesPreselection_iso_chpf(pfc, pv_pos, dz_threshold, dxy_threshold, pt_threshold):
    """Defines preselection used for calculation of isolation-variables for PFs.
    """
    
    if pfc.charge() == 0: return False
    
    if pfc.trackRef().isNull(): return False
    
    if pfc.trackRef().get().numberOfValidHits() == 0: return False
    
    if pfc.trackRef().get().ndof() == 0: return False
    
    if pfc.trackRef().get().charge() == 0: return False
    
    if abs(pfc.trackRef().get().dz(pv_pos)) >= dz_threshold: return False

    if abs(pfc.trackRef().get().dxy(pv_pos)) >= dxy_threshold: return False

    if pfc.pt() <= pt_threshold: return False
    
    return True


def passesPreselection_iso_pf(pfc, pt_threshold):
    """Defines preselection used for calculation of isolation-variables for PFs.
    """

    if pfc.pt() <= pt_threshold: return False

    return True


# def passesPreselection_final_pfc(pfc, pv_pos):
#     """Defines preselection used for BDT for PFs.
#     """
#
#     if not passesPreselection_basic_pfc(pfc): return False
#
#     if pfc.pt() <= 0.7: return False
#
#     return True


def passesPreselection_basic_track(track):
    """Defines basic preselection for tracks.
    """
    
    if track.numberOfValidHits() == 0: return False
    
    if track.ndof() == 0: return False
    
    if track.charge() == 0: return False
    
    # if track.pt() > 5: return False
    
    return True


def passesPreselection_iso_track(track, pv_pos, dz_threshold, dxy_threshold, pt_threshold):
    """Defines preselection used for calculation of isolation-variables for tracks.
    """
        
    if track.numberOfValidHits() == 0: return False
    
    if track.ndof() == 0: return False
    
    if track.charge() == 0: return False
    
    if abs(track.dz(pv_pos)) >= dz_threshold: return False

    if abs(track.dxy(pv_pos)) >= dxy_threshold: return False

    if track.pt() <= pt_threshold: return False
            
    return True


# def passesPreselection_final_track(track, pv_pos):
#     """Defines preselection used for BDT for tracks.
#     """
#
#     if not passesPreselection_basic_track(track): return False
#
#     return True


def passesPreselection_iso_jet(jet, pt_threshold):
    """Defines preselection used for calculation of isolation-variables for jets.
    """
    
    if jet.pt() <= pt_threshold: return False

    return True


###############################################################################################

def vertexFinder(track1, track2, pv_pos):
    """Finds a common vertex for two tracks if the minimin distance is smaller than 0.1 cm
    and returns the vertex' distance to the PV.
    """
    
    mindist, tmin1, tmin2, p1, p2 = minDistanceTrackTrack(track1, track2)
    
    vertex = np.array([0, 0, 0])
    angletoorigin = -1
    angletopv = -1
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


def minDistanceTrackTrack(track1, track2):
    """Gives minimum distance between two tracks and the corresponding points of closest approach on each track.
    """
    
    res = scipy.optimize.minimize(distanceTrackTrack, x0=np.array([0.0, 0.0]), bounds=((-1.57, 1.57), (-1.57, 1.57)), args=(track1, track2))
    
    tmin1 = res.x[0]
    tmin2 = res.x[1]
    
    d, p1, p2 = getDistanceAndPoints(res.x, track1, track2)
    
    return d, tmin1, tmin2, p1, p2


def getDistanceAndPoints(tarray, track1, track2):
    """Distance between a specific point along a track (given by t1) and a specific point along another track (given by t2).
    Additionally returns the points.
    """
    
    t1 = np.array([tarray[0]])
    t2 = np.array([tarray[1]])
    
    p1 = helixOld(t1, track1.phi(), track1.eta(), track1.charge(), track1.pt(), track1.vx(), track1.vy(), track1.vz())
    p2 = helixOld(t2, track2.phi(), track2.eta(), track2.charge(), track2.pt(), track2.vx(), track2.vy(), track2.vz())
    
    d = np.linalg.norm(p1 - p2)
    
    return d, p1, p2


def distanceTrackTrack(tarray, track1, track2):
    """Distance between a specific point along a track (given by t1) and a specific point along another track (given by t2).
    """
    
    t1 = np.array([tarray[0]])
    t2 = np.array([tarray[1]])
    
    d = np.linalg.norm(helixOld(t1, track1.phi(), track1.eta(), track1.charge(), track1.pt(), track1.vx(), track1.vy(), track1.vz())
                     - helixOld(t2, track2.phi(), track2.eta(), track2.charge(), track2.pt(), track2.vx(), track2.vy(), track2.vz()))
    
    return d


@jit(nopython=True)
def helixOld(t, phi, eta, q, pt, vx, vy, vz):
    """Helix parametrization in 3D. Same result as helix(...) but not using CMSSW track parameters.
    """
    
    r = 87.78  # radius [cm] for particle with pT=1GeV in B=3.8T
    
    x = vx + r * q * pt * (np.sin(phi) - np.sin(phi - t))
    
    y = vy + r * q * pt * (-np.cos(phi) + np.cos(phi - t))
    
    z = vz + t * r * q * pt / np.tan(2 * np.arctan(np.exp(-eta)))
    
    return np.concatenate((x, y, z))


@jit(nopython=True)
def distanceOld(t, phi, eta, q, pt, vx, vy, vz, vertex):
    """Distance between a specific point along a track (given by t) and a vertex in 3D.
    """

    d = np.linalg.norm(helixOld(t, phi, eta, q, pt, vx, vy, vz)
                       - np.array(vertex))
    
    return d


def distanceXYOld(t, track, vertex):
    """Distance in XY between a specific point along a track (given by t) and a vertex.
    """
    
    d = np.linalg.norm(helixOld(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())[:2]
                       - np.array([vertex.x(), vertex.y()]))
    
    return d


def distanceZOld(t, track, vertex):
    """Distance in Z between a specific point along a track (given by t) and a vertex.
    """
    
    d = np.linalg.norm(helixOld(t, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())[2:]
                       - np.array([vertex.z()]))
    
    return d


def handmadeDxyDz(track, pv_pos):
    """Impact parameters in XY and Z at PCA to PV in 3D (calculated with helix extrapolation).
    """
    
    res = scipy.optimize.minimize(distanceOld,
                                  x0=np.array([0.0]),
                                  bounds=((-1.57, 1.57),),
                                  args=(track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz(),
                                        (pv_pos.x(), pv_pos.y(), pv_pos.z())))
    
    dxy = distanceXYOld(res.x, track, pv_pos)
    dz = distanceZOld(res.x, track, pv_pos)
    
    return dxy, dz


def handmadeDxyDzTransversePCA(track, pv_pos):
    """Impact parameters in XY and Z at PCA to PV in XY (calculated with helix extrapolation).
    """
    
    res = scipy.optimize.minimize(distanceXYOld, x0=np.array([0.0]), bounds=((-1.57, 1.57),), args=(track, pv_pos))
    
    dxy = distanceXYOld(res.x, track, pv_pos)
    dz = distanceZOld(res.x, track, pv_pos)
    
    return dxy, dz


def handmadeDphiMetPCA(track, pv_pos, met):
    """dPhi between MET and vector from PV to PCA to PV of a track
    """
    
    # gives the helix parameter value for which the distance is minimal
    res = scipy.optimize.minimize(distanceOld,
                                  x0=np.array([0.0]),
                                  bounds=((-1.57, 1.57),),
                                  args=(track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz(),
                                        (pv_pos.x(), pv_pos.y(), pv_pos.z())))
    
    # insert found helix parameter and track parameters to get PCA
    pca = helixOld(res.x, track.phi(), track.eta(), track.charge(), track.pt(), track.vx(), track.vy(), track.vz())
    
    # compute phi angle of vector from PV to PCA
    phipca = np.arctan2(pca[1] - pv_pos.y(), pca[0] - pv_pos.x())
    
    # compute dPhi
    dphi = deltaPhi(phipca, met.phi())
    
    return dphi, phipca, res.x[0]


@jit(nopython=True)
def distance(t, qoverp, lmbd, phi, refX, refY, refZ, vertex):

    return np.linalg.norm(helix(t, qoverp, lmbd, phi, refX, refY, refZ) - np.array(vertex))


# (float64[:](float64[:], float64, float64, float64, float64, float64, float64))
@jit(nopython=True)
def helix(t, qoverp, lmbd, phi, refX, refY, refZ):
    """Helix parametrization using CMSSW track parameters.

    Particles with positive (negative) charge travel along their helix in direction of positive (negative) t parameter.
    To get the azimuthal angle phi for a given t parameter, -t has to be added to phi0 (same for positive and negative charges).
    Apart from phi, the track parameters don't change.

    Mostly taken from Lucas Wiens (Higgs -> tau tau CP analysis: CMS AN-19-192),
    but simplified trigonometric phi term and respecting track charge (differnt direction of rotation).
    """

    # conversion factor to convert from T to GeV/(e * cm) = c [m/s] * 10^-9 * 10^-2
    # 1 T = 1 V s / m^2 = 1 eV s / (e m^2) = 1 s/m * eV / (e m) = c * eV / (e m)
    B = 3.8 * 0.0029979

    x = refX + np.cos(lmbd) / B / qoverp * (np.sin(phi) - np.sin(phi - t))

    y = refY + np.cos(lmbd) / B / qoverp * (-np.cos(phi) + np.cos(phi - t))

    z = refZ + t * np.sin(lmbd) / B / qoverp

    # qoverp = track.parameter(0)
    # lmbd = track.parameter(1)
    # phi = track.parameter(2)
    #
    # x = track.referencePoint().x() + np.cos(lmbd) / B / qoverp * (np.sin(phi) - np.sin(phi - t))
    #
    # y = track.referencePoint().y() + np.cos(lmbd) / B / qoverp * (-np.cos(phi) + np.cos(phi - t))
    #
    # z = track.referencePoint().z() + t * np.sin(lmbd) / B / qoverp

    return np.concatenate((x, y, z))


###############################################################################################
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
            
            leptonTlv = ROOT.TLorentzVector()
            leptonTlv.SetPxPyPzE(lepton.px(),lepton.py(),lepton.pz(),lepton.energy())	
            trkTlv = ROOT.TLorentzVector()
            trkTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*np.cosh(track.eta()))
            
            dr = trkTlv.DeltaR(leptonTlv)
            
            if dr < drmin:
                    
                drmin = dr
                idx = itrack
                matchingTrack = track
                

    if drmin < threshold: match = True
    return match, idx, drmin, matchingTrack
    
"""Foo
"""
def findMinDr_ancestors(svdaughter, gps, ignoreIdx=[None]):

    pdgId = -1
    pdgIdMother = -1
    drmin = 999
    idx = -1

    tlvMother = ROOT.TLorentzVector()	
    genmatchmother = None

    relatives = [None, None]
    charge = 0

    numDaughtersOfMother = 0

    if svdaughter.pt() > 0:

        for igp, gp in enumerate(gps):
                
            if igp in ignoreIdx: continue
            
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

                numDaughtersOfMother = gp.mother().numberOfDaughters()

                
    hasEWancestors = 0

    if not (genmatchmother == None):
        iteration = 0
        
        while genmatchmother.numberOfMothers() > 0:			

            if (genmatchmother.pdgId())>1000000:
                hasEWancestors = 1
                break
            
            genmatchmother = genmatchmother.mother(0)
            iteration += 1
            
    return idx,  pdgId, pdgIdMother, drmin, tlvMother, hasEWancestors, numDaughtersOfMother

	
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
    tlvMother = ROOT.TLorentzVector()	
    if svdaughter.pt() > 0:

        for igp, gp in enumerate(gps):
            
            if igp in ignoreIdx: continue
            
            if not gp.isLastCopy(): continue

            if not gp.charge()== svdaughter.charge(): continue
                                        
            if not abs(gp.eta() - svdaughter.eta()) < 0.1: continue
            
            if not abs(deltaPhi(gp.phi(), svdaughter.phi())) < 1.57: continue
            
            res = scipy.optimize.minimize(distanceXYOld, x0=0.0, bounds=((-1.57, 1.57),), args=(svdaughter, gp.vertex()))  # other minimization method?

            dxyz = distanceXYOld(res.x, svdaughter, gp.vertex())
            if dxyz < dxyzmin:
                    
                dxyzmin = dxyz
                drmin = deltaR(svdaughter.eta(), gp.eta(), gp.phi(), addPhi(svdaughter.phi(), res.x[0]))
                idx = igp
                pdgId = gp.pdgId()
                pdgIdMother = gp.mother().pdgId()
                tlvMother.SetPxPyPzE(gp.mother().px(),gp.mother().py(),gp.mother().pz(),gp.mother().energy())
                genmatchmother = gp.mother()
                charge = gp.charge()

                numDaughtersOfMother  = gp.mother().numberOfDaughters()

                
    hasEWancestors = 0
    if not (genmatchmother == None):
        iteration = 0
        
        while genmatchmother.numberOfMothers() > 0:			

            if (genmatchmother.pdgId())>1000000:
                hasEWancestors = 1
                break
            
            genmatchmother = genmatchmother.mother(0)
            iteration += 1
    return idx,  dxyzmin,pdgId, pdgIdMother, drmin, tlvMother, hasEWancestors, numDaughtersOfMother


def findMatch_track_new(lepton, tracks):
    """Finds matching track for lepton/pion/... via dxyz.
    """

    dxyzmin = 9
    tmin = 9
    drmin = 9
    idx = -1

    if lepton.pt() > 0:

        for itrack, track in enumerate(tracks):

            if not passesPreselection_basic_track(track): continue

            if not track.charge() * lepton.charge() > 0: continue

            if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue

            if not abs(track.eta() - lepton.eta()) < 0.1: continue

            if not abs(deltaPhi(track.phi(), lepton.phi())) < 1.57: continue

            res = scipy.optimize.minimize(distance,
                                          x0=np.array([0.0]),
                                          bounds=((-1.57, 1.57),),
                                          args=(track.parameter(0), track.parameter(1), track.parameter(2),
                                                track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                                                (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z())))

            dxyz = distance(res.x, track.parameter(0), track.parameter(1), track.parameter(2),
                            track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                            (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z()))

            if dxyz < dxyzmin:

                dxyzmin = dxyz
                tmin = res.x[0]
                drmin = deltaR(lepton.eta(), track.eta(), lepton.phi(), addPhi(track.phi(), float(-1. * res.x[0])))
                idx = itrack

    return idx, dxyzmin, tmin, drmin


def findMatch_track_new_random(lepton, tracks):
    """Finds matching track for lepton via dxyz, random with opposite charge.
    """

    dxyzmin = 9
    tmin = 9
    drmin = 9
    idx = -1

    if lepton.pt() > 0:

        for itrack, track in enumerate(tracks):

            if not passesPreselection_basic_track(track): continue

            if not track.charge() * lepton.charge() < 0: continue

            if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue

            if not abs(track.eta() - lepton.eta()) < 0.1: continue

            if not abs(deltaPhi(track.phi(), lepton.phi())) < 1.57: continue

            res = scipy.optimize.minimize(distance,
                                          x0=np.array([0.0]),
                                          bounds=((-1.57, 1.57),),
                                          args=(track.parameter(0), track.parameter(1), track.parameter(2),
                                                track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                                                (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z())))

            dxyz = distance(res.x, track.parameter(0), track.parameter(1), track.parameter(2),
                            track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                            (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z()))

            if dxyz < dxyzmin:

                dxyzmin = dxyz
                tmin = res.x[0]
                drmin = deltaR(lepton.eta(), track.eta(), lepton.phi(), addPhi(track.phi(), float(-1. * res.x[0])))
                idx = itrack

    return idx, dxyzmin, tmin, drmin



"""Finds matching track for track via dxyz.
"""
def findMatch_tracktrack_new(aTrack, tracks):

    dxyzmin = 999
    tmin = 10
    drmin = 999
    idx = -1

    if aTrack.pt() > 0:

        for itrack, track in enumerate(tracks):

            if track == None: continue

            if not passesPreselection_basic_track(track): continue

            if not track.charge() * aTrack.charge() > 0: continue
            ##
            if not abs(track.pt() - aTrack.pt()) / aTrack.pt() < 0.2: continue
                    
            if not abs(track.eta() - aTrack.eta()) < 0.1: continue

            if not abs(deltaPhi(track.phi(), aTrack.phi())) < 1.57: continue
            
            res = scipy.optimize.minimize(distance, x0=np.array([0.0]), bounds=((-1.57, 1.57),),args=(track.parameter(0), track.parameter(1), track.parameter(2),
                                                track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                                                (aTrack.vertex().x(), aTrack.vertex().y(), aTrack.vertex().z())))  # other minimization method?
            dxyz = distance(res.x, track.parameter(0), track.parameter(1), track.parameter(2),
                            track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                            (aTrack.vertex().x(), aTrack.vertex().y(), aTrack.vertex().z()))
            
            if dxyz < dxyzmin:
                    
                dxyzmin = dxyz
                tmin = res.x[0]
                drmin = deltaR(aTrack.eta(), track.eta(), aTrack.phi(), addPhi(track.phi(), res.x[0]))
                idx = itrack
                
    return idx, dxyzmin, tmin, drmin


def findMatch_pfc_new(lepton, pfcands):
    """Finds matching pfc for lepton via dxyz.
    """

    dxyzmin = 9
    tmin = 9
    drmin = 9
    idx = -1

    if lepton.pt() > 0:

        for ipfc, pfc in enumerate(pfcands):

            if not passesPreselection_basic_pfc(pfc): continue

            if not pfc.charge() * lepton.charge() > 0: continue

            if not abs(pfc.pt() - lepton.pt()) / lepton.pt() < 0.2: continue

            if not abs(pfc.eta() - lepton.eta()) < 0.1: continue

            if not abs(deltaPhi(pfc.phi(), lepton.phi())) < 1.57: continue

            track = pfc.trackRef().get()

            res = scipy.optimize.minimize(distance,
                                          x0=np.array([0.0]),
                                          bounds=((-1.57, 1.57),),
                                          args=(track.parameter(0), track.parameter(1), track.parameter(2),
                                                track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                                                (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z())))

            dxyz = distance(res.x, track.parameter(0), track.parameter(1), track.parameter(2),
                            track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                            (lepton.vertex().x(), lepton.vertex().y(), lepton.vertex().z()))

            if dxyz < dxyzmin:

                dxyzmin = dxyz
                tmin = res.x[0]
                drmin = deltaR(lepton.eta(), pfc.eta(), lepton.phi(), addPhi(pfc.phi(), float(-1. * res.x[0])))
                idx = ipfc

    return idx, dxyzmin, tmin, drmin


def findMatch_track_old(lepton, tracks):
    """Finds matching track for lepton via dR.
    """
    
    drmin = 9
    idx = -1
    
    if lepton.pt() > 0:
        
        leptonTlv = ROOT.TLorentzVector()
        leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
    
        for itrack, track in enumerate(tracks):
            
            if track.numberOfValidHits() == 0: continue
            if track.ndof() == 0: continue
            if track.charge() == 0: continue
        
            # if not passesPreselection_basic_track(track): continue
            
            if not track.charge() * lepton.charge() > 0: continue
            
            # if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
                    
            trackTlv = ROOT.TLorentzVector()
            trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*ROOT.TMath.CosH(track.eta()))
            
            dr = trackTlv.DeltaR(leptonTlv)
            
            if dr < drmin:
                
                drmin = dr
                idx = itrack
                
    return idx, drmin


def findMatch_track_old_random(lepton, tracks):
    """Finds matching track for lepton via dR, random with opposite charge.
    """
    
    drmin = 9
    idx = -1
    
    if lepton.pt() > 0:
        
        leptonTlv = ROOT.TLorentzVector()
        leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
    
        for itrack, track in enumerate(tracks):
                        
            if track.numberOfValidHits() == 0: continue
            if track.ndof() == 0: continue
            if track.charge() == 0: continue
            
            # if not passesPreselection_basic_track(track): continue
            
            if not track.charge() * lepton.charge() < 0: continue
            
            # if not abs(track.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
                    
            trackTlv = ROOT.TLorentzVector()
            trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*ROOT.TMath.CosH(track.eta()))
            
            dr = trackTlv.DeltaR(leptonTlv)
            
            if dr < drmin:
                
                drmin = dr
                idx = itrack
                
    return idx, drmin


def findMatch_pfc_old(lepton, pfcands):
    """Finds matching pfc for lepton via deltaR.
    """
    
    drmin = 9
    idx = -1
    
    if lepton.pt() > 0:
        
        leptonTlv = ROOT.TLorentzVector()
        leptonTlv.SetPxPyPzE(lepton.px(), lepton.py(), lepton.pz(), lepton.energy())
        
        for ipfc, pfc in enumerate(pfcands):
            
            # if pfc.trackRef().isNull(): continue
            # if pfc.trackRef().get().numberOfValidHits() == 0: continue
            # if pfc.trackRef().get().ndof() == 0: continue
            # if pfc.trackRef().get().charge() == 0: continue
            
            # if not passesPreselection_basic_pfc(pfc): continue
            
            if not pfc.charge() * lepton.charge() > 0: continue
            
            # if not abs(pfc.pt() - lepton.pt()) / lepton.pt() < 0.2: continue
            
            pfcTlv = ROOT.TLorentzVector()
            pfcTlv.SetPxPyPzE(pfc.px(), pfc.py(), pfc.pz(), pfc.energy())
            
            dr = pfcTlv.DeltaR(leptonTlv)
            
            if dr < drmin:
                
                drmin = dr
                idx = ipfc
                
    return idx, drmin


def findMatch_jet_old(lepton, jets):
    """Finds matching jet for lepton via deltaR.
    """
    
    drmin = 9
    idx = -1
        
    for ijet, jet in enumerate(jets):
        
        dr = deltaR(lepton.eta(), jet.eta(), lepton.phi(), jet.phi())
        
        if dr < drmin:
            
            drmin = dr
            idx = ijet
                
    return idx, drmin


def findMatch_gen_new(track, genparticles):
    """Finds matching gp for track via dxyz.
    """

    dxyzmin = 9
    tmin = 9
    drmin = 9
    idx = -1

    if track.pt() > 0:

        for igp, gp in enumerate(genparticles):

            if not track.charge() * gp.charge() > 0: continue

            if not abs(track.pt() - gp.pt()) / track.pt() < 0.2: continue

            if not abs(track.eta() - gp.eta()) < 0.1: continue

            if not abs(deltaPhi(track.phi(), gp.phi())) < 1.57: continue

            res = scipy.optimize.minimize(distance,
                                          x0=np.array([0.0]),
                                          bounds=((-1.57, 1.57),),
                                          args=(track.parameter(0), track.parameter(1), track.parameter(2),
                                                track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                                                (gp.vertex().x(), gp.vertex().y(), gp.vertex().z())))

            dxyz = distance(res.x, track.parameter(0), track.parameter(1), track.parameter(2),
                            track.referencePoint().x(), track.referencePoint().y(), track.referencePoint().z(),
                            (gp.vertex().x(), gp.vertex().y(), gp.vertex().z()))

            if dxyz < dxyzmin:

                dxyzmin = dxyz
                tmin = res.x[0]
                drmin = deltaR(gp.eta(), track.eta(), gp.phi(), addPhi(track.phi(), float(-1. * res.x[0])))
                idx = igp

    return idx, dxyzmin, tmin, drmin


def findMatch_gen_old(track, genparticles):
    """Finds matching gp for track via dR.
    """
    
    drmin = 9
    idx = -1
    
    if track.pt() > 0:
        
        trackTlv = ROOT.TLorentzVector()
        trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*ROOT.TMath.CosH(track.eta()))
    
        for igp, gp in enumerate(genparticles):
            
            if not track.charge() * gp.charge() > 0: continue
                    
            gpTlv = ROOT.TLorentzVector()
            gpTlv.SetPxPyPzE(gp.px(), gp.py(), gp.pz(), gp.energy())
            
            dr = trackTlv.DeltaR(gpTlv)
            
            if dr < drmin:
                
                drmin = dr
                idx = igp
                
    return idx, drmin


def findMatch_gen_old_easy(track, genparticles):
    """Finds matching gp for track via dR.
    """
    
    drmin = 9
    idx = -1
    
    for igp, gp in enumerate(genparticles):
        
        dr = deltaR(track.eta(), gp.eta(), track.phi(), gp.phi())
        
        if dr < drmin:
                
            drmin = dr
            idx = igp
                
    return idx, drmin

	
"""Finds matching muon for track via dR.
"""
def matchToMuon(track, muons):
	
    drmin = 10
    idx = -1

    if track.pt() > 0:
        
        trackTlv = ROOT.TLorentzVector()
        trackTlv.SetPxPyPzE(track.px(), track.py(), track.pz(), track.pt()*ROOT.TMath.CosH(track.eta()))

        for imuon, muon in enumerate(muons):
            if muon.charge() == 0: continue
            if not muon.charge() * track.charge() > 0: continue
                    
            muonTlv = ROOT.TLorentzVector()
            muonTlv.SetPxPyPzE(muon.px(), muon.py(), muon.pz(), muon.pt()*ROOT.TMath.CosH(muon.eta()))
            
            dr = muonTlv.DeltaR(trackTlv)
            
            if dr < drmin:
                
                drmin = dr
                idx = imuon
                
    return idx, drmin



###############################################################################################

class IPcalculator:

    def __init__(self, track, vertex, verbose=False):

        self.verbose = verbose

        self.track = track

        self.vtx_x = vertex.x()
        self.vtx_y = vertex.y()
        self.vtx_z = vertex.z()
        self.vtx_covariance = np.array([[vertex.covariance(i, j) for j in range(3)] for i in range(3)])

        # 0:qoverp, 1:lambda, 2:phi, 3:dxy, 4:dsz
        self.track_parameters = np.array([track.parameter(i) for i in range(5)])
        self.track_covariance = np.array([[track.covariance(i, j) for j in range(5)] for i in range(5)])

        self.IPvector, self.tmin = self.calculateIPvector()
        self.IPcovariance = self.calculateIPcovariance()

        self.PCA = helix(self.tmin, self.track.parameter(0), self.track.parameter(1), self.track.parameter(2), self.track.referencePoint().x(), self.track.referencePoint().y(), self.track.referencePoint().z())
        self.newPhi = addPhi(self.track.phi(), float(-1.*self.tmin))

        self.dxyz = np.linalg.norm(self.IPvector)
        self.dxy = np.linalg.norm(self.IPvector[:2])
        self.dz = np.linalg.norm(self.IPvector[2:])

        if self.dxyz == 0: self.IPvectorNorm = self.IPvector
        else: self.IPvectorNorm = self.IPvector / self.dxyz

        self.IPerror = np.sqrt(np.dot(self.IPvectorNorm, np.dot(self.IPcovariance, self.IPvectorNorm)))

        self.IPerrorXY = np.sqrt(np.dot(self.IPvectorNorm[:2], np.dot(self.IPcovariance[:2, :2], self.IPvectorNorm[:2])))
        self.IPerrorZ = np.sqrt(np.dot(self.IPvectorNorm[2:], np.dot(self.IPcovariance[2:, 2:], self.IPvectorNorm[2:])))

        if self.verbose:
            print('vertex', (self.vtx_x, self.vtx_y, self.vtx_z))

            print('track parameters')
            print(self.track_parameters)
            print('track covariance matrix')
            print(self.track_covariance)

            print('IPvector', self.IPvector)
            print('IPvectorNorm', self.IPvectorNorm)

            print('IPcovariance', self.IPcovariance)
            print('IPerror', self.IPerror)

            print('IPerrorXY', np.sqrt(np.dot(self.IPvectorNorm[:2], np.dot(self.IPcovariance[:2, :2], self.IPvectorNorm[:2]))))
            print('IPerrorZ', np.sqrt(np.dot(self.IPvectorNorm[2:], np.dot(self.IPcovariance[2:, 2:], self.IPvectorNorm[2:]))))

    def calculateIPvector(self):

        res = scipy.optimize.minimize(distance,
                                      x0=np.array([0.0]),
                                      bounds=((-1.57, 1.57),),
                                      args=(self.track.parameter(0), self.track.parameter(1), self.track.parameter(2),
                                            self.track.referencePoint().x(), self.track.referencePoint().y(), self.track.referencePoint().z(),
                                            (self.vtx_x, self.vtx_y, self.vtx_z)))

        return helix(res.x, self.track.parameter(0), self.track.parameter(1), self.track.parameter(2), self.track.referencePoint().x(), self.track.referencePoint().y(), self.track.referencePoint().z()) - np.array([self.vtx_x, self.vtx_y, self.vtx_z]), res.x

    def calculateIPcovariance(self):

        qoverp = self.track_parameters[0]
        lmbd = self.track_parameters[1]
        phi = self.track_parameters[2]
        dxy = self.track_parameters[3]
        dsz = self.track_parameters[4]

        HelixJacobian = IPcalculator.calculateHelixJacobian(qoverp, lmbd, phi, dxy, dsz, self.tmin)

        HelixCovariance = np.dot(HelixJacobian, np.dot(self.track_covariance, HelixJacobian.T))

        HelixAndPVcovariance = block_diag(HelixCovariance, self.vtx_covariance)

        IPjacobian = np.zeros((3, 6))
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
            print('HelixJacobian')
            print(HelixJacobian)
            print('HelixCovariance')
            print(HelixCovariance)
            print('HelixAndPVcovariance')
            print(HelixAndPVcovariance)
            print('IPjacobian')
            print(IPjacobian)

        return np.dot(IPjacobian, np.dot(HelixAndPVcovariance, IPjacobian.T))

    @staticmethod
    @jit(nopython=True)
    def calculateHelixJacobian(qoverp, lmbd, phi, dxy, dsz, tmin):
        """Moved to stand-alone function in order to optimize with jit
        """

        # conversion factor to convert from T to GeV/(e * cm) = c [m/s] * 10^-9 * 10^-2
        # 1 T = 1 V s / m^2 = 1 eV s / (e m^2) = 1 s/m * eV / (e m) = c * eV / (e m)
        B = 3.8 * 0.0029979
        
        HelixJacobian = np.zeros((3, 5))
        # derivatives of the helix coordinates with respect to
        # qoverp
        HelixJacobian[0][0] = (- (np.cos(lmbd) * (np.sin(phi) - np.sin(phi - tmin))) / (B * qoverp**2.))[0]  # need this indexing for jit
        HelixJacobian[1][0] = (- (np.cos(lmbd) * (-np.cos(phi) + np.cos(phi - tmin))) / (B * qoverp**2.))[0]
        HelixJacobian[2][0] = (- (np.sin(lmbd) * tmin) / (B * qoverp**2.))[0]
        # lambda
        HelixJacobian[0][1] = (- (np.sin(lmbd) * (np.sin(phi) - np.sin(phi - tmin))) / (B * qoverp))[0]
        HelixJacobian[1][1] = (- (np.sin(lmbd) * (-np.cos(phi) + np.cos(phi - tmin))) / (B * qoverp))[0]
        HelixJacobian[2][1] = ((dsz * np.tan(lmbd)) / (np.cos(lmbd)) + (np.cos(lmbd) * tmin) / (B * qoverp))[0]
        # phi
        HelixJacobian[0][2] = (- dxy * np.cos(phi) + (np.cos(lmbd) * (np.cos(phi) - np.cos(phi - tmin))) / (B * qoverp))[0]
        HelixJacobian[1][2] = (- dxy * np.sin(phi) + (np.cos(lmbd) * (np.sin(phi) - np.sin(phi - tmin))) / (B * qoverp))[0]
        HelixJacobian[2][2] = 0.
        # dxy
        HelixJacobian[0][3] = - np.sin(phi)
        HelixJacobian[1][3] = np.cos(phi)
        HelixJacobian[2][3] = 0.
        # dsz
        HelixJacobian[0][4] = 0.
        HelixJacobian[1][4] = 0.
        HelixJacobian[2][4] = 1. / np.cos(lmbd)

        return HelixJacobian

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

    def getIPsignificanceXY(self):
        if self.IPerrorXY == 0: return 0.
        return self.getDxy() / self.IPerrorXY

    def getIPsignificanceZ(self):
        if self.IPerrorZ == 0: return 0.
        return self.getDz() / self.IPerrorZ


###############################################################################################
def computeFastSimBugWeight(genjets, genparticles):
    """event weight correcting FastSim bug (from Sam)
    """

    # weightfile = ROOT.TFile('/data/dust/user/beinsam/pMSSM13TeV/Scan2/FixWeights/rootfiles/fastsim_decay_bug_weights.root')
    weightfile = ROOT.TFile('/data/dust/user/wolfmor/NTupleStuff/fastsim_decay_bug_weights.root')
    weighthist_map = {}
    for iflav in range(1,4):
        weighthist_map[iflav] = weightfile.Get('hRatio_GenJetHadronPtGenJetHadronFlavorLt4')
    weighthist_map[4] = weightfile.Get('hRatio_GenJetHadronPtGenJetHadronFlavorEqEq4')
    weighthist_map[5] = weightfile.Get('hRatio_GenJetHadronPtGenJetHadronFlavorEqEq5')
    weightHistAxis = weighthist_map[5].GetXaxis()
    rangesOfHadFlav = {}
    rangesOfHadFlav[1] = (100,400)
    rangesOfHadFlav[2] = (100,400)
    rangesOfHadFlav[3] = (100,400)
    rangesOfHadFlav[4] = (400,500)
    rangesOfHadFlav[5] = (400,600)

    event_weight = 1.0
    for ig, gjet in enumerate(genjets):
        gjettlv = ROOT.TLorentzVector(gjet.px(),gjet.py(),gjet.pz(),gjet.energy())
        if not gjet.pt()>30: continue
        jetHasOffendingGp = False

        #this is an alternative to miniAOD version which takes flavor from reco jet
        leadingHadronFlavor = 1
        leadHadronPt = 0

        for igp, gp in enumerate(genparticles):
            if not gp.pt()>15: continue
            gptlv = ROOT.TLorentzVector(gp.px(),gp.py(),gp.pz(),gp.energy())
            if not gptlv.DeltaR(gjettlv)<0.4: continue
            pdgid = abs(gp.pdgId())
            if not (pdgid>100 and pdgid<600): continue
            if not (gp.pt()>leadHadronPt): continue
            leadHadronPt = gp.pt()
            leadingHadronFlavor = int(pdgid)/100
        #this is an alternative to miniAOD version which takes flavor from reco jet

        hadronPt = 1
        genJetHadronFlavor = leadingHadronFlavor

        for igp, gp in enumerate(genparticles):
            if not gp.pt()>15: continue
            gptlv = ROOT.TLorentzVector(gp.px(),gp.py(),gp.pz(),gp.energy())
            if not gptlv.DeltaR(gjettlv)<0.4: continue
            pdg = abs(gp.pdgId())
            if not (abs(gp.pdgId())>100 and abs(gp.pdgId())<600): continue
            if pdg>rangesOfHadFlav[genJetHadronFlavor][0] and pdg<rangesOfHadFlav[genJetHadronFlavor][1]:
                if gp.pt()>hadronPt:
                    hadronPt = gp.pt()
            if not gp.numberOfDaughters()>1: continue
            originxy = ROOT.TMath.Sqrt(gp.vx()**2 + gp.vy()**2)
            if not originxy<2.16: continue
            if gp.daughter(0) == None: continue
            decayxy = ROOT.TMath.Sqrt(gp.daughter(0).vx()**2 + gp.daughter(0).vy()**2)
            if (decayxy>2.17 and decayxy<2000):
                jetHasOffendingGp = True
                break

        if jetHasOffendingGp:
            event_weight*=0
        else:
            thept = weightHistAxis.GetBinLowEdge(weightHistAxis.FindBin(min(hadronPt, 498)))
            genjetweight = weighthist_map[genJetHadronFlavor].Interpolate(thept)
            event_weight*=genjetweight
    return event_weight


###############################################################################################
def getCTauCmChi1pmHiggsinoFromMassGeV(dmchipmchi10, mass):

    # https://arxiv.org/pdf/2312.08087
    width_higgsino_gev_2025 = {"mChipm100GeV_dm0p159GeV": 6.8782e-16, "mChipm115GeV_dm0p168GeV": 9.4264e-16,
                               "mChipm140GeV_dm0p179GeV": 1.28457e-15, "mChipm160GeV_dm0p187GeV": 1.55914e-15,
                               "mChipm180GeV_dm0p193GeV": 1.78117e-15, "mChipm200GeV_dm0p198GeV": 1.9772e-15,
                               "mChipm250GeV_dm0p208GeV": 2.4023e-15, "mChipm300GeV_dm0p215GeV": 2.7274e-15,
                               "mChipm500GeV_dm0p23GeV": 3.5047e-15, "mChipm100GeV_dm0p259GeV": 5.35684e-15,
                               "mChipm115GeV_dm0p268GeV": 6.0366e-15, "mChipm140GeV_dm0p279GeV": 6.9396e-15,
                               "mChipm160GeV_dm0p287GeV": 7.6501e-15, "mChipm180GeV_dm0p293GeV": 8.21257e-15,
                               "mChipm200GeV_dm0p298GeV": 8.70202e-15, "mChipm250GeV_dm0p308GeV": 9.73785e-15,
                               "mChipm300GeV_dm0p315GeV": 1.0508e-14, "mChipm500GeV_dm0p33GeV": 1.2302e-14,
                               "mChipm100GeV_dm0p359GeV": 1.63174e-14, "mChipm115GeV_dm0p368GeV": 1.774198e-14,
                               "mChipm140GeV_dm0p379GeV": 1.960104e-14, "mChipm160GeV_dm0p387GeV": 2.10386e-14,
                               "mChipm180GeV_dm0p393GeV": 2.216256e-14, "mChipm200GeV_dm0p398GeV": 2.31328e-14,
                               "mChipm250GeV_dm0p408GeV": 2.51601e-14, "mChipm300GeV_dm0p415GeV": 2.6652e-14,
                               "mChipm500GeV_dm0p43GeV": 3.0067e-14, "mChipm100GeV_dm0p459GeV": 3.74352e-14,
                               "mChipm115GeV_dm0p468GeV": 4.002222e-14, "mChipm140GeV_dm0p479GeV": 4.337052e-14,
                               "mChipm160GeV_dm0p487GeV": 4.593348e-14, "mChipm180GeV_dm0p493GeV": 4.806868e-14,
                               "mChipm200GeV_dm0p498GeV": 5.0053e-14, "mChipm250GeV_dm0p508GeV": 5.40178e-14,
                               "mChipm300GeV_dm0p515GeV": 5.6835e-14, "mChipm500GeV_dm0p53GeV": 6.3145e-14,
                               "mChipm100GeV_dm0p559GeV": 7.62704e-14, "mChipm115GeV_dm0p568GeV": 8.085926e-14,
                               "mChipm140GeV_dm0p579GeV": 8.676776e-14, "mChipm160GeV_dm0p587GeV": 9.126204e-14,
                               "mChipm180GeV_dm0p593GeV": 9.47404e-14, "mChipm200GeV_dm0p598GeV": 9.7725e-14,
                               "mChipm250GeV_dm0p608GeV": 1.038959e-13, "mChipm300GeV_dm0p615GeV": 1.0842e-13,
                               "mChipm500GeV_dm0p63GeV": 1.1853e-13, "mChipm100GeV_dm0p759GeV": 2.43352e-13,
                               "mChipm115GeV_dm0p768GeV": 2.560744e-13, "mChipm140GeV_dm0p779GeV": 2.726332e-13,
                               "mChipm160GeV_dm0p787GeV": 2.853922e-13, "mChipm180GeV_dm0p793GeV": 2.952944e-13,
                               "mChipm200GeV_dm0p798GeV": 3.0382e-13, "mChipm250GeV_dm0p808GeV": 3.21528e-13,
                               "mChipm300GeV_dm0p815GeV": 3.34655e-13, "mChipm500GeV_dm0p83GeV": 3.6414e-13,
                               "mChipm100GeV_dm0p959GeV": 7.17521e-13, "mChipm115GeV_dm0p968GeV": 7.514406e-13,
                               "mChipm140GeV_dm0p979GeV": 7.949828e-13, "mChipm160GeV_dm0p987GeV": 8.278256e-13,
                               "mChipm180GeV_dm0p993GeV": 8.529928e-13, "mChipm200GeV_dm0p998GeV": 8.74514e-13,
                               "mChipm250GeV_dm1p008GeV": 9.18359e-13, "mChipm300GeV_dm1p015GeV": 9.50445e-13,
                               "mChipm500GeV_dm1p03GeV": 1.0212e-12, "mChipm100GeV_dm1p259GeV": 2.6587e-12,
                               "mChipm115GeV_dm1p268GeV": 2.756878e-12, "mChipm140GeV_dm1p279GeV": 2.882924e-12,
                               "mChipm160GeV_dm1p287GeV": 2.976826e-12, "mChipm180GeV_dm1p293GeV": 3.048252e-12,
                               "mChipm200GeV_dm1p298GeV": 3.10932e-12, "mChipm250GeV_dm1p308GeV": 3.23196e-12,
                               "mChipm300GeV_dm1p315GeV": 3.3217e-12, "mChipm500GeV_dm1p33GeV": 3.5173e-12,
                               "mChipm100GeV_dm1p759GeV": 1.139465e-11, "mChipm115GeV_dm1p768GeV": 1.1599754e-11,
                               "mChipm140GeV_dm1p779GeV": 1.1863202e-11, "mChipm160GeV_dm1p787GeV": 1.2050756e-11,
                               "mChipm180GeV_dm1p793GeV": 1.2190312e-11, "mChipm200GeV_dm1p798GeV": 1.23101e-11,
                               "mChipm250GeV_dm1p808GeV": 1.253442e-11, "mChipm300GeV_dm1p815GeV": 1.26987e-11,
                               "mChipm500GeV_dm1p83GeV": 1.30381e-11, "mChipm100GeV_dm2p259GeV": 2.136965e-11,
                               "mChipm115GeV_dm2p268GeV": 2.1598754e-11, "mChipm140GeV_dm2p279GeV": 2.1902202e-11,
                               "mChipm160GeV_dm2p287GeV": 2.2114756e-11, "mChipm180GeV_dm2p293GeV": 2.2272312e-11,
                               "mChipm200GeV_dm2p298GeV": 2.24101e-11, "mChipm250GeV_dm2p308GeV": 2.265442e-11,
                               "mChipm300GeV_dm2p315GeV": 2.28387e-11, "mChipm500GeV_dm2p33GeV": 2.32081e-11,
                               "mChipm100GeV_dm3p259GeV": 4.131965e-11, "mChipm115GeV_dm3p268GeV": 4.1596754e-11,
                               "mChipm140GeV_dm3p279GeV": 4.1980202e-11, "mChipm160GeV_dm3p287GeV": 4.2242756e-11,
                               "mChipm180GeV_dm3p293GeV": 4.2436312e-11, "mChipm200GeV_dm3p298GeV": 4.26101e-11,
                               "mChipm250GeV_dm3p308GeV": 4.289442e-11, "mChipm300GeV_dm3p315GeV": 4.31187e-11,
                               "mChipm500GeV_dm3p33GeV": 4.35481e-11, "mChipm100GeV_dm4p259GeV": 6.126965e-11,
                               "mChipm115GeV_dm4p268GeV": 6.1594754e-11, "mChipm140GeV_dm4p279GeV": 6.2058202e-11,
                               "mChipm160GeV_dm4p287GeV": 6.2370756e-11, "mChipm180GeV_dm4p293GeV": 6.2600312e-11,
                               "mChipm200GeV_dm4p298GeV": 6.28101e-11, "mChipm250GeV_dm4p308GeV": 6.313442e-11,
                               "mChipm300GeV_dm4p315GeV": 6.33987e-11, "mChipm500GeV_dm4p33GeV": 6.38881e-11,
                               "mChipm100GeV_dm5p259GeV": 8.121965e-11, "mChipm115GeV_dm5p268GeV": 8.1592754e-11,
                               "mChipm140GeV_dm5p279GeV": 8.2136202e-11, "mChipm160GeV_dm5p287GeV": 8.2498756e-11,
                               "mChipm180GeV_dm5p293GeV": 8.2764312e-11, "mChipm200GeV_dm5p298GeV": 8.30101e-11,
                               "mChipm250GeV_dm5p308GeV": 8.337442e-11, "mChipm300GeV_dm5p315GeV": 8.36787e-11,
                               "mChipm500GeV_dm5p33GeV": 8.42281e-11, "mChipm225GeV_dm0p203GeV": 2.18435e-15,
                               "mChipm225GeV_dm0p303GeV": 9.21012e-15, "mChipm225GeV_dm0p403GeV": 2.41312e-14,
                               "mChipm225GeV_dm0p503GeV": 5.204045e-14, "mChipm225GeV_dm0p603GeV": 1.0078758e-13,
                               "mChipm225GeV_dm0p803GeV": 3.1257325e-13, "mChipm225GeV_dm1p003GeV": 8.96264e-13,
                               "mChipm225GeV_dm1p303GeV": 3.17023e-12, "mChipm225GeV_dm1p803GeV": 1.242216e-11,
                               "mChipm225GeV_dm2p303GeV": 2.253216e-11, "mChipm225GeV_dm3p303GeV": 4.275216e-11,
                               "mChipm400GeV_dm0p224GeV": 3.1799e-15, "mChipm400GeV_dm0p324GeV": 1.15626e-14,
                               "mChipm400GeV_dm0p424GeV": 2.86678e-14, "mChipm177GeV_dm0p293GeV": 8.21246e-15,
                               "mChipm300GeV_dm0p165GeV": 8.5555e-16, "mChipm400GeV_dm0p174GeV": 1.1239e-15,
                               "mChipm500GeV_dm0p18GeV": 1.3169e-15, "mChipm700GeV_dm0p338GeV": 1.33414e-14,
                               "mChipm900GeV_dm0p342GeV": 1.38818e-14, "mChipm1100GeV_dm0p345GeV": 1.4294e-14,
                               "mChipm700GeV_dm0p238GeV": 3.96778e-15, "mChipm900GeV_dm0p242GeV": 4.21258e-15,
                               "mChipm1100GeV_dm0p245GeV": 4.4025e-15, "mChipm700GeV_dm0p138GeV": 7.6e-18,
                               "mChipm900GeV_dm0p142GeV": 1.9248e-16, "mChipm1100GeV_dm0p145GeV": 2.9766e-16,
                               "mChipm100GeV_dm1p459GeV": 5.44215e-12, "mChipm140GeV_dm1p479GeV": 5.845922e-12,
                               "mChipm180GeV_dm1p493GeV": 6.141112e-12, "mChipm250GeV_dm1p508GeV": 6.46242e-12,
                               "mChipm200GeV_dm1p498GeV": 6.2501e-12, "mChipm500GeV_dm1p53GeV": 6.9361e-12,
                               "mChipm115GeV_dm1p468GeV": 5.6187e-12, "mChipm160GeV_dm1p487GeV": 6.013892e-12,
                               "mChipm300GeV_dm1p515GeV": 6.6147e-12}
    # c in cm / s
    # hbar in GeV s
    ctau = 2.998e10 * 6.582119569509067e-25 / width_higgsino_gev_2025['mChipm' + str(int(mass)) + 'GeV_dm' + str(dmchipmchi10).replace('.', 'p') + 'GeV']

    # outdated:
    # https://arxiv.org/pdf/1410.4549.pdf
    # mpi = 0.135
    # ctau = 1.1*pow(dmchipmchi10*1000./300,-3)*pow(1-pow(mpi/dmchipmchi10,2),-0.5)

    return ctau


def getCTauCmChi1pmWinoFromMassGeV(dmchipmchi10, mass, localpath):
    #https://arxiv.org/pdf/1410.4549.pdf
    fwino = ROOT.TFile(localpath + 'WinoCtauVsDmVsM.root')
    hctau = fwino.Get('hwinoctaucm')
    #winoctau = hctau.GetBinContent(hctau.GetXaxis().FindBin(mass), hctau.GetYaxis().FindBin(dmchipmchi10))
    winoctau = hctau.Interpolate(mass, dmchipmchi10)
    fwino.Close()
    return winoctau


#re-weight ctau from Viktor
def reweight_ctau(ctauIn, ctauOut, list_of_labXY):
    output = 1
    for labXY in list_of_labXY:
        t0 = labXY / 10.0          # convert to cm
        if ctauIn>0 and ctauOut>0:
            output *= ctauIn/ctauOut * ROOT.TMath.Exp(t0/ctauIn - t0/ctauOut)
    return output
