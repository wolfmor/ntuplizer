#! /usr/bin/env python


'''

Runs over AOD files and writes file with histos and tree.

---------------------------------------------------------
python plantTrees.py inputFiles="..." tag="tag1 tag2 ..."
---------------------------------------------------------

tags:
-----

local -> don't open file via xrootd
redirinfn -> use infn redirector instead of global cern redirector
redirfnal -> use fnal redirector instead of global cern redirector
data -> don't load GEN collections, check trigger flags and runnum/lumisec/eventnum
geninfoZ -> check GEN Z boson decay
geninfoW -> check GEN W boson decay
cleanleptons -> perform DY cleaning
signal -> signal GEN variables and FastSim correction for MET
noleptonveto -> don't veto events with leptons
genmatchtracks -> find GEN match to every 10. track
genmatchalltracks -> find GEN match to every track

'''


import os
import sys
import json
import re
import random
from array import array

import ROOT

ROOT.gROOT.SetBatch()        # don't pop up canvases
ROOT.gROOT.SetStyle('Plain') # white background

ROOT.gSystem.Load('libFWCoreFWLite.so')
ROOT.gSystem.Load('libDataFormatsFWLite.so')
ROOT.FWLiteEnabler.enable()

from DataFormats.FWLite import Events, Handle
from FWCore.ParameterSet.VarParsing import VarParsing


from commons import *


'''DY cleaning: replace leptons with "neutrinos" and update collections.
'''
def cleanZllEvent(zl1Idx, zl2Idx, collection, tracks, pfcands, jets, met, hZllLeptonPt, hZllDrTrack, hZllDrPfc, hZllDrJet):

	# clean jets, tracks, pfcands for leptons and adapt MET
	
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


'''Jet energy corrections.
'''
def createJEC(jecSrc, jecLevelList, jetAlgo):
	
	jecParameterList = ROOT.vector('JetCorrectorParameters')()
	
	# Load the different JEC levels (the order matters!)
	for jecLevel in jecLevelList:
		jecParameter = ROOT.JetCorrectorParameters('%s_%s_%s.txt' % (jecSrc, jecLevel, jetAlgo));
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
		for minrun,maxrun,version in inputmap:
			JECMap = {}
			JECMap['jecAK4'] = createJEC('/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_1_7/src/SoftLeptonBDT/JECs/'+version+'/'+version, ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
			self.JECList.append([minrun, maxrun, JECMap])
	
	def GetJECMap(self, run):
		for minrun,maxrun,returnmap in self.JECList:
			if run >= minrun and run <= maxrun:
				return returnmap
		raise Exception('Error! Run '+str(run)+' not found in run ranges')
	
	def jecAK4(self, run):
		JECMap = self.GetJECMap(run)
		return JECMap['jecAK4']


if True:
	
	# https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
	# from https://twiki.cern.ch/twiki/bin/viewauth/CMS/PdmV2016Analysis#Re_reco_datasets_07Aug17
	goldenjson16 = {'273158': [[1, 1279]], '273302': [[1, 459]], '273402': [[100, 292]], '273403': [[1, 53]], '273404': [[1, 18]], '273405': [[2, 25]], '273406': [[1, 112]], '273408': [[1, 6]], '273409': [[1, 309]], '273410': [[1, 90]], '273411': [[1, 29]], '273425': [[62, 352], [354, 733]], '273446': [[1, 33]], '273447': [[1, 113], [115, 412]], '273448': [[1, 391]], '273449': [[1, 214]], '273450': [[1, 214], [219, 647]], '273492': [[71, 71], [73, 282], [284, 325], [327, 338]], '273493': [[1, 233]], '273494': [[1, 192]], '273502': [[73, 256], [258, 318], [320, 813], [815, 1064]], '273503': [[1, 598]], '273554': [[77, 437]], '273555': [[1, 173]], '273725': [[83, 252], [254, 2545]], '273728': [[1, 100]], '273730': [[1, 1814], [1820, 2126]], '274094': [[108, 332]], '274146': [[1, 67]], '274157': [[105, 534]], '274159': [[1, 43]], '274160': [[1, 207]], '274161': [[1, 516]], '274172': [[31, 95]], '274198': [[81, 191]], '274199': [[1, 623]], '274200': [[1, 678]], '274240': [[1, 40], [42, 82]], '274241': [[1, 1152], [1161, 1176]], '274244': [[1, 607]], '274250': [[1, 701]], '274251': [[1, 546]], '274283': [[2, 19]], '274284': [[1, 210]], '274286': [[1, 154]], '274314': [[97, 97], [99, 158]], '274315': [[1, 424]], '274316': [[1, 959]], '274317': [[1, 3]], '274319': [[1, 225]], '274335': [[60, 1003]], '274336': [[1, 14]], '274337': [[3, 17]], '274338': [[1, 698]], '274339': [[1, 29], [31, 31], [33, 33], [35, 93]], '274344': [[1, 632]], '274345': [[1, 170]], '274382': [[94, 144]], '274387': [[88, 439]], '274388': [[1, 1820]], '274420': [[94, 268]], '274421': [[1, 342]], '274422': [[1, 2207]], '274440': [[92, 493]], '274441': [[1, 431]], '274442': [[1, 752]], '274954': [[37, 37], [39, 57]], '274955': [[1, 91]], '274968': [[1, 1192]], '274969': [[1, 1003]], '274970': [[1, 47]], '274971': [[1, 905]], '274998': [[64, 782]], '274999': [[1, 1241]], '275000': [[1, 136]], '275001': [[1, 1781], [1786, 2061]], '275059': [[78, 81], [105, 137]], '275066': [[1, 96]], '275067': [[1, 392]], '275068': [[1, 915]], '275073': [[1, 517]], '275074': [[1, 442], [444, 647]], '275124': [[106, 106], [108, 431]], '275125': [[1, 989]], '275282': [[91, 180]], '275283': [[1, 132]], '275284': [[1, 74]], '275290': [[96, 143]], '275291': [[1, 347]], '275292': [[1, 121]], '275293': [[1, 142], [144, 201]], '275309': [[55, 617]], '275310': [[1, 1929]], '275311': [[1, 1253]], '275319': [[141, 282]], '275337': [[1, 427]], '275338': [[1, 520]], '275344': [[76, 356]], '275345': [[1, 353]], '275370': [[81, 365]], '275371': [[1, 22], [28, 569]], '275375': [[127, 1449]], '275376': [[1, 2667], [2669, 3096]], '275657': [[1, 105]], '275658': [[1, 337]], '275659': [[1, 17]], '275761': [[1, 9]], '275767': [[1, 4]], '275772': [[1, 56]], '275773': [[1, 7]], '275774': [[1, 311], [315, 315]], '275776': [[1, 140]], '275777': [[1, 300]], '275778': [[1, 305]], '275782': [[1, 131], [133, 762]], '275832': [[1, 367]], '275833': [[1, 53], [56, 115], [117, 251]], '275834': [[1, 297]], '275835': [[1, 13]], '275836': [[1, 429], [431, 1163], [1166, 1170], [1184, 1293]], '275837': [[1, 186], [198, 726]], '275847': [[1, 2263]], '275886': [[73, 109]], '275890': [[1, 1393]], '275911': [[62, 298], [300, 354], [356, 440]], '275912': [[1, 289]], '275913': [[1, 475]], '275918': [[1, 318], [348, 361]], '275920': [[5, 463]], '275921': [[1, 2], [4, 5], [17, 20]], '275923': [[3, 53], [63, 64], [66, 126]], '275931': [[1, 14], [19, 89]], '275963': [[82, 139], [141, 172]], '276092': [[74, 149]], '276097': [[1, 507]], '276242': [[1, 7], [18, 61], [72, 1664]], '276243': [[1, 15], [18, 480], [482, 611]], '276244': [[3, 1202]], '276282': [[75, 534], [537, 1142]], '276283': [[3, 1087]], '276315': [[40, 175], [178, 217]], '276317': [[3, 138]], '276318': [[3, 103], [106, 570]], '276355': [[1, 33]], '276361': [[1, 161], [169, 208], [210, 800], [802, 833]], '276363': [[1, 140], [142, 238], [242, 1482]], '276384': [[2, 1117]], '276437': [[63, 224], [227, 1074], [1076, 2190]], '276454': [[1, 527]], '276458': [[1, 341]], '276495': [[87, 268]], '276501': [[4, 221], [223, 2547]], '276502': [[2, 741]], '276525': [[88, 469], [471, 1606], [1626, 2893]], '276527': [[1, 214]], '276528': [[4, 394]], '276542': [[74, 857]], '276543': [[1, 638], [643, 952]], '276544': [[2, 161]], '276545': [[2, 110], [117, 213]], '276581': [[79, 444]], '276582': [[1, 871]], '276583': [[1, 52]], '276584': [[1, 2]], '276585': [[1, 238], [241, 242], [245, 246]], '276586': [[2, 658], [680, 773]], '276587': [[1, 1006]], '276653': [[72, 550]], '276655': [[1, 593], [595, 1106]], '276659': [[1, 127], [129, 252]], '276775': [[96, 1260]], '276776': [[1, 1823]], '276794': [[1, 885]], '276807': [[66, 220]], '276808': [[1, 875]], '276810': [[1, 287]], '276811': [[1, 1270], [1272, 2563]], '276831': [[64, 755], [761, 2702]], '276834': [[1, 720]], '276870': [[78, 1354], [1356, 3108], [3111, 3258], [3260, 3484]], '276935': [[79, 184], [186, 838], [842, 906]], '276940': [[70, 213]], '276946': [[1, 27]], '276947': [[1, 89], [91, 126], [135, 141]], '276948': [[1, 474]], '276950': [[1, 2353]], '277069': [[81, 265], [267, 390]], '277070': [[1, 309], [311, 1059]], '277071': [[1, 82], [90, 178]], '277072': [[1, 253], [256, 466]], '277073': [[1, 90]], '277076': [[1, 3], [5, 7], [9, 35], [38, 1037]], '277087': [[204, 1191]], '277094': [[1, 161], [164, 584]], '277096': [[1, 1309], [1311, 2086]], '277112': [[1, 155]], '277126': [[42, 59]], '277127': [[1, 438], [440, 902]], '277148': [[83, 190], [193, 700]], '277166': [[77, 186], [188, 431]], '277168': [[1, 1708], [1711, 1822], [1824, 2223]], '277180': [[88, 228]], '277194': [[113, 139], [144, 497], [500, 1115], [1117, 1312], [1320, 1749], [1754, 2067], [2070, 2070]], '277305': [[62, 744]], '277420': [[84, 84], [86, 291], [293, 346]], '277981': [[82, 83], [85, 163]], '277991': [[1, 98]], '277992': [[1, 260], [262, 312]], '278017': [[77, 97], [99, 213], [215, 512], [514, 589]], '278018': [[1, 263], [265, 422], [424, 615], [617, 627], [642, 1011], [1020, 1181]], '278167': [[87, 394], [397, 1153], [1155, 1660], [1662, 1707], [1709, 2258]], '278175': [[1, 88]], '278193': [[77, 231]], '278239': [[76, 339], [341, 558], [560, 740]], '278240': [[1, 64], [70, 113], [115, 1121], [1123, 1296], [1299, 1309]], '278273': [[75, 110]], '278274': [[1, 18], [20, 85]], '278288': [[67, 81]], '278289': [[1, 42], [44, 52]], '278290': [[1, 11]], '278308': [[87, 216], [219, 587], [589, 680], [683, 1200], [1217, 1410], [1413, 1848], [1880, 1880]], '278310': [[1, 32], [34, 709]], '278315': [[73, 254], [256, 661], [663, 767]], '278345': [[84, 500], [503, 831]], '278346': [[1, 117]], '278349': [[1, 401], [403, 612], [632, 633]], '278366': [[1, 453]], '278406': [[85, 360], [362, 1682]], '278509': [[91, 1557]], '278769': [[75, 104]], '278770': [[1, 767]], '278801': [[48, 85]], '278802': [[1, 17]], '278803': [[1, 87], [91, 133], [135, 297], [299, 323]], '278804': [[1, 4]], '278805': [[3, 26], [30, 167], [170, 193], [196, 280], [283, 284], [288, 288]], '278808': [[1, 445], [447, 462], [464, 1793]], '278820': [[17, 1533]], '278822': [[1, 1627]], '278873': [[70, 129]], '278874': [[1, 273], [275, 478]], '278875': [[1, 210], [212, 834]], '278923': [[55, 467]], '278957': [[79, 227]], '278962': [[68, 408]], '278963': [[1, 23], [25, 175]], '278969': [[70, 511], [514, 1051], [1053, 1291], [1293, 1397], [1399, 1460]], '278975': [[1, 475], [477, 745], [747, 850]], '278976': [[1, 20]], '278986': [[71, 199]], '279024': [[82, 382]], '279029': [[1, 260]], '279071': [[71, 244]], '279080': [[68, 224]], '279115': [[118, 524]], '279116': [[38, 485]], '279479': [[86, 190]], '279588': [[100, 1259]], '279653': [[77, 77], [82, 261]], '279654': [[1, 108], [110, 1231], [1285, 1299]], '279656': [[1, 43]], '279658': [[1, 689], [691, 713]], '279667': [[68, 1033]], '279681': [[77, 104]], '279682': [[1, 29], [33, 34], [37, 38]], '279683': [[1, 26]], '279684': [[1, 22]], '279685': [[1, 93], [95, 209]], '279691': [[71, 113]], '279694': [[1, 2235]], '279715': [[71, 474], [476, 477], [480, 480], [511, 511], [523, 691]], '279716': [[1, 860], [875, 1528], [1530, 1653]], '279760': [[68, 578], [585, 728]], '279766': [[1, 1689]], '279767': [[1, 776]], '279794': [[77, 1100]], '279823': [[61, 395]], '279841': [[75, 398], [407, 2122]], '279844': [[72, 295]], '279887': [[79, 221], [225, 397]], '279931': [[84, 628], [630, 743], [746, 801], [803, 1043], [1045, 3022]], '279966': [[79, 441]], '279975': [[70, 190], [192, 253], [256, 281], [283, 709], [734, 1121]], '279993': [[85, 156]], '279994': [[1, 47]], '280013': [[1, 25]], '280015': [[1, 39], [41, 56], [59, 554], [560, 580]], '280016': [[1, 149]], '280017': [[1, 608]], '280018': [[1, 1281]], '280020': [[1, 45]], '280024': [[1, 427]], '280187': [[4, 60]], '280188': [[1, 245]], '280191': [[1, 781], [783, 866], [869, 900]], '280194': [[1, 238]], '280242': [[1, 411], [414, 627]], '280249': [[1, 486], [488, 1433]], '280251': [[1, 165], [167, 372]], '280327': [[49, 85]], '280330': [[1, 857]], '280349': [[1, 247], [252, 623], [626, 626]], '280363': [[1, 359]], '280364': [[1, 370], [372, 617], [619, 619], [621, 1090], [1102, 1363]], '280383': [[64, 65]], '280384': [[2, 34]], '280385': [[1, 519], [523, 569], [574, 1187], [1189, 1533], [1536, 2022]], '281613': [[101, 128], [130, 130], [133, 133], [135, 139], [143, 256], [258, 903]], '281639': [[1, 132]], '281641': [[1, 319]], '281693': [[1, 2191]], '281707': [[99, 982], [1000, 1065]], '281726': [[1, 288]], '281727': [[1, 1605]], '281797': [[125, 2176]], '281975': [[1, 215]], '281976': [[1, 2166]], '282033': [[82, 117]], '282034': [[1, 33]], '282035': [[1, 40]], '282037': [[1, 457], [459, 1862]], '282092': [[92, 222], [624, 2276]], '282708': [[1, 8]], '282710': [[1, 2], [8, 8]], '282712': [[1, 1], [10, 68]], '282730': [[89, 164]], '282731': [[1, 172]], '282732': [[1, 69]], '282733': [[1, 177]], '282734': [[1, 327]], '282735': [[1, 642], [645, 1232], [1235, 1823]], '282800': [[1, 377]], '282807': [[1, 326]], '282814': [[1, 1843]], '282842': [[1, 80]], '282917': [[117, 157], [159, 191]], '282918': [[1, 51]], '282919': [[1, 243]], '282922': [[1, 131]], '282923': [[1, 17], [19, 30], [32, 36], [38, 39], [41, 86], [88, 224]], '283042': [[1, 6]], '283043': [[1, 105], [108, 519]], '283049': [[82, 93]], '283050': [[1, 212]], '283052': [[1, 111]], '283059': [[1, 125], [127, 451]], '283270': [[76, 573], [576, 1502], [1504, 1888], [1890, 1912]], '283283': [[4, 1668], [1670, 1748]], '283305': [[79, 85]], '283306': [[1, 289]], '283307': [[1, 153], [156, 456]], '283308': [[1, 547], [549, 571], [573, 895], [897, 948]], '283353': [[80, 822]], '283358': [[1, 243], [245, 981]], '283359': [[1, 428]], '283407': [[82, 114]], '283408': [[1, 27], [29, 2088], [2098, 2125], [2203, 2416], [2528, 2542]], '283416': [[49, 151], [154, 245]], '283453': [[83, 537]], '283469': [[74, 74]], '283478': [[76, 303], [324, 969]], '283548': [[145, 288]], '283680': [[1, 81]], '283681': [[1, 17]], '283682': [[1, 384]], '283685': [[1, 314]], '283820': [[67, 1548]], '283830': [[1, 722]], '283834': [[1, 67], [69, 82]], '283835': [[1, 14], [16, 112]], '283865': [[1, 1177]], '283876': [[65, 211], [215, 724]], '283877': [[1, 1496]], '283884': [[349, 504], [509, 756]], '283885': [[1, 1723]], '283933': [[88, 232]], '283934': [[1, 784], [793, 870], [875, 1245], [1267, 1291]], '283946': [[85, 1448], [1450, 1462]], '283964': [[1, 388]], '284006': [[73, 390]], '284014': [[1, 266]], '284025': [[110, 157]], '284029': [[1, 112]], '284035': [[1, 360]], '284036': [[1, 140], [143, 348]], '284037': [[1, 340]], '284038': [[1, 55]], '284039': [[1, 30]], '284040': [[1, 33]], '284041': [[1, 44]], '284042': [[1, 129]], '284043': [[1, 205], [210, 224]], '284044': [[1, 30]]}

	# from https://twiki.cern.ch/twiki/bin/view/CMS/JECDataMC
	jet_energy_corrections = [
	    [1, 276811, 'Summer16_07Aug2017BCD_V11_DATA'],
	    [276831, 278801, 'Summer16_07Aug2017EF_V11_DATA'],
	    [278802, float('inf'), 'Summer16_07Aug2017GH_V11_DATA']]
	
	
nMaxEventsPerFile = 100000
nMaxTracksPerEvent = 10000

#TODO: check if test
saveOutputFile = True
isTest = True
neventsTest = 100  # number of events to run over in case of test
printevery = 100  # print event number for every Xth event

#TODO: check thresholds for "new" matching
matchingDrThreshold = 0.05
matchingDxyzThreshold = 0.2

#TODO: check met threshold for event selection
metthreshold = 200


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
		
		,('pionptGEN', 'F'), ('pionetaGEN', 'F'), ('pionphiGEN', 'F')
		,('tminmatching', 'F'), ('dxyzmin', 'F'), ('drmin', 'F')
		,('dxyzminrandom', 'F'), ('drminrandom', 'F')
		,('drminold', 'F'), ('drminoldrandom', 'F')
		
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
		,('invmCleaning', 'F')
		,('l1ptCleaning', 'F'), ('l2ptCleaning', 'F')
		,('l1etaCleaning', 'F'), ('l2etaCleaning', 'F')
		,('l1phiCleaning', 'F'), ('l2phiCleaning', 'F')
		,('l1absisodbeta', 'F'), ('l2absisodbeta', 'F')
		,('l1relisodbeta', 'F'), ('l2relisodbeta', 'F')
		]
	event_level_var_names += var_names_cleaning
	
	var_names_event = [
		('cutflow', 'I'), ('randomevent', 'I')
		,('numpvs', 'I'), ('rho', 'F')
		,('metpt', 'F'), ('metphi', 'F')
		,('nofastsimcorrmetpt', 'F'), ('nofastsimcorrmetphi', 'F')
		,('ht', 'F'), ('htmiss', 'F')
		,('numbadjets', 'F'), ('minetaabsbadjets', 'F')
		,('numjets', 'I'), ('numjets30', 'I'), ('numjets50','I'), ('numjets100', 'I'), ('numjets200', 'I')
		,('njetsbtagmedium', 'I'), ('njetsbtagmediumTIGHT', 'I')
		,('mtmetleadingjet', 'F')
		,('mindphimetjets', 'F')
		,('numphotons', 'I'), ('numphotonsiso', 'I')
		,('numpfleptons', 'I'), ('numpfleptonsiso', 'I')
		,('numelectrons', 'I'), ('numelectronsiso', 'I')
		,('nummuons', 'I'), ('nummuonsiso', 'I')
		,('numtaus', 'I'), ('numtausiso', 'I')
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
		
		
	trigger_flags = ['globalSuperTightHalo2016Filter'
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

	
	var_names_jet = [
		('pxjet', 'F'), ('pyjet', 'F'), ('pzjet', 'F')
		,('ptjet', 'F'), ('energyjet', 'F')
		,('etajet', 'F'), ('phijet', 'F')
		,('nconstituentsjet', 'I')
		,('btagjet', 'F')
		]
	
	jet_var_array = {}
	for n in var_names_jet:
		jet_var_array[n[0]] = array('f', 100*[0.])
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
		,('decaymodetau', 'F'), ('mvadiscrtau', 'F')
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
		('randomtrack','I')
		
		,('charge','F')
		,('pxtrack','F'), ('pytrack','F'), ('pztrack','F')
		,('pttrack','F')
        ,('pttrackerror/pttrack','F'), ('log10(pttrackerror/pttrack)','F')
        ,('eta','F'), ('phi','F')
        
		,('dxy','F'), ('dxyhandmade','F'), ('dxyclosestpv','F'), ('dxyclosestpvPU','F')
		,('dz','F'), ('dzhandmade','F'), ('dzclosestpv','F'), ('dzclosestpvPU','F')
		,('dxyerror','F'), ('dzerror','F')
        ,('log10(dxy)','F'), ('log10(dxyhandmade)','F'), ('log10(dxyclosestpv)','F'), ('log10(dxyclosestpvPU)','F')
        ,('log10(dz)','F'), ('log10(dzhandmade)','F'), ('log10(dzclosestpv)','F'), ('log10(dzclosestpvPU)','F')
        ,('log10(dxyerror)','F'), ('log10(dzerror)','F')
        
        ,('trackabsiso','F'), ('trackreliso','F'), ('trackdrmin','F'), ('tracknumneighbours','I')
        ,('pfabsiso','F'), ('pfreliso','F'), ('pfdrmin','F'), ('pfnumneighbours','I')
        ,('chpfabsiso','F'), ('chpfreliso','F'), ('chpfdrmin','F'), ('chpfnumneighbours','I')
        ,('jetiso','F'), ('jetisomulti','F'), ('jetdrmin','F') 
        ,('jetisotight','F'), ('jetisomultitight','F'), ('jetdrmintight','F')
        
        ,('drminphoton', 'F'), ('drminelectron', 'F'), ('drminmuon', 'F')
        ,('drmintau', 'F'), ('closesttaumvadiscr', 'F')
        
        ,('detahighestptjet','F'), ('dphihighestptjet','F')
        ,('dphimet','F'), ('dphimetpca','F')
        
        ,('chi2','F')
        ,('quality','I')
        ,('nvalidhits','I'), ('nlosthits','I')
        
        ,('issignaltrack','I'), ('issusytrack','I'), ('susytrackmother','I'), ('susytrackpdgid','I')
        
        ,('hasGenMatch','I')
        ,('genmatchpdgid','F'), ('genmatchstatus','F'), ('genmatchishardprocess','F'), ('genmatchisfromhardprocess','F')
        ,('genmatchisprompt','F'), ('genmatchisdirecthadrondecayproduct','F'), ('genmatchisdirecttaudecayproduct','F')
        ,('genmatchmotherpdgid','F'), ('genmatchmotherstatus','F'), ('genmatchmotherishardprocess','F')
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

	
'''
###############################################################################################
# get input arguments
###############################################################################################
'''
options = VarParsing('python')
options.parseArguments()

 
# in order to set grid environment on centos7 nodes
if not 'local' in options.tag:
# 	os.system('export XRD_LOGLEVEL=Debug')
# 	print 'xrd debug'
# 	os.system('printenv')
	os.system('xrdfs cms-xrd-global.cern.ch locate /store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6EA5C939-FAC9-E611-86E5-FA163E5F2D05.root')
	os.system('xrdcp -d 1 -f root://cms-xrd-global.cern.ch//store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6EA5C939-FAC9-E611-86E5-FA163E5F2D05.root /dev/null')
# 	os.system('xrd cms-xrd-global.cern.ch locateall /store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6EA5C939-FAC9-E611-86E5-FA163E5F2D05.root')
# 	os.system('export X509_USER_PROXY=/afs/desy.de/user/w/wolfmor/private/x509_voms')
# 	print ''
# 	os.system('source /afs/desy.de/user/w/wolfmor/.bashrc')
# 	print 'bashrc'
# 	print ''
# 	os.system('export XRD_LOGLEVEL=Debug')
# 	print 'xrd debug'
# 	os.system('xrd cms-xrd-global.cern.ch locateall /store/mc/RunIISummer16DR80Premix/ZJetsToNuNu_Zpt-200toInf_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/AODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6-v1/120000/6EA5C939-FAC9-E611-86E5-FA163E5F2D05.root')
# 	os.system('source /cvmfs/grid.cern.ch/etc/profile.d/setup-cvmfs-ui.sh')
# 	os.system('voms-proxy-info')
# if not 'local' in options.tag: os.system('source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh')


'''
###############################################################################################
# define handles and labels
###############################################################################################
'''
if 'data' in options.tag:
	
	handle_trigger = Handle('edm::TriggerResults')
	label_trigger = ('TriggerResults', '', 'RECO')
	
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

handle_taudiscriminatorMVA = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVA = ('hpsPFTauDiscriminationByVLooseIsolationMVArun2v1DBnewDMwLT')

handle_taudiscriminatorMVAraw = Handle('reco::PFTauDiscriminator')
label_taudiscriminatorMVAraw = ('hpsPFTauDiscriminationByIsolationMVArun2v1DBnewDMwLTraw')

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


jettype = 'AK4PFchs'
if 'data' in options.tag: # data
	DataJECs = DataJEC(jet_energy_corrections, jettype)
elif 'signal' in options.tag: # FastSim
	jecAK4 = createJEC('/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_1_7/src/SoftLeptonBDT/JECs/Summer16_FastSimV1_MC/Summer16_FastSimV1_MC', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)
else: # FullSim
	jecAK4 = createJEC('/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_1_7/src/SoftLeptonBDT/JECs/Summer16_07Aug2017_V11_MC/Summer16_07Aug2017_V11_MC', ['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual'], jettype)


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
		
		fin = ROOT.TFile.Open(f.strip())
		
		events = Events(fin)
	
		try:
			nevents = events.size()
		except:
			print 'skipping file ' + f
			continue
		
	else:
		
		redir = 'cms-xrd-global.cern.ch'
		if 'redirinfn' in options.tag: redir = 'xrootd-cms.infn.it'
		if 'redirfnal' in options.tag: redir = 'cmsxrootd.fnal.gov'
		
		fin = ROOT.TFile.Open('root://' + redir + '/' + f.strip())

		retry = 1
		while fin.IsZombie() or not fin.IsOpen():
		
			print 'retry open file ', retry
		
			retry += 1
			if retry > 5: break
		
			fin = ROOT.TFile.Open('root://' + redir + '/' + f.strip())
		
		events = Events(fin)
	
		nevents = events.size()
		
	print '### with ' + str(nevents) + ' events'
	print '### printing every ' + str(printevery) + '. event'

	if saveOutputFile: fout.cd()

	if isTest: nevents = neventsTest
	
	nEventsPerFile = 0

	phifirsttrack = 0
	etafirsttrack = 0
	
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
		
		hCutflow.Fill(0)
		cutflow = 0
		
		random.seed()
		event_level_var_array['randomevent'][0] = random.randrange(10)
		
		'''
		###############################################################################################
		# check runnum/lumisec and add to list for json
		###############################################################################################
		'''
		
		runnum = -1
		lumisec = -1
		eventnum = -1
		if 'data' in options.tag:
			
			runnum = event.eventAuxiliary().run()
			lumisec = event.eventAuxiliary().luminosityBlock()
			eventnum = event.eventAuxiliary().event()
			
			if runnum != lastrun or lumisec != lastlumi:
				if str(runnum) in goldenjson16: goodlumisecs = goldenjson16[str(runnum)]
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
		
		if not 'local' in options.tag:
			if fin.IsZombie() or not fin.IsOpen():
				print 'file not usable'
				sys.exit(1)
		
		if saveOutputFile and ievent % 100 == 0: fout.Write('', ROOT.TObject.kWriteDelete)

		if ievent % printevery == 0: print 'analyzing event %d of %d' % (ievent, min(nevents, nMaxEventsPerFile))

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
		# trigger flags
		###############################################################################################
		'''
		
		trigger_flags_accept = {}
		for tf in trigger_flags:
			trigger_flags_accept[tf] = -1
		
		allfine = True
		
		if 'data' in options.tag:
			
			event.getByLabel(label_trigger, handle_trigger)
			triggerresults = handle_trigger.product()
			
			triggernames = event.object().triggerNames(triggerresults)
			
			for i in range(triggerresults.size()):
				
				tn = triggernames.triggerName(i).replace('Flag_','')
				if tn in trigger_flags:
					if triggerresults.accept(i): trigger_flags_accept[tn] = 1
					else:
						trigger_flags_accept[tn] = 0
						allfine = False
		
		for tf in trigger_flags:
			event_level_var_array[tf][0] = trigger_flags_accept[tf]
			
		if not allfine: continue
		
		###########################################################################################veto
		
		hCutflow.Fill(3)
		cutflow = 3
			
		if not 'data' in options.tag:
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

		if not 'data' in options.tag:
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
		
		event.getByLabel(label_taudiscriminatorMVA, handle_taudiscriminatorMVA)
		taudiscriminatorMVA = handle_taudiscriminatorMVA.product()
		
		event.getByLabel(label_taudiscriminatorMVAraw, handle_taudiscriminatorMVAraw)
		taudiscriminatorMVAraw = handle_taudiscriminatorMVAraw.product()

		tauswithdiscriminators = [(tau, taudiscriminatorDM.value(itau), taudiscriminatorMVA.value(itau), taudiscriminatorMVAraw.value(itau)) for itau, tau in enumerate(taus)]

		
		'''
		###############################################################################################
		# lepton and photon selection
		###############################################################################################
		'''
		
		photons = [p for p in photons if passesPhotID(p, wp='loose')
										and p.pt() > 15 
										and abs(p.eta()) < 2.5]
		
		electrons = [e for e in electrons if passesEleID(e, pv_pos, wp='veto')
											and e.pt() > 10 
											and abs(e.eta()) < 2.5]
		
		muons = [m for m in muons if passesMuonID(m, wp='loose')
									and m.pt() > 10 
									and abs(m.eta()) < 2.4]
		
		tauswithdiscriminators = [t for t in tauswithdiscriminators if t[1] > 0.5
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
			
			Zgammas = [gp for gp in genparticles if gp.isLastCopy() and gp.statusFlags().fromHardProcess() and (abs(gp.pdgId()) == 22 or abs(gp.pdgId()) == 23)]
			
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
						
					if daughter == None: continue
					
					zdaughter_var_array['pdgIdZdaughter'][i] = daughter.pdgId()
					zdaughter_var_array['ptZdaughter'][i] = daughter.pt()
					zdaughter_var_array['etaZdaughter'][i] = daughter.eta()
					zdaughter_var_array['phiZdaughter'][i] = daughter.phi()
					
					if (abs(daughter.pdgId()) == 12 or abs(daughter.pdgId()) == 14 or abs(daughter.pdgId()) == 16):
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
		decayWtau = 0
		if 'geninfoW' in options.tag:
			
			Ws = [gp for gp in genparticles if gp.isLastCopy() and gp.statusFlags().fromHardProcess() and abs(gp.pdgId()) == 24]
			
			numW = len(Ws)
			
			if numW == 1:  # else there's something funny going on...
				
				W = Ws[0]
				
				ptW = W.pt()
				etaW = W.eta()
				phiW = W.phi()
				
				numWDaughters = W.numberOfDaughters()
				
				ptWneutrino = -1
				decayWtau = 0
				
				for i in range(numWDaughters):
					
					daughter = W.daughter(i)
					if not abs(daughter.pdgId()) == 15:
						daughter = getLastCopyStatusOne(daughter)
					else:
						daughter = getLastCopy(daughter)
					
					if daughter == None: continue
					
					wdaughter_var_array['pdgIdWdaughter'][i] = daughter.pdgId()
					wdaughter_var_array['ptWdaughter'][i] = daughter.pt()
					wdaughter_var_array['etaWdaughter'][i] = daughter.eta()
					wdaughter_var_array['phiWdaughter'][i] = daughter.phi()
					
					if abs(daughter.pdgId()) == 15:  # seems to be awfully complicated just to get the tau decay mode
						
						for k in range(daughter.numberOfDaughters()):
							if abs(daughter.daughter(k).pdgId()) == 16: continue
							decayWtau = daughter.daughter(k).pdgId()
							if abs(daughter.daughter(k).pdgId()) == 24:
								decayWtau = daughter.daughter(k).daughter(0).pdgId()
						
					if (abs(daughter.pdgId()) == 12 or abs(daughter.pdgId()) == 14 or abs(daughter.pdgId()) == 16): ptWneutrino = daughter.pt()
			
		event_level_var_array['numW'][0] = numW
		event_level_var_array['ptW'][0] = ptW
		event_level_var_array['etaW'][0] = etaW
		event_level_var_array['phiW'][0] = phiW
		event_level_var_array['numWDaughters'][0] = numWDaughters
		event_level_var_array['ptWneutrino'][0] = ptWneutrino
		event_level_var_array['decayWtau'][0] = decayWtau
		
		
		'''
		###############################################################################################
		# jet ID and JECs
		###############################################################################################
		'''

		jetsP4Raw = []
		jetsP4Corr = []
		jetsIdxGood = []
		badJet = False
		numBadJets = 0
		minetaabsbadjets = 10
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
			nconstituents = jet.numberOfDaughters()
			cm = jet.chargedMultiplicity()
			goodJet = \
				nhf < 0.99 and \
				nef < 0.99 and \
				chf > 0.00 and \
				cef < 0.99 and \
				nconstituents > 1 and \
				cm > 0

			if goodJet:
				jetsIdxGood.append(ijet)
			else:
				numBadJets += 1
				if abs(jet.eta()) < minetaabsbadjets: minetaabsbadjets = abs(jet.eta())
				if abs(jet.eta()) < 2.4 or (jet.pt() > 30 and abs(jet.eta()) < 5):
					badJet = True
					break

		if badJet: continue
		
		###########################################################################################veto
		
		hCutflow.Fill(4)
		cutflow = 4
		
		event_level_var_array['numbadjets'][0] = numBadJets
		event_level_var_array['minetaabsbadjets'][0] = minetaabsbadjets

		jets = [j for ij, j in enumerate(jets) if ij in jetsIdxGood]
		
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
		
		electronsCleaned = 0
		muonsCleaned = 0
		invm = -1
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
			
			if tracks == None: continue
			
			if electronsCleaned: electrons = collection
			else: muons = collection
			
		###########################################################################################veto
		
		hCutflow.Fill(5)
		cutflow = 5	

		event_level_var_array['electronsCleaned'][0] = electronsCleaned
		event_level_var_array['muonsCleaned'][0] = muonsCleaned
		event_level_var_array['invmCleaning'][0] = invm
		event_level_var_array['l1ptCleaning'][0] = l1pt
		event_level_var_array['l2ptCleaning'][0] = l2pt
		event_level_var_array['l1etaCleaning'][0] = l1eta
		event_level_var_array['l2etaCleaning'][0] = l2eta
		event_level_var_array['l1phiCleaning'][0] = l1phi
		event_level_var_array['l2phiCleaning'][0] = l2phi
		event_level_var_array['l1absisodbeta'][0] = l1absisodbeta
		event_level_var_array['l1relisodbeta'][0] = l1relisodbeta
		event_level_var_array['l2absisodbeta'][0] = l2absisodbeta
		event_level_var_array['l2relisodbeta'][0] = l2relisodbeta
			
		'''
		###############################################################################################
		# GEN MET and HT(miss) and FastSim MET correction
		###############################################################################################
		'''
		
		pTneutrinosum = 0
		genmetpt = -1
		genmetphi = -10
		genht = -1
		genhtmiss = -1
		nofastsimcorrmetpt = -1
		nofastsimcorrmetphi = -10
		if not 'data' in options.tag:
			
			allneutrinos = [gp for gp in genparticles if gp.status() == 1 and (abs(gp.pdgId()) == 12 or abs(gp.pdgId()) == 14 or abs(gp.pdgId()) == 16)]
			
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
			
		if 'signal' in options.tag:
			
			nofastsimcorrmetpt = met.pt()
			nofastsimcorrmetphi = met.phi()
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
				
				decaylength3D = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(),2)
					+ pow(gp.vy() - gp.daughter(0).vy(),2)
					+ pow(gp.vz() - gp.daughter(0).vz(),2))

				if not decaylength3D > 0: continue
				
				hasChargino += 1
				
				chipmmGEN = round(gp.mass(), 2)
				
				deltamGEN = round((gp.mass() - gp.daughter(0).mass()), 2)
				
				chipmptGEN = gp.pt()
				chipmetaGEN = gp.eta()
				chipmphiGEN = gp.phi()
				
				decaylengthXY = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(),2)
					+ pow(gp.vy() - gp.daughter(0).vy(),2))

				decaylengthZ = abs(gp.vz() - gp.daughter(0).vz())
				
				c1daughters = findDaughters(gp)
				
				chipmnumdaughters = len(c1daughters)
				
				for c1d in c1daughters:
					
					if abs(c1d.pdgId()) == 211:
				
						pion = c1d
						
						hasPion += 1
						
						chidaughterpdgid = pion.pdgId()
						
						pionptGEN = pion.pt()
						pionetaGEN = pion.eta()
						pionphiGEN = pion.phi()
					
						
						idxold, drminold = findMatch_track_old(pion, tracks)
						_, drminoldrandom = findMatch_track_old_random(pion, tracks)
					
						idx, dxyzmin, tminmatching, drmin = findMatch_track_new(pion, tracks)
						_, dxyzminrandom, _, drminrandom = findMatch_track_new_random(pion, tracks)
	
						if not idx == -1:
							if drmin < matchingDrThreshold and dxyzmin < matchingDxyzThreshold:
								hasMatchedTrack += 1
								if not hasMatchedTrack:
									matchedTrackIdxCharginoPion1 = idx
								else:
									matchedTrackIdxCharginoPion2 = idx
								
# 				else:
# 					
# 					lepton, neutrino = findLeptonNeutrino(gp)
# 					
# 					if not lepton == None:
# 						
# 						chidaughterpdgid = lepton.pdgId()
# 					
# 						idxlepton, dxyzminlepton, _, drminlepton = findMatch_track_new(lepton, tracks)
# 						
# 						if not idxlepton == -1:
# 							if drminlepton < matchingDrThreshold and dxyzminlepton < matchingDxyzThreshold:
# 								susytracks[idxlepton] = 1

		
			for gp in N2s:
				
				decaylength3DN2 = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(),2)
					+ pow(gp.vy() - gp.daughter(0).vy(),2)
					+ pow(gp.vz() - gp.daughter(0).vz(),2))

				if not decaylength3DN2 > 0: continue
				
				chiN2mGEN = round(gp.mass(), 2)
				
				deltamN2GEN = round((gp.mass() - gp.daughter(0).mass()), 2)
				
				chiN2ptGEN = gp.pt()
				chiN2etaGEN = gp.eta()
				chiN2phiGEN = gp.phi()
				
				decaylengthXYN2 = ROOT.TMath.Sqrt(pow(gp.vx() - gp.daughter(0).vx(),2)
					+ pow(gp.vy() - gp.daughter(0).vy(),2))

				decaylengthZN2 = abs(gp.vz() - gp.daughter(0).vz())
				
				n2daughters = findDaughters(gp)
				
				chiN2numdaughters = len(n2daughters)
				
# 				print ''
# 				print 'n2'
# 				
# 				for n2d in n2daughters:
# 					
# 					print n2d.pdgId()
# 					
# 				ancestors = []
# 				
# 				for idaughter in range(gp.numberOfDaughters()):
# 					
# 					daughter = gp.daughter(idaughter)
# 					pdgIdDaughter = abs(daughter.pdgId())
# 					
# 					if pdgIdDaughter == 11 or pdgIdDaughter == 13:
# 						
# 						daughter = getLastCopyStatusOne(daughter)
# 						
# 						if daughter not in ancestors: ancestors.append(daughter)
# 						
# 					else:
# 						
# 						if daughter not in ancestors: ancestors.append(daughter)
# 						
# 						while len([a for a in ancestors if a.status() != 1]) > 0:
# 							
# 							for b in [a for a in ancestors if a.status() != 1]:
# 								
# 								for ibdaughter in range(b.numberOfDaughters()):
# 									
# 									if b.daughter(ibdaughter) not in ancestors: ancestors.append(b.daughter(ibdaughter))
# 									
# 								ancestors.remove(b)
# 								
# 				for a in ancestors:
# 					
# 					if a.charge() == 0: continue
# 					
# 					chiN2daughterpdgid = a.pdgId()
# 					
# 					idxN2, dxyzminN2, _, drminN2 = findMatch_track_new(a, tracks)
# 					
# 					if not idxN2 == -1:
# 						if drminN2 < matchingDrThreshold and dxyzminN2 < matchingDxyzThreshold:
# 							susytracks[idxN2] = abs(a.pdgId())
			
			
			numchidaughters = chipmnumdaughters + chiN2numdaughters
			
			i = 0
			
			for c1d in c1daughters:
				
				chidaughter_var_array['motherchidaughter'][i] = 1
				chidaughter_var_array['pdgIdchidaughter'][i] = c1d.pdgId()
				chidaughter_var_array['ptchidaughter'][i] = c1d.pt()
				chidaughter_var_array['etachidaughter'][i] = c1d.eta()
				chidaughter_var_array['phichidaughter'][i] = c1d.phi()
				
				chidaughter_var_array['hasmatchedtrackchidaughter'][i] = 0
				
				if c1d.charge() != 0:
				
					idxC1, dxyzminC1, _, drminC1 = findMatch_track_new(c1d, tracks)
						
					if not idxC1 == -1:
						if drminC1 < matchingDrThreshold and dxyzminC1 < matchingDxyzThreshold:
							susytracks[idxC1] = (1, c1d.pdgId())
							chidaughter_var_array['hasmatchedtrackchidaughter'][i] = 1
				
				i += 1
				
			for n2d in n2daughters:
				
				chidaughter_var_array['motherchidaughter'][i] = 2
				chidaughter_var_array['pdgIdchidaughter'][i] = n2d.pdgId()
				chidaughter_var_array['ptchidaughter'][i] = n2d.pt()
				chidaughter_var_array['etachidaughter'][i] = n2d.eta()
				chidaughter_var_array['phichidaughter'][i] = n2d.phi()
				
				chidaughter_var_array['hasmatchedtrackchidaughter'][i] = 0
				
				if n2d.charge() != 0: 
				
					idxN2, dxyzminN2, _, drminN2 = findMatch_track_new(n2d, tracks)
						
					if not idxN2 == -1:
						if drminN2 < matchingDrThreshold and dxyzminN2 < matchingDxyzThreshold:
							susytracks[idxN2] = (2, n2d.pdgId())
							chidaughter_var_array['hasmatchedtrackchidaughter'][i] = 1
					
				i += 1
				
		'''
		###############################################################################################
		# get event-level info
		###############################################################################################
		'''
		
		event_level_var_array['numpvs'][0] = len(primaryvertices)
		event_level_var_array['rho'][0] = rho
		
		
		numjets = 0
		numjets30 = 0
		numjets50 = 0
		numjets100 = 0
		numjets200 = 0
		ht = 0
		htmissTlv = ROOT.TLorentzVector()
		idxhighestptjet = 0
		for ijet, jet in enumerate(jets):
			
			if abs(jet.eta()) < 2.4 and jet.pt() > 30: ht += jet.pt()
			if abs(jet.eta()) < 5 and jet.pt() > 30:
				jetTlv = ROOT.TLorentzVector(jet.px(), jet.py(), jet.pz(), jet.energy())
				htmissTlv -= jetTlv

			if abs(jet.eta()) >= 2.4: continue
			
			numjets += 1
			
			jetpt = jet.pt()
			if jetpt > 30: numjets30 += 1
			if jetpt > 50: numjets50 += 1
			if jetpt > 100: numjets100 += 1
			if jetpt > 200: numjets200 += 1

			if jetpt > jets[idxhighestptjet].pt(): idxhighestptjet = ijet
			
		event_level_var_array['numjets'][0] = numjets
		event_level_var_array['numjets30'][0] = numjets30
		event_level_var_array['numjets50'][0] = numjets50
		event_level_var_array['numjets100'][0] = numjets100
		event_level_var_array['numjets200'][0] = numjets200
		event_level_var_array['ht'][0] = ht
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

		if not min(dphimetjets) > 0.5: continue
		
		###########################################################################################veto
		
		hCutflow.Fill(8)
		cutflow = 8
		
		
		chpfcandsforiso = [p for p in pfcands if passesPreselection_iso_pfc(p, pv_pos, dz_threshold=0.1)]
				
		
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
			electron_var_array['pfabsisoelectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, pfcands)
			electron_var_array['pfabsisominielectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, pfcands, isMini=True)
			electron_var_array['chpfabsisoelectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, chpfcandsforiso)
			electron_var_array['chpfabsisominielectron'][ie], _, _, _ = calcIso_pf_or_track_new(e, chpfcandsforiso, isMini=True)
			electron_var_array['jetisoelectron'][ie], electron_var_array['jetisomultielectron'][ie], electron_var_array['jetdrminelectron'][ie] = calcIso_jet_new(e, jets)
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
			muon_var_array['pfabsisomuon'][im], _, _, _ = calcIso_pf_or_track_new(m, pfcands)
			muon_var_array['pfabsisominimuon'][im], _, _, _ = calcIso_pf_or_track_new(m, pfcands, isMini=True)
			muon_var_array['chpfabsisomuon'][im], _, _, _ = calcIso_pf_or_track_new(m, chpfcandsforiso)
			muon_var_array['chpfabsisominimuon'][im], _, _, _ = calcIso_pf_or_track_new(m, chpfcandsforiso, isMini=True)
			muon_var_array['jetisomuon'][im], muon_var_array['jetisomultimuon'][im], muon_var_array['jetdrminmuon'][im] = calcIso_jet_new(m, jets)
			muon_var_array['chhadisomuon'][im] = m.pfIsolationR03().sumChargedHadronPt
			muon_var_array['challisomuon'][im] = m.pfIsolationR03().sumChargedParticlePt
			muon_var_array['neuhadisomuon'][im] = m.pfIsolationR03().sumNeutralHadronEt
			muon_var_array['photisomuon'][im] = m.pfIsolationR03().sumPhotonEt
			muon_var_array['puchhadisomuon'][im] = m.pfIsolationR03().sumPUPt
			muon_var_array['absisodbetamuon'][im] = calcIso_dBeta(m.pfIsolationR03())
			muon_var_array['relisodbetamuon'][im] = calcIso_dBeta(m.pfIsolationR03()) / m.pt()

			if calcIso_dBeta(m.pfIsolationR03()) / m.pt() < 0.2: nummuonsiso += 1
		
		event_level_var_array['nummuonsiso'][0] = nummuonsiso
		
		
		hNumleptons.Fill(numelectronsiso+nummuonsiso)
		
		if not 'noleptonveto' in options.tag:
			if not (numelectronsiso+nummuonsiso) == 0: continue
		
		###########################################################################################veto
		
		hCutflow.Fill(9)
		cutflow = 9
		
				
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
			photon_var_array['pfabsisophoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, pfcands, dontSubtractObject=True)
			photon_var_array['pfabsisominiphoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, pfcands, isMini=True, dontSubtractObject=True)
			photon_var_array['chpfabsisophoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, chpfcandsforiso, dontSubtractObject=True)
			photon_var_array['chpfabsisominiphoton'][ip], _, _, _ = calcIso_pf_or_track_new(p, chpfcandsforiso, isMini=True, dontSubtractObject=True)
			photon_var_array['jetisophoton'][ip], photon_var_array['jetisomultiphoton'][ip], photon_var_array['jetdrminphoton'][ip] = calcIso_jet_new(p, jets)
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
		
		
		pfleptons = [l for l in pfcands if ((abs(l.pdgId()) == 11 or abs(l.pdgId()) == 13)
										and l.pt() > 10
										and abs(l.eta()) < 2.4
										and not (abs(l.eta()) > 1.4442 and abs(l.eta()) < 1.566))]
		
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
			_, pfrelisopflepton, _, _ = calcIso_pf_or_track_new(l, pfcands)
			pflepton_var_array['pfrelisopflepton'][il] = pfrelisopflepton
			_, pflepton_var_array['pfrelisominipflepton'][il], _, _ = calcIso_pf_or_track_new(l, pfcands, isMini=True)
			_, pflepton_var_array['chpfrelisopflepton'][il], _, _ = calcIso_pf_or_track_new(l, chpfcandsforiso)
			_, pflepton_var_array['chpfrelisominipflepton'][il], _, _ = calcIso_pf_or_track_new(l, chpfcandsforiso, isMini=True)
			pflepton_var_array['jetisopflepton'][il], pflepton_var_array['jetisomultipflepton'][il], pflepton_var_array['jetdrminpflepton'][il] = calcIso_jet_new(l, jets)
			
			if pfrelisopflepton < 0.2: numpfleptonsiso += 1
			
		event_level_var_array['numpfleptonsiso'][0] = numpfleptonsiso
		
				
		event_level_var_array['numtaus'][0] = len(tauswithdiscriminators)
		
		numtausiso = 0
		for it, t in enumerate(tauswithdiscriminators):
			
			tau_var_array['chargetau'][it] = t[0].charge()
			tau_var_array['pxtau'][it] = t[0].px()
			tau_var_array['pytau'][it] = t[0].py()
			tau_var_array['pztau'][it] = t[0].pz()
			tau_var_array['pttau'][it] = t[0].pt()
			tau_var_array['energytau'][it] = t[0].energy()
			tau_var_array['etatau'][it] = t[0].eta()
			tau_var_array['phitau'][it] = t[0].phi()
			if (not t[0].leadPFChargedHadrCand() == None) and t[0].leadPFChargedHadrCand().trackRef().isNonnull():
				tau_var_array['dztau'][it] = t[0].leadPFChargedHadrCand().trackRef().get().dz(pv_pos)
				tau_var_array['dxytau'][it] = t[0].leadPFChargedHadrCand().trackRef().get().dxy(pv_pos)
			else:
				tau_var_array['dztau'][it] = -1
				tau_var_array['dxytau'][it] = -1
			tau_var_array['chhadisotau'][it] = t[0].isolationPFChargedHadrCandsPtSum()
			tau_var_array['photisotau'][it] = t[0].isolationPFGammaCandsEtSum()
			tau_var_array['decaymodetau'][it] = t[0].decayMode()
			tau_var_array['mvadiscrtau'][it] = t[3]
			
			if t[2] > 0.5: numtausiso += 1
			
		event_level_var_array['numtausiso'][0] = numtausiso
		hNumtaus.Fill(numtausiso)

		
		btagvalues = []
		i = 0
		for jet in jets:
			
			if abs(jet.eta()) >= 2.4: continue
			
			jet_var_array['ptjet'][i] = jet.pt()
			jet_var_array['etajet'][i] = jet.eta()
			jet_var_array['phijet'][i] = jet.phi()
			jet_var_array['pxjet'][i] = jet.px()
			jet_var_array['pyjet'][i] = jet.py()
			jet_var_array['pzjet'][i] = jet.pz()
			jet_var_array['energyjet'][i] = jet.energy()
			jet_var_array['nconstituentsjet'][i] = jet.numberOfDaughters()
			
			dRmin = 0.1
			thisbtag = -2.
			for ib in range(nbtags):
				dR = deltaR(handle_btag.product().key(ib).get().eta(), jet.eta(), handle_btag.product().key(ib).get().phi(), jet.phi())
				if dR < dRmin:
					dRmin = dR
					thisbtag = max(-1.0, handle_btag.product().value(ib))
			
			jet_var_array['btagjet'][i] = thisbtag
			btagvalues.append((thisbtag, jet.pt(), jet.eta()))
			hBtagjets.Fill(thisbtag)
			
			i += 1

		njetsbtagmedium = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > 0.8484 and jetpt > 30 and abs(jeteta) < 2.4)])
		hNjetsbtagmedium.Fill(njetsbtagmedium)
		event_level_var_array['njetsbtagmedium'][0] = njetsbtagmedium
		event_level_var_array['njetsbtagmediumTIGHT'][0] = len([bt for (bt, jetpt, jeteta) in btagvalues if (bt > 0.8484 and jetpt > 15 and abs(jeteta) < 2.4)])
		
		mtmetleadingjet = ROOT.TMath.Sqrt(2 * met.pt() * jets[idxhighestptjet].pt() 
			* (1 - ROOT.TMath.Cos(deltaPhi(met.phi(), jets[idxhighestptjet].phi()))))
		hMtmetleadingjet.Fill(mtmetleadingjet)
		event_level_var_array['mtmetleadingjet'][0] = mtmetleadingjet

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
		
		event_level_var_array['chidaughterpdgid'][0] = chidaughterpdgid  # TODO: change to nums
		event_level_var_array['chiN2daughterpdgid'][0] = chiN2daughterpdgid
		

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
		
		tracksforiso = [t for t in tracks if passesPreselection_iso_track(t, pv_pos, dz_threshold=0.1)]
		jetsforisotight = [j for j in jets if passesPreselection_iso_jet(j, pt_threshold=30)]
		
		if 'genmatchtracks' in options.tag or 'genmatchalltracks' in options.tag:
			genparticlesformatching = [gp for gp in genparticles if gp.status() == 1]

		numtracksbasicpreselection = 0
		numtracksfinalpreselection = 0
		
		i = 0
		for itrack, track in enumerate(tracks):

			if not passesPreselection_basic_track(track): continue
			
			numtracksbasicpreselection += 1
			
			#TODO: adapt preselection
			if not abs(track.dz(pv_pos)) < 1: continue
			jetisotight, jetisomultitight, jetdrmintight = calcIso_jet_new(track, jetsforisotight, isTrack=True)
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
			
			track_level_var_array['eta'][i] = track.eta()
			track_level_var_array['phi'][i] = track.phi()
			
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
			if abs(track.dz(pv_pos)) >= 0.1: dontSubtractTrackPt = True
			track_level_var_array['trackabsiso'][i], track_level_var_array['trackreliso'][i], track_level_var_array['trackdrmin'][i], track_level_var_array['tracknumneighbours'][i] = calcIso_pf_or_track_new(track, tracksforiso, dontSubtractObject=dontSubtractTrackPt)
			track_level_var_array['pfabsiso'][i], track_level_var_array['pfreliso'][i], track_level_var_array['pfdrmin'][i], track_level_var_array['pfnumneighbours'][i] = calcIso_pf_or_track_new(track, pfcands)
			track_level_var_array['chpfabsiso'][i], track_level_var_array['chpfreliso'][i], track_level_var_array['chpfdrmin'][i], track_level_var_array['chpfnumneighbours'][i] = calcIso_pf_or_track_new(track, chpfcandsforiso, dontSubtractObject=dontSubtractTrackPt)
			track_level_var_array['jetiso'][i], track_level_var_array['jetisomulti'][i], track_level_var_array['jetdrmin'][i] = calcIso_jet_new(track, jets, isTrack=True)
			track_level_var_array['jetisotight'][i] = jetisotight
			track_level_var_array['jetisomultitight'][i] = jetisomultitight
			track_level_var_array['jetdrmintight'][i] = jetdrmintight
			
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
			
			drmintau = 10
			closesttaumvadiscr = -1
			for t in tauswithdiscriminators:
				dr = deltaR(track.eta(), t[0].eta(), track.phi(), t[0].phi())
				if dr < drmintau:
					drmintau = dr
					closesttaumvadiscr = t[3]
			track_level_var_array['drmintau'][i] = drmintau
			track_level_var_array['closesttaumvadiscr'][i] = closesttaumvadiscr
			
			track_level_var_array['detahighestptjet'][i] = abs(track.eta() - jets[idxhighestptjet].eta())
			track_level_var_array['dphihighestptjet'][i] = deltaPhi(track.phi(), jets[idxhighestptjet].phi())
			
			track_level_var_array['dphimet'][i] = deltaPhi(track.phi(), met.phi())
			track_level_var_array['dphimetpca'][i], _, _ = handmadeDphiMetPCA(track, pv_pos, met)
					
			track_level_var_array['chi2'][i] = 	track.normalizedChi2()
			
			quality = 0
			if track.quality(track.qualityByName("loose")): quality = 1
			if track.quality(track.qualityByName("tight")): quality = 2
			if track.quality(track.qualityByName("highPurity")): quality = 3
			track_level_var_array['quality'][i] = quality
			
			track_level_var_array['nvalidhits'][i] = track.numberOfValidHits()
			track_level_var_array['nlosthits'][i] = track.numberOfLostHits()
			
			
			hasGenMatch = -1
			genmatchpdgid = -1
			genmatchmotherpdgid = -1
			genmatchstatus = -1
			genmatchmotherstatus = -1
			genmatchishardprocess = -1
			genmatchmotherishardprocess = -1
			genmatchisfromhardprocess = -1
			genmatchisprompt = -1
			genmatchisdirecthadrondecayproduct = -1
			genmatchisdirecttaudecayproduct = -1
			
			dothegenmatch = False

			if 'genmatchtracks' in options.tag:
				if ievent % 10 == 0 and met.pt() > 250: dothegenmatch = True
				
			if 'genmatchalltracks' in options.tag:
				if met.pt() > 250: dothegenmatch = True
			
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
					while not gp.statusFlags().isFirstCopy():
						for idx in range(gp.numberOfMothers()):
							if gp.pdgId() == gp.mother(idx).pdgId():
								gp = gp.mother(idx)
								break
						mm += 1
						if mm > 100: break
					
					genmatchmother = genmatch.mother(0)
					
					genmatchpdgid = genmatch.pdgId()
					genmatchmotherpdgid = genmatchmother.pdgId()
					
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
						
			track_level_var_array['hasGenMatch'][i] = hasGenMatch
			track_level_var_array['genmatchpdgid'][i] = genmatchpdgid
			track_level_var_array['genmatchmotherpdgid'][i] = genmatchmotherpdgid
			track_level_var_array['genmatchstatus'][i] = genmatchstatus
			track_level_var_array['genmatchmotherstatus'][i] = genmatchmotherstatus
			track_level_var_array['genmatchishardprocess'][i] = genmatchishardprocess
			track_level_var_array['genmatchmotherishardprocess'][i] = genmatchmotherishardprocess
			track_level_var_array['genmatchisfromhardprocess'][i] = genmatchisfromhardprocess
			track_level_var_array['genmatchisprompt'][i] = genmatchisprompt
			track_level_var_array['genmatchisdirecthadrondecayproduct'][i] = genmatchisdirecthadrondecayproduct
			track_level_var_array['genmatchisdirecttaudecayproduct'][i] = genmatchisdirecttaudecayproduct
				
			
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
	
	#Zll cleaning histos
	hZllLeptonPt.Write()
	hZllDrTrack.Write()
	hZllDrPfc.Write()
	hZllDrJet.Write()
	
	#PV histos
	hNPVsPerEvent.Write()
	hPV0x.Write()
	hPV0y.Write()
	hPV0z.Write()
	hPVsx.Write()
	hPVsy.Write()
	hPVsz.Write()
	
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
				fo.write('')
				
			print 'just created empty ' + fout.GetName().replace('.root', '.json')
		
	fout.Close()
