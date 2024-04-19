import os
import sys
from ROOT import *
from array import array
from math import sqrt, pi
import re, array
# import cms_figure
import numpy as np
# from numpy import linspace
# from decimal import *
from glob import glob

import tensorflow as tf
from tensorflow.keras.models import load_model
# from tmva_include import book_tmva, fill_tmva
import uproot
# import gfal2
from itertools import product
import awkward as ak



# gSystem.Load("libFWCoreFWLite.so")
# gSystem.Load("libDataFormatsFWLite.so")
# FWLiteEnabler.enable()
# from DataFormats.FWLite import Events, Handle

# gROOT.SetBatch(True)
# gStyle.SetOptStat(111111)


"""
### , 
### , 
### , 
def addBDTScores(
        fin='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/test.root',
        MVATree='Data/TestTree',
        MVAFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/simpleTMVA.root',
        weightFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/Data/weights/TMVAClassification_DNN.weights.xml',

):
    fout = TFile(fin)
    # fout = TFile(fin, "update")
    # fout = TFile("root://dcache-cms-xrootd.desy.de:1094/"+fin, "update")
    tree = fout.Get(
        'tEvent')  # todo: i think the problem lies here: I should make this tree is still bound th the fut, but gets a new branch: that gives an error
    # tree = TTree()
    # fout.Get('tEvent').Copy(tree)

    print "Retrieve tree"

    max_length = -1

    # bdt_scores = array.array("f", 1000*[-999]) #too Fix
    bdt_scores = std.vector('double')()
    # bdt_branch = tree.Branch("bdt_score", bdt_scores, "bdt_score[n_sv]/F")
    bdt_branch = tree.Branch("bdt_score", bdt_scores)

    print "Retrieve MVA values "
    print "								"

    chain_BDT = TChain('Data/TestTree')
    chain_BDT.Add(MVAFile)

    ### Booking the reader for all variables in the BDT test tree
    reader = TMVA.Reader()
    print chain_BDT
    print chain_BDT.GetListOfBranches()
    var = {}
    BDTvalues = {}
    for b, branch in enumerate(chain_BDT.GetListOfBranches()):
        branch_name = branch.GetName()
        if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
        var[branch_name] = array.array("f", [0])
        reader.AddVariable(branch_name, var[branch_name])

    reader.BookMVA("BDT", weightFile)

    nevents = tree.GetEntries()
    for ievent in range(nevents):
        tree.GetEntry(ievent)
        nsv = tree.n_sv
        BDTvalues = [None] * nsv
        bdt_scores.clear()

        # if nsv==0: continue
        print ievent, nsv
        for isv in range(nsv):
            for b, branch in enumerate(chain_BDT.GetListOfBranches()):
                branch_name = branch.GetName()
                if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70",
                                   "DNN"]: continue
                var[branch_name][0] = getattr(tree, branch_name)[isv]

            mybdt_for_this_sv = reader.EvaluateMVA("BDT")
            print ievent, isv, mybdt_for_this_sv
            # bdt_scores[isv]= mybdt_for_this_sv
            bdt_scores.push_back(mybdt_for_this_sv)

        bdt_branch.Fill()
        if nsv > 0:
            tree.Show(ievent)
            # print " single val"
            # print bdt_scores
            for val in tree.bdt_score: print val
        # exit()
    print "end of the script"

    fout.cd()
    tree.Write("tEvent", TObject.kOverwrite)

    return tree
"""


def joinTreeByFriend(tree, friendFile, friendTree):
    try:
        ff = ROOT.TFile(friendFile)
        tf = ff.Get(friendTree)
        tree.AddFriend(tf)
    except:
        print("failed to add friend", friendFile, "to tree")


### ,
### , 
### , 
# # # >>> import ROOT
# # # >>> f = ROOT.TFile("/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm1p168GeV_Chi20ctau5MM_part5of100_NTuple.root")
# # # >>> t = f.Get("tEvent")
# # # >>> ff = ROOT.TFile("/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm1p168GeV_Chi20ctau5MM_part5of100_NTuple_friend.root")
# # # >>> tf = ff.Get("tFriend")
# # # >>> t.AddFriend(tf)
# # # <ROOT.TFriendElement object ("tFriend") at 0x54748f0>
# # # >>> t.Draw("deltaR")
# # # Info in <TCanvas::MakeDefCanvas>:  created default TCanvas with name c1
# # # >>> t.Draw("bdt_score")

def makeFriendTree(
        fin='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm115GeV_dm1p168GeV_Chi20ctau5MM_part5of100_NTuple.root',
        outpath='/nfs/dust/cms/user/tewsalex/rootfiles/friendTrees_V12/',
        MVATree='Data/TestTree',
        MVAFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/rootfiles/simpleTMVA.root',
        weightFile='/nfs/dust/cms/user/tewsalex/CMSSW_10_5_0/src/Data/weights/TMVAClassification_DNN.weights.xml',
):
    finTree = TFile(fin)
    tree = finTree.Get('tEvent')

    outname = fin.split('/')[-1].split('.root')[0]
    fout = TFile(outpath + outname + '_friend.root', 'recreate')

    friendTree = TTree('tFriend', 'tFriend')

    max_length = -1

    bdt_scores = std.vector('double')()
    is_maxscore = std.vector('int')()
    bdt_branch = friendTree.Branch("bdt_score", bdt_scores)
    bdt_maxbranch = friendTree.Branch("is_maxscoring_sv", is_maxscore)

    print("Retrieve MVA values ")
    print("								")

    chain_BDT = TChain('Data/TestTree')
    chain_BDT.Add(MVAFile)

    ### Booking the reader for all variables in the BDT test tree
    reader = TMVA.Reader()
    print(chain_BDT)
    print(chain_BDT.GetListOfBranches())
    var = {}
    BDTvalues = {}
    for b, branch in enumerate(chain_BDT.GetListOfBranches()):
        branch_name = branch.GetName()
        if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
        var[branch_name] = array.array("f", [0])
        reader.AddVariable(branch_name, var[branch_name])

    reader.BookMVA("BDT", weightFile)

    nevents = tree.GetEntries()
    for ievent in range(nevents):
        tree.GetEntry(ievent)
        nsv = tree.n_sv
        BDTvalues = [None] * nsv
        bdt_scores.clear()
        is_maxscore.clear()

        maxvalue = double(-999)
        # print ievent, nsv
        for isv in range(nsv):
            for b, branch in enumerate(chain_BDT.GetListOfBranches()):
                branch_name = branch.GetName()
                if branch_name in ["classID", "className", "weight", "BDT1000", "BDT150", "BDT40", "BDT70", "DNN"]: continue
                var[branch_name][0] = getattr(tree, branch_name)[isv]

            mybdt_for_this_sv = reader.EvaluateMVA("BDT")
            # print ievent, isv, mybdt_for_this_sv
            if mybdt_for_this_sv > maxvalue: maxvalue = mybdt_for_this_sv

            bdt_scores.push_back(mybdt_for_this_sv)

        if nsv > 0:
            for val in friendTree.bdt_score:
                # print val
                if val == maxvalue:
                    is_maxscore.push_back(1)
                else:
                    is_maxscore.push_back(0)
        else:
            is_maxscore.push_back(-1)

        friendTree.Fill()
    fout.Write('', TObject.kWriteDelete)
    # fout.cd()
    # friendTree.Write("tFriend", TObject.kOverwrite)

    print("added a tFriend to", fin, fout.GetName())


def makeFriendTree_keras(
        fin='/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV11/16UL/ZJetsToNuNu_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/2023_02_22/230222_111558/0000/crab_NTuple_1.root',
        outpath='/nfs/dust/cms/user/wolfmor/FrieNdTuples/Test/',
        my_models=None):

    model_names = [my_model.name for my_model in my_models]
    if len(model_names) > len(list(dict.fromkeys(model_names))):
        print('got models:')
        print(model_names)
        raise Exception('model names cannot be duplicates')

    # my_model.model.summary()
    # print(tf.config.list_physical_devices('GPU'))

    events = uproot.open(fin)['tEvent']

    variable_names = []
    special_names = []
    friendfiles_variables = {}
    for my_model in my_models:
        variable_names += [inp for inp in my_model.inputs if inp not in my_model.parameters.keys() and inp not in my_model.specialinputs.keys()]
        special_names += [my_model.specialinputs[s][0] for s in my_model.specialinputs]
        for ff in my_model.friendfilesvariables.keys():
            if ff in friendfiles_variables.keys():
                friendfiles_variables[ff] += my_model.friendfilesvariables[ff]
            else:
                friendfiles_variables[ff] = my_model.friendfilesvariables[ff]

    # remove duplicates
    variable_names = list(dict.fromkeys(variable_names))
    special_names = list(dict.fromkeys(special_names))
    for ff in friendfiles_variables.keys():
        friendfiles_variables[ff] = list(dict.fromkeys(friendfiles_variables[ff]))

    variables = events.arrays(variable_names, library='np')
    specials = events.arrays(special_names, library='np')

    for my_model in my_models:
        for variable in my_model.specialinputs:
            if variable not in variables:
                variables[variable] = my_model.specialinputs[variable][1](specials[my_model.specialinputs[variable][0]])

    variables.update(events.arrays(['random', 'track_random', 'track_quality',
                                    'track_isSignalTrack', 'track_isSusyTrack', 'track_susyTrackPdgIdMother', 'track_susyTrackPdgId',
                                    'track_hasGenMatch', 'track_genMatchIsPrompt', 'track_genMatchIsFromHardProcess',
                                    'track_genMatchMotherIsTheTau', 'track_genMatchMotherTauDecay'],
                                   library='np'))

    friends = {}
    for ff in friendfiles_variables.keys():
        friends[ff] = uproot.open(fin.replace('/NTuples/', '/FrieNdTuples/').replace('/' + ff[:-(1+ff[::-1].index('p'))] + '/', '/' + ff + '/').replace('.root', '_friend.root'))['tFriend']
        variables.update(friends[ff].arrays(friendfiles_variables[ff], library='np'))

    print('initialize output file')

    outname = fin.split('/')[-1].strip().split('.root')[0]
    foutname = outpath + outname + '_friend.root'

    with uproot.recreate(foutname) as fout:

        # TODO: only one counter for all branches
        branches = {'random_sanitycheck': 'float64', 'track_random_sanitycheck': 'var * float64'}

        parameter_points = {}
        parameter_points_labels = {}

        for my_model in my_models:

            parameter_points[my_model] = [dict(zip(my_model.parameters.keys(), values)) for values in product(*my_model.parameters.values())]
            parameter_points_labels[my_model] = ['_'.join([p + str(pp[p]).replace('.', 'p') for p in pp]) for pp in parameter_points[my_model]]

            for label in parameter_points_labels[my_model]:
                for out in my_model.outputs:

                    branches['track_' + my_model.name + '_' + out + '_' + label] = 'var * float64'

                    if out in my_model.savemaxscoreinfo:

                        branches['track_' + my_model.name + '_' + out + '_' + label + '_ranking'] = 'var * int32'
                        branches['track_' + my_model.name + '_' + out + '_' + label + '_isMaxscore'] = 'var * int32'

                        branches['maxscore_' + my_model.name + '_' + out + '_' + label] = 'float64'

                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_index'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_isSignalTrack'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_isSusyTrack'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_susyTrackPdgIdMother'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_susyTrackPdgId'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_hasGenMatch'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_genMatchIsPrompt'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_genMatchIsFromHardProcess'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_genMatchMotherIsTheTau'] = 'int32'
                        branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_genMatchMotherTauDecay'] = 'int32'

                        for othernode in my_model.outputs:
                            if othernode == out: continue
                            branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_' + othernode] = 'float64'


        fout.mktree('tFriend', branches)
        # fout['tFriend'].show()

        print('evaluate model')

        for event in range(len(variables['track_pt'])):

            n_track = np.array([len(variables['track_pt'][event])], dtype=np.int32)

            branches['random_sanitycheck'] = np.array([variables['random'][event]])  # ak.Array([variables['random'][event]])
            branches['track_random_sanitycheck'] = ak.Array([variables['track_random'][event]])
            branches['ntrack_random_sanitycheck'] = n_track

            for my_model in my_models:

                # TODO: what qualitymask?
                if type(my_model.qualitymask) == str:
                    if my_model.qualitymask == '20240328>0.1':
                        trackmask = (variables['track_PyKeras_V14_20240328_multiclass_Signal_'][event] > 0.1) \
                                    & (variables['track_quality'][event] > 2) \
                                    & (variables['track_abs_eta_'][event] < 2.4)
                    elif my_model.qualitymask == '20240328>0.2':
                        trackmask = (variables['track_PyKeras_V14_20240328_multiclass_Signal_'][event] > 0.2) \
                                    & (variables['track_quality'][event] > 2) \
                                    & (variables['track_abs_eta_'][event] < 2.4)
                    elif my_model.qualitymask == '20240328>0.5':
                        trackmask = (variables['track_PyKeras_V14_20240328_multiclass_Signal_'][event] > 0.5) \
                                    & (variables['track_quality'][event] > 2) \
                                    & (variables['track_abs_eta_'][event] < 2.4)
                    elif my_model.qualitymask == 'dz10':
                        trackmask = (variables['track_quality'][event] > 2) \
                                    & (variables['track_log10_dz_'][event] < 1) \
                                    & (variables['track_abs_eta_'][event] < 2.4)
                    elif my_model.qualitymask == 'dz0p5':
                        trackmask = (variables['track_quality'][event] > 2) \
                                    & (variables['track_log10_dz_'][event] < np.log10(0.5)) \
                                    & (variables['track_abs_eta_'][event] < 2.4)
                    elif my_model.qualitymask == 'dz0p1':
                        trackmask = (variables['track_quality'][event] > 2) \
                                    & (variables['track_log10_dz_'][event] < -1) \
                                    & (variables['track_abs_eta_'][event] < 2.4)
                    elif my_model.qualitymask == 'pt0p5':
                        trackmask = (variables['track_quality'][event] > 2) \
                                    & (variables['track_log10_dz_'][event] < 0) \
                                    & (variables['track_abs_eta_'][event] < 2.4) \
                                    & (variables['track_pt'][event] > 0.5)
                    elif my_model.qualitymask == 'pt1':
                        trackmask = (variables['track_quality'][event] > 2) \
                                    & (variables['track_log10_dz_'][event] < 0) \
                                    & (variables['track_abs_eta_'][event] < 2.4) \
                                    & (variables['track_pt'][event] > 1)
                    else:
                        raise NotImplementedError('cannot understand qualitymask')
                elif my_model.qualitymask:
                    trackmask = (variables['track_quality'][event] > 2) \
                                & (variables['track_log10_dz_'][event] < 0) \
                                & (variables['track_abs_eta_'][event] < 2.4)
                else:
                    trackmask = None

                for pp, ppl in zip(parameter_points[my_model], parameter_points_labels[my_model]):

                    inputs = np.stack([np.full_like(variables['track_pt'][event], pp[inp]) if inp in pp else
                                       (np.full_like(variables['track_pt'][event], variables[inp][event]) if inp in my_model.eventlevelintputs else
                                        variables[inp][event]) for inp in my_model.inputs], axis=1)
                    outputs = my_model.model.predict(inputs)

                    # print(ppl)
                    # print(inputs[0, :])
                    # print(outputs[0, :])

                    if trackmask is not None:
                        outputs[~trackmask, :] *= -1  # TODO: this obviously doesn't work for outputs < 0

                    preds = ak.Array([outputs])

                    # print()
                    # print(preds[:, :, 0])
                    # print(ak.max(preds[:, :, 0]))
                    # print(ak.argmax(preds[:, :, 0]))
                    # print(ak.Array([[int(i == ak.argmax(preds[:, :, 0])) for i in range(len(preds[:, :, 0][0]))]]))
                    # print()
                    # print(preds[:, 1])
                    # print(ak.max(preds[:, 1]))
                    # print(preds[:, 1] == ak.max(preds[:, 1]))
                    #
                    # sys.exit(0)

                    for iout, out in enumerate(my_model.outputs):

                        branches['track_' + my_model.name + '_' + out + '_' + ppl] = preds[:, :, iout]
                        branches['ntrack_' + my_model.name + '_' + out + '_' + ppl] = n_track

                        if out in my_model.savemaxscoreinfo:

                            highscoreidx = ak.argmax(preds[:, :, iout])

                            sorted_preds = sorted(preds[:, :, iout][0], reverse=True)
                            rankings = []
                            for p in preds[:, :, iout][0]:
                                try:
                                    rankings.append(sorted_preds.index(p))
                                except ValueError:  # nan is not in list
                                    rankings.append(-1)

                            branches['track_' + my_model.name + '_' + out + '_' + ppl + '_ranking'] = ak.Array([rankings])
                            branches['ntrack_' + my_model.name + '_' + out + '_' + ppl + '_ranking'] = n_track

                            branches['track_' + my_model.name + '_' + out + '_' + ppl + '_isMaxscore'] = ak.Array([[int(i == highscoreidx) for i in range(len(preds[:, :, iout][0]))]])
                            branches['ntrack_' + my_model.name + '_' + out + '_' + ppl + '_isMaxscore'] = n_track

                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl] = np.array([ak.max(preds[:, :, iout])])  # ak.Array([ak.max(preds[:, :, iout])])

                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_index'] = np.array([highscoreidx])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_isSignalTrack'] = np.array([variables['track_isSignalTrack'][event][highscoreidx]])  # ak.Array([variables['track_isSignalTrack'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_isSusyTrack'] = np.array([variables['track_isSusyTrack'][event][highscoreidx]])  # ak.Array([variables['track_isSusyTrack'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_susyTrackPdgIdMother'] = np.array([variables['track_susyTrackPdgIdMother'][event][highscoreidx]])  # ak.Array([variables['track_susyTrackPdgIdMother'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_susyTrackPdgId'] = np.array([variables['track_susyTrackPdgId'][event][highscoreidx]])  # ak.Array([variables['track_susyTrackPdgId'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_hasGenMatch'] = np.array([variables['track_hasGenMatch'][event][highscoreidx]])  # ak.Array([variables['track_hasGenMatch'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_genMatchIsPrompt'] = np.array([variables['track_genMatchIsPrompt'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_genMatchIsFromHardProcess'] = np.array([variables['track_genMatchIsFromHardProcess'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_genMatchMotherIsTheTau'] = np.array([variables['track_genMatchMotherIsTheTau'][event][highscoreidx]])  # ak.Array([variables['track_genMatchMotherIsTheTau'][event][highscoreidx]])
                            branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_genMatchMotherTauDecay'] = np.array([variables['track_genMatchMotherTauDecay'][event][highscoreidx]])

                            for iothernode, othernode in enumerate(my_model.outputs):
                                if othernode == out: continue
                                branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_' + othernode] = np.array([preds[:, highscoreidx, iothernode][0]])

            # for key in branches:
            #     print(key)
            #     print(branches[key])
            # print(branches)
            # print(ak.zip(branches))
            # print(ak.zip(branches).type)

            # if event == 0:
            #     fout.mktree('tFriend', ak.zip(branches).type)
            #     fout['tFriend'].show()

            fout['tFriend'].extend(branches)

            # if event > 10: break


    print("added a tFriend to", fin, foutname)


class MyModel:
    def __init__(self, name, h5file, inputs, outputs,
                 custom_objects=None, savemaxscoreinfo=None,
                 parameters=None, specialinputs=None, eventlevelintputs=None,
                 qualitymask=False, friendfilesvariables=None):

        self.name = name
        self.model = load_model(h5file, custom_objects=custom_objects)
        self.inputs = inputs
        self.outputs = outputs

        if savemaxscoreinfo is None:
            self.savemaxscoreinfo = outputs
        else:
            self.savemaxscoreinfo = savemaxscoreinfo

        if parameters is None:
            self.parameters = {}
        else:
            self.parameters = parameters

        if specialinputs is None:
            self.specialinputs = {}
        else:
            self.specialinputs = specialinputs

        if eventlevelintputs is None:
            self.eventlevelintputs = []
        else:
            self.eventlevelintputs = eventlevelintputs

        self.qualitymask = qualitymask

        if friendfilesvariables is None:
            self.friendfilesvariables = {}
        else:
            self.friendfilesvariables = friendfilesvariables


class MyVarParsing:

    def __init__(self, args):
        self.inputFiles = ''
        self.tag = ''
        for arg in args:
            unpack = arg.split('=')
            if len(unpack) == 2:
                if unpack[0] == 'inputFiles': self.inputFiles = unpack[1].replace(',', '').split(' ')
                if unpack[0] == 'tag': self.tag = unpack[1].split(' ')



### give the path to the folder that contains the ntuples for datasets on pnfs 
### or give a path to files with wildcards for signal files on dust / local
list_of_datasets = [
    # '/pnfs/desy.de/cms/tier2/store/user/altews/NTuples/NTuplesV12/16UL_preAPV/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/2023_04_28/230428_123855/0000/',
    # '/nfs/dust/cms/user/tewsalex/rootfiles/ntuple_V12/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm*GeV_dm*GeV_Chi20ctau5MM_part*of*_NTuple.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV11/16UL/ZJetsToNuNu_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/2023_02_22/230222_111558/0000/crab_NTuple_1.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p368GeV_part1of25_NTuple_noSVs_job385.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p968GeV_part20of25_NTuple_noSVs_job486.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm*GeV_part*_NTuple_noSVs_job*.root',

    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_500GeV_mChipm400GeV_dm1p0GeV_part16of25_NTuple_noSVs_job155.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_500GeV_mChipm400GeV_dm0p6GeV_part16of25_NTuple_noSVs_job130.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_500GeV_mChipm400GeV_dm0p3GeV_part16of25_NTuple_noSVs_job883.root',

    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV14_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_500GeV_mChipm400GeV_dm1p0GeV_part16of25_NTuple_noSVs_job203.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV14_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_500GeV_mChipm400GeV_dm0p3GeV_part16of25_NTuple_noSVs_job153.root',

    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV14_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_*GeV_mChipm*GeV_dm*GeV_*.root',

    '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV14/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm500GeV_dm0p63GeV_part6of25_NTuple_job3061.root'

    # '/pnfs/desy.de/cms/tier2/store/user/altews/NTuples/NTuplesV13/16UL_preAPV/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/2023_08_12/230812_075400/0000/crab_NTuple_1.root'
]
# outpath = '/nfs/dust/cms/user/tewsalex/rootfiles/friendTrees_V12/'
# outpath = '/nfs/dust/cms/user/wolfmor/FrieNdTuples/TestArthur/'
outpath = '/nfs/dust/cms/user/wolfmor/FrieNdTuples/Test/'


# my_model_Arthur = MyModel(
#     name='Arthur_ckpt_random_nocw',
#     h5file='/nfs/dust/cms/user/tanikulo/NN_keras_checkpoints/ckpt_random_nocw/Arthur_ckpt_random_nocw_Tracks_D64_D64_D64_D64_D3_ep500_b1024_lr0p0001.h5',
#     inputs=['track_pt', 'track_dphiMet',  'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_drminBjetMedium30', 'track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet', 'deltamFile'],
#     outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
#     parameters={'deltamFile': [0.3, 1.0]
# )

# my_model_Arthur = MyModel(
#     name='Arthur_ckpt_random_cw_4_100None',
#     h5file='/nfs/dust/cms/user/tanikulo/NN_keras_checkpoints/ckpt_random_cw_4_100None/Arthur_ckpt_random_cw_4_100None_Tracks_D64_D64_D64_D64_D3_ep5_b10000_lr0p001.h5',
#     inputs=['track_pt', 'track_dphiMet',  'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_drminBjetMedium30', 'track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet', 'deltamFile'],
#     outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
#     parameters={'deltamFile': [0.3, 1.0]}
# )

# my_model_Arthur = MyModel(
#     name='Arthur_best',
#     h5file='/nfs/dust/cms/user/tanikulo/NN_keras_checkpoints/best/Arthur_best_Tracks_D64_D64_D64_D64_D3_ep500_b1024_lr0p0001.h5',
#     inputs=['track_pt', 'track_dphiMet',  'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_drminBjetMedium30', 'track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet', 'deltamFile'],
#     outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
#     parameters={'deltamFile': [0.3, 1.0]}
# )

# my_model_Arthur = MyModel(
#     name='Arthur_ckpt_25',
#     h5file='/nfs/dust/cms/user/tanikulo/NN_keras_checkpoints/ckpt_25/Arthur_ckpt_25_Tracks_3dm_t70v15t15_D64_D64_D64_D64_D3_ep500_b50_lr0p000005_cwbalanced.h5',
#     inputs=['track_pt', 'track_dphiMet',  'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_drminBjetMedium30','track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet', 'deltamFile'],
#     outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
#     parameters={'deltamFile': [0.3, 1.0]}
# )

# my_model_Arthur = MyModel(
#     name='Arthur_ckpt_random_1',
#     h5file='/nfs/dust/cms/user/tanikulo/NN_keras_checkpoints/ckpt_random_1/Arthur_ckpt_random_1_Tracks_D64_D64_D64_D64_D3_ep5_b10000_lr0p001.h5',
#     inputs=['track_pt', 'track_dphiMet',  'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_drminBjetMedium30','track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet', 'deltamFile'],
#     outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
#     parameters={'deltamFile': [0.3, 1.0]}
# )

# my_model_Arthur = MyModel(
#     name='Arthur_final_NN',
#     h5file='/nfs/dust/cms/user/tanikulo/NN_keras_checkpoints/final_NN/Arthur_final_NN_Tracks_5xD32BN_D3_ep500_b1024_lr0p0001.h5',
#     inputs=['track_pt', 'track_dphiMet', 'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_drminBjetMedium30', 'track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet', 'track_pfAbsIso', 'track_numneighboursPf', 'track_drminTrack1', 'track_drminTrack5', 'track_drminTrack10', 'track_drmin2ndTrack1', 'track_drmin2ndTrack5', 'track_drmin2ndTrack10', 'deltamFile'],
#     outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
#     savemaxscoreinfo=['predicted_Signal'],
#     parameters={'deltamFile': [0.3, 0.6, 1.0]},
#     qualitymask=True
# )

# my_model_V11_20230630 = MyModel(
#     name='PyKeras_V11_20230630_multiclass_puretracklevel',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V11_puretracklevel/weights/TrainedModel_PyKeras_V11_20230702_multiclass_puretracklevel.h5',
#     inputs=['deltam', 'track_pt', 'track_abs_eta_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_', 'track_log10_dxy_', 'track_log10_dz_', 'track_log10_dxyPU_', 'track_log10_dzPU_', 'track_log10_dxyError_', 'track_log10_dzError_', 'track_neHadAbsIso0', 'track_tkAbsIso0', 'track_drminTrack10', 'track_drminJet20', 'track_drminJet30'],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     parameters={'deltam': [0.3, 1.0]},
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)]}
# )

# my_model_V14_20231025_1 = MyModel(
#     name='PyKeras_V14_20231025_1_multiclass',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20231025_1_multiclass.h5',
#     inputs=['deltam', 'track_pt', 'track_abs_eta_',
#             'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
#             'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
#             'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
#             'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
#             'track_log10_dxy_', 'track_log10_dz_',
#             'track_log10_dxyPU_', 'track_log10_dzPU_',
#             'track_log10_dxyError_', 'track_log10_dzError_',
#             'track_pfAbsIso',
#             'track_drminTrack5',
#             'track_drmin2ndTrack5',
#             'track_drminJet20',
#             'track_drminJet30',
#             'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
#             'track_abs_dphiMet_',
#             'met_pt'
#             ],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     savemaxscoreinfo=['Signal'],
#     parameters={'deltam': [0.3, 0.6, 1.0]},
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
#                    'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
#                    },
#     eventlevelintputs=['met_pt'],
#     qualitymask=True
# )

# my_model_V14_20231121 = MyModel(
#     name='PyKeras_V14_20231121_multiclass',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20231121_multiclass.h5',
#     inputs=['deltam', 'track_pt', 'track_abs_eta_',
#             'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
#             'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
#             'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
#             'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
#             'track_log10_dxy_', 'track_log10_dz_',
#             'track_log10_dxyPU_', 'track_log10_dzPU_',
#             'track_log10_dxyError_', 'track_log10_dzError_',
#             'track_pfAbsIso',
#             'track_drminTrack5',
#             'track_drmin2ndTrack5',
#             'track_drminJet20',
#             'track_drminJet30',
#             'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
#             'track_abs_dphiMet_',
#             'met_pt'
#             ],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     savemaxscoreinfo=['Signal'],
#     parameters={'deltam': [0.3, 0.6, 1.0]},
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
#                    'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
#                    },
#     eventlevelintputs=['met_pt'],
#     qualitymask=True
# )

# my_model_V14_20240228_1 = MyModel(
#     name='PyKeras_V14_20240228_1_multiclass',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240228_1_multiclass.h5',
#     inputs=['deltam', 'track_pt', 'track_abs_eta_',
#             'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
#             'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
#             'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
#             'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
#             'track_log10_dxy_', 'track_log10_dz_',
#             'track_log10_dxyPU_', 'track_log10_dzPU_',
#             'track_log10_dxyError_', 'track_log10_dzError_',
#             'track_pfAbsIso',
#             'track_drminTrack5',
#             'track_drmin2ndTrack5',
#             'track_drminJet20',
#             'track_drminJet30',
#             'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
#             'track_abs_dphiMet_',
#             'met_pt'
#             ],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     savemaxscoreinfo=['Signal'],
#     parameters={'deltam': [0.3, 0.6, 1.0]},
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
#                    'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
#                    },
#     eventlevelintputs=['met_pt'],
#     qualitymask=True
# )

# my_model_V14_20240319_4 = MyModel(
#     name='PyKeras_V14_20240319_4_multiclass',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240319_4_multiclass.h5',
#     custom_objects={'focal_loss_fn': None},
#     inputs=['deltam', 'track_pt', 'track_abs_eta_',
#             'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
#             'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
#             'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
#             'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
#             'track_log10_dxy_', 'track_log10_dz_',
#             'track_log10_dxyPU_', 'track_log10_dzPU_',
#             'track_log10_dxyError_', 'track_log10_dzError_',
#             'track_pfAbsIso',
#             'track_drminTrack5',
#             'track_drmin2ndTrack5',
#             'track_drminJet20',
#             'track_drminJet30',
#             'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
#             'track_abs_dphiMet_',
#             'met_pt'
#             ],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     savemaxscoreinfo=['Signal'],
#     parameters={'deltam': [0.3, 0.6, 1.0]},
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
#                    'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
#                    },
#     eventlevelintputs=['met_pt'],
#     qualitymask=True
# )

# my_model_V14_20240327 = MyModel(
#     name='PyKeras_V14_20240327_multiclass',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240327_multiclass.h5',
#     custom_objects={'focal_loss_fn': None},
#     inputs=['deltam', 'track_pt', 'track_abs_eta_',
#             'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
#             'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
#             'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
#             'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
#             'track_log10_dxy_', 'track_log10_dz_',
#             'track_log10_dxyPU_', 'track_log10_dzPU_',
#             'track_log10_dxyError_', 'track_log10_dzError_',
#             'track_pfAbsIso',
#             'track_drminTrack5',
#             'track_drmin2ndTrack5',
#             'track_drminJet20',
#             'track_drminJet30',
#             'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
#             'track_abs_dphiMet_',
#             'met_pt'
#             ],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     savemaxscoreinfo=['Signal'],
#     parameters={'deltam': [0.3, 0.6, 1.0]},
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
#                    'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
#                    },
#     eventlevelintputs=['met_pt'],
#     qualitymask=True
# )

# my_model_V14_20240328 = MyModel(
#     name='PyKeras_V14_20240328_multiclass',
#     h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240328_multiclass.h5',
#     custom_objects={'focal_loss_fn': None},
#     inputs=['track_pt', 'track_abs_eta_',
#             'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
#             'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
#             'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
#             'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
#             'track_log10_dxy_', 'track_log10_dz_',
#             'track_log10_dxyPU_', 'track_log10_dzPU_',
#             'track_log10_dxyError_', 'track_log10_dzError_',
#             'track_pfAbsIso',
#             'track_drminTrack5',
#             'track_drmin2ndTrack5',
#             'track_drminJet20',
#             'track_drminJet30',
#             'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
#             'track_abs_dphiMet_',
#             'met_pt'
#             ],
#     outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
#     savemaxscoreinfo=['Signal'],
#     specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
#                    'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
#                    'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
#                    },
#     eventlevelintputs=['met_pt'],
#     qualitymask='dz10'
# )

my_model_V14_20240323 = MyModel(
    name='PyKeras_V14_20240323_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240323_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask=True
)

my_model_V14_20240403 = MyModel(
    name='PyKeras_V14_20240403_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240403_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='20240328>0.1',
    friendfilesvariables={'NTuplesV14p5': ['track_PyKeras_V14_20240328_multiclass_Signal_']}
)

my_model_V14_20240403_1 = MyModel(
    name='PyKeras_V14_20240403_1_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240403_1_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='20240328>0.2',
    friendfilesvariables={'NTuplesV14p5': ['track_PyKeras_V14_20240328_multiclass_Signal_']}
)

my_model_V14_20240403_2 = MyModel(
    name='PyKeras_V14_20240403_2_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240403_2_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='20240328>0.5',
    friendfilesvariables={'NTuplesV14p5': ['track_PyKeras_V14_20240328_multiclass_Signal_']}
)

my_model_V14_20240405_1 = MyModel(
    name='PyKeras_V14_20240405_1_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240405_1_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='20240328>0.2',
    friendfilesvariables={'NTuplesV14p5': ['track_PyKeras_V14_20240328_multiclass_Signal_']}
)

my_model_V14_20240408 = MyModel(
    name='PyKeras_V14_20240408_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240408_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask=True
)

my_model_V14_20240408_1 = MyModel(
    name='PyKeras_V14_20240408_1_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240408_1_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask=True
)

my_model_V14_20240409 = MyModel(
    name='PyKeras_V14_20240409_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240409_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='dz0p5'
)

my_model_V14_20240409_1 = MyModel(
    name='PyKeras_V14_20240409_1_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240409_1_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='dz0p1'
)

my_model_V14_20240411 = MyModel(
    name='PyKeras_V14_20240411_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240411_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='pt0p5'
)

my_model_V14_20240411_1 = MyModel(
    name='PyKeras_V14_20240411_1_multiclass',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V14/weights/TrainedModel_PyKeras_V14_20240411_1_multiclass.h5',
    custom_objects={'focal_loss_fn': None},
    inputs=['deltam', 'track_pt', 'track_abs_eta_',
            'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_',
            'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_',
            'track_log10_IPsigAssPV_', 'track_log10_IPxyAssPV_', 'track_log10_IPzAssPV_',
            'track_log10_IPsigPUAssPV_', 'track_log10_IPxyPUAssPV_', 'track_log10_IPzPUAssPV_',
            'track_log10_dxy_', 'track_log10_dz_',
            'track_log10_dxyPU_', 'track_log10_dzPU_',
            'track_log10_dxyError_', 'track_log10_dzError_',
            'track_pfAbsIso',
            'track_drminTrack5',
            'track_drmin2ndTrack5',
            'track_drminJet20',
            'track_drminJet30',
            'track_abs_detaLeadingJet_', 'track_abs_dphiLeadingJet_',
            'track_abs_dphiMet_',
            'met_pt'
            ],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    savemaxscoreinfo=['Signal'],
    parameters={'deltam': [0.3, 0.6, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)],
                   'track_abs_detaLeadingJet_': ['track_detaLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiLeadingJet_': ['track_dphiLeadingJet', lambda x: abs(x)],
                   'track_abs_dphiMet_': ['track_dphiMet', lambda x: abs(x)]
                   },
    eventlevelintputs=['met_pt'],
    qualitymask='pt1'
)

# TODO: update list
# TODO: add fromtruetau to savemaxscoreinfo?
# TODO: add more deltams?
my_models = [
    my_model_V14_20240411,
    my_model_V14_20240411_1,
]

isjob = False
options = MyVarParsing(sys.argv)
if len(options.inputFiles) > 0:
    list_of_datasets = options.inputFiles
    isjob = True

n_datasets = len(list_of_datasets)
for idataset, dataset in enumerate(list_of_datasets):

    print("----")
    print("adding friends to", dataset)
    print("----")
    # if 'pnfs' in dataset:
    #     inpath = dataset + 'crab_NTuple_*.root'
    #     infiles = glob(inpath)
    #     fout_folder = inpath.split('16UL_preAPV/')[-1].split('/crab_')[0].replace('/', '_') + '/'
    # else:
    #     inpath = dataset
    #     infiles = glob(inpath)
    #     fout_folder = ''

    if dataset.endswith('.root'):
        inpath = dataset.rsplit('/', 1)[0] + '/'
    else:
        inpath = dataset
        dataset += 'crab_NTuple_*.root'

    infiles = glob(dataset)
    n_infiles = len(infiles)

    if not isjob:

        if '/NTuples/' in inpath:
            fout_folder = inpath.split('/NTuples/')[1].rsplit('/', 1)[0] + '/'
        else:
            fout_folder = ''

        if not os.path.exists(outpath + fout_folder):
            # Create the folder
            os.makedirs(outpath + fout_folder)
            print("Folder created successfully.")
        else:
            print("Folder already exists.")

    for ifile, afile in enumerate(infiles):
        print('[' + str(idataset + 1) + '/' + str(n_datasets) + '] [' + str(ifile + 1) + '/' + str(n_infiles) + ']')
        # if isjob:
        #     makeFriendTree(fin=afile, outpath='')
        # else:
        #     makeFriendTree(fin=afile, outpath=outpath + fout_folder)
        if isjob:
            makeFriendTree_keras(fin=afile, outpath='', my_models=my_models)
        else:
            makeFriendTree_keras(fin=afile, outpath=outpath + fout_folder, my_models=my_models)
