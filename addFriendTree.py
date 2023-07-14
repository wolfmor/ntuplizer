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

gROOT.SetBatch(True)
gStyle.SetOptStat(111111)


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
        my_model=None):

    # my_model.model.summary()
    # print(tf.config.list_physical_devices('GPU'))

    events = uproot.open(fin)['tEvent']

    # TODO: handle mix of event-level and track-level variables
    variables = events.arrays([inp for inp in my_model.inputs if inp not in my_model.parameters.keys() and inp not in my_model.specialinputs.keys()], library='np')
    specials = events.arrays([my_model.specialinputs[s][0] for s in my_model.specialinputs], library='np')

    for variable in my_model.specialinputs:
        variables[variable] = my_model.specialinputs[variable][1](specials[my_model.specialinputs[variable][0]])

    variables.update(events.arrays(['random', 'track_random', 'track_isSignalTrack', 'track_isSusyTrack', 'track_susyTrackPdgId', 'track_genMatchMotherIsTheTau'], library='np'))

    outname = fin.split('/')[-1].split('.root')[0]
    foutname = outpath + outname + '_friend.root'

    parameter_points = [dict(zip(my_model.parameters.keys(), values)) for values in product(*my_model.parameters.values())]
    parameter_points_labels = ['_'.join([p + str(pp[p]).replace('.', 'p') for p in pp]) for pp in parameter_points]

    with uproot.recreate(foutname) as fout:

        # TODO: only one counter for all branches
        branches = {'random_sanitycheck': 'float64', 'track_random_sanitycheck': 'var * float64'}
        for label in parameter_points_labels:
            for out in my_model.outputs:

                branches['track_' + my_model.name + '_' + out + '_' + label] = 'var * float64'
                branches['track_' + my_model.name + '_' + out + '_' + label + '_isMaxscore'] = 'var * int32'

                branches['maxscore_' + my_model.name + '_' + out + '_' + label] = 'float64'
                branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_isSignalTrack'] = 'int32'
                branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_isSusyTrack'] = 'int32'
                branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_susyTrackPdgId'] = 'int32'
                branches['maxscore_' + my_model.name + '_' + out + '_' + label + '_genMatchMotherIsTheTau'] = 'int32'

        fout.mktree('tFriend', branches)
        # fout['tFriend'].show()

        for event in range(len(variables['track_pt'])):

            n_track = np.array([len(variables['track_pt'][event])], dtype=np.int32)

            # print(branches)

            for pp, ppl in zip(parameter_points, parameter_points_labels):

                inputs = np.stack([np.full_like(variables['track_pt'][event], pp[inp]) if inp in pp else variables[inp][event] for inp in my_model.inputs], axis=1)
                preds = ak.Array([my_model.model.predict(inputs)])

                # print()
                # print(preds[:, :, 0])
                # print(ak.max(preds[:, :, 0]))
                # print(ak.argmax(preds[:, :, 0]))
                # print(ak.Array([[int(i == ak.argmax(preds[:, :, 0])) for i in range(len(preds[:, :, 0][0]))]]))
                # print()
                # print(preds[:, 1])
                # print(ak.max(preds[:, 1]))
                # print(preds[:, 1] == ak.max(preds[:, 1]))

                # sys.exit(0)

                branches['random_sanitycheck'] = ak.Array([variables['random'][event]])
                branches['track_random_sanitycheck'] = ak.Array([variables['track_random'][event]])
                branches['ntrack_random_sanitycheck'] = n_track
                for iout, out in enumerate(my_model.outputs):

                    highscoreidx = ak.argmax(preds[:, :, iout])

                    branches['track_' + my_model.name + '_' + out + '_' + ppl] = preds[:, :, iout]
                    branches['ntrack_' + my_model.name + '_' + out + '_' + ppl] = n_track

                    branches['track_' + my_model.name + '_' + out + '_' + ppl + '_isMaxscore'] = ak.Array([[int(i == highscoreidx) for i in range(len(preds[:, :, iout][0]))]])
                    branches['ntrack_' + my_model.name + '_' + out + '_' + ppl + '_isMaxscore'] = n_track

                    branches['maxscore_' + my_model.name + '_' + out + '_' + ppl] = ak.Array([ak.max(preds[:, :, iout])])
                    branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_isSignalTrack'] = ak.Array([variables['track_isSignalTrack'][event][highscoreidx]])
                    branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_isSusyTrack'] = ak.Array([variables['track_isSusyTrack'][event][highscoreidx]])
                    branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_susyTrackPdgId'] = ak.Array([variables['track_susyTrackPdgId'][event][highscoreidx]])
                    branches['maxscore_' + my_model.name + '_' + out + '_' + ppl + '_genMatchMotherIsTheTau'] = ak.Array([variables['track_genMatchMotherIsTheTau'][event][highscoreidx]])


            # print(branches)
            # print(ak.zip(branches))
            # print(ak.zip(branches).type)

            # if event == 0:
            #     fout.mktree('tFriend', ak.zip(branches).type)
            #     fout['tFriend'].show()

            fout['tFriend'].extend(branches)


    print("added a tFriend to", fin, foutname)


class MyModel:
    def __init__(self, name, h5file, inputs, outputs, parameters=None, specialinputs=None):

        self.name = name
        self.model = load_model(h5file)
        self.inputs = inputs
        self.outputs = outputs
        if parameters is None:
            self.parameters = {}
        else:
            self.parameters = parameters
        if specialinputs is None:
            self.specialinputs = {}
        else:
            self.specialinputs = specialinputs


### give the path to the folder that contains the ntuples for datasets on pnfs 
### or give a path to files with wildcards for signal files on dust / local
list_of_datasets = [
    # '/pnfs/desy.de/cms/tier2/store/user/altews/NTuples/NTuplesV12/16UL_preAPV/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/2023_04_28/230428_123855/0000/',
    # '/nfs/dust/cms/user/tewsalex/rootfiles/ntuple_V12/step3_higgsino_RunIISpring21UL16FS_susyall_mChipm*GeV_dm*GeV_Chi20ctau5MM_part*of*_NTuple.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV11/16UL/ZJetsToNuNu_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/2023_02_22/230222_111558/0000/crab_NTuple_1.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p368GeV_part1of25_NTuple_noSVs_job385.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm0p968GeV_part20of25_NTuple_noSVs_job486.root',
    # '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalV4_16/step3_higgsinoDm0eDmpm_RunIISpring21UL16FS_susyall_mChipm115GeV_dm*GeV_part*_NTuple_noSVs_job*.root',
    '/pnfs/desy.de/cms/tier2/store/user/mowolf/NTuples/NTuplesV12_noSVs/16UL/SignalStopV4_16/step3_higgsino_RunIISpring21UL16FS_stopstop_500GeV_mChipm400GeV_dm1p0GeV_part16of25_NTuple_noSVs_job155.root',
]
# outpath = '/nfs/dust/cms/user/tewsalex/rootfiles/friendTrees_V12/'
outpath = '/nfs/dust/cms/user/wolfmor/FrieNdTuples/TestArthur/'
# outpath = '/nfs/dust/cms/user/wolfmor/FrieNdTuples/Test/'


my_model_V11_20230630 = MyModel(
    name='PyKeras_V11_20230630_multiclass_puretracklevel',
    h5file='/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V11_puretracklevel/weights/TrainedModel_PyKeras_V11_20230702_multiclass_puretracklevel.h5',
    inputs=['deltam', 'track_pt', 'track_abs_eta_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_', 'track_log10_IPsigPU_', 'track_log10_IPxyPU_', 'track_log10_IPzPU_', 'track_log10_dxy_', 'track_log10_dz_', 'track_log10_dxyPU_', 'track_log10_dzPU_', 'track_log10_dxyError_', 'track_log10_dzError_', 'track_neHadAbsIso0', 'track_tkAbsIso0', 'track_drminTrack10', 'track_drminJet20', 'track_drminJet30'],
    outputs=['Signal', 'Background_nogenmatch', 'Background_prompt', 'Background_secondary', 'Background_fromtruetau'],
    parameters={'deltam': [0.3, 1.0]},
    specialinputs={'track_abs_eta_': ['track_eta', lambda x: abs(x)]}
)

my_model_Arthur = MyModel(
    name='Arthur_ckpt_24fast',
    h5file='/nfs/dust/cms/user/tanikulo/Arthur_ckpt_24fast_Tracks_3dm_t70v15t15_D64_BN_D64_BN_D64_BN_D64_BN_D3_ep10_b40_lr0p00001_weights.h5',
    inputs=['deltamFile', 'track_pt', 'track_dphiMet',  'track_eta', 'track_log10_IPsigPU_', 'track_log10_IPsig_', 'track_log10_IPxy_', 'track_log10_IPz_','track_drminBjetMedium30','track_drminJet30', 'track_log10_dxy_', 'track_log10_dxyError_', 'track_log10_dz_', 'track_log10_dzError_', 'track_detaLeadingJet'],
    outputs=['predicted_Background_PU', 'predicted_Background_nonPU', 'predicted_Signal'],
    parameters={'deltamFile': [0.3, 1.0]}  # TODO: add more dMs
)

for dataset in list_of_datasets:

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
        print('[' + str(ifile + 1) + '/' + str(n_infiles) + ']')
        # makeFriendTree(fin=afile, outpath=outpath + fout_folder)
        makeFriendTree_keras(fin=afile, outpath=outpath + fout_folder, my_model=my_model_Arthur)
