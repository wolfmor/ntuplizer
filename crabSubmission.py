#! /usr/bin/env python

import os, sys, shutil
from glob import glob
from random import shuffle 
from time import sleep
import datetime
now = datetime.datetime.now()


# TODO: add low-mass DY?
datasets_MC = {}

datasets_MC['era16_UL_APV'] = [
    # '/ZJetsToNuNu_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM'
    # '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    #
    # '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    # '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    # '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    # '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',

    # '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8_ext1-v1/AODSIM',

    '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    # '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',

    '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    # '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v3/AODSIM',
    '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',

    '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',
    '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v2/AODSIM',

    # '/DYJetsToLL_M-50_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECOAPV-106X_mcRun2_asymptotic_preVFP_v8-v1/AODSIM',
]

datasets_MC['era16_UL'] = [
    # '/QCD_HT100to200_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/QCD_HT200to300_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/QCD_HT300to500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/QCD_HT500to700_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    # '/QCD_HT700to1000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/QCD_HT1000to1500_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/QCD_HT1500to2000_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/QCD_HT2000toInf_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',

    # '/ZJetsToNuNu_Zpt-100to200_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',

    # '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    # '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    # '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    # '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',

    # '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v1/AODSIM',
    # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13_ext1-v1/AODSIM',

    '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',

    # '/TTZToLLNuNu_M-10_TuneCP5_13TeV-amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    # '/TTZToQQ_TuneCP5_13TeV_amcatnlo-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    # '/TTWJetsToLNu_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    # '/TTWJetsToQQ_TuneCP5_13TeV-amcatnloFXFX-madspin-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',

    '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v3/AODSIM',
    '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',

    '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',
    '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',

    # '/DYJetsToLL_M-50_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v2/AODSIM',

    # '/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL16RECO-106X_mcRun2_asymptotic_v13-v1/AODSIM',
]

datasets_MC['era17_UL'] = [
    # '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    
    # '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',

    # '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v3/AODSIM',

    # # '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v3/AODSIM',
    # # '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v3/AODSIM',
    # # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6_ext2-v3/AODSIM',

    # '/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',

    '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',

    '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',
    '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v1/AODSIM',

    '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
    '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL17RECO-106X_mc2017_realistic_v6-v2/AODSIM',
]

datasets_MC['era18_UL'] = [
    # '/ZJetsToNuNu_Zpt-200toInf_BPSFilter_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    
    # '/ZJetsToNuNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ZJetsToNuNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ZJetsToNuNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ZJetsToNuNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ZJetsToNuNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ZJetsToNuNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',

    # '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',

    # # '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    # # '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext1 has more events
    # '/WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v3/AODSIM',
    # # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',
    # # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',  # ext2 has more events
    # # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext1-v1/AODSIM',  # ext2 has more events
    # '/WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1_ext2-v3/AODSIM',

    # '/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-200to400_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-600to800_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-800to1200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-1200to2500_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',

    # '/WW_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    # '/WZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    # '/ZZ_TuneCP5_13TeV-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

    # '/ST_t-channel_antitop_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    '/ST_t-channel_top_4f_InclusiveDecays_TuneCP5_13TeV-powheg-madspin-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ST_tW_antitop_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',
    # '/ST_tW_top_5f_inclusiveDecays_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v1/AODSIM',

    # '/TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    # '/TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',
    '/TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18RECO-106X_upgrade2018_realistic_v11_L1v1-v2/AODSIM',

]


datasets_MET = {}

datasets_MET['era16_UL_APV'] = [
    '/MET/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/MET/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD',
]

datasets_MET['era16_UL'] = [
    '/MET/Run2016F-21Feb2020_UL2016-v1/AOD',
    '/MET/Run2016G-21Feb2020_UL2016-v1/AOD',
    '/MET/Run2016H-21Feb2020_UL2016-v2/AOD',
]

datasets_MET['era17_UL'] = [
    '/MET/Run2017B-09Aug2019_UL2017_rsb-v1/AOD',
    '/MET/Run2017C-09Aug2019_UL2017_rsb-v1/AOD',
    '/MET/Run2017D-09Aug2019_UL2017_rsb-v1/AOD',
    '/MET/Run2017E-09Aug2019_UL2017_rsb-v1/AOD',
    '/MET/Run2017F-09Aug2019_UL2017_rsb-v1/AOD',
]

datasets_MET['era18_UL'] = [
    # '/MET/Run2018A-15Feb2022_UL2018-v1/AOD',
    '/MET/Run2018B-15Feb2022_UL2018-v1/AOD',
    # '/MET/Run2018C-15Feb2022_UL2018-v1/AOD',
    # '/MET/Run2018D-15Feb2022_UL2018-v1/AOD',
]


datasets_SingleMuon = {}

datasets_SingleMuon['era16_UL_APV'] = [
    # '/SingleMuon/Run2016B-21Feb2020_ver1_UL2016_HIPM-v1/AOD',
    # '/SingleMuon/Run2016B-21Feb2020_ver2_UL2016_HIPM-v1/AOD',
    # '/SingleMuon/Run2016C-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/SingleMuon/Run2016D-21Feb2020_UL2016_HIPM-v1/AOD',
    # '/SingleMuon/Run2016E-21Feb2020_UL2016_HIPM-v1/AOD',
    '/SingleMuon/Run2016F-21Feb2020_UL2016_HIPM-v1/AOD',
]

datasets_SingleMuon['era16_UL'] = [
    # '/SingleMuon/Run2016F-21Feb2020_UL2016-v1/AOD',
    '/SingleMuon/Run2016G-21Feb2020_UL2016-v1/AOD',
    # '/SingleMuon/Run2016H-21Feb2020_UL2016-v1/AOD',
]

datasets_SingleMuon['era17_UL'] = [
    # '/SingleMuon/Run2017B-15Feb2022_UL2017-v1/AOD',
    # '/SingleMuon/Run2017C-15Feb2022_UL2017-v1/AOD',
    # '/SingleMuon/Run2017D-15Feb2022_UL2017-v1/AOD',
    # '/SingleMuon/Run2017E-15Feb2022_UL2017-v1/AOD',
    '/SingleMuon/Run2017F-15Feb2022_UL2017-v1/AOD',
]

datasets_SingleMuon['era18_UL'] = [
    # '/SingleMuon/Run2018A-15Feb2022_UL2018-v1/AOD',
    # '/SingleMuon/Run2018B-15Feb2022_UL2018-v1/AOD',
    # '/SingleMuon/Run2018C-15Feb2022_UL2018-v1/AOD',
    '/SingleMuon/Run2018D-15Feb2022_UL2018-v1/AOD',
]


era = 'era18_UL'  # era16_UL_APV era16_UL era17_UL era18_UL

test = False
limitNJobs = False
moritz = True

data = True
sameSign = False
cleanDY = False
skipSVs = True

if data:
    datasets = datasets_MET[era]  # datasets_SingleMuon[era]  # # TODO: adapt MET or SingleMuon
else:
    datasets = datasets_MC[era]

ntuplizer_tag = "crab,genmatchalltracks," + era
if data:
    ntuplizer_tag = ntuplizer_tag.replace('genmatchalltracks', 'data')
if sameSign:
    ntuplizer_tag += ",samesign"
if cleanDY:
    ntuplizer_tag += ",cleanleptons"
if skipSVs:
    ntuplizer_tag += ",skipSVs"

if moritz:
    crabprojectpath = "/nfs/dust/cms/user/wolfmor/"
    inputpath = "/afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_6_34/src/SoftDisplacedPion/ntuplizer/"
    outputpath = "/store/user/mowolf/NTuples/NTuplesV15" + "/" if skipSVs else "_withSVs/"  # TODO: adapt ntuple version
    outputpath += {'era16_UL_APV': '16UL_preAPV/', 'era16_UL': '16UL/', 'era17_UL': '17UL/', 'era18_UL': '18UL/'}[era]
else:
    crabprojectpath = ""
    inputpath = "/nfs/dust/cms/user/tewsalex/CMSSW_10_2_18/src/ntuplizer/"
    outputpath = "/store/user/altews/NTuples/NTuplesV12/16UL_preAPV/"


def submitToCrab():
    nJobs = 0
    for ijob, dataset in enumerate(datasets):
        # if data:
        #     task_tag = ((dataset.split('/')[1]).split("/")[0])+"_Run"+(dataset.split('/Run')[1]).split("/")[0]
        # else:
        #     task_tag = ((dataset.split('/')[1]).split("_")[0])+((dataset.split('/')[1]).split("_")[1])+(dataset.split('/')[2]).split("-")[0]
        if data:
            task_tag = dataset.split('/')[2]
        else:
            task_tag = dataset.split('/')[2].split('-')[0]

        task_tag += dataset.replace('/', '_')[:50]

        if sameSign:
            task_tag += "_sameSign"
        if cleanDY:
            task_tag += "_cleaned"
        if skipSVs:
            task_tag += "_skipSVs"

        task_tag = task_tag[:100]  # make sure the tag is not longer than 100 characters

        print("use tag:", task_tag)

        if os.path.exists(crabprojectpath + "crab_projects/crab_" + task_tag):
            shutil.rmtree(crabprojectpath + "crab_projects/crab_" + task_tag)
        if os.path.exists("crabSubmissionScript_"+str(ijob)+".py"): os.remove("crabSubmissionScript_"+str(ijob)+".py")

        f = open("crabSubmissionScript_"+str(ijob)+".py", "w")

        f.write("from CRABClient.UserUtilities import config\n")
        f.write("import datetime\n")
        f.write("import glob\n")
        f.write("now = datetime.datetime.now()\n")
        f.write("config = config()\n")
        f.write("config.General.requestName = \"" + task_tag + "\"\n")
        f.write("config.General.workArea = \"" + crabprojectpath + "crab_projects\"\n")
        f.write("config.General.transferOutputs = True\n")
        f.write("config.General.transferLogs = True\n")
        f.write("config.JobType.pluginName = \"Analysis\"\n")
        if skipSVs:
            f.write("config.JobType.psetName = \"construct_secondary_vertices_dummy_cfg.py\"\n")
        else:
            if sameSign:
                f.write("config.JobType.psetName = \"construct_secondary_vertices_ss_cfg.py\"\n")
            else:
                f.write("config.JobType.psetName = \"construct_secondary_vertices_cfg.py\"\n")
        f.write("config.JobType.scriptArgs = [\"tag=" + ntuplizer_tag + "\" , \"dataset=" + dataset + "\"]\n")
        f.write("config.JobType.inputFiles = [\"" + inputpath + "plantTrees.py\",\n")
        f.write("\"" + inputpath + "commons.py\",\n")
        if skipSVs:
            f.write("\"" + inputpath + "FrameworkJobReport.xml\",\n")
        f.write("\"/nfs/dust/cms/user/wolfmor/NTupleStuff/\"]\n")
        f.write("config.JobType.allowUndistributedCMSSW = True\n")
        if skipSVs:
            f.write("config.JobType.scriptExe = \"run_plantTrees_python3.sh\"\n")  # TODO: _python3?
        else:
            f.write("config.JobType.scriptExe = \"allinonejob.sh\"\n")
        if sameSign:
            if data: f.write("config.JobType.outputFiles = [\"crab_ss_NTuple.root\", \"crab_ss_NTuple.json\"]\n")
            else: f.write("config.JobType.outputFiles = [\"crab_ss_NTuple.root\"]\n")
        elif data: f.write("config.JobType.outputFiles = [\"crab_NTuple.root\", \"crab_NTuple.json\"]\n")
        else: f.write("config.JobType.outputFiles = [\"crab_NTuple.root\"]\n")
        # f.write("config.Data.partialDataset = True\n")  # TODO: run over incomplete dataset?
        f.write("config.Data.inputDataset = \"" + dataset + "\"\n")
        f.write("config.Data.inputDBS = \"global\"\n")
        f.write("config.Data.splitting = \"FileBased\"\n")
        f.write("config.Data.unitsPerJob = 1\n")  # TODO: change to 2 to decrease number of jobs?
        outputpathappend = ''
        if sameSign:
            outputpathappend += "samesign/"
        if cleanDY:
            outputpathappend += "cleaned/"
        f.write("config.Data.outLFNDirBase = \"" + outputpath + outputpathappend + "\"\n")
        if limitNJobs:
            f.write("NJOBS = 1 \n")
            f.write("config.Data.totalUnits = config.Data.unitsPerJob * NJOBS \n")
        f.write("config.Data.publication = False \n")
        if data:
            f.write("config.Data.outputDatasetTag = \"" + dataset.split('/')[2].split('-')[0] + "\" + \"_\"  + now.strftime(\"%Y_%m_%d\")\n")
        else:
            f.write("config.Data.outputDatasetTag = now.strftime(\"%Y_%m_%d\")\n")
        f.write("config.Site.storageSite = \"T2_DE_DESY\"\n")
        # f.write("config.Site.ignoreGlobalBlacklist = True\n")  # TODO: remove
        # f.write("config.Site.blacklist = [\"T2_IN_TIFR\", \"T2_BR_SPRACE\", \"T1_RU_JINR\"]\n")  # TODO: specify blacklist?
        f.write("config.Site.whitelist = [\"T2_DE_DESY\"]\n")  # TODO: specify whitelist?
        f.write("config.Data.ignoreLocality = True\n")  # TODO: to only run on whitelist # TODO: remove!
        f.write("config.JobType.maxJobRuntimeMin = 2750\n")  # TODO: ?
        f.write("config.JobType.maxMemoryMB = 5000\n")  # TODO: 5000 ? 2500 ?

        f.close()

        f = open("crabSubmissionScript_" +str(ijob)+ ".py", "r")
        command = "crab submit -c crabSubmissionScript_" +str(ijob)+ ".py & "  # TODO:  --proxy /afs/desy.de/user/w/wolfmor/k5-ca-proxy.pem ?
        print("command", command)
        if not test: os.system(command)

        sleep(120)
        nJobs += 1

    if not test: print("submitted " + str(nJobs) + " tasks")


submitToCrab()
