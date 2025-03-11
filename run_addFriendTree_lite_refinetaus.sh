#!/bin/bash
# source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
cp /afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_2_18/src/SoftDisplacedTrack/training/NNmulticlass_V15/weights/TrainedModel_PyKeras_V15_20240711_multiclass.h5 .
cp /data/dust/user/wolfmor/Refinement/SoftTracks/model_refinement_regression_20241221_2_cpu.pt .
cp /data/dust/user/wolfmor/Refinement/SoftTracks/model_refinement_regression_20250102_cpu.pt .
cp /data/dust/user/wolfmor/Refinement/SoftTracks/model_refinement_regression_20250102_1_cpu.pt .
cp /data/dust/user/wolfmor/Refinement/SoftTracks/model_refinement_regression_20250102_2_cpu.pt .
cp /data/dust/user/wolfmor/Refinement/SoftTracks/model_refinement_regression_20250102_3_cpu.pt .
python /afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_6_34/src/SoftDisplacedPion/ntuplizer/addFriendTree_lite_refinetaus.py "$@"
