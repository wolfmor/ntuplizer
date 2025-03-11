#!/bin/bash
# source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos9-gcc11-opt/setup.sh
python /afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_6_34/src/SoftDisplacedPion/ntuplizer/addFriendTree_lite.py "$@"
