import os
from glob import glob
from mydatasets import mydatasets_V15

nameout = '/nfs/dust/cms/user/wolfmor/gridcontrol/workdir/friends_FRIENDSVERSION_ERA.conf'

datasets = mydatasets_V15
friends_version = 'V15p15'
eras = [
    'era16_UL_APV',
    'era16_UL',
    'era17_UL',
    'era18_UL',
]

# TODO: restrict backgrounds?
include = [
    # 'DataMuon',
    # 'WJetsToLNu_HT-100To200',
    # 'WW',
    # 'WZ',
    # 'ZZ',
    # 'ST_t-channel_antitop',
    # 'ST_t-channel_top',
    # 'ST_tW_antitop',
    # 'ST_tW_top',
    # 'TTJets_DiLept',
    # 'TTJets_SingleLeptFromT',
    # 'TTJets_SingleLeptFromTbar',
    'DYJetsToLL_M-50_HT-100to200',
]
exclude = [
    # 'DYJetsToLL',
    # 'DataMuon',
]

nameout = nameout.replace('.conf', '_onlycleaneddymc_HT-100to200.conf')


# TODO: change "files per job" to something larger?
S = '''
[global]
task 				= CMSSW 				; Select grid-control CMSSW task
backend 			= local					; host = local machine, local = condor batch

[jobs]
;jobs				= 1
in queue     	  	= 1000
wall time         	= 2:59
memory            	= 2048
max retry 	  		= 0
defect tries		= 0

[CMSSW]
scram project 	  	= CMSSW CMSSW_13_0_3  ; CMSSW CMSSW_10_2_18
scram arch	  		= el9_amd64_gcc11  ; slc7_amd64_gcc700
se runtime        	= True 					; Large project areas need to be transferred via the SE

epilog executable 	= /afs/desy.de/user/w/wolfmor/cmssw/CMSSW_10_6_34/src/SoftDisplacedPion/ntuplizer/run_addFriendTree.sh

epilog arguments  	= inputFiles="@FILE_NAMES@"

dataset provider  	= scan

dataset				=
    DATASET

dataset splitter 	= FileBoundarySplitter
files per job    	= 2

[storage]
se path           	= srm://dcache-se-cms.desy.de:8443//pnfs/desy.de/cms/tier2/store/user/mowolf/FrieNdTuples/NTuplesFRIENDSVERSION/ERA
se output files   	= *.root	; Name of the CMSSW output file
se output pattern 	= @NICK@/@XBASE@.root  ; @XBASE@_job@MY_JOBID@.root

[backend]
proxy = VomsProxy

[task]
depends += glite

'''

for era in eras:

    with open(nameout.replace('FRIENDSVERSION', friends_version).replace('ERA', era), 'w') as f:

        s = S

        s = s.replace('FRIENDSVERSION', friends_version)
        s = s.replace('ERA', era.replace('era', '').replace('_', '').replace('APV', '_preAPV'))

        datasetlist = []
        for d in sorted(datasets[era].keys()):

            if len(include) > 0 and not any([i in d for i in include]): continue
            if len(exclude) > 0 and any([e in d for e in exclude]): continue

            paths = []
            if type(datasets[era][d]) == list:
                paths += datasets[era][d]
            else:
                paths = [datasets[era][d]]

            for p in paths:
                if len(p) > 0:
                    if '/*/' in p:
                        datasetlist += [p.replace('/*/', '/' + folder + '/') for folder in os.listdir('/'.join(p.split('/')[:-2]))]
                    else:
                        datasetlist += [p]

        s = s.replace('DATASET', '\n    '.join(['/'.join(ds.split('/')[11:-1]) + ' : ' + ds for ds in datasetlist]))

        f.write(s)
