# ntuplizer
A simple FWLite ntuplizer taking AOD files as inputs, running on the NAF.

Various event selection cuts are applied and individual tracks are stored for further analysis. The considered signal process is chargino/chargino or chargino/neutralino pair production. Depends on helper functions in `commons.py`.
## How to use
Make sure CMSSW is set up (tested with CMSSW_10_2_18).
```
git clone https://github.com/wolfmor/ntuplizer.git
cd ntuplizer/
python plantTrees.py inputFiles="file1, file2, ..." tag="tag1 tag2 ..."
```
By default, files are opened via xrootd, so be sure to have a valid grid proxy. To use local files (e.g. on pnfs) include the tag `local`.
## Tags
Tag | Meaning
------------ | -------------
`local` | don't open file via xrootd
`redirinfn` | use infn redirector for xrootd
`redirfnal` | use fnal redirector for xrootd
`data` | check trigger flags and runnum/lumisec, save json file
`geninfoZ` | save GEN info for DY process
`geninfoW` | save GEN info for W boson production process
`cleanleptons` | perform DY cleaning
`signal` | signal GEN variables and FastSim correction for MET/JEC
`noleptonveto` | don't veto events with leptons
`genmatchtracks` | try to find GEN match to every 10. track
`genmatchalltracks` | try to find GEN match to every track
`era16_07Aug17` | use corresponding golden json and JECs
`era16_UL` | ...
`era16_UL_APV` | ...
## Input files
Can be a single file or list of files (either local or starting with /store/...). Also possible to use with grid-control (see examples and https://grid-control.github.io/index.html):
```
grid-control -cG example.conf
```
