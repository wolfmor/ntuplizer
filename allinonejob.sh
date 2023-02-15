

#echo "================= getting the input files ===================="
#python -c "import PSet; print(PSet.process.source.fileNames.value())"
#echo "step 1"  > simpleoutput.txt
##ls -la >> simpleoutput.txt
## Ok, let's stop fooling around and execute the job:
#cmsRun -j FrameworkJobReport.xml -p PSet.py #config option crab 

#echo "step 2"  >> simpleoutput.txt
##ls -la >> simpleoutput.txt

##echo "================= producing the ntuples ===================="
#CRABFILES=$(python -c "import PSet; print(','.join(PSet.process.source.fileNames.value()))" 2>/dev/null| tail -1)
#python -c "import PSet; print(','.join(PSet.process.source.fileNames.value()))" 2>/dev/null| tail -1 >> simpleoutput.txt

#echo "I am a simple output for job "$1 >> simpleoutput.txt
#echo "Skript excecuted from" ${PWD} >> simpleoutput.txt
#BASEDIR=$(dirname $0)
#echo "Skript location" $BASEDIR >> simpleoutput.txt
#WHEREAMI=${PWD}

#python plantTrees.py inputFiles="$CRABFILES" tag="era16_07Aug17, crab, cleanleptons"
##python plantTrees.py inputFiles="$CRABFILES" tag="era16_07Aug17, crab"

##python plantTrees.py inputFiles="$CRABFILES" tag="era16_07Aug17, test, crab, debug"

#echo "step 3"  >> simpleoutput.txt
#ls -la >> simpleoutput.txt



echo "================= getting the input files ===================="
python -c "import PSet; print(PSet.process.source.fileNames.value())"
python -c "import PSet; print(PSet.process.dumpPython())" 
cmsRun -j FrameworkJobReport.xml -p PSet.py #config option crab 

echo "================= producing the ntuples ===================="

exitCode=1
exitMessage="The Ntupelizer Failed"
errorType="ntupelizer crash"

CRABFILES=$(python -c "import PSet; print(','.join(PSet.process.source.fileNames.value()))" 2>/dev/null| tail -1)
echo $2
python plantTrees.py inputFiles="$CRABFILES" $2 $3

# At the end of the script modify the FJR
ret=$?
if [ $ret -ne 0 ] && [ -e FrameworkJobReport.xml ]
then
	cat << EOF > FrameworkJobReport.xml.tmp
<FrameworkJobReport>
<FrameworkError ExitStatus="$exitCode" Type="$errorType" >
$exitMessage
</FrameworkError>
EOF

	tail -n+2 FrameworkJobReport.xml >> FrameworkJobReport.xml.tmp
	mv FrameworkJobReport.xml.tmp FrameworkJobReport.xml
elif [[ $ret -ne 0 ]]
then
	cat << EOF > FrameworkJobReport.xml
<FrameworkJobReport>
<FrameworkError ExitStatus="$exitCode" Type="$errorType" >
$exitMessage
</FrameworkError>
</FrameworkJobReport>
EOF
fi
