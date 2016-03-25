#! /bin/tcsh
#
# Sep. 22, 2014: modify the output directory so that the output goes to the
#                lpclljj disk
#
# submit all jobs based on a dataset, 
# GOOD RUN/LUMI sections is given in the json_file
# merge condor/crab submission together 
#
######################################################################
#output directories
if ($#argv <1) then 

echo "Usage: ${0} prodtag AOD data[/MC] user True[False on grid] runconfig dataset [runrange] [jsonfile]";
exit(0);
endif

echo $#argv;

set prodtag    = ${1}
set dataformat = ${2}
set dataflag   = ${3}
set username   = ${4}
set atfnalt1 = "${5}"
set runconfig  = ${6}
set dataset    = ${7}

set runrange   = "0-999999"
set jsonfile   = "Cert_13TeV_16Dec2015ReReco_Collisions15_25ns_JSON_Silver_v2.txt"


set useaod        = "True"
set isdata        = "False"
set subdir        = "noreplica/Summer2015/${prodtag}"
set temp          = `echo ${dataset} | awk -F/ '{print $2"_"$3}'`

set user          = "group/lpclljj"
if (${username} == "${USER}" ) then 

set user          = "user/${USER}"

endif


#data or MC
if (${dataflag} == "data") then

set isdata          = "True"
set subdir          = "noreplica/Run2015/${prodtag}"
endif

# output datatag
set outputdatatag   = "${temp}"
if ($#argv >= 8) then

set runrange        = "${8}"
set outputdatatag   = "${temp}_${8}"
endif

set crabfile_template = "./crab3_config.template"
set config_py         = "${prodtag}_${outputdatatag}.py"
set output_rootfile   = "${prodtag}_${outputdatatag}.root"
set crab_file         = "${prodtag}_${outputdatatag}_crab3.py"


#json file
if ($#argv >= 9) then

set jsonfile   = "${9}"
endif

echo "production tag = ${prodtag}";
echo "useaod  = ${useaod}";
echo "isdata = ${isdata}";
echo "subdir = ${subdir}";
echo "user   = ${user}";
echo "at T1    = ${atfnalt1}";
echo "runconfig  = $runconfig";
echo "dataset    = ${dataset}";
echo "config     = ${config_py}";
echo "runrange   = [$runrange]";
echo "ouputrootfile = ${output_rootfile}";
echo "jsonfile      = ${jsonfile}";

# 1) configuration file
sed -e "s:?useaod:${useaod}:g" \
    -e "s:?outputrootfile:${output_rootfile}:g" \
    -e "s:?isdata:${isdata}:g" \
        ${runconfig} > ${config_py}

sed -e "s:?isdata:${isdata}:g" \
    -e "s:?subdir:${subdir}:g" \
    -e "s:?atfnalt1:${atfnalt1}:g" \
    -e "s:?runconfig:${config_py}:g" \
    -e "s:?dataset:${dataset}:g" \
    -e "s:?outputdatatag:${outputdatatag}:g" \
    -e "s:?user:${user}:g" \
    -e "s:?jsonfile:${jsonfile}:g" \
    -e "s:?runrange:${runrange}:g" \
	${crabfile_template} > ${crab_file}

source /cvmfs/cms.cern.ch/crab3/crab.csh; eval `scramv1 runtime -csh`; crab submit ${crab_file}


cp  ${crab_file} ${subdir}

