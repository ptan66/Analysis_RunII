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

echo "Usage: ${0} prodtag AOD data[/MC] user True[False on grid] runconfig dataset ehlt muehlt muhlt [runrange] [jsonfile]";
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
set jsonfile   = "Cert_271036-280385_13TeV_PromptReco_Collisions16_JSON_NoL1T.txt"


set useaod        = "True"
set isdata        = "False"
set subdir        = "noreplica/Summer2016/${prodtag}"
set temp          = `echo ${dataset} | awk -F/ '{print $2"_"$3}'`

set user          = "group/lpclljj"
if (${username} == "${USER}" ) then 

set user          = "user/${USER}"

endif


#data or MC
if (${dataflag} == "data") then

set isdata          = "True"
set subdir          = "noreplica/Run2016/${prodtag}"
endif



set hlt_nonisoelectron  = `echo ${8} | awk -F: '{print $1}'`
set hlt_electron        = `echo ${8} | awk -F: '{print $2}'`
set hlt_electronl       = `echo ${8} | awk -F: '{print $3}'`

set hlt_muonelectron  = `echo ${9} | awk -F: '{print $1}'`
set hlt_electronmuon  = `echo ${9} | awk -F: '{print $2}'`

set hlt_nonisomuon  = `echo ${10} | awk -F: '{print $1}'`
set hlt_muon        = `echo ${10} | awk -F: '{print $2}'`
set hlt_muonl       = `echo ${10} | awk -F: '{print $3}'`
set hlt_muon2       = `echo ${10} | awk -F: '{print $4}'`



# output datatag
set outputdatatag   = "${temp}"
if ($#argv >= 11) then

set runrange        = "${11}"
set outputdatatag   = "${temp}_${11}"
endif

set crabfile_template = "./crab3_config.template"
set config_py         = "${prodtag}_${outputdatatag}.py"
set output_rootfile   = "${prodtag}_${outputdatatag}.root"
set crab_file         = "${prodtag}_${outputdatatag}_crab3.py"


#json file
if ($#argv >= 12) then

set jsonfile   = "${12}"
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
    -e "s:?hltSINGLEELE1:${hlt_nonisoelectron}:g" \
    -e "s:?hltSINGLEELE2:${hlt_electron}:g" \
    -e "s:?hltDOUBLEELE1:${hlt_electronl}:g" \
    -e "s:?hltMUE:${hlt_muonelectron}:g" \
    -e "s:?hltEMU:${hlt_electronmuon}:g" \
    -e "s:?hltSINGLEMU1:${hlt_nonisomuon}:g" \
    -e "s:?hltSINGLEMU2:${hlt_muon}:g" \
    -e "s:?hltDOUBLEMU1:${hlt_muonl}:g" \
    -e "s:?hltDOUBLEMU2:${hlt_muon2}:g" \
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

source /cvmfs/cms.cern.ch/crab3/crab.csh; eval `scramv1 runtime -csh`; 
#crab submit ${crab_file}


cp  ${crab_file} ${subdir}

