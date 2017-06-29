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

echo "Usage: ${0} prodtag cfg_py input_dataset dbs dataset_tag [eosfiles]";
exit(0);
endif

echo $#argv;

set prodtag        = ${1}
set cfg_py         = ${2}
set input_dataset  = ${3}
set dbs            = ${4}
set dataset_tag    = ${5}
set output_tag     = ${1}
set crab_file      = "${1}_${5}_crab3.py"
set config_py      = "${1}_${5}.py"

set output_rootfile = "${1}_${5}.root" 
set isdata = "False"
set useaod = "True"

set subdir        = "noreplica/Private_Ntuple/${prodtag}"


set crabfile_template = "crab3_template.py";

set eos_files = "";
if ($#argv >5) then 

set eos_files = ${6}
set crabfile_template = "crab3_template_eosfiles.py";


endif     


echo $crabfile_template;


set hlt_nonisoelectron = "dummy"
set hlt_electron = "dummy"
set hlt_electronl = "dummy"
set hlt_muonelectron = "dummy"
set hlt_electronmuon = "dummy"
set hlt_nonisomuon = "dummy"
set hlt_muon = "dummy"
set hlt_muonl = "dummy"
set hlt_muon2 = "dummy"


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
        ${cfg_py} > ${config_py}



sed -e "s:?config_py:${config_py}:g" \
    -e "s:?input_dataset:${input_dataset}:g" \
    -e "s:?dbs:${dbs}:g" \
    -e "s:?subdir:${subdir}:g" \
    -e "s:?eos_files:${eos_files}:g" \
    -e "s:?output_tag:${output_tag}:g" \
    -e "s:?dataset_tag:${dataset_tag}:g" \
	${crabfile_template} > ${crab_file}

source /cvmfs/cms.cern.ch/crab3/crab.csh; eval `scramv1 runtime -csh`; crab submit ${crab_file}


cp  ${crab_file} ${subdir}

