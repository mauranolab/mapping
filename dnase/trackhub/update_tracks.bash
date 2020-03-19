#!/bin/bash
############################################################################
# Usage:
# This is the main driver script for the construction of a UCSC track hub.
# This script can be called from anywhere.
#
# update_tracks.bash hub_type hub_target short_label long_label URLbase
# hub_type:  One of: CEGS, MAURANOLAB, or SARS
# hub_target:  The full path to the new hub location, which must already exist and be empty.
# short_label, long_label:  The short and long labels that appear on the UCSC browser Track Data Hubs page. If you have multiple hubs you should make sure these are unique, so you know what you are clicking on.
# URLbase:      Complete URL (e.g. https://<login:pwd>@cascade.isg.med.nyu.edu/cegs) to the directory which contains the hub directory (e.g. "trackhub", trackhub_dev", "trackhub_sars"), used for hubCheck
#
#
# Standard execution of this script:
#/vol/cegs/src/trackhub/src_dev/update_tracks.bash CEGS /vol/cegs/public_html/trackhub_dev "CEGS Dev" "CEGS Development Hub" https://cascade.isg.med.nyu.edu/cegs
#/vol/cegs/src/trackhub/src_prod/update_tracks.bash MAURANOLAB /home/cadlej01/public_html/trackhub "Maurano Lab" "Maurano Lab Hub" https://cascade.isg.med.nyu.edu/~cadlej01
#/vol/mauranolab/mapped/src/dnase/trackhub/update_tracks.bash SARS /vol/mauranolab/sars/public_html/trackhub "SARS-CoV2" "SARS-CoV2 Hub" https://cascade.isg.med.nyu.edu/sars
#
# There are two sets of code, defined as development and production.
# They are located here:
#     /vol/cegs/src/trackhub/src_dev
#     /vol/cegs/src/trackhub/src_prod
#
# Standard execution of this script (for CEGS prod):
# To avoid UCSC errors appearing in the production hub, the production version of the CEGS hub is generated via two scripts:
# /vol/cegs/src/trackhub/make_prod_hub_alpha.bash calls update_tracks.bash as above, except the hub is placed in:
#     /home/cadlej01/public_html/prod_test_trackhub
#
# /vol/cegs/src/trackhub/make_prod_hub_final.bash clears out the old production hub, then copies the contents
# of prod_test_trackhub into /vol/cegs/public_html/trackhub
#
# The intent of the two-stage process is to allow for an opportunity to check the output of the
# production code prior to clearing out the old production directory.
#
############################################################################
#
# Outline of the script architecture: 
#
# update_tracks.bash       - The main script, which calls all the others.
#
# makeAssemblyTracks.bash  - Updates the bigBED files for the assemblies.
#                          - Creates track files for these bigBED files.
#                          - Is only executed with the CEGS arguement set.
#
# make_tracks.bash         - Updates the "flowcell" and "aggregates" tracks.
#
# makeDescFiles.bash       - When possible, creates description HTML files for tracks with readcounts files.
#
# samplesforTrackhub.R     - Called by make_tracks.bash
#                          - Scans directory trees, and builds an input file for MakeTrackhub.py
#
# MakeTrackhub.py          - Called by make_tracks.bash
#                          - Uses the Daler trackhub package to generate track files.
#
# makeDescHtml.R           - Called by makeDescFiles.bash
#                          - Creates html files from simple data tables.
#
############################################################################
module load ucsckentutils/379
module load bedops/2.4.37
module load R/3.5.0
module load python/3.8.1


############################################################################

hub_type=$1
hub_target=$2
short_label=$3
long_label=$4
URLbase=$5

# Check the inputs:
if [[ "${hub_type}" == "CEGS" ]]; then
    customGenomeAssembly="cegsvectors"
    assemblyBaseDir="/vol/cegs"
elif [[ "${hub_type}" == "MAURANOLAB" ]]; then
    customGenomeAssembly="mauranolab"
    assemblyBaseDir="/vol/mauranolab"
elif [[ "${hub_type}" == "SARS" ]]; then
    customGenomeAssembly="NA"
    assemblyBaseDir="/vol/mauranolab/sars"
else
    echo "ERROR You need to enter a valid hub type. Either: CEGS or MAURANOLAB or SARS. Exiting..."
    exit 1
fi

# We only accept an absolute path, not a relative path.
if [[ ! "${hub_target}" = /* ]]; then
    echo "ERROR The target directory path is relative, not absolute. Exiting..."
    exit 2
fi

# Make sure the target directory has already been built.
if [ ! -d "${hub_target}" ]; then
    echo "ERROR The target directory does not exist.  Exiting..."
    exit 3
fi

# Make sure the target directory is empty.
if [ "$(ls -A ${hub_target})" ]; then
    echo "ERROR The target directory is not empty.  Exiting..."
    exit 4
fi

echo "Building new ${hub_type} track hub in ${hub_target}"

# Make soft links in hub directory, so we can have multiple hubs with their own data sources in public_html.
for subdir in mapped aggregations publicdata; do
    if [[ -d "${assemblyBaseDir}/${subdir}" ]]; then
        ln -s ${assemblyBaseDir}/${subdir} ${hub_target}/
    fi
done

############################################################################
# We need a tmp directory to store intermediate files.

echo Creating a temporary directory
TMPDIR=`mktemp -d`   # TMPDIR has no trailing slash
echo TMPDIR is: ${TMPDIR}

# For testing (remember to comment out the rm at the bottom of this script)
# TMPDIR=/vol/cegs/src/trackhub/src_dev/myTMP
# rm -rf ${TMPDIR}
# mkdir ${TMPDIR}
# echo TMPDIR is: ${TMPDIR}

############################################################################

# Where is this file located? We use this info to find other required resources.
# The path will have no trailing slash.
path_to_main_driver_script=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${path_to_main_driver_script}

# What genomes will we be working with?
# These are defined in the CEGS_genomes, MAURANOLAB_genomes, and SARS_genomes files.
# Read the applicable one here:
readarray -t genome_array < "${path_to_main_driver_script}/${hub_type}_genomes"
echo " "
echo The genomes in this hub will be:
for i in "${genome_array[@]}"; do
    # Make a subdirectory for each genome we are working with.
    mkdir "${hub_target}/${i}"
    echo ${i}
done
echo " "

# We also need to make tracks for the various assemblies. Do it here:
echo "Starting makeAssemblyTracks.bash"
./makeAssemblyTracks.bash ${path_to_main_driver_script} ${hub_target} ${TMPDIR} ${hub_type} ${customGenomeAssembly} ${assemblyBaseDir} "${genome_array[@]}"

# Now construct the "flowcell" and "aggregation" tracks in TMPDIR.
hub_dir=`basename ${hub_target}`
echo "Starting make_tracks.bash"
./make_tracks.bash ${TMPDIR} ${hub_type} ${path_to_main_driver_script} ${assemblyBaseDir} ${hub_target} "${genome_array[@]}"

######################################################################################
# Now copy the track information to the hub location.
echo
echo "Updating track files"

make_track_include () {
    local tracks=$1
    if [ -f ${tracks} ]; then
       echo -e "include ${tracks}\n" >> trackDb_001.txt
    fi
}

update_genome () {
    genome=$1
    cd "${hub_target}/${genome}"
    
    if [ -f "${TMPDIR}/assembly_tracks/trackDb_assemblies_${genome}.txt" ]; then
        cp "${TMPDIR}/assembly_tracks/trackDb_assemblies_${genome}.txt" trackDb_001.txt
    fi
    
    if [[ "${genome}" == "${customGenomeAssembly}" ]]; then
        cp "${TMPDIR}/assembly_tracks/cytoBandIdeo.bigBed" data
    fi
    
    make_track_include trackDb.Flowcells.txt
    make_track_include trackDb.Aggregations.txt
    make_track_include trackDb.Public_Data.txt
    make_track_include trackDb.By_Locus.txt
    make_track_include trackDb.noSupertrack.txt
}

for i in "${genome_array[@]}"; do
    update_genome $i 
done

# Move the GC percentage file:
if [ -f "${TMPDIR}/assembly_tracks/${customGenomeAssembly}.gc.bw" ]; then
    cp ${TMPDIR}/assembly_tracks/${customGenomeAssembly}.gc.bw ${hub_target}/${customGenomeAssembly}/data
fi

######################################################################################

echo
echo "Making description files"
cd ${path_to_main_driver_script}
./makeDescFiles.bash ${path_to_main_driver_script} ${hub_type} ${hub_target} ${TMPDIR} "${genome_array[@]}"
######################################################################################
# Make the hub.txt and genomes.txt files, and populate the structure with other fixed, hand made assets.


echo
echo "Finalizing trackhub"
hub_id=${short_label// /_}
echo "hub hub_id_${hub_id}" > "${hub_target}/hub.txt"
echo "shortLabel ${short_label}" >> "${hub_target}/hub.txt"
echo "longLabel ${long_label}" >> "${hub_target}/hub.txt"
echo "genomesFile genomes.txt" >> "${hub_target}/hub.txt"
echo "email cadley.mauranolab@gmail.com" >> "${hub_target}/hub.txt"
echo "descriptionUrl description.html" >> "${hub_target}/hub.txt"


echo
echo "Deploying trackhub"
cd ${path_to_main_driver_script}

cp -R assets/${hub_type}/. ${hub_target}

for i in "${genome_array[@]}"; do
    [ "${i}" = "${customGenomeAssembly}" ] && continue
    echo "genome ${i}" >> "${hub_target}/genomes.txt"
    echo "trackDb ${i}/trackDb_001.txt" >> "${hub_target}/genomes.txt"
    echo " " >> "${hub_target}/genomes.txt"
done

if [ -f "${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit" ]; then
    cp ${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit ${hub_target}/${customGenomeAssembly}/data/${customGenomeAssembly}.2bit
fi
#########################################################

time_stamp=$(date +"%m-%d-%y %T")
echo "<pre>Hub was constructed at: ${time_stamp} </pre>" >> "${hub_target}/description.html"


######################################################################################
# Check for hub errors.
# Make sure some version of the ucsckentutils module has already been loaded.

echo
echo "Running hubCheck"
hubCheck -noTracks -udcDir=${TMPDIR} "${URLbase}/${hub_dir}/hub.txt"
ierr=$?
echo "hubCheck exit code is: ${ierr}"

echo

if [ ${ierr} -eq 0 ]; then
   echo "Removing TMP directory."
   rm -rf ${TMPDIR}
   echo "Done!"
else
   echo "Retaining TMPDIR."
   echo "It is: ${TMPDIR}"
fi
######################################################################################

