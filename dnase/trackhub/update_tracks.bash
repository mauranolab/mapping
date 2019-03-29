#!/bin/bash
############################################################################
# This is the main driver script for the construction of a UCSC track hub.
# This script can be called from anywhere, as it figures out where it is located, 
# then moves there to continue working.
#
# For inputs, it needs to know if we are building a CEGS hub or a MAURANOLAB hub.
#
# It also needs to know a destination directory for the new hub. Full path plase!
# In practice, we expect this will be a CEGS development destination, a CEGS
# production destination, or a standard MAURANOLAB destination. The destination 
# can be any valid path though. It requires the directory to already exist and be empty.
#
# Next, it needs a "short label" and a "long label" for the hub.txt file.
# These labels appear in the browser Track Data Hubs page, and should be
# unique, so you know what hubs you are clicking on.
#
# Lastly, it needs the url and its corresponding full path to the relevant public_html 
# directory for the hub being built. It needs both of these so it can construct the proper 
# hub address, which is needed to bring up the hub in the UCSC browser.
############################################################################
# There are two sets of code, defined as development and production.
# They are located here:
#     /vol/cegs/src/trackhub/src_dev
#     /vol/cegs/src/trackhub/src_prod
#
# The main driver script, which exists in each directory, is: update_tracks.bash
# It can be called from anywhere. For example, if the working directory 
# is /vol/cegs/src/trackhub , then the development code could be launched as follows:
#
#
#     src_dev/update_tracks.bash  CEGS   /vol/cegs/public_html/trackhub_dev  "CEGS Dev"  "CEGS Development Hub"   \
#     ==========================  =====  ==================================  =========   ======================
#       Call the main script      Arg 1               Arg 2                    Arg 3             Arg 4
#
#                                 https://***REMOVED***@cascade.isg.med.nyu.edu/cegs/trackhub_dev
#                                 =============================================================
#                                                     Arg 5
#
# Arg 1:  Either CEGS or MAURANOLAB 
#
# Arg 2:  The full path to the new hub location, which must already exist and be empty.
#
# Args 3 & 4:  The short and long labels that appear on the UCSC browser Track Data Hubs page.
#              If you have multiple hubs you should make sure these are unique, so you know what you are clicking on.
#
# Args 5:      The complete url to the new hub location
#
############################################################################
# Standard execution of this script (for CEGS dev and MAURANOLAB):
#
# First, make sure the target directory exists and is empty.
# Then:
#
# /vol/cegs/src/trackhub/src_dev/update_tracks.bash CEGS /vol/cegs/public_html/trackhub_dev "CEGS Dev" "CEGS Development Hub" \
#                                                https://***REMOVED***@cascade.isg.med.nyu.edu/cegs/trackhub_dev
#
#
# /vol/cegs/src/trackhub/src_prod/update_tracks.bash MAURANOLAB /home/cadlej01/public_html/trackhub "Maurano Lab" "Maurano Lab Hub" \
#                                   https://***REMOVED***@cascade.isg.med.nyu.edu/~cadlej01/trackhub
#
############################################################################
# Standard execution of this script (for CEGS prod):
#
# To avoid UCSC errors appearing in the production hub, the production version of the CEGS hub is generated via two scripts:
#     /vol/cegs/src/trackhub/make_prod_hub_alpha.bash
#     /vol/cegs/src/trackhub/make_prod_hub_final.bash
#
# make_prod_hub_alpha.bash calls update_tracks.bash as above, except the hub is placed in:
#     /home/cadlej01/public_html/prod_test_trackhub
#
# make_prod_hub_final.bash clears out the old production hub, then copies the contents
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

hub_type=$1
hub_target=$2
short_label=$3
long_label=$4
URL_path=$5


# Check the inputs:

if [ "${hub_type}" != "CEGS" ] && [ "${hub_type}" != "MAURANOLAB" ]; then
    echo You need to enter a hub type. Either: CEGS or MAURANOLAB. Exiting...
    exit 1
fi

# We only accept an absolute path, not a relative path.
if [[ ! "${hub_target}" = /* ]]; then
    echo The target directory path is relative, not absolute. Exiting...
    exit 2
fi

# Make sure the target directory has already been built.
if [ ! -d "${hub_target}" ]; then
    echo The target directory does not exist.  Exiting...
    exit 3
fi

# Make sure the target directory is empty.
if [ "$(ls -A ${hub_target})" ]; then
    echo The target directory is not empty.  Exiting...
    exit 4
fi

echo Building new ${hub_type} track hub in ${hub_target}
echo Working...

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
# These are defined in the CEGS_genomes and MAURANOLAB_genomes files.
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

# For a CEGS hub, we also need to make tracks for the various assemblies. Do it here:
if [ ${hub_type} = "CEGS" ]; then
    echo Starting makeAssemblyTracks.bash
    ./makeAssemblyTracks.bash ${path_to_main_driver_script} ${hub_target} ${TMPDIR} "${genome_array[@]}"
fi

# Now construct the "flowcell" and "aggregation" tracks in TMPDIR.
echo Starting make_tracks.bash
./make_tracks.bash ${TMPDIR} ${hub_type} ${path_to_main_driver_script} "${genome_array[@]}"

######################################################################################
# Now copy the track information to the hub location.
echo Updating track files

update_genome () {
    genome=$1
    cd "${hub_target}/${genome}"

    if [ ${hub_type} = "CEGS" ]; then
        # Only CEGS has assembly tracks.
        cp "${TMPDIR}/assembly_tracks/trackDb_assemblies_${genome}.txt" trackDb_001.txt
    fi

    # Process the "flowcell" tracks.
    cat "${TMPDIR}/MakeTrackhub_${genome}_consolidated.out" >> trackDb_001.txt

    # Process the "aggregation" tracks.
    num_line=`(wc -l < "${TMPDIR}/samplesforTrackhub_${genome}_consolidated_agg.tsv")`
    if [ ${num_line} -gt 1 ]; then
       # If num_line == 1, then there are no aggregation tracks for this genome.
       cat "${TMPDIR}/MakeTrackhub_${genome}_consolidated_agg.out" >> trackDb_001.txt
    fi

    # Process the "publicdata" tracks.
    num_line=`(wc -l < "${TMPDIR}/samplesforTrackhub_${genome}_consolidated_pub.tsv")`
    if [ ${num_line} -gt 1 ]; then
       # If num_line == 1, then there are no publicdata tracks for this genome.
       cat "${TMPDIR}/MakeTrackhub_${genome}_consolidated_pub.out" >> trackDb_001.txt
    fi
}

for i in "${genome_array[@]}"; do
    update_genome $i 
done

######################################################################################

echo Making description files
cd ${path_to_main_driver_script}
./makeDescFiles.bash ${path_to_main_driver_script} ${hub_type} ${hub_target} \
                     ${TMPDIR} "${genome_array[@]}"

######################################################################################
# Make the hub.txt and genomes.txt files, and populate the structure with 
# other fixed, hand made assets.

hub_id=${short_label// /_}
echo "hub hub_id_${hub_id}" > "${hub_target}/hub.txt"
echo "shortLabel ${short_label}" >> "${hub_target}/hub.txt"
echo "longLabel ${long_label}" >> "${hub_target}/hub.txt"
echo "genomesFile genomes.txt" >> "${hub_target}/hub.txt"
echo "email cadley.nyulangone@gmail.com" >> "${hub_target}/hub.txt"
echo "descriptionUrl hub_description.html" >> "${hub_target}/hub.txt"

cd ${path_to_main_driver_script}

if [ "${hub_type}" = "CEGS" ]; then
    cp -R assets/CEGS/. ${hub_target}
else
    cp -R assets/MAURANOLAB/. ${hub_target}
fi

for i in "${genome_array[@]}"; do
    [ "${i}" = "cegsvectors" ] && continue
    echo "genome ${i}" >> "${hub_target}/genomes.txt"
    echo "trackDb ${i}/trackDb_001.txt" >> "${hub_target}/genomes.txt"
    echo " " >> "${hub_target}/genomes.txt"
done

time_stamp=$(date +"%m-%d-%y %T")
echo "<pre>Hub was constructed at: ${time_stamp} </pre>" >> "${hub_target}/genomes_description.html"
echo "<pre>Hub was constructed at: ${time_stamp} </pre>" >> "${hub_target}/hub_description.html"

######################################################################################
# Provide a README with the hub link.
echo "The UCSC browser url to this hub is:" | tee "${hub_target}/README"

# Construct url for the new hub.
echo "${URL_path}/hub.txt" | tee -a "${hub_target}/README"

######################################################################################
# Check for hub errors. Load kent module here in case makeAssemblyTracks.bash is ever not previously executed.
module load ucsckentutils/12152017

hubCheck -noTracks -udcDir=${TMPDIR} "${URL_path}/hub.txt"
ierr=$?
echo hubCheck error is: ${ierr}

if [ ${ierr} -eq 0 ]; then
   echo Removing TMP directory.
   rm -rf ${TMPDIR}
   echo Hub construction complete.
else
   echo Retaining TMPDIR.
   echo It is: ${TMPDIR}
fi
######################################################################################

