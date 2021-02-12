#!/bin/bash
set -eu -o pipefail
############################################################################
#
# This is the main driver script for the construction of a UCSC track hub.
# This script can be called from anywhere.
#
# Usage:
# update_tracks.bash <hub_type> <hub_target> <short_label> <long_label>
#     hub_type:   One of: CEGS, MAURANOLAB, or SARS
#     hub_target:  The full path to the new hub location, which must already exist.
#     short_label, long_label:  The short and long labels that appear on the UCSC browser Track Data Hubs page. If you have multiple hubs you should make sure these are unique, so you know what you are clicking on.
#
# Standard execution of this script:
#    For development:
#        /vol/cegs/src/trackhub/src_dev/update_tracks.bash CEGS /vol/cegs/public_html/trackhub_dev "CEGS Dev" "CEGS Development Hub"
#        /vol/cegs/src/trackhub/src_dev/update_tracks.bash SARS /home/cadlej01/public_html/trackhub_sars "SARS" "SARS Hub"
#
#    For production:
#        /vol/cegs/src/trackhub/src_prod/update_tracks.bash CEGS /vol/cegs/public_html/trackhub "CEGS" "CEGS Hub"
#        /vol/cegs/src/trackhub/src_prod/update_tracks.bash MAURANOLAB /home/cadlej01/public_html/trackhub "Maurano Lab" "Maurano Lab Hub"
#        /vol/mauranolab/mapped/src/dnase/trackhub/update_tracks.bash SARS /vol/sars/public_html/trackhub "SARS-CoV2" "SARS-CoV2 Hub"
#
# There are two sets of code, defined as development and production.
# They are located here:
#     /vol/cegs/src/trackhub/src_dev
#     /vol/cegs/src/trackhub/src_prod
#
############################################################################
#
# Outline of the script architecture: 
#
# update_tracks.bash       - The main script, which calls all the others.
#
# makeAssemblyTracks.bash  - Updates the bigBED files for the assemblies.
#                          - Creates supertracks for these bigBED files.
#
# make_tracks.bash         - Updates the Flowcells, Aggregations, Public_Data, and By_Locus supertracks.
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
module load miller/5.4.0
############################################################################

hub_type=$1
hub_target_final=$2
short_label=$3
long_label=$4

# Check the inputs:
if [[ "${hub_type}" == "CEGS" ]]; then
    customGenomeAssembly="cegsvectors"
    assemblyBaseDir="/vol/cegs"
elif [[ "${hub_type}" == "MAURANOLAB" ]]; then
    customGenomeAssembly="mauranolab"
    assemblyBaseDir="/vol/mauranolab"
elif [[ "${hub_type}" == "SARS" ]]; then
    customGenomeAssembly="NA"
    assemblyBaseDir="/vol/sars"
else
    echo "ERROR You need to enter a valid hub type. Either: CEGS or MAURANOLAB or SARS. Exiting..."
    exit 1
fi

# We only accept an absolute path, not a relative path.
if [[ ! "${hub_target_final}" = /* ]]; then
    echo "ERROR The target directory path is relative, not absolute. Exiting..."
    exit 2
fi

# Make sure the target directory has already been built.
if [ ! -d "${hub_target_final}" ]; then
    echo "ERROR The target directory does not exist.  Exiting..."
    exit 3
fi

echo "Building new ${hub_type} track hub in ${hub_target_final}"

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

hub_target="/${TMPDIR}/trackhub"
mkdir ${hub_target}

# Make soft links in hub directory, so we can have multiple hubs with their own data sources in public_html.
for subdir in mapped aggregations publicdata; do
    if [[ -d "${assemblyBaseDir}/${subdir}" ]]; then
        ln -s ${assemblyBaseDir}/${subdir} ${hub_target}/
    fi
done

############################################################################

# Where is this file located? We use this info to find other required resources.
# The path will have no trailing slash.
src=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )
cd ${src}

# What genomes will we be working with?
# These are defined in the CEGS_genomes, MAURANOLAB_genomes, and SARS_genomes files.
# Read the applicable one here:
readarray -t genome_array < "${src}/assets/${hub_type}_genomes"
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
./makeAssemblyTracks.bash ${src} ${hub_target} ${TMPDIR} ${hub_type} ${customGenomeAssembly} ${assemblyBaseDir} "${genome_array[@]}"

# Now construct the "flowcell" and "aggregation" tracks in TMPDIR.
hub_dir=`basename ${hub_target}`
echo "Starting make_tracks.bash"
./make_tracks.bash ${TMPDIR} ${hub_type} ${src} ${assemblyBaseDir} ${hub_target} "${genome_array[@]}"

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
cd ${src}
./makeDescFiles.bash ${src} ${assemblyBaseDir} ${hub_type} ${hub_target} ${TMPDIR} "${genome_array[@]}"
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
cd ${src}

cp -R assets/${hub_type}/. ${hub_target}

# The values held in genome_array are the genome names, like hg38, mm10, rn6, etc.
for cur_genome in "${genome_array[@]}"; do
    # The genomes native to the UCSC browser have simple, standard, two line stanzas in the genomes.txt file.
    # They are written to genomes.txt in this loop. The custom assemblies have more complicated, non-standard stanzas.
    # Those stanza lines are stored in the stub genomes.txt file in the assets/CEGS subdirectory, so we skip building them in this loop.
    [ "${cur_genome}" = "${customGenomeAssembly}" ] || [ "${cur_genome}" = "t2t" ] && continue

    echo "genome ${cur_genome}" >> "${hub_target}/genomes.txt"
    echo "trackDb ${cur_genome}/trackDb_001.txt" >> "${hub_target}/genomes.txt"
    echo " " >> "${hub_target}/genomes.txt"
done

if [ -f "${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit" ]; then
    cp ${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit ${hub_target}/${customGenomeAssembly}/data/${customGenomeAssembly}.2bit
fi

if [[ "${hub_type}" == "CEGS" ]]; then
    mkdir "${hub_target}/t2t/data"
    cp -p /vol/isg/annotation/fasta/t2t/t2t.2bit ${hub_target}/t2t/data/t2t.2bit
fi
#########################################################

time_stamp=$(date +"%m-%d-%y %T")
echo "<pre>Hub was constructed at: ${time_stamp} </pre>" >> "${hub_target}/description.html"

######################################################################################
# Make a reference vs chromosomes table for the cegsvectors description page.

chrom_sizes_file_list=$(find /vol/cegs/sequences/cegsvectors_* -mindepth 1 -maxdepth 1 -type f -name "*.chrom.sizes")

for chrom_sizes_file in ${chrom_sizes_file_list}; do
    IFS=/ read -a get_assmbly <<< ${chrom_sizes_file}
    outputLine="${get_assmbly[4]#cegsvectors_}|"
    
    while read chromSizes_line_in; do
        read chrom all_other <<< ${chromSizes_line_in}
    outputLine="${outputLine}<a href=hgTracks?genome=cegsvectors&position=${chrom}>${chrom}</a><br>"
    done < ${chrom_sizes_file}
    
    outputLine="${outputLine%<br>}"  # Delete the last instance of <br> from $outputLine
    echo "${outputLine}"
done | sort -k1,1 > "${TMPDIR}/chroms_per_cegsvector.txt"

Rscript makeChroms_per_cegsvectorHtml.R "${TMPDIR}/chroms_per_cegsvector.txt" "${TMPDIR}/chroms_per_cegsvector.html"

cat "${TMPDIR}/chroms_per_cegsvector.html" >> "${hub_target}/description.html"

######################################################################################
# Check for hub errors.
# Make sure some version of the ucsckentutils module has already been loaded.

echo
echo "Running hubCheck"
hubCheck -noTracks -udcDir=${TMPDIR} "${hub_target}/hub.txt"
echo

# Clean out old hub:
for i in "${genome_array[@]}"; do
    rm -rf "${hub_target_final}/${i}"
done
rm -f  ${hub_target_final}/genomes.txt
rm -f  ${hub_target_final}/description.html
rm -f  ${hub_target_final}/hub.txt
rm -f  ${hub_target_final}/publicdata
rm -f  ${hub_target_final}/aggregations
rm -f  ${hub_target_final}/mapped

# Make sure there is no junk left in there:
if [ "$(ls -A ${hub_target_final})" ]; then
    echo "ERROR The final target directory is not empty.  Exiting..."
    echo "Not updating hub. Retaining TMPDIR."
    echo "It is: ${TMPDIR}"
    exit 4
fi

# Replace old hub with new hub:
echo "Copying the hub."
cp -rpd ${hub_target}/* ${hub_target_final}

echo "Removing TMPDIR directory."
rm -rf ${TMPDIR}
echo "Done!"
######################################################################################

# Enable BLAT for custom assemblies...
if [[ "${hub_type}" == "CEGS" ]]; then
    blatport=17779
elif [[ "${hub_type}" == "MAURANOLAB" ]]; then
    blatport=17778
else
    blatport=0
fi

if [ "${blatport}" != "0" ]; then
    echo "Copying ${customGenomeAssembly}.2bit to cadlej01_shared..."
    cp -p ${hub_target_final}/${customGenomeAssembly}/data/${customGenomeAssembly}.2bit /vol/isg/cadlej01_shared
    twoBitFiles="${customGenomeAssembly}.2bit"

    if [[ "${hub_type}" == "CEGS" ]]; then
        echo "Copying t2t.2bit to cadlej01_shared..."
        cp -p /vol/isg/annotation/fasta/t2t/t2t.2bit /vol/isg/cadlej01_shared
        twoBitFiles="${twoBitFiles} t2t.2bit"
    fi

    echo "Restarting gfServer..."
    ssh isglcdcpvm001.nyumc.org "/usr/local/bin/blat/gfServer stop localhost ${blatport} -log=/vol/isg/cadlej01_shared/stop_gfServer_VM_${hub_type}.log; cd /vol/isg/cadlej01_shared; /usr/local/bin/blat/gfServer start localhost ${blatport} ${twoBitFiles} -canStop -log=/vol/isg/cadlej01_shared/start_gfServer_VM_${hub_type}.log -stepSize=5 > /dev/null &"
    echo "Finished restart of gfServer."
fi

