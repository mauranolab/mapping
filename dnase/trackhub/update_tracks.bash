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
elif [[ "${hub_type}" == "HOLTLAB" ]]; then
    customGenomeAssembly="NA"
    assemblyBaseDir="/vol/mauranolab/flowcells/public_html/holtlab"
else
    echo "ERROR You need to enter a valid hub type: CEGS, MAURANOLAB, HOLTLAB, SARS. Exiting..."
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
hub_target="${TMPDIR}/trackhub"
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
${src}/makeAssemblyTracks.bash ${src} ${hub_target} ${TMPDIR} ${hub_type} ${customGenomeAssembly} ${assemblyBaseDir} "${genome_array[@]}"

# Now construct the "flowcell" and "aggregation" tracks in TMPDIR.
echo
echo "Starting make_tracks.bash"
#
#make_tracks.bash is called only here
${src}/make_tracks.bash ${TMPDIR} ${hub_type} ${src} ${assemblyBaseDir} ${hub_target} "${genome_array[@]}"

######################################################################################
# Now copy the track information to the hub location.
echo
echo "Updating track files"
for genome in "${genome_array[@]}"; do
    if [ -f "${TMPDIR}/assembly_tracks/trackDb_assemblies_${genome}.txt" ]; then
        cp "${TMPDIR}/assembly_tracks/trackDb_assemblies_${genome}.txt" ${hub_target}/${genome}/trackDb.txt
    fi
    
    for curSection in Flowcells Aggregations Public_Data By_Locus; do
        if [ -f "${hub_target}/${genome}/trackDb.${curSection}.txt" ]; then
            #split into chunks to avoid "maxSize to udcFileReadAll is 16777216" error on hubCheck
            #like split -l 175000 -a 3 -d but takes care not to split stanzas which causes an error
            awk -v max=160000 -v outbase="${hub_target}/${genome}/trackDb.${curSection}." 'BEGIN {wantToSplit=0; curChunk=1} {print > outbase "_" curChunk ".txt"} NR % max == 0 {wantToSplit=1} wantToSplit && $0~/^[ \t]*$/ {curChunk+=1; wantToSplit=0}' ${hub_target}/${genome}/trackDb.${curSection}.txt
            
            rm -f ${hub_target}/${genome}/trackDb.${curSection}.txt
            for chunk in ${hub_target}/${genome}/trackDb.${curSection}.*.txt; do
                chunkFilename=`basename ${chunk}`
                echo -e "include ${chunkFilename}\n" >> ${hub_target}/${genome}/trackDb.txt
            done
        fi
    done
done


######################################################################################
echo
echo "Making description files"
#TODO this seems to be pretty slow, looks like the R code is the bottleneck
#makeDescFiles.bash is called only here
${src}/makeDescFiles.bash ${src} ${assemblyBaseDir} ${hub_type} ${hub_target} ${TMPDIR} "${genome_array[@]}"
######################################################################################
# Make the hub.txt and genomes.txt files, and populate the structure with other fixed, hand made assets.
echo
echo "Finalizing trackhub"

echo "hub hub_id_${short_label// /_}" > ${hub_target}/hub.txt
echo "shortLabel ${short_label}" >> ${hub_target}/hub.txt
echo "longLabel ${long_label}" >> ${hub_target}/hub.txt
echo "genomesFile genomes.txt" >> ${hub_target}/hub.txt
echo "email cadley.mauranolab@gmail.com" >> ${hub_target}/hub.txt
echo "descriptionUrl description.html" >> ${hub_target}/hub.txt


echo
cp -R --preserve=timestamps ${src}/assets/${hub_type}/. ${hub_target}

#Write out two line stanzas for genomes native to the UCSC browser
for cur_genome in "${genome_array[@]}"; do
    # The custom assemblies have more complicated, non-standard stanzas stored in the stub genomes.txt file in the assets/CEGS subdirectory, so skip building them in this loop.
    [ "${cur_genome}" = "${customGenomeAssembly}" ] || [ "${cur_genome}" = "t2t" ] && continue
    
    echo "genome ${cur_genome}" >> ${hub_target}/genomes.txt
    echo "trackDb ${cur_genome}/trackDb.txt" >> ${hub_target}/genomes.txt
    echo >> ${hub_target}/genomes.txt
done

if [ -f "${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit" ]; then
    cp --preserve=timestamps ${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit ${hub_target}/${customGenomeAssembly}/data/${customGenomeAssembly}.2bit
fi

if [[ "${hub_type}" == "CEGS" ]]; then
    mkdir "${hub_target}/t2t/data"
    cp --preserve=timestamps /vol/isg/annotation/fasta/t2t/t2t.2bit ${hub_target}/t2t/data/t2t.2bit
fi


######################################################################################
#Finalize description

time_stamp=$(date +"%m-%d-%y %T")
echo "<pre>Hub was constructed at: ${time_stamp} </pre>" >> ${hub_target}/description.html


if [[ "${customGenomeAssembly}" != "NA" ]]; then
    echo "Adding chromosomes table for the description page"
    for chrom_sizes_file in `find ${assemblyBaseDir}/sequences/${customGenomeAssembly}_* -mindepth 1 -maxdepth 1 -type f -name "*.chrom.sizes"`; do
        IFS=/ read -a get_assmbly <<< ${chrom_sizes_file}
        
        outputLine="${get_assmbly[4]#${customGenomeAssembly}_}|"
        while read chromSizes_line_in; do
            read chrom all_other <<< ${chromSizes_line_in}
        outputLine="${outputLine}<a href=hgTracks?genome=${customGenomeAssembly}&position=${chrom}>${chrom}</a><br>"
        done < ${chrom_sizes_file}
        
        outputLine="${outputLine%<br>}"  # Delete the last instance of <br> from $outputLine
        echo "${outputLine}"
    done | sort -k1,1 > ${TMPDIR}/chroms_per_cegsvector.txt
    
    Rscript ${src}/makeChroms_per_cegsvectorHtml.R ${TMPDIR}/chroms_per_cegsvector.txt ${TMPDIR}/chroms_per_cegsvector.html
    cat ${TMPDIR}/chroms_per_cegsvector.html >> ${hub_target}/description.html
fi

######################################################################################
# Check for hub errors.
echo
echo "Running hubCheck"
hubCheck -noTracks -udcDir=${TMPDIR} ${hub_target}/hub.txt
echo


######################################################################################
# Enable BLAT for custom assemblies.

if [[ "${hub_type}" == "CEGS" ]]; then
    blatport=17779
elif [[ "${hub_type}" == "MAURANOLAB" ]]; then
    blatport=17778
else
    blatport=0
fi

if [ "${blatport}" != "0" ]; then
    echo "Copying ${customGenomeAssembly}.2bit to shared directory for blat..."
    cp --preserve=timestamps ${hub_target_final}/${customGenomeAssembly}/data/${customGenomeAssembly}.2bit /vol/isg/blat_data
    twoBitFiles="${customGenomeAssembly}.2bit"
    
    if [[ "${hub_type}" == "CEGS" ]]; then
        echo "Copying t2t.2bit to shared directory for blat..."
        cp --preserve=timestamps ${hub_target}/t2t/data/t2t.2bit /vol/isg/blat_data
        twoBitFiles="${twoBitFiles} t2t.2bit"
    fi
    
    #NB requires password-less access via ssh
    echo "Restarting gfServer..."
    ssh isglcdcpvm001.nyumc.org "/usr/local/bin/blat/gfServer stop localhost ${blatport} -log=/vol/isg/blat_data/stop_gfServer_VM_${hub_type}.log; cd /vol/isg/blat_data; /usr/local/bin/blat/gfServer start localhost ${blatport} ${twoBitFiles} -canStop -log=/vol/isg/blat_data/start_gfServer_VM_${hub_type}.log -stepSize=5 > /dev/null &"
    echo "Finished restart of gfServer."
fi


######################################################################################
# Clean out old hub and deploy

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
echo "Deploying trackhub."
cp -rpd ${hub_target}/* ${hub_target_final}

echo "Removing TMPDIR directory."
rm -rf ${TMPDIR}
echo "Done!"

