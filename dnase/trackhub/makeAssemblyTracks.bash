#!/bin/bash
set -eu -o pipefail
####################################################################
# Get input parameters
src=$1
hub_target=$2
TMPDIR=$3
#unused
hub_type=$4
customGenomeAssembly=$5
assemblyBaseDir=$6
shift 6
genome_array=("$@")


####################################################################
# We need to make bigBED files from bed files. To do this we need 
# chrom.sizes files. So move to the right place for downloading them.
mkdir ${TMPDIR}/assembly_tracks

####################################################################
# make_bigBED  will make ".bb" files from ".bed" files.
# The output filename will be the same as the input filename, except with a "bb" suffix.
#
# bed_N is a utility function used by make_bigBED. It returns the "N" of the bedN suffix.
bed_N () {
    local filename=$1
    # filename is in the form:  <stuff>.bed or <stuff>.bedN or <stuff>.bedNN
    
    # This retrieves everything after the last "." in filename.
    # So bedN has the form of either:  "bed" or "bedN" or "bedNN".
    local bedN=${filename##*.}
    
    local N
    if [ "${bedN}" = "bed" ]; then
        N="4"
    else
        # This cuts the text "bed" out of the bedN variable.
        # So it leaves just the numeric part of bedN.
        N=${bedN/bed/}
    fi
    
    echo ${N}  # Returns N via $()
}

make_bigBED () {
    local myBEDfile=$1                # This is a full path name to a bed file.
    local expected_N=$2              #Integer number of bed fields
    local chrom_sizes=$3            # This is a full path name of one of the NYU chrom.sizes files.
    local output_file=$4
    
#    echo "[makeAssemblyTracks] ${myBEDfile} bed${expected_N}" > /dev/stderr
    
    # This scans the lines of myBEDfile, and looks for non-comment, non-empty lines.
    # It then checks if the number of fields in each of these lines is equal to expected_N.
    # It returns the number of these lines for which the number of fields is not equal to expected_N.
    local num_bad_lines=$(grep -v '^#' ${myBEDfile} | awk -F "\t" 'BEGIN {OFS="\t"} NF!=0 {print NF}' | grep -v -c ${expected_N} | sort | uniq)
    if [ ${num_bad_lines} -ne 0 ]; then
        # Issue a warning when a line in myBEDfile contains an unexpected number of fields.
        echo "[makeAssemblyTracks] WARNING Unexpected number of fields in ${myBEDfile}" > /dev/stderr
    fi
    
    grep -v '^#' ${myBEDfile} | cut -f1-${expected_N} | sort-bed - | awk -F "\t" 'BEGIN {OFS="\t"} \
        NF >= 5 { \
           field5 = int($5 + 0.5);
           $5 = (field5 <= 1000) ? field5 : 1000; \
        } {print}' > $TMPDIR/assembly_tracks/assembly_bed_tmp_sorted.bed
    
    # drop bigBed chatter from stderr ("pass1...", "pass2...")
    bedToBigBed -tab -type=bed${expected_N} $TMPDIR/assembly_tracks/assembly_bed_tmp_sorted.bed ${chrom_sizes} ${output_file} 2>&1 | grep -v "^pass[12] \- " > /dev/stderr || true  # Avoid pipefail crash since grep usually returns nothing.
}


###################################################################################
# Initialize an output file to hold paths to new bigBED files.
OUTFILE="${TMPDIR}/assembly_tracks/assembly_bb_files"
touch ${OUTFILE}

# Initialize the error log file where bedToBigBed problems will show up.
echo "[makeAssemblyTracks] Making annotation tracks..."

# This loop generates the ".bb" files which are found in the trackhub "<genome>/data" directories.
# genome_array just contains genomes to be included, e.g. hg38, mm10, rn6, and cegsvectors/mauranolab.
for genome in "${genome_array[@]}"; do
    if [ ! -d "${assemblyBaseDir}/sequences/${genome}" ]; then
        echo "[makeAssemblyTracks] Skipping ${genome} - no annotation tracks found"
        continue
    fi
    
    genome="${genome}/"
    if [ "${genome}" = "${customGenomeAssembly}/" ]; then
        chrom_sizes="${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.chrom.sizes"
        
        cat <<- HEREDOC_02 > ${TMPDIR}/assembly_tracks/trackDb_assemblies_${customGenomeAssembly}.txt
track GC_percent
shortLabel GC Percent
longLabel GC Percent
group ${customGenomeAssembly}
visibility full
maxHeightPixels 36
viewLimits 0:100
windowingFunction mean
type bigWig
bigDataUrl data/${customGenomeAssembly}.gc.bw
priority 1

HEREDOC_02
    else
        chrom_sizes="/vol/isg/annotation/fasta/${genome/\//}/${genome/\//}.chrom.sizes"
    fi
    
    for f in $(find ${assemblyBaseDir}/sequences/${genome} -mindepth 1 -maxdepth 1 -type d -not -path "./bak*" -not -path "./trash*" -printf '%f/\n'); do
        # assmbly looks like:  HPRT1
        assmbly=`basename ${f}`
        
        # Create an html file for this assembly. Later it will get moved to "descriptions" subdirectory.
        echo "<pre>" > ${TMPDIR}/descriptions/${genome%/}_assembly_${assmbly}_${genome%/}.html
        echo "Source data located in: ${assemblyBaseDir}/sequences/${genome}${assmbly}" >> ${TMPDIR}/descriptions/${genome%/}_assembly_${assmbly}_${genome%/}.html
        echo "</pre>" >> ${TMPDIR}/descriptions/${genome%/}_assembly_${assmbly}_${genome%/}.html
        
        
        # Note that we ignore emacs backups in the next line via the [0-9]*$
        for bed_file in $(find "${assemblyBaseDir}/sequences/${genome}${assmbly}" -mindepth 1 -maxdepth 1 -type f -regex '.*[.]bed[0-9]*$'); do
            # Make directories for the .bb files, as needed.
            mkdir -p ${hub_target}/${genome}data
            mkdir -p ${hub_target}/${genome}data/${assmbly}
            
            # Make the .bb files
            expected_N=$(bed_N ${bed_file})
            out_filename=`basename ${bed_file} | perl -pe 's/\.bed\d?$//g;'`".bb"
            make_bigBED ${bed_file} ${expected_N} ${chrom_sizes} ${hub_target}/${genome}data/${assmbly}/${out_filename}
            echo "${genome}data/${assmbly}/${out_filename} ${expected_N}" >> ${OUTFILE}
        done
    done
done


#Now that the directory exists, can create G+C and cytoBandIdeo tracks
if [ -f "${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.chrom.sizes" ] ; then
    cat ${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.chrom.sizes | LC_COLLATE=C sort -k1,1 -k2,2n | awk '{print $1,0,$2,$1,"gneg"}' > ${TMPDIR}/assembly_tracks/cytoBandIdeo.bed
    bedToBigBed -type=bed4 ${TMPDIR}/assembly_tracks/cytoBandIdeo.bed -as=${src}/cytoband.as ${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.chrom.sizes ${hub_target}/${customGenomeAssembly}/data/cytoBandIdeo.bb
    
    if [ -f "${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.2bit" ]; then
        # Divert again to make the GC percentage file:
        hgGcPercent -wigOut -doGaps -file=stdout -win=5 -verbose=0 ${customGenomeAssembly} ${assemblyBaseDir}/sequences/${customGenomeAssembly} > ${TMPDIR}/assembly_tracks/${customGenomeAssembly}.wig
        wigToBigWig ${TMPDIR}/assembly_tracks/${customGenomeAssembly}.wig ${assemblyBaseDir}/sequences/${customGenomeAssembly}/${customGenomeAssembly}.chrom.sizes ${hub_target}/${customGenomeAssembly}/data/${customGenomeAssembly}.gc.bw
    fi
fi


# Read in the lines from OUTFILE which contain the paths to the bigBED files we just made.
#  sort OUTFILE, so that the lines are all grouped by:
#      First:  genome
#      Second: assmbly
#      Last:   track names
old_genome=""

sort ${OUTFILE} | while read line_in; do
    # Get path to bb file and the bed_type of the file
    read -r line bed_type <<< "${line_in}"
    
    # Divide "line" into standard components. These will always have the same structure.
    # array: <genome>/data/<assmbly>/<track>.bb
    #           0      1       2        3        
    IFS='/' read -r -a array <<< "$line"
    genome="${array[0]}"
    assmbly="${array[2]}"
    
    if [ "${genome}" != "${old_genome}" ]; then
        # We found a new genome. Create a new file for it, and initialize it with supertrack lines.
        out_file="${TMPDIR}/assembly_tracks/trackDb_assemblies_${genome}.txt"
        old_genome=${genome}
        
        if [[ "${genome}" == "${customGenomeAssembly}" ]]; then
            # Add the cytoband track:
            echo "track cytoBandIdeo" >> ${out_file}
            echo "type bigBed" >> ${out_file}
            echo "shortLabel cytoBandIdeo" >> ${out_file}
            echo "longLabel Chromosome ideogram with cytogenetic bands" >> ${out_file}
            echo "bigDataUrl data/cytoBandIdeo.bb" >> ${out_file}
            echo >> ${out_file}
        fi
        
        echo "track ${genome}_Assemblies" >> ${out_file}
        echo "shortLabel Assemblies" >> ${out_file}
        echo "longLabel Assemblies" >> ${out_file}
        
        # We need this to avoid having the ${customGenomeAssembly} Assemblies being shown in the "Other" control group.
        if [[ "${genome}" == "${customGenomeAssembly}" ]]; then
            echo "group ${customGenomeAssembly}" >> ${out_file}
        fi
        echo "superTrack on show" >> ${out_file}
        
        echo "priority 10" >> ${out_file}
        echo >> ${out_file}
        
        # Since we have a new genome, we also need to give it new composite tracks.
        # The line below forces execution of the next if block, which makes composite tracks.
        old_assmbly=""
    fi
    
    if [ "${assmbly}" != "${old_assmbly}" ]; then
        # A new assembly was found. So give it new composite tracks.
        old_assmbly=${assmbly}
        echo "    track ${genome}_${assmbly}" >> ${out_file}
        echo "    compositeTrack on" >> ${out_file}
        echo "    visibility pack" >> ${out_file}
        echo "    type bigBed" >> ${out_file}
        echo "    html descriptions/${genome}_assembly_${assmbly}.html" >> ${out_file}
        echo "    shortLabel ${assmbly}" >> ${out_file}
        echo "    longLabel ${assmbly}" >> ${out_file}
        
        # Display assembly tracks by default if genome=${customGenomeAssembly}. Otherwise not.
        if [[ "${genome}" == "${customGenomeAssembly}" ]]; then
            # This does not work with the supertrack:  echo "    parent ${genome}_Assemblies on" >> ${out_file}
            echo "    parent ${genome}_Assemblies" >> ${out_file}
        else
            echo "    parent ${genome}_Assemblies off" >> ${out_file}
        fi
        
        # This keeps the assmbly composite on top of the other tracks in the browser (except for GC Content).
        echo "    priority 25" >> ${out_file}
        echo >> ${out_file}
    fi
    
    # Now we can print the required track information.
    BASE_bb=${array[3]}
    BASE=${BASE_bb/.bb/}         # Chop off trailing ".bb"
    
    tr_name=${BASE//./_}         # Replace periods with underscores, for track names.
                                 # Should no longer be needed, but doing it just in case.
    
    lbl_name=${tr_name//_/ }     # Replace underscores with spaces, for track labels.
    
    tr_name="${genome}_sequences_${tr_name}_bed"
    
    echo "        track ${tr_name}" >> ${out_file}
    echo "        parent ${genome}_${assmbly} on" >> ${out_file}
    echo "        visibility pack" >> ${out_file}
    echo "        longLabel ${lbl_name}" >> ${out_file}
    echo "        shortLabel ${lbl_name}" >> ${out_file}
    echo "        bigDataUrl ${array[1]}/${array[2]}/${array[3]}" >> ${out_file}
    
    echo "        type bigBed ${bed_type}" >> ${out_file}
    if [ "${bed_type}" -ge 9 ]; then
        echo "        itemRgb On" >> ${out_file}
    elif [ "${bed_type}" -ge 6 ]; then
            echo "        colorByStrand 255,0,0 0,0,255" >> ${out_file}
    fi
    echo >> ${out_file}
done
