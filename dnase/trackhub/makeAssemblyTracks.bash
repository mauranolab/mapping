#!/bin/bash
####################################################################
# Get input parameters
path_to_main_driver_script=$1
hub_target=$2
TMPDIR=$3
shift 3
genome_array=("$@")

####################################################################
# We use the Kent utilities.
module load ucsckentutils/12152017
####################################################################
# We need to make bigBED files from bed files. To do this we need 
# chrom.sizes files. So move to the right place for downloading them.
mkdir "${TMPDIR}/assembly_tracks"
cd "${TMPDIR}/assembly_tracks"

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
    local myBEDfile=$1                 # This is a full path name to a bed file.
    local chrom_sizes=$2               # This is a full path name of one of the NYU chrom.sizes files.

    local out_tmp=${myBEDfile%.bed*}        # Chops off the trailing ".bed" or ".bedN" from myBEDfile.

    local output_file=${out_tmp##*/}".bb"   # This retrieves everything after the last / in out_tmp.
                                            # Since the previous statment chopped off the bed extension,
                                            # this will leave only the short file name. We then add the ".bb".
                                            # So output_file has a short form:  <name>.bb
                                            # Note that output_file does not contain path elements.

    local expected_N=$(bed_N ${myBEDfile})  # myBEDfile is in the form <stuff>.bed or <stuff>.bedN or <stuff>.bedNN
                                            # The function bed_N returns the N or NN from the myBEDfile name.
                                            # If myBEDfile is in the form <stuff>.bed, then bed_N returns 3.

    # This scans the lines of myBEDfile, and looks for non-comment, non-empty lines.
    # It then checks if the number of fields in each of these lines is equal to expected_N.
    # It returns the number of these lines for which the number of fields is not equal to expected_N.
    local num_bad_lines=$(grep -v '^#' ${myBEDfile} | awk -F "\t" 'BEGIN {OFS="\t"} {print NF}' | \
         uniq | sort | uniq | grep -v '^0$' | grep -v -c ${expected_N})

    if [ ${num_bad_lines} -ne 0 ]; then
        # Issue a warning when a line in myBEDfile contains an unexpected number of fields.
        echo "Unexpected number of fields in ${myBEDfile}"  >> make_bigBED.log
        (>&2 echo "Unexpected number of fields in ${myBEDfile}")
    fi

    sort -k1,1 -k2,2n ${myBEDfile} | cut -d $'\t' -f1-${expected_N} > myBEDfile_sorted.bed

    if [ "${expected_N}" -ge 5 ]; then
        # Round column 5 of the bed file to an integer.
        # It also needs to be less than or equal to 1,000.

awk -f <(cat << "AWK_HEREDOC_01"
BEGIN{OFS="\t"}
{
   if(NF == 0) next;

   field5 = int($5 + 0.5)
   $5 = (field5 <= 1000) ? field5 : 1000
   print $0
}
END{}
AWK_HEREDOC_01
) < myBEDfile_sorted.bed > myBEDfile_sorted_awk.bed

        mv myBEDfile_sorted_awk.bed myBEDfile_sorted.bed
    fi

    # For debugging
    echo ${myBEDfile} "bed${expected_N}" >> make_bigBED.log

    # Make the bigBed file, and redirect stderr to the log file.
    bedToBigBed -tab -type="bed${expected_N}" \
                myBEDfile_sorted.bed ${chrom_sizes} ${output_file} 2>> make_bigBED.log

    # Clean up
    rm myBEDfile_sorted.bed

    # Return the output filename. Acquired by calling routine via $()
    echo ${output_file}
}

############################################################################
# This function makes the hub genome data directories as required, if they do not already exist.
make_directory () {
    local genome
    genome=$1

    if [ ! -d "${hub_target}/${genome}/data" ]; then
        mkdir "${hub_target}/${genome}/data"
    fi

    # Locus subdirectory...
    if [ ! -d $2 ]; then
        mkdir $2
    fi
}

###################################################################################
# Initialize an output file to hold paths to new bigBED files. We'll need the full path below. Clean up, just in case.
OUTFILE="${TMPDIR}/assembly_tracks/output_full_paths"

# Initialize the error log file where bedToBigBed problems will show up.
echo Starting calls to bedToBigBed... > make_bigBED.log

# Move to where the assembly sequences are kept, and get the names of the genome directories.
BASE="/vol/cegs/sequences/"
cd ${BASE}
genome_dirs=($(ls -d */))   # Elements will look like:  hg38/

# All the elements of "genome_dirs" should be in "genome_array" as well.
# genome_array could be bigger than genome_dirs though.

for genome in "${genome_dirs[@]}"; do
    # cegsvectors assembly tracks are hand made at this time.
    if [ "${genome}" = "cegsvectors/" ]; then
        cp "${path_to_main_driver_script}/assets/CEGS/trackDb_assemblies_cegsvectors.txt" "${TMPDIR}/assembly_tracks"
        continue
    fi

    cd ${BASE}${genome}
    assmbly_dirs=($(ls -d */))   # Elements will look like:  HPRT1/

    for assmbly in "${assmbly_dirs[@]}"; do
        # Note we that ignore emacs backups in the next line via the [/d]?$
        bed_files=($(ls "${BASE}${genome}${assmbly}"* | egrep *[.]bed[0-9]*$))
    
        cd "${TMPDIR}/assembly_tracks"
        for bed_file in "${bed_files[@]}"; do
            # 'bed_file' is in the form: <BASE><genome><assmbly><bed name>.bedN  (or ".bed" or ".bedNN")

            # Get just the bedN part
            bed_type=${bed_file##*.}                 # This retrieves everything after the last "." in bed_file
            if [ "${bed_type}" = "bed" ]; then
                bed_type="bed4"
            fi
        
            # Make directories for the .bb files, as needed.
            make_directory ${genome} "${hub_target}/${genome}data/${assmbly}"

            # Make the .bb files
            chrom_sizes="/vol/isg/annotation/fasta/${genome/\//}/${genome/\//}.chrom.sizes"
            out_file=$(make_bigBED ${bed_file} ${chrom_sizes})
        
            # Move the .bb files into the appropriate directories
            mv ${out_file} "${hub_target}/${genome}data/${assmbly}"
            echo "${hub_target}/${genome}data/${assmbly}${out_file}" "${bed_type}" >> ${OUTFILE}
        done
    done
done
# At the end of the above for loops, we're now in ${TMPDIR}/assembly_tracks

# Sort OUTFILE (made in the above for loops), so that the lines are all grouped by:
#      First:  genome
#      Second: assmbly
#      Last:   track names
# This will do that:
sort ${OUTFILE} > ${OUTFILE}_tmp
mv ${OUTFILE}_tmp ${OUTFILE}

# Clear out old track files, just in case.
for i in "${genome_array[@]}"; do
    [ "${i}" = "cegsvectors" ] && continue
    rm -f "trackDb_assemblies_${i}.txt"
done

# This is used in the below while loop.
old_genome=""

# Read in the lines from OUTFILE.
# These lines contain the paths to the bigBED files we just made.
# Since we sorted OUTFILE above, they are already organized by genome/assmbly/trackname.
# So we know that when we find a new genome, we must already be done with the old genome, etc.
while read -r line_in ; do
    # Get path to bb file and the bed_type of the file
    read -r line_in2 bed_type <<< "${line_in}"

    # Strip out the base path, which could be anything.
    line=${line_in2/${hub_target}\//}

    # Divide "line" into standard components. These will always have the same structure.
    # array: <genome>/data/<assmbly>/<track>.bb
    #           0      1       2        3        
    IFS='/' read -r -a array <<< "$line"
    genome="${array[0]}"
    assmbly="${array[2]}"

    if [ "${genome}" != "${old_genome}" ]; then
        # We found a new genome. Create a new file for it, and initialize it with supertrack lines.
        out_file="trackDb_assemblies_"${genome}".txt"
        old_genome=${genome}

        echo track Assemblies >> ${out_file}
        echo shortLabel Assemblies >> ${out_file}
        echo longLabel Assemblies >> ${out_file}
        echo superTrack on >> ${out_file}
        echo " " >> ${out_file}

        # Since we have a new genome, we also need to give it new composite tracks.
        # The line below forces execution of the next if block, which makes composite tracks.
        old_assmbly=""
    fi

    if [ "${assmbly}" != "${old_assmbly}" ]; then
        # A new assmbly was found. So give it new composite tracks.
        old_assmbly=${assmbly}
        echo "    track ${assmbly}" >> ${out_file}
        echo "    compositeTrack on" >> ${out_file}
        echo "    type bigBed" >> ${out_file}
        echo "    shortLabel ${assmbly}" >> ${out_file}
        echo "    longLabel ${assmbly}" >> ${out_file}
        echo "    parent Assemblies on" >> ${out_file}
        echo " " >> ${out_file}
    fi

    # Now we can print the required track information.
    URL=${array[1]}"/"${array[2]}"/"${array[3]}
    BASE_bb=${array[3]}
    BASE=${BASE_bb/.bb/}         # Chop off trailing ".bb"

    tr_name=${BASE//./_}         # Replace periods with underscores, for track names.
                                 # Should no longer be needed, but doing it just in case.

    lbl_name=${tr_name//_/ }     # Replace underscores with spaces, for track labels.

    echo "        track ${tr_name}" >> ${out_file}
    echo "        parent ${assmbly} off" >> ${out_file}
    echo "        visibility pack" >> ${out_file}
    echo "        longLabel ${lbl_name}" >> ${out_file}
    echo "        shortLabel ${lbl_name}" >> ${out_file}
    echo "        bigDataUrl ${URL}" >> ${out_file}

    # bed_type is in the form "bedN".
    # This cuts the text "bed" out of the bed_type variable.
    # So it leaves just the numeric part of bedN.
    N=${bed_type/bed/}
    echo "        type bigBed ${N}" >> ${out_file}
    echo " " >> ${out_file}
done < ${OUTFILE}

# Clean up
rm -f ${OUTFILE}

#########################################################
