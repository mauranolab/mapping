#!/bin/bash
##############################################################################
# Pass in the input parameters:

path_to_main_driver_script=$1
hub_type=$2
hub_target=$3

# We need a tmp directory to store intermediate files.
# One was established in the calling script.
TMPDIR=$4

shift 4
genome_array=("$@")

##############################################################################
module load R/3.5.0
##############################################################################
# This function calls the R code that makes the data table html code
# It is called by the below "find_readcounts" function.

make_html () {
    local infile=$1
    local desc_file=$2
    local target=$3
    
    # If only the header is there now, that means we have no data for that genome.
    local n_lines=$(wc -l < ${infile})
    
    if (( n_lines == 1)); then
        rm ${infile}
    else
        # Make the html version of the genome's readcounts data.
        Rscript ${path_to_main_driver_script}/makeDescHtml.R ${infile} "${infile}.html"
    fi
    
    local tmpdir=$(dirname "${desc_file}")
    echo "<pre>" > "${tmpdir}/make_html_tmp"

    if [ -f ${desc_file} ]; then
        cat "${desc_file}" >> "${tmpdir}/make_html_tmp"
        echo " " >> "${tmpdir}/make_html_tmp"
    fi

    echo "Source data located in: ${target%/readcounts.summary.txt}" >> "${tmpdir}/make_html_tmp"
    echo " " >> "${tmpdir}/make_html_tmp"

    echo "</pre>" >> "${tmpdir}/make_html_tmp"
    echo "<hr>" >> "${tmpdir}/make_html_tmp"
    
    if (( n_lines != 1)); then
        cat "${infile}.html" >> "${tmpdir}/make_html_tmp"
    fi
    
    mv "${tmpdir}/make_html_tmp" "${infile}.html"
}
##############################################################################
# These functions makes html versions of the readcounts files, by genome.

find_targets () {
    # This function loads target_vec for the flowcell directories.
    
    local subdir_names=("$@")
    local target
    local BASE
    local BASE2
    
    for i in "${subdir_names[@]}"; do
        BASE2=${i%/}       # Kill trailing slash
        BASE=${BASE2##*/}  # Get subdir name
        
        if [ "${BASE}" = "trash" ] || [ "${BASE}" = "bak" ] || \
           [ "${BASE}" = "trash2" ] || [ "${BASE}" = "trash.oops" ] ; then
             continue
        fi
        
        target="${i}readcounts.summary.txt"
        target_vec+=( ${target} )
    done
}

find_readcounts () {
    local dir_base=$1
    local suffix=$2
    shift 2
    local dir_names=("$@")
    
    local flowcell
    local ymd
    local BASE
    local BASE2
    local target
    local desc_file
    local tmp_flowcell
    local d_out
    local yyyy
    local mm
    local dd
    local info_file
    local infile

    # Each flowcell directory may be associated with a readcounts file.
    # Get the data from the readcounts file if it is there, and make a table from the data.
    for flowcell in "${dir_names[@]}"; do
        # No readcounts files exist in these two directories.
        # If encountered, the subsequent $target test will fail, and we move on to the next flowcell.
        #     /vol/mauranolab/mapped/src/
        #     /vol/mauranolab/mapped/trackhub/
        
        # Fill up target_vec:
        if [ ${dir_base} = "/vol/cegs/mapped/" ] || [ ${dir_base} = "/vol/mauranolab/mapped/" ] ; then
            subdir_names=($(ls -d ${dir_base}${flowcell}*/ 2> /dev/null))    # These are full paths, with trailing slashes.
            numElements="${#subdir_names[@]}"
            if [ "${numElements}" = "0" ]; then
                echo "[makeDescFiles.bash] WARNING No subdirectories in ${dir_base}${flowcell}"
                continue
            fi
            target_vec=()
            find_targets  "${subdir_names[@]}"
        else
            # For non-flowcell directories, there is only one place for the readcounts file to be:
            target=${dir_base}${flowcell}${suffix}
            target_vec=()
            target_vec+=( ${target} )
        fi
        
        for target in ${target_vec[@]}; do
            # Find a date to associate with the flowcell from the info.txt file
            info_file="/vol/mauranolab/flowcells/data/"${flowcell}"info.txt"
            
            # Make sure the file is there:
            if [ ! -f ${info_file} ]; then
                # It is not.
                ierr=0
            else
                # It is, but is the date there?
                ierr=$(grep -c '#Load date' ${info_file})
            fi
            
            if (( ierr == 0)); then
                ymd="00000000"
            else
                d_out=$(grep '#Load date' ${info_file})
                yyyy=$(echo $d_out | cut -d' ' -f3 | cut -d'-' -f1)
                mm=$(echo $d_out | cut -d' ' -f3 | cut -d'-' -f2)
                dd=$(echo $d_out | cut -d' ' -f3 | cut -d'-' -f3)
                ymd=${yyyy}${mm}${dd}
            fi
            
            # We now have a date. Next, construct a genome related filename for the output.
            if [ "${ymd}" = "00000000" ]; then
                tmp_flowcell="/"${flowcell%/}
            else
                tmp_flowcell="/"${ymd}"_"${flowcell%/}
            fi
            
            if [ ${dir_base} = "/vol/cegs/mapped/" ] || [ ${dir_base} = "/vol/mauranolab/mapped/" ] ; then
                # The name of the flowcell subdirectory that the readcounts.summary.txt file lives in is also the name of the relevant assay type.
                # Use that assay type to help make a unique html name.
                BASE=${target%/readcounts.summary.txt}       # Kill trailing /readcounts.summary.txt
                BASE2=${BASE##*/}  # Get subdir name
                tmp_flowcell="${tmp_flowcell}_${BASE2}"
            fi
            # tmp_flowcell will be used in the for loop below.
            
            if [ -f ${target} ]; then
                # There is a readcounts file in the flowcell directory.
                echo ${target}
                
                # Initialize the output files with header lines.
                head -n1 ${target} > "${TMPDIR}/tmp_header"
                
                # Use the header we just made to initialize the genome specific output files.
                # Break up the readcounts files into their genomic subsections.
                for i in "${genome_array[@]}"; do
                    declare "${i}"="${TMPDIR}${tmp_flowcell}_${i}"
                    ref="${i}"
                    cat "${TMPDIR}/tmp_header" > ${!ref}
                    grep ${i} ${target} >> ${!ref}
                    
                    BASE=${target%readcounts.summary.txt}       # Kill trailing readcounts.summary.txt
                    desc_file="${BASE}description.html"
                    
                    declare "${i}_desc"="${TMPDIR}${tmp_flowcell}_${i}_desc"
                    ref2="${i}_desc"
                    
                    if [ -f ${desc_file} ]; then
                        cat ${desc_file} >  ${!ref2}
                    else
                        # Make sure ref2 does not exist. Maybe unnecessary ?
                        rm -f ${!ref2}
                    fi
                done
                
                for i in "${genome_array[@]}"; do
                    ref="${i}"
                    ref2="${i}_desc"
                    make_html ${!ref} ${!ref2} ${target}
                done
            else
                # No readcounts file exists. UCSC browser description area is just
                # blank if there is no file, or a bad url. No error msg pops up.
                echo ${target} does not exist
                
                BASE=${target%readcounts.summary.txt}       # Kill trailing readcounts.summary.txt
                desc_file="${BASE}description.html"
                
                if [ -f ${desc_file} ]; then
                    for i in "${genome_array[@]}"; do
                        infile="${TMPDIR}${tmp_flowcell}_${i}"
                        echo "<pre>" > "${infile}.html"
                        cat "${desc_file}" >> "${infile}.html"
                        echo "</pre>" >> "${infile}.html"
                        echo "<hr>" >> "${infile}.html"
                    done
                fi
            fi
        done
    done
}
##############################################################################
# Make the html files.
if [ "${hub_type}" = "CEGS" ]; then
    # CEGS flowcells
    dir_base="/vol/cegs/mapped/"
    cd ${dir_base}
    dir_names=($(ls -d */))    # These have trailing slashes.
    suffix="NA"
    find_readcounts ${dir_base} ${suffix} "${dir_names[@]}"
    
    # CEGS aggregations
    dir_base="/vol/cegs/aggregations/"
    cd ${dir_base}
    dir_names=($(ls -d */))     # These have trailing slashes.
    suffix="readcounts.summary.txt"
    find_readcounts ${dir_base} ${suffix} "${dir_names[@]}"
    
    # CEGS publicdata
    dir_base="/vol/cegs/publicdata/"
    cd ${dir_base}
    dir_names=($(ls -d */))     # These have trailing slashes.
    suffix="readcounts.summary.txt"
    find_readcounts ${dir_base} ${suffix} "${dir_names[@]}"
else
    # Maurano flowcells
    dir_base="/vol/mauranolab/mapped/"
    cd ${dir_base}
    dir_names=($(ls -d */))     # These have trailing slashes.
    suffix="NA"
    find_readcounts ${dir_base} ${suffix} "${dir_names[@]}"
    
    # Maurano aggregations
    dir_base="/vol/mauranolab/aggregations/"
    cd ${dir_base}
    dir_names=($(ls -d */))     # These have trailing slashes.
    suffix="readcounts.summary.txt"
    find_readcounts ${dir_base} ${suffix} "${dir_names[@]}"
    
    # Maurano publicdata
    dir_base="/vol/mauranolab/publicdata/"
    cd ${dir_base}
    dir_names=($(ls -d */))     # These have trailing slashes.
    suffix="readcounts.summary.txt"
    find_readcounts ${dir_base} ${suffix} "${dir_names[@]}"
fi

###########################################################################
# Move the html files to trackhub_dev, and rename them.
# Ultimate file names are in the form:  <date>_<flowcell>.html
# Files reside in:  .../trackhub_dev/<genome>_descriptions/

# This function does the moving and renaming.  It is called once for each genome.
move_html () {
    local genome=$1
    local suffix="*_"${genome}".html"
    local fnames=($(ls *${suffix}))
    
    # mkdir /vol/cegs/public_html/trackhub_dev/${genome}/descriptions
    mkdir ${hub_target}/${genome}/descriptions
    
    local i
    for i in "${fnames[@]}"; do
        # Put a line break at the top of the HTML file, so the table does not 
        # overlap the horizontal line placed by the browser.
        echo "<br>" > ${i}".tmp"
        
        # Add the HTML file
        cat ${i} >> ${i}".tmp"
        
        local j=${i%${suffix}}
        mv ${i}".tmp" ${hub_target}/${genome}/descriptions/${j}.html
    done
}
###########################################################################

# Move to where the description tables have been stored.
cd ${TMPDIR}

# Now move the files into the hub.
for genome in "${genome_array[@]}"; do
    move_html ${genome}
done

###########################################################################
echo "Done processing Description files."

