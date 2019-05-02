#!/bin/bash
##############################################################################
# Pass in the input parameters:

path_to_main_driver_script=$1
hub_type=$2
hub_target=$3

# We need a tmp directory to store intermediate files.
# One was establihed in the calling script.
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

    # If only the header is there now, that means we have no data for that genome.
    local n_lines=$(wc -l < ${infile})

    if (( n_lines == 1)); then
        rm ${infile}
    else
        # Make the html version of the genome's readcounts data.
        Rscript ${path_to_main_driver_script}/makeDescHtml.R ${infile} "${infile}.html"
    fi
}
##############################################################################
# This function makes html versions of the readcounts files, by genome.

find_readcounts () {
    local dir_base=$1
    local suffix=$2
    shift 2
    local dir_names=("$@")

    local flowcell
    local ymd

    # Each flowcell directory may be associated with a readcounts file.
    # Get the data from the readcounts file if it is there, and make a table from the data.
    for flowcell in "${dir_names[@]}"; do
        # No readcounts files exist in these two directories.
        # If encountered, the subsequent $target test will fail, and we move on to the next flowcell.
        #     /vol/mauranolab/mapped/src/
        #     /vol/mauranolab/mapped/trackhub/

        # This is the name of the flowcell's readcounts file, if it exists.
        local target=${dir_base}${flowcell}${suffix}

        if [ -f ${target} ]; then
            # There is a readcounts file in the flowcell directory.
            echo ${target}

            # Find a date to associate with the flowcell from the info.txt file
            local info_file="/vol/mauranolab/flowcells/data/"${flowcell}"info.txt"

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
                local d_out=$(grep '#Load date' ${info_file})
                local yyyy=$(echo $d_out | cut -d' ' -f3 | cut -d'-' -f1)
                local mm=$(echo $d_out | cut -d' ' -f3 | cut -d'-' -f2)
                local dd=$(echo $d_out | cut -d' ' -f3 | cut -d'-' -f3)
                ymd=${yyyy}${mm}${dd}
            fi

            # We now have a date. Next, construct a genome related filename for the output.
            local tmp_flowcell
            if [ "${ymd}" = "00000000" ]; then
                tmp_flowcell="/"${flowcell%/}
            else
                tmp_flowcell="/"${ymd}"_"${flowcell%/}
            fi
            # tmp_flowcell will be used after the HEREDOC sections below.


            # Initialize the output files with header lines.
               # Also replace header underscores with a space. But R's read.table function bashes the space with a . !!!
               # We use this later, to replace the .'s with a space, but not touch the underscores in the html table tags.
               # See the move_html function below, where this is done.
               #
               # This HEREDOC also surrounds field values with quotes, to retain the spaces we are making.
            awk -f <(cat << "AWK_HEREDOC_01"
            BEGIN { FS = "\t" ; q="\"";  OFS = q FS q }
            {
                if(NR == 2) exit;

                $1 = $1;  # Causes awk to remove excess white space from all the fields. Necessary!
                gsub("_"," ", $0)
                print q $0 q;
            }
            END{}
AWK_HEREDOC_01
) < ${target} > "${TMPDIR}/tmp_header"

            # Use the header we just made to initialize the genome specific output files.
            # Break up the readcounts files into their genomic subsections.
            for i in "${genome_array[@]}"; do
                declare "${i}"="${TMPDIR}${tmp_flowcell}_${i}"
                ref="${i}"
                cat "${TMPDIR}/tmp_header" > ${!ref}
                grep ${i} ${target} >> ${!ref}
            done

            for i in "${genome_array[@]}"; do
                ref="${i}"
                make_html ${!ref}
            done
        else
            # No readcounts file exists. UCSC browser description area is just
            # blank if there is no file, or a bad url. No error msg pops up.
            echo ${target} does not exist
        fi
    done
}
##############################################################################
# Make the html files.
if [ "${hub_type}" = "CEGS" ]; then
    # CEGS flowcells
    dir_base="/vol/cegs/mapped/"
    cd ${dir_base}
    dir_names=($(ls -d */))    # These have trailing slashes.
    suffix="Project_CEGS/readcounts.summary.txt"
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
    # suffix="Project_Maurano/readcounts.summary.txt"
    suffix="dnase/readcounts.summary.txt"
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
echo Done processing Description files.

