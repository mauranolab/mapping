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

find_readcounts () {
    local dir_base=$1
    local suffix=$2
    
    if [ ! -d "${dir_base}" ] ; then
        return 0
    fi
    
    local flowcell
    local BASE
    local BASE2
    local BASEwoSlash
    local target
    local desc_file
    local partialHtmlName
    local HtmlName
    local info_file
    local infile
    
    
    # Each flowcell directory may be associated with a readcounts file.
    # Get the data from the readcounts file if it is there, and make a table from the data.
    for flowcelldir in `find ${dir_base} -mindepth 1 -maxdepth 1 -type d`; do
        flowcell=`basename ${flowcelldir}`
        # Fill up target_vec with directory names that house readcounts.summary.txt:
        if [[ "${dir_base}" = */mapped/ ]]; then
            mapfile -t target_vec < <( find ${dir_base}${flowcell} -maxdepth 2 \( -name "trash*" -o -name "bak" \) -prune -o -type f -name "readcounts.summary.txt" -print )
        else
            # For non-flowcell directories, there is only one place for the readcounts file to be:
            mapfile -t target_vec < <( find ${dir_base}${flowcell} -maxdepth 1 -type f -name "readcounts.summary.txt" )
        fi
        
        # Find a date to associate with the flowcell from the info.txt file
        info_file="/vol/mauranolab/flowcells/data/${flowcell}/info.txt"
        
        # Make sure the file is there:
        local yyyymmdd="00000000"
        partialHtmlName="/${flowcell%/}"             # partialHtmlName will be used in the for loop below.
        if [ -f ${info_file} ]; then
            # It is, but is the date there?
            ierr=$(grep -c '#Load date' ${info_file})
            if (( ierr != 0 )); then
                # Find a date to associate with the flowcell from the info.txt file
                yyyymmdd=$(grep '#Load date' ${info_file} | cut -f2 | perl -pe 's/\-//g;')  # Relies on tab delimited field split.
                partialHtmlName="/${yyyymmdd}_${flowcell%/}"
            fi
        fi

        numElements="${#target_vec[@]}"
        if [ "${numElements}" = "0" ]; then
            echo "[makeDescFiles.bash] WARNING No subdirectories with a readcounts file in ${dir_base}${flowcell}"
            # But print the data directory path:
            
            subdir_names=(`find ${flowcelldir} -mindepth 1 -maxdepth 1 -type d`)
            for j in "${subdir_names[@]}"; do
                BASEwoSlash=${j%/}        # Strip trailing slash
                BASE=${BASEwoSlash##*/}   # Get subdir name
                
                for i in "${genome_array[@]}"; do
                    echo "<pre>" > "${TMPDIR}${partialHtmlName}_${BASE}_${i}.html"
                    echo "Source data located in: ${dir_base}${flowcell}" >> "${TMPDIR}${partialHtmlName}_${BASE}_${i}.html"
                    echo " " >> "${TMPDIR}${partialHtmlName}_${BASE}_${i}.html"
                    echo "</pre>" >> "${TMPDIR}${partialHtmlName}_${BASE}_${i}.html"
                    echo "<hr>" >> "${TMPDIR}${partialHtmlName}_${BASE}_${i}.html"
                done
            done
            continue
        fi
        

        for target in ${target_vec[@]}; do
            htmlName=${partialHtmlName}

            if [[ "${dir_base}" = */mapped/ ]]; then
                # The name of the flowcell subdirectory that the readcounts.summary.txt file lives in is also the name of the relevant assay type.
                # Use that assay type to help make a unique html name.
                BASE=${target%/readcounts.summary.txt}       # Kill trailing /readcounts.summary.txt
                BASE2=${BASE##*/}  # Get subdir name
                htmlName="${htmlName}_${BASE2}"
            fi
            
            # Initialize the output files with header lines.
            head -n1 ${target} > "${TMPDIR}/tmp_header"
            
            # Use the header we just made to initialize the genome specific output files.
            # Break up the readcounts files into their genomic subsections.
            for i in "${genome_array[@]}"; do
                declare "${i}"="${TMPDIR}${htmlName}_${i}"
                ref="${i}"
                cat "${TMPDIR}/tmp_header" > ${!ref}
                grep ${i} ${target} >> ${!ref}
                
                BASE=${target%readcounts.summary.txt}       # Kill trailing readcounts.summary.txt
                desc_file="${BASE}description.html"
                
                declare "${i}_desc"="${TMPDIR}${htmlName}_${i}_desc"
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
        done
    done
}
##############################################################################
# Make the html files.
if [[ "${hub_type}" == "CEGS" ]]; then
    # CEGS flowcells
    hub_basedir="/vol/cegs"
elif [[ "${hub_type}" == "MAURANOLAB" ]]; then
    # Maurano
    hub_basedir="/vol/mauranolab"
elif [[ "${hub_type}" == "SARS" ]]; then
    # SARS
    hub_basedir="/vol/mauranolab/sars"
fi

# flowcells
find_readcounts ${hub_basedir}/mapped/ NA

# aggregations
find_readcounts ${hub_basedir}/aggregations/ readcounts.summary.txt

# publicdata
find_readcounts ${hub_basedir}/publicdata/ readcounts.summary.txt

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
