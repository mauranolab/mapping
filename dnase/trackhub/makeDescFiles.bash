#!/bin/bash
set -eu -o pipefail

# Pass in the input parameters:
src=$1
assemblyBaseDir=$2
hub_target=$3
TMPDIR=$4
shift 4
#shouldn't this be able to infer genomes from what's in readcounts.summary.txt?
genome_array=("$@")


# These functions makes html versions of the readcounts files, by genome.
find_readcounts () {
    local dir_base=$1
    
    echo "[makeDescFiles] Looking for readcounts.summary.txt files in ${dir_base}"
    
    if [ ! -d "${dir_base}" ] ; then
        return 0
    fi
    
    local flowcell
    local BASE
    local readcountsFile
    local desc_file
    local partialHtmlName
    local HtmlName
    local info_file
    
    #Has already been made by makeAssemblyTracks.bash
    #mkdir -p ${TMPDIR}/descriptions/
    
    # Each flowcell directory may be associated with a readcounts file.
    # Get the data from the readcounts file if it is there, and make a table from the data.
    for flowcelldir in `find ${dir_base} -mindepth 1 -maxdepth 1 -type d`; do
        # Fill up target_vec with directory names that house readcounts.summary.txt:
        if [[ `basename ${dir_base}` == "mapped" ]]; then
            target_vec=`find ${flowcelldir} -maxdepth 2 \( -name "trash*" -o -name "bak" \) -prune -o -type f -name "readcounts.summary.txt" -print`
        else
            # For non-flowcell directories, there is only one place for the readcounts file to be:
            target_vec=`find ${flowcelldir} -maxdepth 1 -type f -name "readcounts.summary.txt"`
        fi
        
        
        # Find a date to associate with the flowcell from the info.txt file
        flowcell=`basename ${flowcelldir}`
        info_file="/vol/mauranolab/flowcells/data/${flowcell}/info.txt"
        # Make sure the file is there:
        local yyyymmdd="00000000"
        partialHtmlName="/${flowcell}"             # partialHtmlName will be used in the for loop below.
        if [ -f ${info_file} ]; then
            # It is, but is the date there?
            ierr=$(grep -c '#Load date' ${info_file})
            if (( ierr != 0 )); then
                # Find a date to associate with the flowcell from the info.txt file
                yyyymmdd=$(grep '#Load date' ${info_file} | cut -f2 | perl -pe 's/\-//g;')  # Relies on tab delimited field split.
                partialHtmlName="/${yyyymmdd}_${flowcell}"
            fi
        fi
        
#        # Just print the data directory path if there is no readcounts file
#        #BUGBUG are the results ever moved out of TMPDIR?
#        if [ "${target_vec}" = "" ]; then
#            echo "[makeDescFiles.bash] WARNING No subdirectories with a readcounts file in ${flowcelldir}"
#            for BASE in `find ${flowcelldir} -mindepth 1 -maxdepth 1 -type d | xargs -I {} basename {}`; do
#                for curGenome in "${genome_array[@]}"; do
#                    echo "<pre>" > ${TMPDIR}${partialHtmlName}_${BASE}_${curGenome}.html
#                    echo "Source data located in: ${flowcelldir}" >> ${TMPDIR}${partialHtmlName}_${BASE}_${curGenome}.html
#                    echo " " >> ${TMPDIR}${partialHtmlName}_${BASE}_${curGenome}.html
#                    echo "</pre>" >> ${TMPDIR}${partialHtmlName}_${BASE}_${curGenome}.html
#                    echo "<hr>" >> ${TMPDIR}${partialHtmlName}_${BASE}_${curGenome}.html
#                done
#            done
#            continue
#        fi
        
        for readcountsFile in ${target_vec}; do
            BASE=${readcountsFile%readcounts.summary.txt}       # Kill trailing readcounts.summary.txt
            htmlName="${partialHtmlName}"
            if [[ `basename ${dir_base}` == "mapped" ]]; then
                # The name of the flowcell subdirectory that the readcounts.summary.txt file lives in is also the name of the relevant assay type.
                # Use that assay type to help make a unique html name.
                # Use subdir name from BASE
                htmlName="${htmlName}_"`basename ${BASE}`
            fi
            
            # Break up the readcounts files into their genomic subsections.
            for curGenome in "${genome_array[@]}"; do
                # Initialize the output files with header lines
                head -n1 ${readcountsFile} > ${TMPDIR}/descriptions/${htmlName}_${curGenome}
                grep ${curGenome} ${readcountsFile} >> ${TMPDIR}/descriptions/${htmlName}_${curGenome} || true  # Avoid pipefail crash when there are no readcounts lines for this genome.
                
                #BUGBUG makes placeholder file for every genome in genome_array, even if there is no data mapped there
                echo "<pre>" > ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                if [ -f ${BASE}/description.html ]; then
                    cat ${BASE}/description.html >> ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                    echo " " >> ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                fi
                echo "Source data located in: ${readcountsFile%/readcounts.summary.txt}" >> ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                echo " " >> ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                echo "</pre>" >> ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                echo "<hr>" >> ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                
                # If only the header is there now, that means we have no data for that genome.
                if [[ `wc -l < ${TMPDIR}/descriptions/${htmlName}_${curGenome}` > 1 ]]; then
                    # Make the html version of the genome's readcounts data.
                    #Since this gets called frequently, use --vanilla to try to speed up
                    Rscript --vanilla ${src}/makeDescHtml.R ${TMPDIR}/descriptions/${htmlName}_${curGenome} ${TMPDIR}/descriptions/make_html_readcounts_tmp.html
                    cat ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html ${TMPDIR}/descriptions/make_html_readcounts_tmp.html > ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html.new
                    mv ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html.new ${TMPDIR}/descriptions/${htmlName}_${curGenome}.html
                fi
            done
        done
    done
}

for dir in mapped aggregations publicdata; do
    find_readcounts ${assemblyBaseDir}/${dir}
done


###########################################################################
# Move the html files to trackhub_dev, and rename them.
# Ultimate file names are in the form:  <date>_<flowcell>.html
# Files reside in:  .../trackhub_dev/<genome>_descriptions/
echo "[makeDescFiles] Combining descriptions by genome"
for curGenome in "${genome_array[@]}"; do
    mkdir ${hub_target}/${curGenome}/descriptions
    #Errors if nothing exists in ${TMPDIR}/*_${curGenome}
    for descfile in `find ${TMPDIR}/descriptions -name "*_${curGenome}.html"`; do
        outfile=`basename ${descfile} _${curGenome}.html`
        # Put a line break at the top of the HTML file, so the table does not overlap the horizontal line placed by the browser
        echo "<br>" > ${hub_target}/${curGenome}/descriptions/${outfile}.html
        cat ${descfile} >> ${hub_target}/${curGenome}/descriptions/${outfile}.html
    done
done

###########################################################################
echo "[makeDescFiles] Done"
