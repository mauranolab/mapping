#!/bin/bash
##############################################################################
# Pass in the input parameters:

# The temporary output file directory. It's used for the various tmp working files.
TMP_OUT=$1

# hub_type: CEGS has additional features
hub_type=$2

path_to_main_driver_script=$3

assemblyBaseDir=$4

URLbase=$5

shift 5
genome_array=("$@")

##############################################################################
# This is the samplesforTrackhub.R section for "flowcell" samples.
# It finds the samples within each flowcell directory, and dumps relevant 
# sample information into an output file.

# Here we define the output file for samplesforTrackhub.R
# "outfile_base" is used later to construct file names for extracts of "outfile".
outfile_base=${TMP_OUT}"/samplesforTrackhub"
outfile=${outfile_base}"_all.tsv"


echo
Rscript --vanilla ${path_to_main_driver_script}/samplesforTrackhub.R \
        --out ${outfile} \
        --workingDir ${assemblyBaseDir}/mapped \
        --descend \
        --project byFC \
        --quiet

# Split up the samplesforTrackhub.R output into separate files for each genome.
for i in "${genome_array[@]}"; do
    declare "outfile_${i}"="${outfile_base}_${i}_consolidated.tsv"
    ref="outfile_${i}"
    head -n 1 ${outfile} > ${!ref}
    grep ${i} ${outfile} >> ${!ref}
done

# Lastly, a "header" is needed below for various aggregation output files.
# Save it now from the Rscript outfile.
head -n 1 ${outfile} > "${TMP_OUT}/header"

# This is the end of the samplesforTrackhub.R section for "flowcell" samples.
#
###########################################################################
# CEGS_byLocus section:

if [ ${hub_type} = "CEGS" ]; then
    echo
    Rscript --vanilla ${path_to_main_driver_script}/samplesforTrackhub.R \
            --out ${outfile} \
            --workingDir ${assemblyBaseDir}/mapped \
            --descend \
            --project CEGS_byLocus \
            --quiet
    
    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for i in "${genome_array[@]}"; do
        declare "outfile_${i}"="${outfile_base}_${i}_consolidated_locus.tsv"
        ref="outfile_${i}"
        head -n 1 ${outfile} > ${!ref}
        grep ${i} ${outfile} >> ${!ref}
    done
fi

# This is the end of the samplesforTrackhub.R section for making the CEGS_byLocus files.
#
###########################################################################
###########################################################################
# This begins the samplesforTrackhub.R section for aggregation & publicdata
# samples. It finds the samples within each directory, and dumps relevant
# sample information into an output file.


# First, obtain two parameters: "filebase_col" and "Group_col".
# They are needed in agg_pub_loop function below.

# Get number of the header column containing "filebase"
filebase_col=$(awk -f <(cat << "AWK_HEREDOC_01"
BEGIN{FS="\t"}
{
    for (i=1; i<=NF; i++)
    {
        if($i == "filebase") {print i; exit}
    }
}
AWK_HEREDOC_01
) < ${TMP_OUT}"/header")

echo
echo "filebase_col is ${filebase_col}"


# Get number of the header column containing "Group"
Group_col=$(awk -f <(cat << "AWK_HEREDOC_02"
BEGIN{FS="\t"}
{
    for (i=1; i<=NF; i++)
    {
        if($i == "Group") {print i; exit}
    }
}
AWK_HEREDOC_02
) < ${TMP_OUT}"/header")

echo "Group_col is ${Group_col}"
echo

###########################################################################
# Before moving on, below we define a few functions. There will be a comment
# when we return to the main line of code.
###########################################################################
# This function is called once for each aggregation/publicdata directory.
# It calls samplesforTrackhub.R for that directory, then splits the output by genome.

agg_pub_loop () {
    local dir_loop_name=$1
    
    # aggregations or publicdata
    local loop_type=$2
    
    local workingDir="${assemblyBaseDir}/${loop_type}/${dir_loop_name}"
    
    #Pull in annotation if available
    local inputfile=""
    if [ -s "${workingDir}/sampleannotation.txt" ]; then
        inputfile="--inputfile ${workingDir}/sampleannotation.txt"
    fi
    
    Rscript --vanilla ${path_to_main_driver_script}/samplesforTrackhub.R \
            --out ${outfile} \
            --workingDir ${workingDir} \
            ${inputfile} \
            --quiet

###
# Make some adjustments to the "outfile" columns, as the aggregations are
# structured somewhat differently than the flowcells, and the CEGS version 
# of samplesforTrackhub.R was written with the flowcells in mind.
awk -v f_col=${filebase_col} \
    -v G_col=${Group_col} \
    -v dir_name=${dir_loop_name} \
    -f <(cat << "AWK_HEREDOC_03"
BEGIN{FS="\t"; OFS="\t"; dir_name2=sprintf("%s%s", dir_name, "/")}
{
    if(NR == 1) next;

    sub(/^/, dir_name2, $f_col)
    sub(/NA/, dir_name, $G_col)
    print
}
AWK_HEREDOC_03
) < ${outfile} > "${TMP_OUT}/tmp"

    cp "${TMP_OUT}/tmp" ${outfile}
    rm "${TMP_OUT}/tmp"
# Adjustment of outfile columns complete.
###

    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for i in "${genome_array[@]}"; do
        # Note that prior to entering this function, "outfile" was set to be: "${outfile_base}_all_agg.tsv".
        # This never changes, so each call to this function over-writes the previous "outfile".
        # However, the output from the below call to grep is APPENDED to previous output from calls
        # to this function. So the genome specific output files get bigger as we do more aggregation directories,
        # and call this function for each one of them.

        ref="outfile_${i}"
        grep ${i} ${outfile} >> ${!ref}
    done
    echo
}
# End of the agg_pub_loop function section.

###########################################################################
# Back in the main line of code.
###########################################################################
# aggregations:

# Initialize aggregation output files
for i in "${genome_array[@]}"; do
    declare "outfile_${i}"="${outfile_base}_${i}_consolidated_agg.tsv"
    ref="outfile_${i}"
    cat "${TMP_OUT}/header" > "${!ref}"
done
outfile="${outfile_base}_all_agg.tsv"

# Get the names of the aggregation directories
cd ${assemblyBaseDir}/aggregations/
dir_names=($(ls -d */))

# Get rid of trailing slashes 
agg_names=()
for i in "${dir_names[@]}"; do
    x=$i
    y=${x%/}
    agg_names+=($y)
done

# Call samplesforTrackhub.R for each aggregation directory.
for i in "${agg_names[@]}"; do
    agg_pub_loop $i aggregations
done

# End of the samplesforTrackhub.R section for "aggregation" samples.

###########################################################################
# publicdata:
if [ "${hub_type}" != "SARS" ]; then
    # Initialize publicdata output files
    for i in "${genome_array[@]}"; do
        declare "outfile_${i}"="${outfile_base}_${i}_consolidated_pub.tsv"
        ref="outfile_${i}"
        cat "${TMP_OUT}/header" > "${!ref}"
    done
    outfile="${outfile_base}_all_pub.tsv"
    
    # Get the names of the publicdata directories
    cd ${assemblyBaseDir}/publicdata/
    dir_names=($(ls -d */))
    
    # Get rid of trailing slashes 
    pub_names=()
    for i in "${dir_names[@]}"; do
        x=$i
        y=${x%/}
        pub_names+=($y)
    done
    
    # Call samplesforTrackhub.R for each publicdata directory.
    for i in "${pub_names[@]}"; do
        agg_pub_loop $i publicdata
    done
fi
# End of the samplesforTrackhub.R section for "publicdata" samples.

###########################################################################
###########################################################################
# This is the start of the MakeTrackhub.py section.
# It calls the Daler track construction code to build the track hub.
# All output is placed in the temporary directory.

# This function will be called once for each genome.
make_tracks () {
    local mappedgenome=$1
    local consol_suffix=$2
    local URLbase=$3
    local supertrack=$4
    
    local mappedgenome_consol=${mappedgenome}${consol_suffix}
    local infile=${TMP_OUT}"/samplesforTrackhub_"${mappedgenome_consol}".tsv"
    local outfile=${TMP_OUT}"/MakeTrackhub_"${mappedgenome_consol}".out"
    
    local includeSampleIDinSampleCol=""
    if [ ${supertrack} != "Aggregations" ] && [ ${supertrack} != "Public_Data" ] && [ ${supertrack} != "By_Locus" ]; then
        includeSampleIDinSampleCol="--includeSampleIDinSampleCol"
    fi
    
    local tracknameprefix=""
    local generateHTMLdescription="--generateHTMLdescription"
    if [ ${supertrack} = "By_Locus" ]; then
        tracknameprefix="--tracknameprefix byLocus"
        subgroupprefix="--subgroupnames Project,Assembly,Type"
        generateHTMLdescription=""
    else
        subgroupprefix="--subgroupnames Replicate"
    fi
    
    if [ ${supertrack} = "By_Locus" ]; then
        supertrackPriority=20
    elif [ ${supertrack} = "Aggregations" ]; then
        supertrackPriority=30
    elif [ ${supertrack} = "Public_Data" ]; then
        supertrackPriority=40
    elif [ ${supertrack} = "Flowcells" ]; then
        supertrackPriority=50
    else
        supertrackPriority=99
    fi
    
    python ${path_to_main_driver_script}/MakeTrackhub.py ${infile} \
           ${generateHTMLdescription} \
           ${includeSampleIDinSampleCol} \
           ${tracknameprefix} \
           ${subgroupprefix} \
           --supertrack ${supertrack} \
           --supertrackPriority ${supertrackPriority} \
           --genome ${mappedgenome} \
           --checksamples \
           --URLbase ${URLbase} > ${outfile}
}
###############################
consol_suffix_in="_consolidated"
supertrack_in="Flowcells"

for i in "${genome_array[@]}"; do
    make_tracks ${i} ${consol_suffix_in} ${URLbase}/mapped/ ${supertrack_in}
done

###############################
consol_suffix_in="_consolidated_agg"
supertrack_in="Aggregations"

for i in "${genome_array[@]}"; do
    make_tracks ${i} ${consol_suffix_in} ${URLbase}/aggregations/ ${supertrack_in}
done

###############################
if [ "${hub_type}" != "SARS" ]; then
    consol_suffix_in="_consolidated_pub"
    supertrack_in="Public_Data"
    
    for i in "${genome_array[@]}"; do
        make_tracks ${i} ${consol_suffix_in} ${URLbase}/publicdata/ ${supertrack_in}
    done
fi

###############################
if [ ${hub_type} = "CEGS" ]; then
    consol_suffix_in="_consolidated_locus"
    supertrack_in="By_Locus"
    
    for i in "${genome_array[@]}"; do
        make_tracks ${i} ${consol_suffix_in} ${URLbase}/mapped/ ${supertrack_in}
    done
fi
###############################
echo
echo "Done with Daler python code."

