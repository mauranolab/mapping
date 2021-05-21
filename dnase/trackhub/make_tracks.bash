#!/bin/bash
set -eu -o pipefail
##############################################################################
# Pass in the input parameters:

# The temporary output file directory. It's used for the various tmp working files.
TMP_OUT=$1

# hub_type: CEGS has additional features
hub_type=$2

src=$3

assemblyBaseDir=$4

hub_target=$5

shift 5
genome_array=("$@")

##############################################################################
# This is the samplesforTrackhub.R section for "flowcell" samples.
# It finds the samples within each flowcell directory, and dumps relevant 
# sample information into an output file.

# Here we define the output file for samplesforTrackhub.R
# "outfile_base" is used later to construct file names for extracts of "outfile".
outfile_base="${TMP_OUT}/samplesforTrackhub"
outfile="${outfile_base}_all.tsv"


echo
Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile} --workingDir ${assemblyBaseDir}/mapped --descend --project byFC --quiet

# Split up the samplesforTrackhub.R output into separate files for each genome.
# BUGBUG It is important that the cegsvectors assemblies do not unintentionally include standard genome names in their own names.  For example, if a Mapped_genome is cegsvectors_Myspecialmm10gene, then the track will appear in both the cegsvectors genome files and in the mm10 files.
for genome in "${genome_array[@]}"; do
    mlr --tsv filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} > ${outfile_base}_${genome}_consolidated.tsv
done

# Lastly, a "header" is needed below for various aggregation output files.
# Save it now from the Rscript outfile.
head -n 1 ${outfile} > "${TMP_OUT}/header"

# This is the end of the samplesforTrackhub.R section for "flowcell" samples.

###########################################################################
# CEGS_byLocus section:

if [ ${hub_type} = "CEGS" ]; then
    echo
    Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile} --workingDir ${assemblyBaseDir}/mapped --descend --project CEGS_byLocus --quiet
    
    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for genome in "${genome_array[@]}"; do
        mlr --tsv filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} > ${outfile_base}_${genome}_consolidated_locus.tsv
    done
fi

# This is the end of the samplesforTrackhub.R section for making the CEGS_byLocus files.

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
    
    Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile} --workingDir ${workingDir} ${inputfile} --quiet
    
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
    
    cat "${TMP_OUT}/header" "${TMP_OUT}/tmp" > ${outfile}
    rm "${TMP_OUT}/tmp"
    # Adjustment of outfile columns complete.
    ###
    
    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for genome in "${genome_array[@]}"; do
        # Note that prior to entering this function, "outfile" was set to be: "${outfile_base}_all_agg.tsv".
        # This never changes, so each call to this function over-writes the previous "outfile".
        # However, the output from the below call to mlr is APPENDED to previous output from calls
        # to this function. So the genome specific output files get bigger as we do more aggregation directories,
        # and call this function for each one of them.
        
        if [ "${loop_type}" = "aggregations" ]; then
            # aggregations
            mlr --tsv --headerless-csv-output filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} >> ${outfile_base}_${genome}_consolidated_agg.tsv
        else
            # publicdata
            mlr --tsv --headerless-csv-output filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} >> ${outfile_base}_${genome}_consolidated_pub.tsv
        fi
    done
    echo
}
# End of the agg_pub_loop function section.

###########################################################################
# Back in the main line of code.
###########################################################################
# aggregations:

# Initialize aggregation output files
for genome in "${genome_array[@]}"; do
    cat "${TMP_OUT}/header" > "${outfile_base}_${genome}_consolidated_agg.tsv"
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
    for genome in "${genome_array[@]}"; do
        cat "${TMP_OUT}/header" > "${outfile_base}_${genome}_consolidated_pub.tsv"
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
    local consolidated_suffix=$2
    local supertrack=$3
    
    local infile="${TMP_OUT}/samplesforTrackhub_${mappedgenome}${consolidated_suffix}.tsv"
    
    num_line=`(wc -l < ${infile})`
    if [ ${num_line} -le 1 ]; then
        # There are no tracks of this type
        return 0
    fi
    
    local tracknameprefix=""
    local generateHTMLdescription="--generateHTMLdescription"
    local includeSampleIDinSampleCol=""
    subgroupprefix="--subgroupnames Replicate"
    if [ ${supertrack} = "By_Locus" ]; then
        tracknameprefix="--tracknameprefix byLocus"
        subgroupprefix="--subgroupnames Project,Assembly,Type"
        generateHTMLdescription=""
        supertrackPriority=20
        URLbase="../mapped/"
    elif [ ${supertrack} = "Aggregations" ]; then
        supertrackPriority=30
        URLbase="../aggregations/"
    elif [ ${supertrack} = "Public_Data" ]; then
        supertrackPriority=40
        URLbase="../publicdata/"
    elif [ ${supertrack} = "Flowcells" ]; then
        supertrackPriority=50
        URLbase="../mapped/"
        includeSampleIDinSampleCol="--includeSampleIDinSampleCol"
    else
        echo "ERROR impossible"
        exit 1
    fi
    
    python ${src}/MakeTrackhub.py ${infile} ${generateHTMLdescription} ${includeSampleIDinSampleCol} ${tracknameprefix} ${subgroupprefix} --supertrack ${supertrack} --supertrackPriority ${supertrackPriority} --genome ${mappedgenome} --checksamples --URLbase ${URLbase} > ${hub_target}/${mappedgenome}/trackDb.${supertrack}.txt
}

for genome in "${genome_array[@]}"; do
    make_tracks ${genome} "_consolidated" "Flowcells"
done

for genome in "${genome_array[@]}"; do
    make_tracks ${genome} "_consolidated_agg" "Aggregations"
done

if [ "${hub_type}" != "SARS" ]; then
    for genome in "${genome_array[@]}"; do
        make_tracks ${genome} "_consolidated_pub" "Public_Data"
    done
fi

if [ ${hub_type} = "CEGS" ]; then
    for genome in "${genome_array[@]}"; do
        make_tracks ${genome} "_consolidated_locus" "By_Locus"
    done
fi


echo
echo "[make_tracks] Done"

