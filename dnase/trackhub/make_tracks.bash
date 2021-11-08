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
echo "[make_tracks] generating sample list"
Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile} --workingDir ${assemblyBaseDir}/mapped --descend --project byFC --quiet

# a "header" is needed below for various aggregation output files.
head -n 1 ${outfile} > ${TMP_OUT}/header

# Split up the samplesforTrackhub.R output into separate files for each genome.
# BUGBUG It is important that the cegsvectors assemblies do not unintentionally include standard genome names in their own names.  For example, if a Mapped_genome is cegsvectors_Myspecialmm10gene, then the track will appear in both the cegsvectors genome files and in the mm10 files.
for genome in "${genome_array[@]}"; do
    mlr --tsv filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} > ${outfile_base}_${genome}_consolidated.tsv
done



echo
echo "[make_tracks] generating sample list for CEGS_byLocus"
if [ ${hub_type} = "CEGS" ]; then
    Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile} --workingDir ${assemblyBaseDir}/mapped --descend --project CEGS_byLocus --quiet
    
    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for genome in "${genome_array[@]}"; do
        mlr --tsv filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} > ${outfile_base}_${genome}_consolidated_locus.tsv
    done
fi



echo
echo "[make_tracks] generating sample list for aggregations and public_data"
# First, obtain number of the header column for: "filebase_col" and "group_col" for agg_pub_loop function below.
filebase_col=$(awk -v colname="filebase" -F "\t" '{for (i=1; i<=NF; i++) {if($i == colname) {print i; exit}}}' < ${TMP_OUT}/header)
group_col=$(awk -v colname="Group" -F "\t" '{for (i=1; i<=NF; i++) {if($i == colname) {print i; exit}}}' < ${TMP_OUT}/header)
echo "[make_tracks] filebase_col is ${filebase_col}; group_col is ${group_col}"

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
    cat ${outfile} | awk -v filebase_col=${filebase_col} -v group_col=${group_col} -v dir_name=${dir_loop_name} -F "\t" 'BEGIN {OFS="\t"; dir_name2=sprintf("%s%s", dir_name, "/")} \
    NR>1 { \
        sub(/^/, dir_name2, $filebase_col); \
        sub(/NA/, dir_name, $group_col); \
        print; \
    }' | cat ${TMP_OUT}/header - > ${outfile}.new
    mv ${outfile}.new ${outfile}
    
    # Adjustment of outfile columns complete.
    ###
    
    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for genome in "${genome_array[@]}"; do
        # Note that prior to entering this function, "outfile" was set to be: "${outfile_base}_all_agg.tsv".
        # This never changes, so each call to this function over-writes the previous "outfile".
        # However, the output from the below call to mlr is APPENDED to previous output from calls
        # to this function. So the genome specific output files get bigger as we do more aggregation directories,
        # and call this function for each one of them.
        
        mlr --tsv --headerless-csv-output filter -S "'\$Mapped_Genome =~ \".*${genome}.*\"'" ${outfile} >> ${outfile_base}_${genome}_consolidated_${loop_type}.tsv
    done
    echo
}


for curbase in aggregations publicdata; do
    if [ -d "${assemblyBaseDir}/${curbase}/" ]; then
        # Initialize aggregation output files
        for genome in "${genome_array[@]}"; do
            cat ${TMP_OUT}/header > ${outfile_base}_${genome}_consolidated_${curbase}.tsv
        done
        outfile="${outfile_base}_all_agg.tsv"
        
        # Call samplesforTrackhub.R for each directory.
        for i in `find ${assemblyBaseDir}/${curbase} -mindepth 1 -maxdepth 1 -type d | xargs -I {} basename {}`; do
            agg_pub_loop $i ${curbase}
        done
    fi
done


echo "[make_tracks] running MakeTrackhub"
make_tracks () {
    local mappedgenome=$1
    local consolidated_suffix=$2
    local supertrack=$3
    
    local infile="${TMP_OUT}/samplesforTrackhub_${mappedgenome}${consolidated_suffix}.tsv"
    
    if [ ! -f "${infile}" ] ||  [ `wc -l < ${infile}` -le 1 ]; then
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
    make_tracks ${genome} "_consolidated_agggregations" "Aggregations"
done

for genome in "${genome_array[@]}"; do
    make_tracks ${genome} "_consolidated_publicdata" "Public_Data"
done

for genome in "${genome_array[@]}"; do
    make_tracks ${genome} "_consolidated_locus" "By_Locus"
done


echo
echo "[make_tracks] Done"

