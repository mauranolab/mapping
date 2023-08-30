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


echo
echo "[make_tracks] generating sample list for Flowcells"
outfile_Flowcells="${TMP_OUT}/samplesforTrackhub_Flowcells.tsv"
Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile_Flowcells} --workingDir ${assemblyBaseDir}/mapped --descend --project byFC --quiet


# Split up the samplesforTrackhub.R output into separate files for each genome.
# BUGBUG It is important that the cegsvectors assemblies do not unintentionally include standard genome names in their own names.  For example, if a Mapped_genome is cegsvectors_Myspecialmm10gene, then the track will appear in both the cegsvectors genome files and in the mm10 files.
for curGenome in "${genome_array[@]}"; do
    mlr --tsv filter -S "'\$Mapped_Genome =~ \".*${curGenome}.*\"'" ${outfile_Flowcells} > ${TMP_OUT}/samplesforTrackhub_${curGenome}_Flowcells.tsv
done



if [[ "${hub_type}" == "CEGS" ]]; then
    echo
    echo "[make_tracks] generating sample list for By_Locus"
    outfile_byLocus="${TMP_OUT}/samplesforTrackhub_By_Locus.tsv"
    Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile_byLocus} --workingDir ${assemblyBaseDir}/mapped --descend --project CEGS_byLocus --quiet
    
    # Split up the samplesforTrackhub.R output into separate files for each genome.
    for curGenome in "${genome_array[@]}"; do
        mlr --tsv filter -S "'\$Mapped_Genome =~ \".*${curGenome}.*\"'" ${outfile_byLocus} > ${TMP_OUT}/samplesforTrackhub_${curGenome}_By_Locus.tsv
    done
fi


# a "header" is needed below for various aggregation output files. Take it from the outfile_Flowcells output from above.
head -n 1 ${outfile_Flowcells} > ${TMP_OUT}/header_Flowcells
# obtain number of the header column for: "filebase_col" and "group_col" for loop below.
filebase_col=$(awk -v colname="filebase" -F "\t" '{for (i=1; i<=NF; i++) {if($i == colname) {print i; exit}}}' < ${TMP_OUT}/header_Flowcells)
group_col=$(awk -v colname="Group" -F "\t" '{for (i=1; i<=NF; i++) {if($i == colname) {print i; exit}}}' < ${TMP_OUT}/header_Flowcells)
echo
echo "[make_tracks] filebase_col is ${filebase_col}; group_col is ${group_col}"

# Call samplesforTrackhub.R for aggregations or publicdata, and split the output by genome.
for curbase in aggregations publicdata; do
    if [ -d "${assemblyBaseDir}/${curbase}/" ]; then
        echo
        echo "[make_tracks] generating sample list for ${curbase}"
        # Initialize aggregation output files
        for curGenome in "${genome_array[@]}"; do
            cat ${TMP_OUT}/header_Flowcells > ${TMP_OUT}/samplesforTrackhub_${curGenome}_${curbase}.tsv
        done
        outfile_aggpub="${TMP_OUT}/samplesforTrackhub_all_${curbase}.tsv"
        
        # Call samplesforTrackhub.R for each directory.
        for dir_loop_name in `find ${assemblyBaseDir}/${curbase} -mindepth 1 -maxdepth 1 -type d | xargs -I {} basename {}`; do
            workingDir="${assemblyBaseDir}/${curbase}/${dir_loop_name}"
            
            #Pull in annotation if available
            inputfile=""
            if [ -s "${workingDir}/sampleannotation.txt" ]; then
                inputfile="--inputfile ${workingDir}/sampleannotation.txt"
            fi
            
            Rscript --vanilla ${src}/samplesforTrackhub.R --out ${outfile_aggpub} --workingDir ${workingDir} ${inputfile} --quiet
            
            # Make some adjustments to the "outfile_aggpub" columns, as the aggregations are structured somewhat differently than the flowcell versions
            cat ${outfile_aggpub} | awk -v filebase_col=${filebase_col} -v group_col=${group_col} -v dir_name=${dir_loop_name} -F "\t" 'BEGIN {OFS="\t"; dir_name2=sprintf("%s%s", dir_name, "/")} \
            NR>1 { \
                sub(/^/, dir_name2, $filebase_col); \
                sub(/NA/, dir_name, $group_col); \
                print; \
            }' | cat ${TMP_OUT}/header_Flowcells - > ${outfile_aggpub}.new
            mv ${outfile_aggpub}.new ${outfile_aggpub}
            
            # Split up the samplesforTrackhub.R output into separate files for each genome and append output for this dir_loop_name.
            for curGenome in "${genome_array[@]}"; do
                mlr --tsv --headerless-csv-output filter -S "'\$Mapped_Genome =~ \".*${curGenome}.*\"'" ${outfile_aggpub} >> ${TMP_OUT}/samplesforTrackhub_${curGenome}_${curbase}.tsv
            done
            echo
        done
    fi
done


echo "[make_tracks] running MakeTrackhub"
make_tracks () {
    local mappedgenome=$1
    local track_type=$2
    local supertrack_name=$3
    
    local infile="${TMP_OUT}/samplesforTrackhub_${mappedgenome}${track_type}.tsv"
    
    if [ ! -f "${infile}" ] ||  [ `wc -l < ${infile}` -le 1 ]; then
        # There are no tracks of this type
        echo "[make_tracks] no tracks found for ${infile}"
        return 0
    fi
    
    local tracknameprefix=""
    local generateHTMLdescription="--generateHTMLdescription"
    local includeSampleIDinSampleCol=""
    subgroupprefix="--subgroupnames Replicate"
    if [ ${supertrack_name} = "By_Locus" ]; then
        tracknameprefix="--tracknameprefix Locus"
        subgroupprefix="--subgroupnames Project,Assembly,Type"
        generateHTMLdescription=""
        supertrackPriority=20
        URLbase="../mapped/"
    elif [ ${supertrack_name} = "Aggregations" ]; then
        tracknameprefix="--tracknameprefix Agg"
        supertrackPriority=30
        URLbase="../aggregations/"
    elif [ ${supertrack_name} = "Public_Data" ]; then
        supertrackPriority=40
        URLbase="../publicdata/"
    elif [ ${supertrack_name} = "Flowcells" ]; then
        supertrackPriority=50
        URLbase="../mapped/"
        includeSampleIDinSampleCol="--includeSampleIDinSampleCol"
    else
        echo "ERROR impossible"
        exit 1
    fi
    
    #BUGBUG I think this outputs in the order of ${infile}, which is only sorted by sample ID. Perhaps it should be sorted by data type too?
    python ${src}/MakeTrackhub.py ${infile} ${generateHTMLdescription} ${includeSampleIDinSampleCol} ${tracknameprefix} ${subgroupprefix} --supertrack ${supertrack_name} --supertrackPriority ${supertrackPriority} --genome ${mappedgenome} --checksamples --URLbase ${URLbase} > ${hub_target}/${mappedgenome}/trackDb.${supertrack_name}.txt
}

for curGenome in "${genome_array[@]}"; do
    make_tracks ${curGenome} "_Flowcells" "Flowcells"
done

for curGenome in "${genome_array[@]}"; do
    make_tracks ${curGenome} "_aggregations" "Aggregations"
done

for curGenome in "${genome_array[@]}"; do
    make_tracks ${curGenome} "_publicdata" "Public_Data"
done

for curGenome in "${genome_array[@]}"; do
    make_tracks ${curGenome} "_By_Locus" "By_Locus"
done


echo
echo "[make_tracks] Done"

