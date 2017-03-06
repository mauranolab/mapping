#!/bin/bash
set -e # -o pipefail

#qsub -S /bin/bash -terse -j y -N Job -o ${sample} -b y "~/scratch/transposon/src/extractBCcounts.sh ${sample} $amplicon $bclen $bcread $chunksize" | perl -pe 's/[^\d].+$//g;' > sgeid.${sample}
#
python -m cProfile -s time ../src/DeDup_adjacencyOutput.py --col 1 -o output/test_$1 $1 >cProfile/test_$1
#qsub name J
#
#qsub -S /bin/bash -t 1-${numjobs} -terse -j y -N extract.${sample} -o ${sample} -b y "~/scratch/transposon/src/extractBCcounts.sh ${sample} $amplicon $bclen $bcread $chunksize" | perl -pe 's/[^\d].+$//g;' > sgeid.${sample}
#
for i in `ls Data/`; do  qsub -S /bin/bash -terse -j y -N JobIds/$i -o cProfile/$i -b y 'python -m cProfile -s time ../src/DeDup_adjacencyOutput.py --col 1 -o output/test_'$i' Data/'$i'' ;done

#for i in {1000,10000,100000}; do  qsub -S /bin/bash -terse -j y -N JobIds/test_${i} -o cProfile/test_${i} -b y head -$i BS00067A-5xBGlo_K562d4_2hDpn_iPCR.barcodes.txt| py -m cProfile -s time ../src/DeDup_adjacencyOutput.py --col 1 -o output/test${i} -;done

#edit_distance 
for i in `ls Data/`; do  qsub -S /bin/bash -t 8 -terse -j y -N JobIds/Edit_$i -o cProfile/Edit_$i -b y 'python -m cProfile -s time ../src/DeDup_adjacencyOutput_edit_distance.py --col 1 -o output/Edit_'$i' Data/'$i'' ;done
for i in `ls Data/`; do  qsub -S /bin/bash -t 8 -terse -j y -N Edit2_$i -o cProfile/Edit2_$i -b y 'python -m cProfile -s time ../src/DeDup_adjacencyOutput_edit_distance.py --col 1 -o output/Edit2_'$i' Data/'$i'' ;done
#edit_distance + swap
for i in `ls Data/`; do  qsub -S /bin/bash -t 8 -terse -j y -N JobIds/EditSwap_$i -o cProfile/EditSwap_$i -b y 'python -m cProfile -s time ../src/DeDup_adjacencyOutput_edit_distance_Swap.py --col 1 -o output/EditSwap_'$i' Data/'$i'' ;done
#Levenstein applied directly 
for i in `ls Data/`; do  qsub -S /bin/bash -t 8 -terse -j y -N JobIds/Leven_$i -o cProfile/Leven_$i -b y 'python -m cProfile -s time ../src/DeDup_adjacencyOutput_SpeedUP.py --col 1 -o output/Leven_'$i' Data/'$i'' ;done