######
## this script is to run all steps of pipeline to quantify barcodes and read counts for sgRNAs
######
## load mudules 
module load bedtools
module load samtools
module load fastx-toolkit/0.0.13;
module load bowtie/0.12.5
mkdir -p $PWD/logs

## define index, input and output
indexedGuides="30kguides" ## sgRNA index
RAW="CB5LRANXX_7_20170825B_20170827.bam";
#OUT="OUT_temp"
#mkdir -p $OUT

## file name
Sam="30k_lane8_bw1MM";
out_merge_bowtie_shortsam="out_0_30kHagar.txt";
out_make_masterdic="MD_Hagar_i0MM_sg1MM.txt"
out_IncAbu="Hagar_1MM_10X_2min";

guide_file="30kMs_unique.tsv"
exp_index="Hagar_expDesign.csv";

## parameters


## alignment with bowtie1
if [ ! -e ${Sam}.sam ]; then
    echo "alignment started ..."
    bowtie -m 1 --best --strata -S -p 16 -v 1 $indexedGuides --max ${Sam}.notUniq.fastq <( bamToFastq -i ${RAW} -fq /dev/stdout | fastx_trimmer -Q33 -l 20 ) | samtools view -S -F 0x0004 - > ${Sam}.sam;
    echo "alignment finished ..."
fi

## make shorter sam file to save information required for the following processing
if [ ! -e ${Sam}_short.sam ]; then
    echo "shorter sam started..."
    samtools view -h $RAW | cut -f 1,10,12,14 > ${Sam}_short.sam
    echo "shorter sam finished..."
fi

## extract aligned and sam files
if [ ! -e ${out_merge_bowtie_shortsam} ]; then 
    echo "extract from sam files started... "
    echo python Merge_bowtie_shortsam.py ${Sam}.sam ${Sam}_short.sam $out_merge_bowtie_shortsam;
    python Merge_bowtie_shortsam.py ${Sam}.sam ${Sam}_short.sam $out_merge_bowtie_shortsam;
    echo "extracting from bam files finished"
fi

## make master dict
if [ ! -e ${out_make_masterdic} ]; then 
    echo "master dict started... "
    echo python make_Masterdict.py $out_merge_bowtie_shortsam $guide_file $exp_index $out_make_masterdic;
    python make_Masterdict.py $out_merge_bowtie_shortsam $guide_file $exp_index $out_make_masterdic;
    echo "master dict finished..."
fi

## generate IncAb table
if [ ! -e "out_5_${out_IncAbu}_IncAbu.txt" ]; then
    echo "IncAb table started..."
    python IncidenceAbundance.py $out_make_masterdic $exp_index $out_IncAbu;
    echo "IncAb table finished..."
fi