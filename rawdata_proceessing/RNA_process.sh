#!/usr/bin/bash
set -x

samtools='samtools'
hisat2='hisat2'
bcftools='bcftools'
htseq='htseq-count'
bedtools='bedtools'
python='python'
trim_galore='trim_galore'

SARS_CoV_2_ht2_idx=/gshare/xielab/jianfc/COVID/pseudovirus_mutation/reference/SARS-CoV-2-wt/SARS-CoV-2-wt-genome
#SARS_Cov_2_ht2_idx=/share/home/baiyl/database/SARS_Cov_2/hisat_index/SARS_Cov_2
green_monkey_idx=/share/home/baiyl/database/chlSab2/Index/hisat_index/chlSab2
SARS_Tor2_ht2_idx=/share/home/baiyl/database/SARS_Tor2/hisat_index/SARS_Tor2
rVSVdel21_idx=/gshare/xielab/jianfc/COVID/pseudovirus_mutation/reference/rVSVdel21/rVSVdel21
rVSVdel21_beta_idx=/gshare/xielab/jianfc/COVID/pseudovirus_mutation/reference/rVSVdel21_beta/rVSVdel21_beta
rVSVdel21_delta_idx=/gshare/xielab/jianfc/COVID/pseudovirus_mutation/reference/rVSVdel21_delta/rVSVdel21_delta
rVSVOmicron_idx=/gshare/xielab/jianfc/COVID/pseudovirus_mutation/reference/rVSVOmicron/rVSVOmicron
SARS_CoV_2_spike_ht2_idx=/gshare/xielab/jianfc/COVID/pseudovirus_mutation/reference/SARS-CoV-2-spike/SARS-CoV-2-spike

declare -A ref_ht2
ref_ht2["SARS-CoV-2"]=$SARS_CoV_2_ht2_idx
ref_ht2["SARS-CoV-2-spike"]=$SARS_CoV_2_spike_ht2_idx
ref_ht2["SARS_Tor2"]=$SARS_Tor2_ht2_idx
ref_ht2["green_monkey"]=$green_monkey_idx
ref_ht2["rVSVdel21"]=$rVSVdel21_idx
ref_ht2["rVSVOmicron"]=$rVSVOmicron_idx
ref_ht2["rVSVdel21_beta"]=$rVSVdel21_beta_idx
ref_ht2["rVSVdel21_delta"]=$rVSVdel21_delta_idx

id=$1
sp=$2
R1=$3
R2=$4
OUT=$5

ht2_idx=${ref_ht2[$sp]}
function trim_pre_se {
mkdir -p ${OUT}/trim/$id
$trim_galore --fastqc --output_dir $OUT/trim/$id --cores $SLURM_CPUS_PER_TASK $R1
}

function trim_pre {
Read1=$R1
Read2=$R2
read_len=20
mkdir -p ${OUT}/trim/${id}

$trim_galore --fastqc --paired \
        --phred33 \
        --length $read_len \
        --retain_unpaired \
        --output_dir ${OUT}/trim/${id} \
        --cores $SLURM_CPUS_PER_TASK \
        $Read1 $Read2

}

function hisat2_align_se {
mkdir -p $OUT/hisat2
Read1=$OUT/trim/$id/$(ls ${OUT}/trim/${id} | grep _trimmed.fq.gz)

$hisat2 -p $SLURM_CPUS_PER_TASK \
	-x $ht2_idx \
	-U $Read1 \
	-S ${OUT}/hisat2/${id}.${sp}.sam

$samtools view -h -q 2 -@ $SLURM_CPUS_PER_TASK -O BAM ${OUT}/hisat2/${id}.${sp}.sam|$samtools sort -@ $SLURM_CPUS_PER_TASK > ${OUT}/hisat2/${id}.${sp}.bam
$samtools index ${OUT}/hisat2/${id}.${sp}.bam
}

function hisat2_align {
mkdir -p ${OUT}/hisat2
Read1=$OUT/trim/$id/$(ls ${OUT}/trim/${id} | grep _val_1.fq.gz)
Read2=$OUT/trim/$id/$(ls ${OUT}/trim/${id} | grep _val_2.fq.gz)

$hisat2 -p $SLURM_CPUS_PER_TASK \
	-x $ht2_idx \
	-1 $Read1 \
	-2 $Read2 \
	-S ${OUT}/hisat2/${id}.${sp}.sam

$samtools view -h -q 2 -@ $SLURM_CPUS_PER_TASK -O BAM ${OUT}/hisat2/${id}.${sp}.sam|$samtools sort -@ $SLURM_CPUS_PER_TASK > ${OUT}/hisat2/${id}.${sp}.bam
$samtools index ${OUT}/hisat2/${id}.${sp}.bam

}
function Dedup {
picard='picard'
$picard -Xmx20g MarkDuplicates \
        I=$OUT/hisat2/${id}.${sp}.bam \
        O=$OUT/hisat2/${id}.${sp}.rmdup.bam \
        M=$OUT/hisat2/${id}.${sp}.rmdup.txt \
        REMOVE_DUPLICATES=true
$samtools index $OUT/hisat2/${id}.${sp}.rmdup.bam
}

function call_snps {
mkdir -p ${OUT}/Snps
genome=$(ls ${ht2_idx}* | grep .f | grep a$)

$bcftools mpileup -d 900000 -C 50 -Q 30 -q 5 -E ${OUT}/hisat2/${id}.${sp}.rmdup.bam -f $genome --output ${OUT}/Snps/${id}.${sp}.snp_dup.vcf

$bcftools view  ${OUT}/Snps/${id}.${sp}.snp_dup.vcf |sed 's/;I16.*QS=/\t/g'|sed 's/;.*//g'|sed 's/DP=//g'|awk '$0!~"^#"' | grep -v INDEL > ${OUT}/Snps/${id}.${sp}.variants_dup.txt
$python Allele_Snp_rate.py ${OUT}/Snps/${id}.${sp}.variants_dup.txt ${OUT}/Snps/${id}.${sp}.variants_dup_clean.txt
grep -v Chrom ${OUT}/Snps/${id}.${sp}.variants_dup_clean.txt |awk 'BEGIN{depths=0;Alts=0}{if($3=="A"){alt=1-$5}else if($3=="T"){alt=1-$6}else if($3=="C"){alt=1-$7}else if($3=="G"){alt=1-$8};depths+=$4;Alts+=alt;aa=$4*alt;bb=sprintf("%.0f",aa);print $0,alt,bb}' > ${OUT}/Snps/${id}.${sp}.variants_dup_clean.mut.txt
# $bedtools intersect -a <(awk '{$2=$2" "$2;print $0}' ${OUT}/Snps/${id}.${sp}.variants_dup_clean.mut.txt|sed 's/ /\t/g') -b <(awk '{print "NC_045512.2",$3,$4,$2}' /share/home/baiyl/database/SARS_Cov_2/Annotation.txt |sed 's/ /\t/g') -wao |cut -f -11,15|awk 'BEGIN{print "Chrom Pos Ref Depth A T C G MutRate MutNum Gene"}{print $0}'|sed 's/ /\t/g' > ${OUT}/Snps/${id}.${sp}.variants_dup_clean.mut.tsv
}

function stat {
mkdir -p $OUT/stat
set -x
Raw_reads=$(($(less $OUT/trim/$(ls ${OUT}/trim | grep ${id}_ | grep trimming_report) | grep "Total reads processed"|awk '{print $4}'|sed 's/,//g')*2))
Trimed_reads=$(($(zcat $OUT/trim/$(ls ${OUT}/trim | grep ${id}_ | grep _val_1.fq.gz) | wc -l)/2))
Filter_reads_r=`echo "scale=4;100*(1-$Trimed_reads/$Raw_reads)"|bc`

$samtools stats ${OUT}/hisat2/${id}.${sp}.bam > ${OUT}/hisat2/${id}.${sp}.bam.stat
reads_aligned=`grep "reads mapped:" ${OUT}/hisat2/${id}.${sp}.bam.stat|cut -f 3|xargs`
Align_rate=`echo "scale=2;100*$reads_aligned/$Trimed_reads"|bc|xargs`
Dup_ratio=`cat $OUT/hisat2/${id}.${sp}.rmdup.txt | grep "Unknown" | awk '{print $10}'`

echo -e "Raw_reads\tTrimed_reads\tFilter_reads_percent\treads_aligned\tAlign_percent\tDup_ratio" > $OUT/stat/Align_stat.txt
echo -e "${id}\t${sp}\t" \
"${Raw_reads}\t${Trimed_reads}\t${Filter_reads_r}\t" \
"${reads_aligned}\t${Align_rate}\t${Dup_ratio}" >> $OUT/stat/Align_stat.txt

set +x
}

if [ "$R2" == "NONE" ];then
trim_pre_se
hisat2_align_se
else
trim_pre
hisat2_align
fi

stat

Dedup

call_snps

set +x
