#job submission: qsub -cwd -l vf=40G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log [bash.sh]
python='/hwfswh2/BC_PUB/Software/07.User-defined/02.Research/xieyeming/conda/anaconda3/bin/python'
bedtools="/ifswh1/BC_PUB/biosoft/pipeline/10x/RNA/10x_RNA_ATAC/cellranger-atac-1.2.0/cellranger-atac-1.2.0/miniconda-atac-cs/4.3.21-miniconda-atac-cs-c10/bin/bedtools"
sample='SCActrl_hg19mt'
chunk_dir='/zfswh1/BC_RD_P0/P19Z12200N0089_xym/project_2021/yaning_SCActrl_20220314/01.bwasw_nanopolish_hg19mt/chunk/'
juicer_tools="/research/liwenqing/Software/juicer/juicer_tools.1.9.9_jcuda.0.8.jar"
Rscript="/mnt/Software/anaconda3/envs/r/bin/Rscript"

date
chunk=($(ls ${chunk_dir}))
echo ${chunk[@]}
for i in ${chunk[@]}; do
  ${bedtools} bamtobed -bed12 -i ${chunk_dir}/${i}/${sample}.sorted.${i}.bam |cut -f1,2,3,4,6 >> ${sample}_combine.bed
done
mkdir -p tmp
sort -T tmp -k4,4 -k1,1 -k2,2n ${sample}_combine.bed > ${sample}_resort.bed
${python} scripts/merge_overlap_range_no_index.py -i ${sample}_resort.bed -r 150 -o ${sample}_merged.bed
cut -f1-5 ${sample}_merged.bed > ${sample}_merged.bed.final

${python} scripts/add_segment_id.py --input_bed ${sample}_merged.bed.final --output_bed_seg ${sample}_merged.bed_seg
cat ${sample}_merged.bed_seg |cut -f1-4|sed 's/\t/_/'|sed 's/\t/_/' > ${sample}.hic_tmp1
${Rscript} scripts/cartesian_join.R ${sample}.hic_tmp1 HiC_merge.csv

cat HiC_merge.csv|sed 's/\_/\t/g' > HiC_merge.csv.tmp1
${python} scripts/change_file_format_v3.py HiC_merge.csv.tmp1 HiC_merge.csv.tmp1.changed
awk '{if ($3 > $7){ print $1, $6, $7, $8, $9, $2, $3, $4, $5, $11,$10}else {print}}' HiC_merge.csv.tmp1.changed > sort_changed_HiC_merge.csv.tmp1
sort -T tmp -k3,3d -k7,7d sort_changed_HiC_merge.csv.tmp1 > sort_changed_HiC_merge.csv.tmp2
java -jar ${juicer_tools} pre sort_changed_HiC_merge.csv.tmp2 ${sample}.hic hg19 -t tmp

date
