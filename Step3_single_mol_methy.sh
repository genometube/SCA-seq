#job submission: qsub -cwd -l vf=20G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log [bash.sh]
python='/hwfswh2/BC_PUB/Software/07.User-defined/02.Research/xieyeming/conda/anaconda3/bin/python'
sample='SCActrl_hg19mt'
chunk_dir='/zfswh1/BC_RD_P0/P19Z12200N0089_xym/project_2021/yaning_SCActrl_20220314/01.bwasw_nanopolish_hg19mt/chunk/'

date
chunk=($(ls ${chunk_dir}))
echo ${chunk[@]}
for i in ${chunk[@]}; do
  cat ${chunk_dir}/${i}/${sample}_cpggpc_methy_calls.${i}.xls |tail -n +2 >> ${sample}_cpggpc_methy_calls.xls
done
mkdir -p tmp
echo "methy table"
 ${python} scripts/single_mol_site_v2.py -i ${sample}_cpggpc_methy_calls.xls -o ${sample}_cpggpc_single_mol_methy.xls -O ${sample}_cpggpc_single_mol_site.xls
sort -k3,3 -k1,1 -k2,2n -T tmp ${sample}_cpggpc_single_mol_methy.xls > sorted_${sample}_cpggpc_single_mol_methy.xls
sort -k3,3 -k1,1 -k2,2n -T tmp ${sample}_cpggpc_single_mol_site.xls > sorted_${sample}_cpggpc_single_mol_site.xls

grep GpC sorted_${sample}_cpggpc_single_mol_methy.xls > sorted_${sample}_gpc_single_mol_methy.xls
grep CpG sorted_${sample}_cpggpc_single_mol_methy.xls > sorted_${sample}_cpg_single_mol_methy.xls
grep GpC sorted_${sample}_cpggpc_single_mol_site.xls > sorted_${sample}_gpc_single_mol_site.xls
grep CpG sorted_${sample}_cpggpc_single_mol_site.xls > sorted_${sample}_cpg_single_mol_site.xls

bed_seg="../02.hic_hg19mt/${sample}_merged.bed_seg"

${python} scripts/link_bed_methy_v1.py -b ${bed_seg} -m sorted_${sample}_gpc_single_mol_methy.xls -o ${sample}_chr7.bed_seg_gpc_m -z 0
${python} scripts/link_bed_methy_v1.py -b ${sample}_chr7.bed_seg_gpc_m -m sorted_${sample}_gpc_single_mol_site.xls -o ${sample}_chr7.bed_seg_gpc_m_um -z 0
${python} scripts/link_bed_methy_v1.py -b ${sample}_chr7.bed_seg_gpc_m_um -m sorted_${sample}_cpg_single_mol_methy.xls -o ${sample}_chr7.bed_seg_gpc_m_um_cpg_m -z 0
${python} scripts/link_bed_methy_v1.py -b ${sample}_chr7.bed_seg_gpc_m_um_cpg_m -m sorted_${sample}_cpg_single_mol_site.xls -o ${sample}_chr7.bed_seg_gpc_m_um_cpg_m_um -z 0

echo "combine single_mol_methy_seg table"
echo "chr_name segment_start segment_end read_id strand segment_num GpC_methy GpC_unmethy CpG_methy CpG_unmethy"|tr [:blank:] \\t > ${sample}_single_mol_methy_seg.xls
cat ${sample}_chr7.bed_seg_gpc_m_um_cpg_m_um |awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8-$7"\t"$9"\t"$10-$9}' >>  ${sample}_single_mol_methy_seg.xls
gzip -c ${sample}_single_mol_methy_seg.xls > ${sample}_single_mol_methy_seg.xls.gz

echo "combine single_mol_methy_site table"
echo "chr_name methy_site read_id strand likelihood_ratio methy_type"|tr [:blank:] \\t > ${sample}_single_mol_methy_site.xls
cat sorted_${sample}_cpggpc_single_mol_site.xls >> ${sample}_single_mol_methy_site.xls
gzip -c ${sample}_single_mol_methy_site.xls > ${sample}_single_mol_methy_site.xls.gz

date

