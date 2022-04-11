#qsub -cwd -l vf=1G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log
export PATH=/zfswh1/BC_RD_P0/zhangchen3/software/bin/bzip2/bin:$PATH
export LD_LIBRARY_PATH=/zfswh1/BC_RD_P0/zhangchen3/software/gcc/lib64:/zfswh1/BC_RD_P0/zhangchen3/software/gcc/lib:/zfswh1/BC_RD_P0/zhangchen3/software/bin/bzip2/lib:$LD_LIBRARY_PATH
export CFLAGS="-fPIC -I/zfswh1/BC_RD_P0/zhangchen3/software/bin/bzip2/include"
export LDFLAGS="-fPIC -L/zfswh1/BC_RD_P0/zhangchen3/software/bin/bzip2/lib"
export PATH=$PATH:/zfswh1/BC_RD_P0/zhangchen3/software/gcc/bin
export C_INCLUDE_PATH=$C_INCLUDE_PATH:/zfswh1/BC_RD_P0/zhangchen3/software/gcc/include
export CPLUS_INCLUDE_PATH=$CPLUS_INCLUDE_PATH:/zfswh1/BC_RD_P0/zhangchen3/software/gcc/include
export LIBRARY_PATH=$LIBRARY_PATH:/zfswh1/BC_RD_P0/zhangchen3/software/gcc/lib

bwa="/share/app/bwa-0.7.12/bwa"
samtools="/share/app/samtools-1.2/bin/samtools"
nanopolish='/hwfswh2/BC_PUB/Software/07.User-defined/02.Research/zhangchen/software/nanopolish/nanopolish'
export HDF5_PLUGIN_PATH=/zfswh1/BC_RD_P0/P19Z12200N0089_xym/archive/plugin
thread='1'

sample='SCActrl_hg19mt'
ref='/zfswh1/BC_RD_P0/P19Z12200N0089_xym/db/genome/hg19/Sequence/bwa_index_hg19mt/genome.fa'

chip_id="PAK52156_pass_795ebebf"
fast5_pass="/ifswh1/rawdata/Nanopore/PC48A014/WHONT20220301-3/WHpooling20220301-3/20220302_0122_6F_PAK52156_0f915e9a/fast5_pass"
fastq_pass="/ifswh1/rawdata/Nanopore/PC48A014/WHONT20220301-3/WHpooling20220301-3/20220302_0122_6F_PAK52156_0f915e9a/fastq_pass"

num="$1"
interval="$2"
start=$(((num-1)*interval))
end=$((num*interval - 1))

date
#list=($(find ${fastq_pass} -name "*.fastq" -type f|sed 's/\_/\t/g'|rev|cut -f1|rev|cut -d '.' -f1|sort -n))
#echo {0..500}
echo "$start $end"
for k in $(seq $start $end); do
  i=chunk_${k}
  mkdir -p chunk/${i}
  if [ ! -f chunk/${i}/${sample}.sorted.${i}.bam ];then
    echo ${i}
    date
    mkdir -p chunk/${i}/fast5 chunk/${i}/tmp
    ${bwa} bwasw -t ${thread} -b 5 -q 2 -r 1 -T 15 -z 10 ${ref} ${fastq_pass}/${chip_id}_${k}.fastq |\
    awk '$6 !~ /H/ { print }'|awk '$5 > 29'|awk 'length($10) > 49' >> chunk/${i}/${i}.sam
    ${samtools} view -@ ${thread} -bT ${ref} chunk/${i}/${i}.sam > chunk/${i}/${i}.bam
    ${samtools} sort -T chunk/${i}/tmp -@ ${thread} -o chunk/${i}/${sample}.sorted.${i}.bam chunk/${i}/${i}.bam

    cp ${fastq_pass}/${chip_id}_${k}.fastq chunk/${i}/${i}.fq
    cp ${fast5_pass}/${chip_id}_${k}.fast5 chunk/${i}/fast5/${chip_id}_${k}.fast5

    ${nanopolish} index -d chunk/${i}/fast5 chunk/${i}/${i}.fq
    ${samtools} index chunk/${i}/${sample}.sorted.${i}.bam

    echo "call gpc cpg"
    ${nanopolish} call-methylation -t ${thread} -r chunk/${i}/${i}.fq -b chunk/${i}/${sample}.sorted.${i}.bam -g ${ref} --methylation cpggpc > chunk/${i}/${sample}_cpggpc_methy_calls.${i}.xls
    rm -rf chunk/${i}/${i}.fq* chunk/${i}/fast5/ chunk/${i}/${i}.sam chunk/${i}/${i}.bam chunk/${i}/tmp
    date
  fi

done

date

