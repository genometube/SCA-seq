#qsub -cwd -l vf=3G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log
date

interval="80"
start=1
end=100

for i in $(seq $start $end); do
  echo "qsub -cwd -l vf=10G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log bwasw_nanopolish_single_run.sh $i ${interval}"
  qsub -cwd -l vf=10G,num_proc=1 -P RNAProj -q bc_rd.q,bc.q -e log -o log bwasw_nanopolish_single_run.sh $i ${interval}
done
date