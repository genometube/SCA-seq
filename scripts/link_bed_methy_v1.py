import argparse
import gzip
from operator import itemgetter
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
parser = argparse.ArgumentParser(description='python /home/yeming/PycharmProjects/pythonProject/yaning_poreC/link_bed_methy.py -b nome75_chr7_resort.bed_seg -m sorted_nome75_gpc_single_mol_methy.xls -o nome75_chr7_resort.bed_m')
parser.add_argument('-b', '--input_bed', required=True, help='')
parser.add_argument('-m', '--input_methy', required=True, help='')
parser.add_argument('-o', '--output_bed_methy', required=True, help='')
parser.add_argument('-z', '--zip_bed', required=True, help='[0,1]')
args = parser.parse_args()

if __name__ == '__main__':
    # write
    o = open(args.output_bed_methy, 'w')
    # read
    id_to_methy = {}
    if args.zip_bed == '1':
        with gzip.open(args.input_methy, 'rt') as f:
            for line in f:
                tmp = line.strip().split()
                pos,read_id = int(tmp[1]),tmp[2]
                if read_id in id_to_methy:
                    id_to_methy[read_id].append(pos)
                else:
                    id_to_methy[read_id]=[pos]
    if args.zip_bed == '0':
        with open(args.input_methy, 'r') as f:
            for line in f:
                tmp = line.strip().split()
                pos,read_id = int(tmp[1]),tmp[2]
                if read_id in id_to_methy:
                    id_to_methy[read_id].append(pos)
                else:
                    id_to_methy[read_id]=[pos]

    with open(args.input_bed, 'r') as f:
        for line in f:
            tmp = line.strip().split()
            bed_start,bed_end,read_id = int(tmp[1]),int(tmp[2]),tmp[3]
            if read_id in id_to_methy:
                pos_list = id_to_methy[read_id]
                count = 0
                for i in range(len(pos_list)):
                    if pos_list[i] >= bed_start and pos_list[i] <= bed_end:
                        count += 1
                o.write('\t'.join(tmp) + '\t' + str(count) +'\n')
            else:
                o.write('\t'.join(tmp) + '\t' + '0' + '\n')

    o.close()
