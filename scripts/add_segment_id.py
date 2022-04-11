import argparse
from operator import itemgetter
# import numpy as np
# import matplotlib.pyplot as plt
# import seaborn as sns
parser = argparse.ArgumentParser(description='python /home/yeming/PycharmProjects/pythonProject/yaning_poreC/add_segment_id.py -i nome75_chr7_resort.bed -o nome75_chr7_resort.bed_seg')
parser.add_argument('-i', '--input_bed', required=True, help='')
parser.add_argument('-o', '--output_bed_seg', required=True, help='')
args = parser.parse_args()

if __name__ == '__main__':
    # write
    o = open(args.output_bed_seg, 'w')
    # read
    with open(args.input_bed, 'r') as f:
        line = f.read().strip().split('\n')
    o.write(line[0] + '\t' + '1'+'\n')
    count = 1
    for i in range(1,len(line)):
        if line[i-1].split()[3] == line[i].split()[3]:
            count += 1
            o.write(line[i] + '\t' + str(count) + '\n')
        else:
            count = 1
            o.write(line[i] + '\t' + '1' + '\n')
    o.close()

