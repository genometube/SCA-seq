#YX Sun Sep 20 20:08:08 CST 2020
"""input example
/share/app/python-3.6.9/bin/python3.6 /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_xieyeming1/proj_2021/zhichao_4Dhic_20210419/script/py_script/merge_overlap_4k.py -i /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_xieyeming1/proj_2021/zhichao_4Dhic_20210419/output/_04_slide_window/barcode_index_sorted.bed \
-o /zfssz5/BC_RD_P1/PROJECT/P18Z12200N0396_xieyeming1/proj_2021/zhichao_4Dhic_20210419/output/_04_slide_window/barcode_index_merge_4k.bed
cat barcode_index_merge_4k.bed|awk '{print $3 - $2}'|awk '{A[$1]++}END{for(i in A)print i,A[i]}'|sort -n -k1|tr [:blank:] \\t|awk '{b[$1]=$2;sum=sum+$2} END{for (i in b) print i,b[i],(b[i]/sum)*100}'|tr [:blank:] \\t|sort -k1 -n > barcode_index_merge_4k.bed.list.pct
cat barcode_index_merge_4k.bed|awk '{print $3 - $2}'> barcode_index_merge_4k.bed.list
"""

import argparse
from operator import itemgetter
import subprocess
import re
import os
import datetime

parser = argparse.ArgumentParser(description='')
parser.add_argument('-i', '--input_file', required=True, help = '')
parser.add_argument('-r', '--range', required=True, help = '')
parser.add_argument('-o', '--output', required=True, help = 'output')
args = parser.parse_args()

def intersection(lst1, lst2):
    lst3 = [value for value in lst1 if value in lst2]
    return lst3

def mergeIntervals(arr):
    # Sorting based on the increasing order
    # of the start intervals
    arr.sort(key=lambda x: x[0])
    # array to hold the merged intervals
    # 'max' value gives the last point of
    # that particular interval
    # 's' gives the starting point of that interval
    # 'm' array contains the list of all merged intervals
    m = []
    s = -10000
    max = -100000
    merge_count = 0
    merge_count_list = []
    for i in range(len(arr)):
        #merge_count = 0
        a = arr[i]
        #print(a)
        if a[0] > max:
            if i != 0:
                m.append([s, max])
            merge_count_list.append(1)
            max = a[1]
            s = a[0]
        else:
            if a[1] >= max:
                max = a[1]
            merge_count += 1
            merge_count_list[i - merge_count] += 1

        #print(m)
        #print(merge_count_list)
        #print(i)
        #print(merge_count)
        #print(s,max)
    if max != -100000 and [s, max] not in m:
        m.append([s, max])
        merge_count_list.append(1)

    return m,merge_count,merge_count_list

if __name__ == '__main__':
    extend_len = int(int(args.range)/2)
    chr_id_to_start_end = {}
    with open(args.input_file, 'r') as f1:
        for line in f1:
#        line = f1.read().strip().split('\n')
#        for i in range(len(line)):
            tmp = line.strip().split()
            chr,start,end,barcode,strand = tmp[0],int(tmp[1])-extend_len,int(tmp[2])+extend_len,tmp[3],tmp[4]
            #strand = '.'
            chr_id = chr + '%' + barcode + '%' + strand
            if chr_id in chr_id_to_start_end:
                chr_id_to_start_end[chr_id].append([start,end])
            else:
                chr_id_to_start_end[chr_id] = [[start,end]]
            #print(tmp)
    f_out = open(args.output, 'w')
    for k in chr_id_to_start_end:
        #print(chr_id_to_start_end[k])
        m,merge_count,merge_count_list = mergeIntervals(chr_id_to_start_end[k])
        #print(m)
        #print(k + '\t' + str(merge_count))
        tmp = k.split('%')
        chr,barcode,strand = tmp[0],tmp[1],tmp[2]
        for i in range(len(m)):
            f_out.write('\t'.join([chr, str(m[i][0] + extend_len), str(m[i][1] - extend_len), barcode, strand, str(merge_count_list[i])]) + '\n')

    f_out.close()