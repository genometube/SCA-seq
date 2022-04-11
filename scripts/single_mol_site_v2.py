import argparse
from operator import itemgetter
import numpy as np
#import matplotlib.pyplot as plt
#import seaborn as sns
parser = argparse.ArgumentParser(description='${python} single_mol_site_v2.py -i test_cpggpc.tsv -o test_cpggpc.tsv.out -O test_cpggpc.tsv.out.site')
parser.add_argument('-i', '--input', required=True, help = '')
parser.add_argument('-o', '--out_methy', required=True, help = '')
parser.add_argument('-O', '--out_site', required=True, help = '')
args = parser.parse_args()

if __name__ == '__main__':
    # write
    o = open(args.out_methy, 'w')
    O = open(args.out_site, 'w')
    # read
    with open(args.input, 'r') as f:
        GCG_num = 0
        count = 0
        prev_line = ''
        while f:
            line = f.readline().strip()
            if chr == 'chromosome':
                continue
            if prev_line == '':
                prev_line = line
                continue

            tmp_p = prev_line.split('\t')
            chr_p, strand_p, start_p, end_p, read_id_p, log_lik_ratio_p, num_motifs_p, motif_type_p, seq_p = tmp_p[0], tmp_p[1], tmp_p[2], tmp_p[3], tmp_p[4],tmp_p[5], tmp_p[9], tmp_p[10], tmp_p[11]
            if line == "" and count == 1:
                break
            if line == "" and count == 0:
                # print(prev_line)
                methy_pos_all = np.where(np.asarray(list(seq_p)) == 'M')
                for i in range(len(methy_pos_all[0])):
                    if seq_p[methy_pos_all[0][i] + 1] != 'G' and motif_type_p == 'GC':  # exclude GCG
                        O.write('\t'.join([chr_p, str(2+int(start_p) + int(methy_pos_all[0][i])-int(methy_pos_all[0][0])), read_id_p, strand_p, log_lik_ratio_p,'GpC']) + '\n')
                        if float(log_lik_ratio_p) > 1:
                            o.write('\t'.join([chr_p, str(2 + int(start_p) + int(methy_pos_all[0][i]) - int(methy_pos_all[0][0])),read_id_p, strand_p, log_lik_ratio_p, 'GpC']) + '\n')
                    elif seq_p[methy_pos_all[0][i] - 1] != 'G' and motif_type_p == 'CG':  # exclude GCG
                        O.write('\t'.join([chr_p, str(1+int(start_p) + int(methy_pos_all[0][i])-int(methy_pos_all[0][0])), read_id_p, strand_p, log_lik_ratio_p,'CpG']) + '\n')
                        if float(log_lik_ratio_p) > 1.5:
                            o.write('\t'.join([chr_p, str(1 + int(start_p) + int(methy_pos_all[0][i]) - int(methy_pos_all[0][0])),read_id_p, strand_p, log_lik_ratio_p, 'CpG']) + '\n')
                    elif seq_p[methy_pos_all[0][i] - 1] == 'G' and motif_type_p == 'CG':
                        GCG_num += 1
                    elif seq_p[methy_pos_all[0][i] + 1] == 'G' and motif_type_p == 'GC':
                        GCG_num += 1
                break
            tmp = line.split('\t')
            chr, strand, start, end, read_id, log_lik_ratio, num_motifs, motif_type, seq = tmp[0], tmp[1], tmp[2], tmp[3], tmp[4], tmp[5], tmp[9], tmp[10], tmp[11]
            ### prev_line locus == line locus
            if chr == chr_p and start == start_p and end == end_p:
                m_pos_CG = np.where(np.asarray(list(seq_p)) == 'M')
                m_pos_GC = np.where(np.asarray(list(seq)) == 'M')
                m_pos_all = np.unique(np.concatenate((m_pos_CG[0],m_pos_GC[0]),0))
                m_pos_CG_clean, m_pos_GC_clean = np.in1d(m_pos_CG, m_pos_GC), np.in1d(m_pos_GC, m_pos_CG)

                for i in range(len(m_pos_all)):
                    if m_pos_all[i] in m_pos_CG[0] and m_pos_all[i] in m_pos_GC[0]:
                        GCG_num += 1
                    elif m_pos_all[i] in m_pos_CG[0] and m_pos_all[i] not in m_pos_GC[0]:
                        if m_pos_CG[0][0] >= m_pos_GC[0][0]:
                            O.write('\t'.join([chr_p, str(2 + int(start_p) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id_p,strand_p, log_lik_ratio_p, 'CpG']) + '\n')
                            if float(log_lik_ratio_p) > 1.5:
                                o.write('\t'.join([chr_p, str(2 + int(start_p) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id_p,strand_p, log_lik_ratio_p, 'CpG']) + '\n')
                        elif m_pos_CG[0][0] < m_pos_GC[0][0]:
                            O.write('\t'.join([chr_p, str(1 + int(start_p) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id_p,strand_p, log_lik_ratio_p, 'CpG']) + '\n')
                            if float(log_lik_ratio_p) > 1.5:
                                o.write('\t'.join([chr_p, str(1 + int(start_p) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id_p,strand_p, log_lik_ratio_p, 'CpG']) + '\n')
                    elif m_pos_all[i] not in m_pos_CG[0] and m_pos_all[i] in m_pos_GC[0]:
                        if m_pos_CG[0][0] >= m_pos_GC[0][0]:
                            O.write('\t'.join([chr, str(2 + int(start) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id,strand, log_lik_ratio, 'GpC']) + '\n')
                            if float(log_lik_ratio) > 1:
                                o.write('\t'.join([chr, str(2 + int(start) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id, strand,log_lik_ratio, 'GpC']) + '\n')
                        elif m_pos_CG[0][0] < m_pos_GC[0][0]:
                            O.write('\t'.join([chr, str(1 + int(start) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id, strand,log_lik_ratio, 'GpC']) + '\n')
                            if float(log_lik_ratio) > 1:
                                o.write('\t'.join([chr, str(1 + int(start) + int(m_pos_all[i]) - int(m_pos_all[0])), read_id, strand,log_lik_ratio, 'GpC']) + '\n')

                count = 1
                prev_line = line
            elif count == 1:    ### line after prev_line locus == line locus
                count = 0
                prev_line = line
                continue
            else:    ### prev_line locus != line locus
                methy_pos_all = np.where(np.asarray(list(seq_p)) == 'M')
                for i in range(len(methy_pos_all[0])):
                    if seq_p[methy_pos_all[0][i] + 1] != 'G' and motif_type_p == 'GC':  # exclude GCG
                        O.write('\t'.join([chr_p, str(2+int(start_p) + int(methy_pos_all[0][i])-int(methy_pos_all[0][0])), read_id_p, strand_p, log_lik_ratio_p,'GpC']) + '\n')
                        if float(log_lik_ratio_p) > 1:
                            o.write('\t'.join([chr_p, str(2 + int(start_p) + int(methy_pos_all[0][i]) - int(methy_pos_all[0][0])),read_id_p, strand_p, log_lik_ratio_p, 'GpC']) + '\n')
                    elif seq_p[methy_pos_all[0][i] - 1] != 'G' and motif_type_p == 'CG':  # exclude GCG
                        O.write('\t'.join([chr_p, str(1+int(start_p) + int(methy_pos_all[0][i])-int(methy_pos_all[0][0])), read_id_p, strand_p, log_lik_ratio_p,'CpG']) + '\n')
                        if float(log_lik_ratio_p) > 1.5:
                            o.write('\t'.join([chr_p, str(1 + int(start_p) + int(methy_pos_all[0][i]) - int(methy_pos_all[0][0])),read_id_p, strand_p, log_lik_ratio_p, 'CpG']) + '\n')
                    elif seq_p[methy_pos_all[0][i] - 1] == 'G' and motif_type_p == 'CG':
                        GCG_num += 1
                    elif seq_p[methy_pos_all[0][i] + 1] == 'G' and motif_type_p == 'GC':
                        GCG_num += 1

                prev_line = line

        print('GCG_num: ' + str(GCG_num) + '\n')
    o.close()
