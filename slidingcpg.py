# slidingcpg.py a simple python script that calculates odds ratio of CpG over a designated window size -w and a step size -s of an input fasta file.
# v0.1 Rhys Parry r.parry@uq.edu.au July 2024
#Usage:  python slidingcpg.py input.fa -w 200 -s 100
import argparse
from Bio import SeqIO

def calc_cpg_ratio(sequence, window_size, step_size):
    results = []
    for i in range(0, len(sequence) - window_size + 1, step_size):
        window = sequence[i:i+window_size]
        c_count = window.count('C')
        g_count = window.count('G')
        cg_count = window.count('CG')
        try:
            o_e_ratio = (cg_count * len(window)) / (c_count * g_count)
        except ZeroDivisionError:
            o_e_ratio = 0
        results.append((i, o_e_ratio))
    return results

def main():
    parser = argparse.ArgumentParser(description='Perform sliding window analysis of CpG content.')
    parser.add_argument('input_file', type=str, help='Input FASTA file.')
    parser.add_argument('-w', type=int, default=500, help='Window size.')
    parser.add_argument('-s', type=int, default=250, help='Step size.')
    args = parser.parse_args()

    with open(args.input_file, 'r') as f:
        for record in SeqIO.parse(f, 'fasta'):
            results = calc_cpg_ratio(str(record.seq).upper(), args.w, args.s)
            for location, ratio in results:
                print(f'{location}\t{ratio}')

if __name__ == '__main__':
    main()

