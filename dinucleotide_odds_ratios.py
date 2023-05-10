# dinucleotide_odds_ratios.py a simple python script that calculates odds ratio of the 16 dinucleotide sequences.
# v0.1 Rhys Parry r.parry@uq.edu.au
# Usage python dinucleotide_odds_ratios.py inputCDSmultifastafile.fa > output.txt

import sys
from collections import defaultdict

# Define a function to calculate the odds ratio
def calc_odds_ratio(seq):
    # Initialize dictionary to count dinucleotide frequencies
    counts = defaultdict(int)
    for i in range(len(seq) - 1):
        dinuc = seq[i:i+2]
        counts[dinuc] += 1

    # Calculate odds ratio for each dinucleotide
    odds_ratios = {}
    for dinuc in counts:
        X = dinuc[0]
        Y = dinuc[1]
        f_XY = counts[dinuc] / (len(seq) - 1)
        f_X = seq.count(X) / len(seq)
        f_Y = seq.count(Y) / len(seq)
        odds_ratios[dinuc] = f_XY / (f_X * f_Y)

    return odds_ratios


# Open multifasta file
fasta_file = sys.argv[1]

# Parse sequences from multifasta file
sequences = {}
with open(fasta_file) as f:
    name = None
    seq = ''
    for line in f:
        if line.startswith('>'):
            if name is not None:
                sequences[name] = seq.upper()
                seq = ''
            name = line[1:].strip()
        else:
            seq += line.strip()
    if name is not None:
        sequences[name] = seq.upper()

# Calculate odds ratios for each sequence
dinucs = ['AA', 'AT', 'AC', 'AG', 'TA', 'TT', 'TC', 'TG', 'CA', 'CT', 'CC', 'CG', 'GA', 'GT', 'GC', 'GG']
header = '\t'.join(dinucs)
print('Sequence\t' + header)
for name, seq in sequences.items():
    row = []
    row.append(name)
    odds_ratios = calc_odds_ratio(seq)
    for dinuc in dinucs:
        odds_ratio = odds_ratios.get(dinuc, 0)
        row.append(str(round(odds_ratio, 2)))
    print('\t'.join(row))
