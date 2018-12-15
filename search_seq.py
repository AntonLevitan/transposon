from Bio.SeqUtils import nt_search
from Bio import SeqIO

query = 'NTNNNNAN'
# file = '/home/user/Desktop/transposon/SCTrainingSet/WT_1/FreadI1_trimmed.fastq'
file = '/home/user/Desktop/transposon/S288C_reference_sequence_R64-2-1_20150113.fsa'

# find positions of a query sequence in a genome:
def find_positions(file, query):
    genome = SeqIO.parse(file, 'fasta')
    positions = nt_search(str(genome), query)
    return positions


# counts number of appearances of a query sequence in a file 
def count_appearance(file, query, filetype='fasta'):
    count = 0
    total_length = 0
    for record in SeqIO.parse(file, filetype):
        x = nt_search(str(record.seq), query)
        count += len(x[1:])
        total_length += len(record.seq)
    print(total_length)
    print(count)
    print(total_length / count)


# counts lines in a txt file
def count_lines(file):
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

count_appearance(file, query)