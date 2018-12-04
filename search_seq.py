from Bio.SeqUtils import nt_search
from Bio import SeqIO

query = 'NTNNNNAN'
file = '/home/user/Desktop/transposon/SCTrainingSet/WT_1/FreadI1_trimmed.fastq'

# counts number of appearances of a query sequence in a .fastq file 
def count_appearance(file, query):
    count = 0
    for record in SeqIO.parse('/home/user/Desktop/transposon/SCTrainingSet/WT_1/FreadI1_trimmed.fastq', 'fastq'):
        x = nt_search(str(record), query)
        count += len(x[1:])

    return count

# x = count_appearance(file, query)
# print(x)

# counts lines in a txt file
def count_lines(file):
    with open(file) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

# print(count_lines('/home/user/Desktop/transposon/SCTrainingSet/mapped_sorted_reads.txt'))
# print(count_lines('/home/user/Desktop/transposon/SCTrainingSet/unmapped_sorted_reads.txt'))
