from Bio.SeqUtils import nt_search
from Bio import SeqIO

query = 'NTNNNNAN'
file = '/home/user/Desktop/transposon/SCTrainingSet/WT_1/FreadI1_trimmed.fastq'
# file = '/home/user/Desktop/transposon/S288C_reference_sequence_R64-2-1_20150113.fsa'


def find_positions(seqfile, queryseq, filetype='fasta'):

    """ Finds positions of a query sequence in a file """

    genome = SeqIO.parse(seqfile, filetype)
    positions = nt_search(str(genome), queryseq)
    return positions


def freq_appearance(file, query, filetype='fasta'):

    """ Frequency of occurrence of a query sequence in a file

    :param file: sequence file
    :param query: sequence query (string)
    :param filetype: parameter for SeqIO.parse (default='fasta')

    """

    count = 0
    total_length = 0
    for record in SeqIO.parse(file, filetype):
        x = nt_search(str(record.seq), query)
        count += len(x[1:])
        total_length += len(record.seq)
    print(total_length)
    print(count)
    print(count / total_length)


def count_lines(txtfile):

    """ Counts lines in a txt file """

    with open(txtfile) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


freq_appearance(file, query, 'fastq')
