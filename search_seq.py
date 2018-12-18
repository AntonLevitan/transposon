from Bio.SeqUtils import nt_search
from Bio import SeqIO
import gffutils
import pandas as pd
import csv
import numpy as np

query = 'NTNNNNAN'
# seq_file = '/home/user/Desktop/transposon/SCTrainingSet/WT_1/FreadI1_trimmed.fastq'
seq_file = 'S288C_reference_sequence_R64-2-1_20150113.fsa'

# gffutils creating database
# db = gffutils.create_db('saccharomyces_cerevisiae_R64-2-1_20150113.gff',
#                         'sc_gffutils_database')

db = gffutils.FeatureDB('sc_gffutils_database')
sc_genes = db.all_features(featuretype='gene')


def count_lines(txtfile):

    """ Counts lines in a txt file """

    with open(txtfile) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def find_positions(seqfile, queryseq, filetype='fasta'):

    """ Finds positions of a query sequence in a file """

    genome = SeqIO.parse(seqfile, filetype)
    positions = nt_search(str(genome), queryseq)
    return positions


def get_features(gffutils_db):

    """ Retrieves coordinates of all genes from a gffutils database """

    gene_chroms = []
    gene_ids = []
    gene_starts = []
    gene_ends = []
    for gene in gffutils_db:
        gene_chroms.append(gene.chrom)
        gene_ids.append(gene.id)
        gene_starts.append(gene.start)
        gene_ends.append(gene.end)

    gene_coords_df = pd.DataFrame(list(zip(gene_ids, gene_chroms, gene_starts, gene_ends)), columns=['id', 'chrom',
                                                                                                     'start', 'end'])
    gene_coords_df.to_csv('sc_features_coords.csv')

    return gene_coords_df, list(pd.unique(gene_chroms))


def freq_appearance(file, query, gffutils_db, filetype='fasta'):
    """ Frequency of occurrence of a query sequence in a file

    :param file: sequence file
    :param query: sequence query (string)
    :param gffutils_db: gffutils all_features database
    :param filetype: parameter for SeqIO.parse (default='fasta')

    """
    sc_features, chromes = get_features(gffutils_db)
    chroms_locs = pd.DataFrame(np.nan, index=range(250000), columns=chromes)
    count = 0
    total_length = 0
    i = 0
    for record in SeqIO.parse(file, filetype):
        x = nt_search(str(record.seq), query)
        count += len(x[1:])
        total_length += len(record.seq)
        chroms_locs[chromes[i]] = pd.Series(x[1:])
        i += 1

    print(total_length)
    print(count)
    print(float(count) / total_length)
    chroms_locs = chroms_locs.dropna(how='all')
    chroms_locs.to_csv('seq_locations_on_chromosomes.csv')

    return chroms_locs, sc_features


def count_target_seq():
    # chr_loc = pd.read_csv('seq_locations_on_chromosomes.csv').drop('Unnamed: 0', axis=1)
    # genes_coords = pd.read_csv('sc_features_coords.csv').drop('Unnamed: 0', axis=1)
    chr_loc, genes_coords = freq_appearance(seq_file, query, sc_genes)
    chroms = list(pd.unique(chr_loc.columns))
    genes_coords['target_seq_counts'] = np.zeros(len(genes_coords.index))
    genes_coords['500bp_up_target_seq_counts'] = np.zeros(len(genes_coords.index))
    genes_coords['500bp_down_target_seq_counts'] = np.zeros(len(genes_coords.index))
    target_counts = []
    up_target_counts = []
    down_target_counts = []

    for chrom in chroms:

        chr_coords = genes_coords[genes_coords['chrom'] == chrom]
        seq_intervals = []
        up_intervals = []
        down_intervals = []

        for j in chr_coords.id:
            seq_intervals.append(frozenset(range(chr_coords[chr_coords.id == j].start,
                                                 chr_coords[chr_coords.id == j].end + 1)))

        seq_counts = [0] * len(seq_intervals)

        for n in chr_loc[chrom]:
            for i, inter in enumerate(seq_intervals):
                if n in inter:
                    seq_counts[i] += 1

        target_counts.extend(seq_counts)

        for j in chr_coords.id:
            up_intervals.append(frozenset(range(chr_coords[chr_coords.id == j].start - 500,
                                                chr_coords[chr_coords.id == j].start + 1)))

        up_counts = [0] * len(up_intervals)

        for n in chr_loc[chrom]:
            for i, inter in enumerate(up_intervals):
                if n in inter:
                    up_counts[i] += 1

        up_target_counts.extend(up_counts)

        for j in chr_coords.id:
            down_intervals.append(frozenset(range(chr_coords[chr_coords.id == j].end,
                                                  chr_coords[chr_coords.id == j].end + 500 + 1)))

        down_counts = [0] * len(down_intervals)

        for n in chr_loc[chrom]:
            for i, inter in enumerate(down_intervals):
                if n in inter:
                    down_counts[i] += 1

        down_target_counts.extend(down_counts)

    target_counts = np.array(target_counts)
    up_target_counts = np.array(up_target_counts)
    down_target_counts = np.array(down_target_counts)

    genes_coords['target_seq_counts'] = target_counts
    genes_coords['500bp_up_target_seq_counts'] = up_target_counts
    genes_coords['500bp_down_target_seq_counts'] = down_target_counts

    genes_coords['length'] = genes_coords['end'] + 1 - genes_coords['start']
    genes_coords['target_per_bp'] = genes_coords['target_seq_counts'] / genes_coords['length']
    genes_coords['up_target_per_bp'] = genes_coords['500bp_up_target_seq_counts'] / 500
    genes_coords['down_target_per_bp'] = genes_coords['500bp_down_target_seq_counts'] / 500

    print(genes_coords.head())
    genes_coords.to_csv('target_seq_counts.csv')


# count_target_seq()

target_seq = pd.read_csv('target_seq_counts.csv')
translation = pd.read_csv('nomenclature_translation_SGD.csv')
tn_features = pd.read_csv('fr1_2_sorted__analysis.csv')

translated = target_seq.merge(translation, on='id')
final_data = tn_features.merge(translated, on='Standard name')
final_data.to_csv('fr1_2_target_sequence_with_all_features.csv')
