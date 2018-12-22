from Bio.SeqUtils import nt_search
from Bio import SeqIO
import gffutils
import pandas as pd
import csv
import numpy as np

query = 'NTNNNNAN'
# seq_file = '/home/user/Desktop/transposon/SCTrainingSet/WT_1/FreadI1_trimmed.fastq'
seq_file = 'WtReadsSet1_180906/C_glabrata_CBS138_current_chromosomes.fasta'

# gffutils creating database
# db = gffutils.create_db('/home/anton/Desktop/transposon/GlabrataFeatures.gff3',
#                         'cg_gffutils_database')

db = gffutils.FeatureDB('cg_gffutils_database')
cg_genes = db.all_features(featuretype='ORF')


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
        if 'C_glabrata' in gene.chrom:
            gene_ids.append(str(gene.id).split(':', 2)[1])
        else:
            gene_ids.append(gene.id)
        gene_starts.append(gene.start)
        gene_ends.append(gene.end)

    gene_coords_df = pd.DataFrame(list(zip(gene_ids, gene_chroms, gene_starts, gene_ends)), columns=['id', 'chrom',
                                                                                                     'start', 'end'])
    gene_coords_df.to_csv('cg_features_coords.csv')

    return gene_coords_df


def freq_appearance(file, query, gffutils_db, filetype='fasta'):
    """ Frequency of occurrence of a query sequence in a file

    :param file: sequence file
    :param query: sequence query (string)
    :param gffutils_db: gffutils all_features database
    :param filetype: parameter for SeqIO.parse (default='fasta')

    """
    chromes = list(pd.read_csv('cg_features_coords.csv')['chrom'].unique())
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

    chroms_locs.to_csv('cg_seq_locations_on_chromosomes.csv')

    return chroms_locs


# freq_appearance(seq_file, query, cg_genes)

def count_target_seq(file, query, genes):
    # chr_loc = pd.read_csv('seq_locations_on_chromosomes.csv').drop('Unnamed: 0', axis=1)
    # genes_coords = pd.read_csv('sc_features_coords.csv').drop('Unnamed: 0', axis=1)
    # chr_loc = freq_appearance(file, query, genes)

    chr_loc = pd.read_csv('cg_seq_locations_on_chromosomes.csv', index_col=0)
    genes_coords = pd.read_csv('cg_features_coords.csv', index_col=0)
    hit_map = pd.read_csv('FreadsWT_sorted_Hits.csv')
    chroms = list(genes_coords['chrom'].unique())
    offset = 600
    neigh_offset = 10000
    target_counts = []
    up_target_counts = []
    down_target_counts = []
    hits = []
    hits_up = []
    hits_down = []
    neighbourhood = []
    reads = []
    max_free = []

    for chrom in chroms:

        chr_coords = genes_coords[genes_coords['chrom'] == chrom]
        max_free.extend((max_free_in_range(chr_coords, hit_map[hit_map['Chromosome'] == chrom]['Position'], 'start',
                                           'end')))
        reads.extend(count_reads_in_range(chr_coords, hit_map[hit_map['Chromosome'] == chrom], 'start', 'end'))
        neighbourhood.extend(count_in_range(chr_coords, hit_map[hit_map['Chromosome'] == chrom]['Position'], 'start',
                                            'end', s_offset=-neigh_offset, e_offset=neigh_offset))
        hits.extend(count_in_range(chr_coords, hit_map[hit_map['Chromosome'] == chrom]['Position'], 'start', 'end'))

        hits_up.extend(count_in_range(chr_coords, hit_map[hit_map['Chromosome'] == chrom]['Position'], 'start',
                                      'start', s_offset=-offset))
        hits_down.extend(count_in_range(chr_coords, hit_map[hit_map['Chromosome'] == chrom]['Position'], 'end',
                                        'end', e_offset=offset))

        target_counts.extend(count_in_range(chr_coords, chr_loc[chrom], 'start', 'end'))
        up_target_counts.extend(count_in_range(chr_coords, chr_loc[chrom], 'start', 'start', s_offset=-offset))
        down_target_counts.extend(count_in_range(chr_coords, chr_loc[chrom], 'end', 'end', e_offset=offset))

    genes_coords['Length'] = genes_coords['end'] + 1 - genes_coords['start']
    genes_coords['Hits'] = hits
    genes_coords['Reads'] = reads
    genes_coords['Hits_' + str(offset) + '_bp_upstream'] = hits_up
    genes_coords['Hits_' + str(offset) + '_bp_downstream'] = hits_down
    genes_coords['Neighbourhood Hits'] = neighbourhood
    genes_coords['Neighbourhood Index'] = (genes_coords['Hits'] / genes_coords['Length']) / \
                                          (genes_coords['Neighbourhood Hits'] / neigh_offset * 2)
    genes_coords['Max Free Region'] = max_free
    genes_coords['Freedom Index'] = genes_coords['Max Free Region'] / genes_coords['Length']
    genes_coords['target_seq_counts'] = target_counts
    genes_coords[str(offset) + '_bp_up_target_seq_counts'] = up_target_counts
    genes_coords[str(offset) + '_bp_down_target_seq_counts'] = down_target_counts
    genes_coords['target_per_bp'] = genes_coords['target_seq_counts'] / genes_coords['Length']
    genes_coords[str(offset) + '_up_target_per_bp'] = genes_coords[str(offset) +
                                                                   '_bp_up_target_seq_counts'] / offset
    genes_coords[str(offset) + '_down_target_per_bp'] = genes_coords[str(offset) +
                                                                     '_bp_down_target_seq_counts'] / offset

    genes_coords = genes_coords.round({'Neighbourhood Index': 3, 'Freedom Index': 3, 'target_per_bp': 4,
                                       str(offset) + '_up_target_per_bp': 4,
                                       str(offset) + '_down_target_per_bp': 4})
    print(genes_coords.head())
    genes_coords.to_csv('cg_target_seq_counts_test.csv')


def count_in_range(coords, segment, start, end, s_offset=0, e_offset=1):

    intervals = []
    for j in coords.id:
        intervals.append(frozenset(range(coords[coords.id == j][start] + s_offset,
                                         coords[coords.id == j][end] + e_offset)))

    seq_counts = [0] * len(intervals)
    for n in segment:
        for i, inter in enumerate(intervals):
            if n in inter:
                seq_counts[i] += 1
    return seq_counts


def count_reads_in_range(coords, segment, start, end, s_offset=0, e_offset=1):
    reads = segment.set_index(['Position', 'Strand'])
    segment = segment['Position']
    intervals = []

    for j in coords.id:
        intervals.append(frozenset(range(coords[coords.id == j][start] + s_offset,
                                         coords[coords.id == j][end] + e_offset)))

    seq_counts = [0] * len(intervals)
    for n in segment:
        for i, inter in enumerate(intervals):
            if n in inter:
                x = reads.loc[n, 'Reads'].sum()
                seq_counts[i] += x
    return seq_counts


def max_free_in_range(coords, segment, start, end, s_offset=0, e_offset=1):

    intervals = []
    hit_locs = [[] for x in range(len(coords.id))]
    dist = [[] for x in range(len(coords.id))]
    k = 0
    for j in coords.id:
        intervals.append(frozenset(range(coords[coords.id == j][start] + s_offset,
                                         coords[coords.id == j][end] + e_offset)))
        hit_locs[k].append(int(coords[coords.id == j][start].to_string(index=False)))
        hit_locs[k].append(int(coords[coords.id == j][end].to_string(index=False)) + 1)
        k += 1

    for n in segment:
        for i, inter in enumerate(intervals):
            if n in inter:
                hit_locs[i].append(int(n))

    for i in range(len(hit_locs)):
        hit_locs[i] = sorted(hit_locs[i])

    n = 0
    max_dist = []
    for i in hit_locs:
        for j in range(1, len(i)):
            dist[n].append(abs(i[j] - i[j-1]))
        max_dist.append(max(dist[n]))
        n += 1
    return max_dist


count_target_seq(seq_file, query, cg_genes)

# target_seq = pd.read_csv('target_seq_counts.csv')
# translation = pd.read_csv('nomenclature_translation_SGD.csv')
# tn_features = pd.read_csv('fr1_2_sorted__analysis.csv')
#
# translated = target_seq.merge(translation, on='id')
# final_data = tn_features.merge(translated, on='Standard name')
# final_data.to_csv('fr1_2_target_sequence_with_all_features.csv')
