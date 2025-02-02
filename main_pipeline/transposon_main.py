import subprocess

import csv
import pysam
import gffutils
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import nt_search

from transposon_constants import *

# todo: add reference indexing optionality for different species
# wrap reference index in subprocess
# for now stored in dependencies directory

# todo: add cutudapt functionality
# remove tn sequences
# remove adapters

offset = 100
class_features = data_directory + file_prefix + '_' + str(offset) + '.csv'


def bash_pipeline():

    if arguments().sp == 'sc':
        ref = sc_ref
    elif arguments().sp == 'cg':
        ref = cg_ref
    elif arguments().sp == 'ca':
        ref = ca_ref
    elif arguments().sp == 'sp':
        ref = sp_ref
    else:
        raise ValueError('Unknown species flag specified')

    if not os.path.exists(data_directory):
        os.makedirs(data_directory)

    if os.path.exists(data_directory + sam_file):
        print('sam file found... moving on')

    if not os.path.exists(data_directory + sam_file):
        print('aligning reads...')
        if arguments().pe:
            subprocess.call(['bowtie2', '-p', '8', '--end-to-end', '-x', dependencies_dir + ref,
                                    '--fr', '-1', data_directory + arguments().fastq,
                                    '-2', data_directory + arguments().fastq2,
                                    '-S', data_directory + sam_file])
        else:
            subprocess.call(['bowtie2', '-p', '8', '--end-to-end', '-x', dependencies_dir + ref,
                                    '-U', data_directory + arguments().fastq,
                                    '-S', data_directory + sam_file])

    if os.path.exists(data_directory + bam_file):
        print('bam file found... moving on')

    if not os.path.exists(data_directory + bam_file):
        print('converting sam to bam...')
        bam_logfile = open(data_directory + bam_file, 'w')
        proc = subprocess.Popen(['samtools', 'view', '-@', '8', '-Sb', data_directory + sam_file], stdout=bam_logfile)
        proc.wait()

    if os.path.exists(data_directory + sorted_bam):
        print('sorted_bam file found... moving on')

    if not os.path.exists(data_directory + sorted_bam):
        print('sorting reads...')
        subprocess.call(['samtools', 'sort', '-@', '8', '-o', data_directory + sorted_bam,
                                data_directory + bam_file])

    if os.path.exists(data_directory + bam_index):
        print('bam index found... moving on')

    if not os.path.exists(data_directory + bam_index):
        print('indexing bam...')
        pysam.index(data_directory + sorted_bam)

    print('---- done! ----')


def get_features():

    """ Retrieves coordinates of all genes from a gffutils database """
    if arguments().sp == 'sc':
        database = sc_database
        features = sc_features
    elif arguments().sp == 'cg':
        database = cg_database
        features = cg_features
    elif arguments().sp == 'ca':
        database = ca_database
        features = ca_features
    elif arguments().sp == 'sp':
        database = sp_database
        features = sp_features
    else:
        raise ValueError('Unknown species flag specified')

    if not os.path.exists(dependencies_dir + database):
        print('creating database')
        gffutils.create_db(features, dependencies_dir + database)

    db = gffutils.FeatureDB(dependencies_dir + database)

    if arguments().sp == 'cg':
        genes = db.all_features(featuretype='ORF')
    elif arguments().sp == 'sc':
        genes = db.all_features(featuretype='gene')
    elif arguments().sp == 'ca':
        genes = db.all_features(featuretype='gene')
    elif arguments().sp == 'sp':
        genes = db.all_features(featuretype='gene')
    else:
        raise ValueError('Unknown species flag specified')

    print('extracting genomic features')
    gene_chroms = []
    gene_ids = []
    gene_starts = []
    gene_ends = []
    for gene in genes:
        gene_chroms.append(gene.chrom)
        if arguments().sp == 'cg':
            gene_ids.append(str(gene.id).split(':', 2)[1])
        else:
            gene_ids.append(gene.id)
        gene_starts.append(gene.start)
        gene_ends.append(gene.end)

    gene_coords_df = pd.DataFrame(list(zip(gene_ids, gene_chroms, gene_starts, gene_ends)), columns=['Standard name',
                                                                                                     'chrom',
                                                                                                     'start', 'end'])
    if arguments().sp == 'cg':
        gene_coords_df.to_csv(cg_features_file)
    elif arguments().sp == 'sc':
        gene_coords_df.to_csv(sc_features_file)
    elif arguments().sp == 'ca':
        gene_coords_df.to_csv(ca_features_file)
    elif arguments().sp == 'sp':
        gene_coords_df.to_csv(sp_features_file)
    else:
        raise ValueError('Unknown species flag specified')

    return gene_coords_df


def glabrata_hitmap():

    """ Generates hitmap on chromosomes from a bam file"""

    min_mapq = arguments().q
    bam = pysam.AlignmentFile(data_directory + sorted_bam, "rb")

    if arguments().sp == 'cg':
        features_file = cg_features_file
    elif arguments().sp == 'sc':
        features_file = sc_features_file
    elif arguments().sp == 'ca':
        features_file = ca_features_file
    elif arguments().sp == 'sp':
        features_file = sp_features_file
    else:
        raise ValueError('Unknown species flag specified')

    if not os.path.exists(features_file):
        get_features()

    features = pd.read_csv(features_file)

    if arguments().sp == 'cg':
        gene_chroms = features['chrom'].unique()
        hit_map = {chrom: {'W': {}, 'C': {}} for chrom in gene_chroms}
    elif arguments().sp == 'sc':
        gene_chroms = ['ref|NC_001133|', 'ref|NC_001134|', 'ref|NC_001135|', 'ref|NC_001136|', 'ref|NC_001137|',
                       'ref|NC_001138|', 'ref|NC_001139|', 'ref|NC_001140|', 'ref|NC_001141|', 'ref|NC_001142|',
                       'ref|NC_001143|', 'ref|NC_001144|', 'ref|NC_001145|', 'ref|NC_001146|', 'ref|NC_001147|',
                       'ref|NC_001148|', 'ref|NC_001224|']
        hit_map = {chrom: {'W': {}, 'C': {}} for chrom in gene_chroms}
    elif arguments().sp == 'ca':
        gene_chroms = features['chrom'].unique()
        hit_map = {chrom: {'W': {}, 'C': {}} for chrom in gene_chroms}
    elif arguments().sp == 'sp':
        gene_chroms = features['chrom'].unique()
        hit_map = {chrom: {'W': {}, 'C': {}} for chrom in gene_chroms}
    else:
        raise ValueError('Unknown species flag specified')

    print('mapping hits')

    for line in bam:
        if line.mapq < min_mapq:
            continue

        raw_chrom = bam.getrname(line.reference_id)

        if raw_chrom not in hit_map:
            continue

        if raw_chrom not in gene_chroms:
            continue

        # Since start < end always, in alignments which are reversed (along the
        # Crick strand) the start of the fragment is actually at the 'end' point.
        if line.is_reverse:
            pos = line.reference_end
            strand = 'C'
        else:
            # BAM files use 0-based indexing, and we work in 1-based indexing,
            # so we have to add one.
            pos = line.reference_start + 1
            strand = 'W'

        hit_map[raw_chrom][strand][pos] = hit_map[raw_chrom][strand].get(pos, 0) + 1

    with open(os.path.splitext(data_directory + sorted_bam)[0] + "_Hits.csv", "wb") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Chromosome", "Strand", "Position", "Reads"])
        for chrom in sorted(hit_map.keys()):
            for strand in hit_map[chrom].keys():
                for pos in sorted(hit_map[chrom][strand].keys()):
                    writer.writerow([chrom, strand, pos, hit_map[chrom][strand][pos]])

    sc_trans()


def target_seq(genome, query, filetype='fasta'):

    """ Finds a target sequence on a genome """

    print('finding target sequences')

    if arguments().sp == 'cg':
        features_file = cg_features_file
    elif arguments().sp == 'sc':
        features_file = sc_features_file
    elif arguments().sp == 'ca':
        features_file = ca_features_file
    elif arguments().sp == 'sp':
        features_file = sp_features_file
    else:
        raise ValueError('Unknown species flag specified')

    chromes = list(pd.read_csv(features_file)['chrom'].unique())
    if arguments().sp == 'ca':
        chromes = [str(chrom)[:9] for chrom in chromes if str(chrom)[8] == 'A']
    if arguments().sp == 'sp':
        chromes = ['I', 'II', 'III']

    chroms_locs = pd.DataFrame(np.nan, index=range(250000), columns=chromes)
    count = 0
    total_length = 0
    i = 0
    for record in SeqIO.parse(genome, filetype):
        if arguments().sp == 'sp':
            if record.id in chromes:
                x = nt_search(str(record.seq), query)
                count += len(x[1:])
                total_length += len(record.seq)
                chroms_locs[chromes[i]] = pd.Series(x[1:])
                i += 1
            else:
                pass
        else:
            x = nt_search(str(record.seq), query)
            count += len(x[1:])
            total_length += len(record.seq)
            chroms_locs[chromes[i]] = pd.Series(x[1:])
            i += 1

    chroms_locs = chroms_locs.dropna(how='all')

    if arguments().sp == 'cg':
        chroms_locs.to_csv(dependencies_dir + cg_hermes_on_chr)
    elif arguments().sp == 'sc':
        chroms_locs.to_csv(dependencies_dir + sc_hermes_on_chr)
    elif arguments().sp == 'ca':
        chroms_locs.to_csv(dependencies_dir + ca_hermes_on_chr)
    elif arguments().sp == 'sp':
        chroms_locs.to_csv(dependencies_dir + sp_hermes_on_chr)
    else:
        raise ValueError('Unknown species flag specified')

    return chroms_locs


def classifier_features(offset):

    """ Generates features for the classifier """

    if arguments().sp == 'cg':
        chr_loc = pd.read_csv(dependencies_dir + cg_hermes_on_chr, index_col=0)
        genes_coords = pd.read_csv(cg_features_file, index_col=0)
    elif arguments().sp == 'sc':
        chr_loc = pd.read_csv(dependencies_dir + sc_hermes_on_chr, index_col=0)
        genes_coords = pd.read_csv(sc_features_file, index_col=0)
    elif arguments().sp == 'ca':
        chr_loc = pd.read_csv(dependencies_dir + ca_hermes_on_chr, index_col=0)
        genes_coords = pd.read_csv(ca_features_file, index_col=0)
        chromes = list(genes_coords['chrom'])
        genes_coords['chrom'] = [str(chrom)[:9] for chrom in chromes]
        genes_coords = genes_coords[genes_coords['chrom'].str.contains('A')]
    elif arguments().sp == 'sp':
        chr_loc = pd.read_csv(dependencies_dir + sp_hermes_on_chr, index_col=0)
        genes_coords = pd.read_csv(sp_features_file, index_col=0).rename(columns={'Standard name': 'id'})
        genes_coords['id'] = genes_coords['id'].str[5:]
        genes_coords = genes_coords[(genes_coords['chrom'] == 'I') | (genes_coords['chrom'] == 'II') |
                                    (genes_coords['chrom'] == 'III')]
    else:
        raise ValueError('Unknown species flag specified')

    hit_map2 = pd.read_csv(os.path.splitext(data_directory + sorted_bam)[0] + '_Hits.csv', index_col=0)
    if arguments().sp == 'sp':
        hit_map2['Chromosome'] = hit_map2['Chromosome'].str[3:]

    if arguments().sp == 'ca':
        hit_map_chr = list(hit_map2['Chromosome'].unique())
        chroms = list(genes_coords['chrom'].unique())
        chroms_tra = dict(zip(hit_map_chr, chroms))
        hit_map2['chrom'] = hit_map2['Chromosome'].map(chroms_tra)
        hit_map2['Chromosome'] = hit_map2['chrom']
        hit_map2.pop('chrom')

    chroms = list(genes_coords['chrom'].unique())
    if arguments().sp == 'sp':
        hit_map2 = hit_map2[(hit_map2['Chromosome'] == 'I') | (hit_map2['Chromosome'] == 'II') |
                            (hit_map2['Chromosome'] == 'III')]

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

    print('generating features')

    for chrom in chroms:

        chr_coords = genes_coords[genes_coords['chrom'] == chrom]
        max_free.extend((max_free_in_range(chr_coords, hit_map2[hit_map2['Chromosome'] == chrom]['Position'], 'start',
                                           'end')))
        reads.extend(count_reads_in_range(chr_coords, hit_map2[hit_map2['Chromosome'] == chrom], 'start', 'end'))
        neighbourhood.extend(count_in_range(chr_coords, hit_map2[hit_map2['Chromosome'] == chrom]['Position'], 'start',
                                            'end', s_offset=-neigh_offset, e_offset=neigh_offset))
        hits.extend(count_in_range(chr_coords, hit_map2[hit_map2['Chromosome'] == chrom]['Position'], 'start', 'end'))

        hits_up.extend(count_in_range(chr_coords, hit_map2[hit_map2['Chromosome'] == chrom]['Position'], 'start',
                                      'start', s_offset=-offset))
        hits_down.extend(count_in_range(chr_coords, hit_map2[hit_map2['Chromosome'] == chrom]['Position'], 'end',
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
    genes_coords['Neighbourhood Index 2'] = ((genes_coords['Hits_' + str(offset) + '_bp_upstream']
                                              + genes_coords['Hits']
                                              + genes_coords['Hits_' + str(offset) + '_bp_downstream'])
                                              / genes_coords['Length']) / \
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
    print('----------Done!---------')
    genes_coords.to_csv(class_features, index=False)


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


def sc_trans():

    data = pd.read_csv(os.path.splitext(data_directory + sorted_bam)[0] + "_Hits.csv")
    chroms = pd.read_csv(dependencies_dir + sc_hermes_on_chr)
    chroms = chroms.columns.values[1:]
    alct_chroms = data['Chromosome'].unique()
    chrom_dict = dict(zip(alct_chroms, chroms))

    for key, value in chrom_dict.items():
        data['Chromosome'] = data['Chromosome'].replace(key, value)

    data.to_csv(os.path.splitext(data_directory + sorted_bam)[0] + "_Hits.csv")


if __name__ == '__main__':
    if arguments().sp == 'cg':
        features_file = cg_features_file
    elif arguments().sp == 'sc':
        features_file = sc_features_file
    elif arguments().sp == 'ca':
        features_file = ca_features_file
    elif arguments().sp == 'sp':
        features_file = sp_features_file

    if not os.path.exists(data_directory + bam_index):
        bash_pipeline()

    if not os.path.exists(features_file):
        get_features()

    if not os.path.exists(os.path.splitext(data_directory + sorted_bam)[0] + '_Hits.csv'):
        glabrata_hitmap()

    if arguments().sp == 'cg':
        if not os.path.exists(dependencies_dir + cg_hermes_on_chr):
            target_seq(cg_ref_file, target_Hermes)
    if arguments().sp == 'sc':
        if not os.path.exists(dependencies_dir + sc_hermes_on_chr):
            target_seq(sc_ref_file, target_Hermes)
    if arguments().sp == 'ca':
        if not os.path.exists(dependencies_dir + ca_hermes_on_chr):
            target_seq(ca_ref_file, target_PiggyBac)
    if arguments().sp == 'sp':
        if not os.path.exists(dependencies_dir + sp_hermes_on_chr):
            target_seq(sp_ref_file, target_PiggyBac)

    if not os.path.exists(class_features):
        classifier_features(offset)
