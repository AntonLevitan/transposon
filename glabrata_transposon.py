import pysam
import argparse
import os
import csv
# import GenomicFeatures
import gffutils
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.SeqUtils import nt_search

cg_database = 'cg_gffutils_database'

query = 'NTNNNNAN'
seq_file = 'C_glabrata_CBS138_current_chromosomes.fasta'


def get_features(gffutils_db):

    """ Retrieves coordinates of all genes from a gffutils database """

    gene_chroms = []
    gene_ids = []
    gene_starts = []
    gene_ends = []
    for gene in gffutils_db:
        gene_chroms.append(str(gene.chrom).split('_', 1)[0])
        gene_ids.append(str(gene.id).split(':', 2)[1])
        gene_starts.append(gene.start + 1)
        gene_ends.append(gene.end)

    gene_coords_df = pd.DataFrame(list(zip(gene_ids, gene_chroms, gene_starts, gene_ends)), columns=['id', 'chrom',
                                                                                                     'start', 'end'])
    gene_coords_df.to_csv('cg_features_coords.csv')

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
    chroms_locs.to_csv('gc_seq_locations_on_chromosomes.csv')

    return chroms_locs, sc_features


# freq_appearance(seq_file, query, cg_all_features)


def glabrata_hitmap():

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--bam", required=True)
    parser.add_argument("-q", type=int, default=20, help="The minimum mapping quality to consider.")

    args = parser.parse_args()
    min_mapq = args.q
    bam = pysam.AlignmentFile(args.bam, "rb")

    if not os.path.exists(cg_database):
        gffutils.create_db('GlabrataFeatures.gff3', cg_database)

    cg_db = gffutils.FeatureDB(cg_database)
    cg_genes = cg_db.all_features(featuretype='ORF')

    gene_chroms = []
    gene_ids = []
    gene_starts = []
    gene_ends = []
    for gene in cg_genes:
        gene_chroms.append(str(gene.chrom).split('_', 1)[0])
        gene_ids.append(str(gene.id).split(':', 2)[1])
        gene_starts.append(gene.start + 1)
        gene_ends.append(gene.end)

    gene_coords_df = pd.DataFrame(list(zip(gene_ids, gene_chroms, gene_starts, gene_ends)), columns=['id', 'chrom',
                                                                                                     'start', 'end'])
    gene_coords_df.to_csv('cg_features_coords.csv')

    for line in bam:
        if line.mapq < min_mapq:
            continue

        raw_chrom = bam.getrname(line.reference_id)

        #  TODO: chrom_names? maybe chrom directly from database? check back after creating database.
        # if implement database here, change it.

        if raw_chrom not in set(gene_chroms):
            continue

        chrom = str(gene.chrom).split('_', 1)[0]
        if chrom not in set(gene_chroms):
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

        hit_map[chrom][strand][pos] = hit_map[chrom][strand].get(pos, 0) + 1

    with open(os.path.splitext(args.bam)[0] + "_Hits.csv", "wb") as out_file:
        writer = csv.writer(out_file)
        writer.writerow(["Chromosome", "Strand", "Position", "Reads", "Gene"])
        for chrom in sorted(hit_map.keys()):
            for strand in hit_map[chrom].keys():
                for pos in sorted(hit_map[chrom][strand].keys()):
                    # features = pom_db.get_features_at_location(chrom, pos)
                    features = []  # It appears finding the hit gene is not necessary at this point

                    if len(features) > 2:
                        print "More than 1 feature at position", chrom, pos

                    writer.writerow([chrom, strand, pos,
                                     hit_map[chrom][strand][pos],
                                     "nan" if not features else features[0].standard_name])

glabrata_hitmap()