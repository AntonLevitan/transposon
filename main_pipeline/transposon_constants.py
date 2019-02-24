import os
import argparse

species_flag = '--sp'
paired_end_flag = '--pe'
fastq_flag = '--fastq'
fastq2_flag = '--fastq2'
quality_flag = '-q'

data_directory = 'data' + os.sep
dependencies_dir = 'dependencies' + os.sep

fastq_suffix = '.fastq'
sam_suffix = '.sam'
bam_suffix = '.bam'
bai_suffix = '.bai'

sc_ref = 'sc'
cg_ref = 'cg'
ca_ref = 'ca'
sp_ref = 'sp'

cg_database = 'cg_gffutils_database'
ca_database = 'C_albicans_SC5314_version_A22-s07-m01-r08_features.gff.gffutils_db.sqlite'
sc_database = 'sc_gffutils_database'
sp_database = 'schizosaccharomyces_pombe.chr.gff3.gffutils_db.sqlite'

cg_features = dependencies_dir + 'GlabrataFeatures.gff3'
ca_features = dependencies_dir + 'C_albicans_SC5314_A22_current_features.gff'
sc_features = dependencies_dir + 'saccharomyces_cerevisiae_R64-2-1_20150113.gff'
sp_features = dependencies_dir + 'schizosaccharomyces_pombe.chr.gff3'

cg_reference = 'C_glabrata_CBS138_current_chromosomes.fasta'
sc_reference = 'S288C_reference_sequence_R64-2-1_20150113.fsa'
ca_reference = 'C_albicans_SC5314_version_A22-s07-m01-r08_chromosomes_HapA.fasta'
sp_reference = 'Schizosaccharomyces_pombe.ASM294v2.30.dna.genome.fa'

target_Hermes = 'NTNNNNAN'
target_PiggyBac = 'TTAA'

cg_ref_file = dependencies_dir + cg_reference
ca_ref_file = dependencies_dir + ca_reference
sc_ref_file = dependencies_dir + sc_reference
sp_ref_file = dependencies_dir + sp_reference

cg_features_file = dependencies_dir + cg_ref + '_features_coords.csv'
ca_features_file = dependencies_dir + ca_ref + '_features_coords.csv'
sc_features_file = dependencies_dir + sc_ref + '_features_coords.csv'
sp_features_file = dependencies_dir + sp_ref + '_features_coords.csv'

cg_hermes_on_chr = 'cg_hermes_seq_locations_on_chromosomes.csv'
ca_hermes_on_chr = 'ca_piggybac_seq_locations_on_chromosomes.csv'
sc_hermes_on_chr = 'sc_hermes_seq_locations_on_chromosomes.csv'
sp_hermes_on_chr = 'sp_piggybac_seq_locations_on_chromosomes.csv'


def arguments():

    parser = argparse.ArgumentParser()
    parser.add_argument(species_flag, dest='sp', required=True)
    parser.add_argument(paired_end_flag, dest='pe')
    parser.add_argument(fastq_flag, dest='fastq', required=True)
    parser.add_argument(fastq2_flag, dest='fastq2')
    parser.add_argument(quality_flag, type=int, default=20, help='The minimum mapping quality to consider.')
    return parser.parse_args()


file_prefix = str(arguments().fastq)[:-6]
sam_file = file_prefix + sam_suffix
bam_file = file_prefix + bam_suffix
sorted_bam = file_prefix + '_sorted' + bam_suffix
bam_index = sorted_bam + bai_suffix
