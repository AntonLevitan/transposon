import os
import subprocess
import argparse
import pysam


PAIRED_END_FLAG = '--pe'

DATA_DIRECTORY = 'WtReadsSet1_180906' + os.sep
FASTQ_SUFFIX = '.fastq'
SAM_SUFFIX = '.sam'
BAM_SUFFIX = '.bam'
BAI_SUFFIX = '.bai'

FILE_PREFIX = 'FreadsWT'

SAM_FILE = FILE_PREFIX + SAM_SUFFIX
BAM_FILE = FILE_PREFIX + BAM_SUFFIX
SORTED_BAM = FILE_PREFIX + '_sorted' + BAM_SUFFIX
BAM_INDEX = SORTED_BAM + BAI_SUFFIX

# TODO: add reference indexing optionality for different species
# wrap reference index in subprocess

SC_REF = 'sc'
CG_REF = 'cg'


def arguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(PAIRED_END_FLAG, dest='pe', action='store_true')
    parser.add_argument('--fastq', required=True)

    return parser.parse_args()


if not os.path.exists(DATA_DIRECTORY):
    os.makedirs(DATA_DIRECTORY)

if not os.path.exists(DATA_DIRECTORY + SAM_FILE):
    print('Aligning reads...')
    if arguments().pe:
        FASTQ_2 = input('Enter reverse fastq: ') + FASTQ_SUFFIX
        proc = subprocess.call(['bowtie2', '-p', '8', '--end-to-end', '-x', DATA_DIRECTORY + CG_REF,
                                '--fr', '-1', DATA_DIRECTORY + arguments().fastq, '-2', DATA_DIRECTORY + FASTQ_2,
                                '-S', DATA_DIRECTORY + SAM_FILE])
    else:
        proc = subprocess.call(['bowtie2', '-p', '8', '--end-to-end', '-x', DATA_DIRECTORY + CG_REF,
                                '-U', DATA_DIRECTORY + arguments().fastq,
                                '-S', DATA_DIRECTORY + SAM_FILE])

if not os.path.exists(DATA_DIRECTORY + BAM_FILE):    
    print('Converting sam to bam...')
    bam_logfile = open(DATA_DIRECTORY + BAM_FILE, 'w')
    proc = subprocess.Popen(['samtools', 'view', '-@', '8', '-Sb', DATA_DIRECTORY + SAM_FILE], stdout=bam_logfile)
    proc.wait()
    proc.kill()

if not os.path.exists(DATA_DIRECTORY + SORTED_BAM):
    print('Sorting reads...')
    proc = subprocess.call(['samtools', 'sort', '-@', '8', '-o', DATA_DIRECTORY + SORTED_BAM, DATA_DIRECTORY + BAM_FILE])

if not os.path.exists(DATA_DIRECTORY + BAM_INDEX):
    print('Indexing bam...')
    pysam.index(DATA_DIRECTORY + SORTED_BAM)

print('---- Done! ----')
