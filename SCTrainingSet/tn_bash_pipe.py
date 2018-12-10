import os
import subprocess
import argparse

PAIRED_END_FLAG = '--pe'

FASTQ_SUFFIX = '.fastq'
SAM_SUFFIX = '.sam'
BAM_SUFFIX = '.bam'
BAI_SUFFIX = '.bai'

FASTQ_1 = input('Enter forward fastq: ') + FASTQ_SUFFIX

FILE_PREFIX = input('Enter output filename: ')

SAM_FILE = FILE_PREFIX + SAM_SUFFIX
BAM_FILE = FILE_PREFIX + BAM_SUFFIX
SORTED_BAM = FILE_PREFIX + '_sorted_'+ BAM_SUFFIX
BAM_INDEX = SORTED_BAM + BAI_SUFFIX

# TODO: add reference indexing optionality for different species
# wrap reference index in subprocess

SC_REF = 'sc'


def arguments():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(PAIRED_END_FLAG, dest='pe', action='store_true')

    return parser.parse_args()


if not os.path.exists(SAM_FILE):
    print('Aligning reads...')
    if arguments().pe:
        FASTQ_2 = input('Enter reverse fastq: ') + FASTQ_SUFFIX
        proc = subprocess.call(['bowtie2', '-p', '8', '--end-to-end', '-x', SC_REF, '--fr', '-1', FASTQ_1, '-2', FASTQ_2, '-S', SAM_FILE])
    else:
        proc = subprocess.call(['bowtie2', '-p', '8', '--end-to-end', '-x', SC_REF, '-U', FASTQ_1, '-S', SAM_FILE])

if not os.path.exists(BAM_FILE):    
    print('Converting sam to bam...')
    bam_logfile = open(BAM_FILE, 'w')
    proc = subprocess.Popen(['samtools', 'view', '-Sb', SAM_FILE], stdout=bam_logfile)
    proc.wait()
    proc.kill()

if not os.path.exists(SORTED_BAM):
    print('Sorting reads...')
    proc = subprocess.call(['samtools', 'sort', '-o', SORTED_BAM, BAM_FILE])

if not os.path.exists(BAM_INDEX):
    print('Indexing bam...')
    proc = subprocess.call(['samtools', 'index', SORTED_BAM])
