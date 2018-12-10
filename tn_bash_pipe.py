import os
import subprocess

FASTQ_SUFFIX = '.fastq'
SAM_SUFFIX = '.sam'
BAM_SUFFIX = '.bam'

file_preffix = input('Enter .fastq filename: ')

fastq_file = file_preffix + FASTQ_SUFFIX
sam_file = file_preffix + SAM_SUFFIX
bam_file = file_preffix + BAM_SUFFIX

if not os.path.exists(sam_file):
    count_logfile = open(sam_file, 'w')
    proc = subprocess.Popen(['bowtie2', '-p 8', '--end-to-end', '-x sc', '--fr -1 FreadI1_trimmed.fastq -2 RreadI1.fastq', '-S Both_readI1_trimmed.sam'], stdout=count_logfile)
    proc.wait()
    proc.kill()