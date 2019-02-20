from Bio import SeqIO
cg_reference = 'Schizosaccharomyces_pombe_all_chromosomes.fa'

length = 0
for record in SeqIO.parse(cg_reference, 'fasta'):
    length += len(record)
    # print len(record)
    print record.id
# print length
