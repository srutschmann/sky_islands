#! /usr/bin/env python

''''
- online blast sequences with blastx program and nr database against insect database
- filter with e-value < 0.01

- input file: sequence.fasta # needs to be specified in line 19
- call script: ./blast_sequence.py

last modified: 2016/08
'''

# load modules
import Bio
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# open fasta file
fasta_string = open("sequence.fasta").read()

# blast search
result_handle = NCBIWWW.qblast("blastx", "nr", fasta_string, ncbi_gi="TRUE", format_type = "XML", entrez_query = "txid50557[Organism:exp]") # use blastx, maybe add description (default is 500)
blast_records = NCBIXML.parse(result_handle) # get xml report from blast search

f = open('loci.blast.out', 'a')

# filter trough xml output file and get hits with e-values threshold into loci.blast.out file
E_VALUE_THRESH = 0.01
for record in blast_records:
  for alignment in record.alignments:
    for hsp in alignment.hsps:
      if hsp.expect < E_VALUE_THRESH:
        f.write("%s\t%d\t%g\t%s\n" % (record.query, record.alignments[0].length, record.alignments[0].hsps[0].expect, record.alignments[0].title.split('>')[0]))
        
f.close()

        
# filter trough xml output file and print hits with e-values threshold
# E_VALUE_THRESH = 0.01
# for record in blast_records:
#   for alignment in record.alignments:
#     for hsp in alignment.hsps:
#       if hsp.expect < E_VALUE_THRESH:
#         print('***Alignment***')
#         print('sequence:', alignment.title)
#         print('length:', alignment.length)
#         print('e value:', hsp.expect)
#         print(hsp.query)
#         print(hsp.match)
#         print(hsp.sbjct)
#         
