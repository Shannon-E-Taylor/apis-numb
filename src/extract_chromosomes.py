#!/usr/bin/env python3

from Bio import SeqIO

###########
# GLOBALS #
###########

fa = snakemake.input['fa']
output_fa = snakemake.output['fa']
#target_chr = ['Group1.4', 'Group3.4']

########
# MAIN #
########

scaffolds = [x for x in SeqIO.parse(fa, 'fasta')]
target_chr = [x.id for x in scaffolds if x.id.startswith('Group3.') 
              or x.id == 'Group1.4']
kept_scaffolds = [x for x in scaffolds if x.id in target_chr]
#kept_scaffolds = scaffolds
SeqIO.write(kept_scaffolds, output_fa, 'fasta')
