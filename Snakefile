#!/usr/bin/env python3

import pathlib

#############
# FUNCTIONS # 
#############

def resolve_path(x):
    return str(pathlib.Path(x).resolve())

#########
# RULES #
#########

rule target:
    input:
        'output/040_process/results.out.bed' 

rule process:
	input: 
		'data/ncbi-genomes-2018-04-30/GCF_000002195.4_Amel_4.5_genomic.fna.nsd',
		scrm = 'output/030_scrmshaw/hits/hexmcd',
		db = 'data/ncbi-genomes-2018-04-30/GCF_000002195.4_Amel_4.5_genomic.fna'
	output: 
		hits = 'output/040_process/all_hits.txt', 
		top = 'output/040_process/top_hits.txt',
		csv = 'output/040_process/results.out.csv',
		bed = 'output/040_process/results.out.bed'

	log: 
		'output/logs/040_process'
	run: 
		shell(
			#grab all the hit files and put them into one file 
			'cat {input.scrm}/*/*.hits.ranked > {output.hits} ; '
			#filter top hits only 
			'grep -A 1 "lr=[123]," {output.hits} --no-group-separator '
			'| tr -d "[[:blank:]]" > {output.top} ; '
			#blast search 
			'blastn -db {input.db} -query {output.top} -out {output.csv} -outfmt "10 sacc sstart send" ; '
			#make it a .bed file for IGV 
			'<{output.csv} tr "," "\t" > {output.bed} '
		)


rule scrmshaw:
    input:
        fa = 'output/020_remove_repeats/masked_chromosomes.fa',
        genes = 'output/010_ref/genes.txt',
        exons = 'output/010_ref/exons.txt',
        traindirlst = 'data/dros/trainSet' 
        #file containing the paths (from apis-numb) to the training sets 
    output:
        'output/030_scrmshaw/hits/hexmcd/mapping0.ap'
    params:
        outdir = 'output/030_scrmshaw'
    log:
        'output/logs/030_scrmshaw.log'
    shell:
        'SCRMshaw/code/scrm.pl '
        '--thitg 300 --imm --hexmcd --pac '
        '--genome {input.fa} '
        '--exon {input.exons} '
        '--gene {input.genes} '
        '--outdir {params.outdir} '
        '--traindirlst {input.traindirlst} '
        '&> {log}'

rule rename_output: #rename Group1.4 to Group1p4 to get around SCRMshaw bug
   input:
       fa = 'output/020_remove_repeats/chromosomes.fa.2.7.7.80.10.50.500.mask',
       genes = 'output/010_ref/genes.txt',
       exons = 'output/010_ref/exons.txt'
   output:
       fa = 'output/020_remove_repeats/masked_chromosomes.fa'
   shell: #change dots to p to make SCRMshaw happy
       "sed -e 's/\./p/g'  {input.fa} > {output.fa}; " 
       "sed -i -e 's/\./p/g' {input.genes}; "
       "sed -i -e 's/\./p/g' {input.exons}; "

rule remove_repeats:
    input:
        fa = 'output/010_ref/chromosomes.fa'
    output: 
        fa = temp('output/020_remove_repeats/'
                  'chromosomes.fa.2.7.7.80.10.50.500.mask')
    params: 
        wd = 'output/020_remove_repeats'
    run: 
        my_fasta = resolve_path(input.fa)
        #exit 0 because trf409.linux64 returns error 2
        shell(
              'cd {params.wd} || exit 1 ; '
              '../../trf409.linux64 '
              '{my_fasta} '
              '2 7 7 80 10 50 500 -m -h || exit 0 '
            )


rule extract_chromosomes:
    input:
        fa = 'data/Amel_4.5_scaffolds.fa'
    output:
        fa = 'output/010_ref/chromosomes.fa'
    threads: 
        1
    script:
        'src/extract_chromosomes.py'


rule extract_genes_and_exons:
    input:
        gff = 'data/amel_OGSv3.2.gff3'
    output:
        genes = 'output/010_ref/genes.txt',
        exons = 'output/010_ref/exons.txt'
    log:
        log = 'output/logs/010_extract_genes_and_exons.log'
    threads: 
        1
    script:
        'src/extract_genes_and_exons.R'


rule make_blast_db: #this should probably be done upstream but idk how. Also, outputs to /data folder. 
	input: #from ncbi- if the output is being fed into IGV the genomes need to be the same
		db = 'data/ncbi-genomes-2018-04-30/GCF_000002195.4_Amel_4.5_genomic.fna' 
	output: 
		db = 'data/ncbi-genomes-2018-04-30/GCF_000002195.4_Amel_4.5_genomic.fna.nsd'
	shell:
		'makeblastdb -in {input.db} -parse_seqids -dbtype nucl '
