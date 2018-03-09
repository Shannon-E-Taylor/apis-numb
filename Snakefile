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
        'output/030_scrmshaw/hits'

rule scrmshaw:
    input:
        fa = 'output/020_remove_repeats/masked_chromosomes.fa',
        genes = 'output/010_ref/genes.txt',
        exons = 'output/010_ref/exons.txt',
        traindirlst = ('data/data2generateGBEresults/'
                       'data2generateGBEresults/'
                       'data/CRM.train/trainSet')
    output:
        'output/030_scrmshaw/hits'
    params:
        outdir = 'output/030_scrmshaw',
        universeMap = ('data/data2generateGBEresults/'
                       'data2generateGBEresults/data/universalGeneSet')

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
        '--universeMap {params.universeMap} '
        '&> {log}'

rule rename_output:
    input: 
        fa = 'output/020_remove_repeats/chromosomes.fa.2.7.7.80.10.50.500.mask'
    output: 
        fa = 'output/020_remove_repeats/masked_chromosomes.fa'
    shell: 
        'cp {input.fa} {output.fa}'


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
              'trf409.linux64 '
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



