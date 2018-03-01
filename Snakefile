#!/usr/bin/env python3

#########
# RULES #
#########

rule target:
    input:
        'output/020_scrmshaw/hits/hexmcd/mapping0.ap.hits'

rule scrmshaw:
    input:
        fa = 'output/010_ref/chromosomes.fa',
        genes = 'output/010_ref/genes.txt',
        exons = 'output/010_ref/exons.txt',
        traindirlst = ('data/data2generateGBEresults/'
                       'data2generateGBEresults/'
                       'data/CRM.train/trainSet')
    output:
        'output/020_scrmshaw/hits/hexmcd/mapping0.ap.hits'
    params:
        outdir = 'output/020_scrmshaw'
    log:
        'output/logs/020_scrmshaw.log'
    shell:
        'SCRMshaw/code/scrm.pl '
        '--thitg 300 --imm --hexmcd --pac '
        '--genome {input.fa} '
        '--exon {input.exons} '
        '--gene {input.genes} '
        '--outdir {params.outdir} '
        '--traindirlst {input.traindirlst}'
        '&> {log}'

rule extract_chromosomes:
    input:
        fa = 'data/Amel_4.5_scaffolds.fa'
    output:
        fa = 'output/010_ref/chromosomes.fa'
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
    script:
        'src/extract_genes_and_exons.R'
