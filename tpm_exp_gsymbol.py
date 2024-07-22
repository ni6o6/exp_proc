#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 13:12:23 2023

@author: Naoko Iida

when multiple ids for a gene, choose larger count.
"""

import gzip
import sys

count_file = sys.argv[1]
summary_file = sys.argv[2]
gtf_file = sys.argv[3]
output_file1 = sys.argv[4]
output_file2 = sys.argv[5]

gene_id2gene_name = {}
with open(gtf_file, 'r') as hin:
    for line in hin:
        if line.startswith('#'): continue
        F = line.rstrip('\n').split('\t')
        
        if F[2] == "gene":
            infos = F[8].split(';')
            gene_id, gene_name = '', ''
            for elm in infos:
                
                if elm.startswith("gene_id"):
                    gene_id = elm.replace("gene_id ", '').strip('"')
                if elm.startswith(" gene_name"):
                    gene_name = elm.replace(" gene_name ", '').strip('"')
            if gene_name.startswith("ENSG"): continue
            else:
                gene_id2gene_name[gene_id] = gene_name

with open(summary_file, 'r') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0] == "Assigned": 
            total_read_num = int(F[1])

total_read_ratio = float(total_read_num) / 1000000

gene2fpkm = {}
gene2count = {}
with gzip.open(count_file, 'rt') as hin:
    for line in hin:
        F = line.rstrip('\n').split('\t')
        if F[0].startswith('#'): continue
        if F[0] == "Geneid": continue
        gene_id = F[0]
        if gene_id not in gene_id2gene_name: continue
        gene = gene_id2gene_name[F[0]]
        gene_length = float(F[5])
        count = float(F[6])

        fpkm = count / ((gene_length / 1000) * total_read_ratio)
        if gene not in gene2fpkm: 
            gene2fpkm[gene] = fpkm
            gene2count[gene] = count
        else: 
            if fpkm > gene2fpkm[gene]: 
                gene2fpkm[gene] = fpkm
                gene2count[gene] = count

with open(output_file1, "w") as hout, open(output_file2, "w") as kout:
    for gene in gene2fpkm:
        hout.write(gene + '\t' + str(round(gene2fpkm[gene], 4)) + '\n')
    for gene in gene2count:
        kout.write(gene + '\t' + str(int(gene2count[gene])) + '\n')


