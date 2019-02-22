
# coding: utf-8

import sys
import os
import numpy as np
import re
import collections

modelID = sys.argv[1] # query genome name '1035191.3.fas'
#swissprotv = int(sys.argv[2]) # 1 (old) or 2 (newer) 

num_genomes = 71 # sys.argv[2]
#if len(sys.argv)>2:
#    num_genomes = sys.argv[2] #sys.argv[0] = xxx.py

Dir1 = '/data/results/' + modelID + '/'
profile_name = Dir1 + 'profile/' + modelID + '.profile.txt'

#if swissprotv == 1:
#    gene_ec_link_only_fname = Dir1 + modelID + '.gene_ec_link_only_sw1'
#    out_fname = Dir1 + modelID + '.pc_sw1'
#else:
gene_ec_link_only_fname = Dir1 + modelID + '.gene_ec_link_only'
out_fname = Dir1 + modelID + '.pc' 

gene_in_order = [] # list of unique gene names. Gene whose ec # is in the global network

with open(gene_ec_link_only_fname,'r') as f: # 0=index, 1=gene, 2=# of ec annotated to the gene, 3~=ec numbers
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        gene_in_order.append(line_entry[1])
         
network_gene_list = gene_in_order # include =1, =2..
#multi_domain_genes = [item for item, count in collections.Counter(gene_names).items() if count > 1]
num_network_genes = len(network_gene_list)


# read profile.txt genes x 71 length vec
phylo_corr_table = {} # key = gene
with open(profile_name,'r') as f: 
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        if num_genomes != (len(line_entry)-1):
            print 'num_genomes is ' + str(len(line_entry)-1) + '. not ' + str(num_genomes) + '\n'
            
        phylo_corr_table[line_entry[0]] = [-1]*num_genomes
        for i in np.arange(num_genomes):
            phylo_corr_table[line_entry[0]][i] = int(line_entry[i+1])
            
gene_name_w_pc = [] # =1, =2
gene_w_pc = []

#GP TO FIX FOR ALL GENES IN NETWORK
gene_names_no_multi = []

for i in np.arange(0,num_network_genes): # a gene in the network
    if re.search('=', gene_in_order[i]):
        multi_d = re.search('=', gene_in_order[i])
        gene = gene_in_order[i][0:multi_d.start()]
    else:
        gene = gene_in_order[i]
    
    gene_names_no_multi.append(gene) # GP
    
    if phylo_corr_table.has_key(gene):
        gene_name_w_pc.append(gene_in_order[i])
        gene_w_pc.append(gene)
        

def pearson_corr_BN(BNlist1, BNlist2):
    # BNvec1 and BNvec2 are binary integer vectors of the same length
    # if (len(BNvec1) != len(BNvec2)):
    #    print 'pearson_corr vector size'
    #    break
    #print 'diff' if np.sum(BNvec1 - BNvec2)!=0 else ''
    BNvec1 = np.array(BNlist1)
    BNvec2 = np.array(BNlist2)
    
    N = float(len(BNvec1))
    in_1 = float(np.sum(BNvec1))
    #print in_1
    in_2 = float(np.sum(BNvec2))
    #print in_2
    intersect_1_2 = float(np.sum(BNvec1 * BNvec2))
    #print intersect_1_2
    denominator = np.sqrt( (N*in_1 - in_1*in_1)*(N*in_2 - in_2*in_2) )
    #print denominator
    if denominator > 0:
        corr = (N*intersect_1_2 - in_1*in_2) / denominator
        #print corr
    else:
        corr = np.nan
    return corr

fout = open(out_fname,'wb')
#[fout.write(gene + '\t') for gene in gene_name_w_pc]# genes in the network having pc
[fout.write(gene + '\t') for gene in gene_in_order]# GP ALL GENES
fout.write('\n')
num_genes = len(gene_name_w_pc)

#for i in np.arange(0,num_genes): # a gene in the network having pc
#    x_gene_profile = phylo_corr_table[gene_w_pc[i]]
#    for j in np.arange(i+1,num_genes):
#        y_gene_profile = phylo_corr_table[gene_w_pc[j]]
#        
#        corr = pearson_corr_BN(x_gene_profile, y_gene_profile)
#        if np.isnan(corr) == False:
#            fout.write(str(corr) + '\t')
#        else:
#            fout.write('nan\t')
#    fout.write('\n')
#fout.close()

#GP MUST REPLACE SO IT PRINTS CORR FOR EVERY GENE

for i in np.arange(0,num_network_genes):
    if phylo_corr_table.has_key(gene_names_no_multi[i]):
      x_gene_profile = phylo_corr_table[gene_names_no_multi[i]]
      for j in np.arange(i+1,num_network_genes):
        if phylo_corr_table.has_key(gene_names_no_multi[j]):
          y_gene_profile = phylo_corr_table[gene_names_no_multi[j]]
          corr = pearson_corr_BN(x_gene_profile, y_gene_profile)
          if np.isnan(corr) == False:
            fout.write(str(corr) + '\t')
          else:
            fout.write('nan\t')
        else:
          fout.write('nan\t')
    else:
      for j in np.arange(i+1,num_network_genes):
        fout.write('nan\t')
    fout.write('\n')
fout.close()