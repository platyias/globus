import os
import sys #only needed to determine Python version number
import numpy as np
import pandas as pd
import re

modelID = sys.argv[1] # .fas
#yeast = int(sys.argv[2]) # 0 false or 1 true
#swissprotv = int(sys.argv[3]) # 1 (old) or 2 (newer) 
#global_netv = int(sys.argv[2]) # 1 (old) or 2 (newer) 

Dir1 = '/data/results/' + modelID + '/'
Dir2 = '/data/database/'

#if swissprotv == 1:
#    gene_hit_fname = Dir1 + modelID + '_M_sw1.txt' # from Parse_blastp_outfmt6-All_models.ipynb
#    out_filename1 = Dir1 + modelID + '.gene_sw1'
#    out_filename2 = Dir1 + modelID + '.gene_ec_link_only_sw1'
#else:
gene_hit_fname = Dir1 + modelID + '_M.txt' # from Parse_blastp_outfmt6-All_models.ipynb
out_filename1 = Dir1 + modelID + '.gene'
out_filename2 = Dir1 + modelID + '.gene_ec_link_only'

#if yeast == 1:
#    yeast_g_ec_fname = Dir2 + 'ans_gene_ec_dict.npy' 
#else:
#    yeast_g_ec_fname = Dir2 + 'ans_gene_ec_dict_ibsu1103.npy' 

#if global_netv == 2:
global_network_fname = Dir2 + 'GRGENEECnet_2014.txt'
#else:
#    global_network_fname = Dir2 + 'GRGENEECnet.org.txt'


ec_index_dict = {} # global network
# key = EC number, val = index
with open(global_network_fname,'r') as f: 
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        if ec_index_dict.has_key(line_entry[1]):
            print 'error: same ec numbers: ' + line_entry[1]
        else:
            ec_index_dict[line_entry[1]] = line_entry[0]

# **f[0] = q_df['qseqid'].iloc[ind] (HMPREF0185_00018)
# **f[1] = ec (2.1.1.113)
# **f[3] = str('%s' %(q_df['pident'].iloc[ind])) (26.7)

# In[68]:

# Read BlastP results ( bacteria genome e.g. 1035191.3.fas to Swist-Prot enzyme sequence )
# for each gene, among all ec # hits, Keep only ones that are in the global network
gene_ecIndex_dict = {} # key = a gene in genome, val = list of ec index in the global network

# for each gene, among all ec # hits, list sequence identity of ones that are in the global network
gene_pident_dict = {}

gene_ec_by_pident_dict = {} # key = (gene+ec), value = pident

gene_in_order = [] # list of unique gene names. Gene whose ec # is in the global network
with open(gene_hit_fname,'r') as f:
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        if ec_index_dict.has_key(line_entry[1]): # ec_index_dict = global network: key = EC number, val = index
            
            gene_ec_by_pident_dict[line_entry[0] + ' ' + line_entry[1]] = line_entry[4]
            
            if gene_ecIndex_dict.has_key(line_entry[0]):
                gene_ecIndex_dict[line_entry[0]].append(ec_index_dict[line_entry[1]])
                gene_pident_dict[line_entry[0]].append(line_entry[4])
            else: 
                gene_in_order.append(line_entry[0])
                gene_ecIndex_dict[line_entry[0]] = [ec_index_dict[line_entry[1]]]
                gene_pident_dict[line_entry[0]] = [line_entry[4]]                


#### iLL672 of Saccharomyces cerevisiae S288c
#yeast_gene_ec_index_dict = {} # key = gene, val = ec index 
#yeast_gene_ec_dict = np.load(yeast_g_ec_fname).item()  # key = gene, val = a list of ecs
#for gene in yeast_gene_ec_dict.keys():
#    ecs = yeast_gene_ec_dict[gene]
#    ec_indices = []
#    for ec in ecs:
#        if ec_index_dict.has_key(ec):
#            ec_indices.append(ec_index_dict[ec])

#    yeast_gene_ec_index_dict[gene] = ec_indices

fout1 = open(out_filename1,'wb')
fout2 = open(out_filename2,'wb') # gene_ec_link only used for get_orthology.py
for i in np.arange(len(gene_in_order)): # gene_in_order=list of unique gene names. Gene whose ec # is in the global network
    ecIndices = gene_ecIndex_dict[gene_in_order[i]]
    fout1.write(str(i) + '\t' + gene_in_order[i] + '\t' + str(len(ecIndices)) + '\t')
    fout2.write(str(i) + '\t' + gene_in_order[i] + '\t' + str(len(ecIndices)) + '\t')
    [fout1.write(ecInd+'\t') for ecInd in ecIndices]
    fout1.write('///\t')
    [fout2.write(ecInd+'\t') for ecInd in ecIndices]
    fout2.write('///\n')
    pidents = gene_pident_dict[gene_in_order[i]]
    [fout1.write(pid+'\t') for pid in pidents]
    fout1.write('///\t')
    fout1.write('0\t///\n')
    
#    if re.search('=',gene_in_order[i]):
#        org_name = gene_in_order[i][0:-2]
#    else:
#        org_name = gene_in_order[i]
#    if yeast_gene_ec_index_dict.has_key(org_name):
#        yeast_ecIndices = yeast_gene_ec_index_dict[org_name]
#        fout1.write(str(len(yeast_ecIndices)) + '\t')
#        [fout1.write(ecInd+'\t') for ecInd in yeast_ecIndices]
#        fout1.write('///\n')
#    else:
     
    
fout1.close()
fout2.close()

