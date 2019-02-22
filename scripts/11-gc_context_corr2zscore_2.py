import sys
import os
import numpy as np


modelID = sys.argv[1] # query genome name
#swissprotv = int(sys.argv[3]) # 1 (old) or 2 (newer) 
#swissprotv = int(sys.argv[2]) # 1 (old) or 2 (newer) 

Dir1 = '/data/results/' + modelID + '/'
#if swissprotv == 1: 
#    gc_corr_fname = Dir1 + modelID + '.gc_sw1'
#    SwProt_gene_fname = Dir1 + modelID + '.gene_ec_link_only_sw1'
#    out_fname = Dir1 + modelID + '.' + 'gc.Z_sw1'
#else:
gc_corr_fname = Dir1 + modelID + '.gc'
SwProt_gene_fname = Dir1 + modelID + '.gene_ec_link_only'
out_fname = Dir1 + modelID + '.' + 'gc.Z'

gene_ecIndex_dict = {} # key = a gene index in genome, val = list of ec index in the global network
gene_in_order = {} # list of unique gene names. Gene whose ec # is in the global network
# gene_in_order['xxx=1'] = xxx
#SwProt_gene_fname = Dir1 + modelID_1 + '.gene_ec_link_only'
with open(SwProt_gene_fname,'r') as f: # 0=index, 1=gene, 2=# of ec annotated to the gene, 3~=ec numbers
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        gene_in_order[line_entry[1]] = line_entry[1].split('=')[0]
        for i in np.arange(int(line_entry[2])):
            if gene_ecIndex_dict.has_key(line_entry[0]):
                gene_ecIndex_dict[line_entry[0]].append(line_entry[3+i])
            else:
                gene_ecIndex_dict[line_entry[0]] = [line_entry[3+i]]


'''
Total number of genomes: 110 + 1
YAL001C YAL002W 2.73750514
YAL001C YAL003W 2.47955202
YAL001C YAL005C 2.69146405
YAL001C YAL007C 2.26612897
'''
num_genes = 1
curr_gene = ''
gene_name_vec = {} # name -> index from 0, ... num_genes-1
with open(gc_corr_fname,'r') as f:  
    line = f.readline()
    for line in f:
        line = line.rstrip() # YAL001C YAL002W 2.73750514
        line_entry = line.split()
        if curr_gene == '':
            curr_gene = line_entry[0]
            #gene_name_vec.append(curr_gene)
            gene_name_vec[curr_gene] = num_genes-1
        else:
            if curr_gene != line_entry[0]:
                break
        num_genes = num_genes + 1
        #gene_name_vec.append(line_entry[1])
        gene_name_vec[line_entry[1]] = num_genes-1


Gene_Gene_context_mat = np.zeros([num_genes,num_genes])

with open(gc_corr_fname,'r') as f:  
    line = f.readline()
    
    curr_gene = ''
    row=-1

    for line in f:
        line = line.rstrip() # YAL001C YAL002W 2.73750514
        line_entry = line.split()
        if curr_gene != line_entry[0]:
            if curr_gene != '':
                # row begins with 0
                for j in np.arange(row+1,num_genes):
                    Gene_Gene_context_mat[row][j] = vec[j-row-1]
                    Gene_Gene_context_mat[j][row] = vec[j-row-1]
            vec = []
            row += 1
            curr_gene = line_entry[0]
        vec.append(float(line_entry[-1]))    
        #vec.append(line_entry[-1])
    for j in np.arange(row+1,num_genes):
        Gene_Gene_context_mat[row][j] = vec[j-row-1]
        Gene_Gene_context_mat[j][row] = vec[j-row-1]


sub_index_vec = []
geneName2index = {} # index of sub_index_vec and 

for gene in gene_in_order.keys(): # gene name in .gene file
    name = gene.split('=')[0]
    if name in gene_name_vec.keys():
        if gene_name_vec[name] not in sub_index_vec:
            sub_index_vec.append(gene_name_vec[name])
            geneName2index[name] = len(sub_index_vec)-1
    else:
        print 'err'
        
sub_index_vec_dup = [] # index of sub_index_vec and 
for gene in gene_in_order.keys(): # gene name in .gene file
    name = gene.split('=')[0]
    sub_index_vec_dup.append(geneName2index[name])


sub_Gene_Gene_context_mat = Gene_Gene_context_mat[sub_index_vec,:][:,sub_index_vec]


sub_num_genes = len(sub_index_vec)
gene_context_mean_vec = np.zeros(sub_num_genes)
gene_context_std_vec = np.zeros(sub_num_genes)
for i in np.arange(0,sub_num_genes): # a gene in the network having pc
    x_gene_corr = sub_Gene_Gene_context_mat[i]
    x_gene_corr = x_gene_corr[ ~np.isnan(x_gene_corr) ]
    gene_context_mean_vec[i] = np.mean(x_gene_corr)
    gene_context_std_vec[i] = np.std(x_gene_corr)


Gene_Gene_Zscore_mat = np.zeros([sub_num_genes,sub_num_genes]) 

for i in np.arange(0,sub_num_genes): # a gene in the network having pc
    Gene_Gene_Zscore_mat[i][i] = -3
        
    for j in np.arange(i+1,sub_num_genes):
        if ~np.isnan(sub_Gene_Gene_context_mat[i][j]) and ~np.isnan(gene_context_mean_vec[i]) and ~np.isnan(gene_context_std_vec[i]) and ~np.isnan(gene_context_mean_vec[j]) and ~np.isnan(gene_context_std_vec[j]) and gene_context_std_vec[i] <> 0 and gene_context_std_vec[j] <> 0:
            # warning is from z1 and z2 calculation
            z1 = (sub_Gene_Gene_context_mat[i][j] - gene_context_mean_vec[i])/gene_context_std_vec[i]
            z2 = (sub_Gene_Gene_context_mat[i][j] - gene_context_mean_vec[j])/gene_context_std_vec[j]
            z_score = np.sqrt(z1*z1 + z2*z2) ##/2.0 GP Error in calculating Z-score
            if np.isnan(z_score):
                z_score = -2
        else:
            z_score = -2
        Gene_Gene_Zscore_mat[i][j] = z_score
        Gene_Gene_Zscore_mat[j][i] = z_score


Gene_Gene_Zscore_mat_dup = Gene_Gene_Zscore_mat[sub_index_vec_dup,:][:,sub_index_vec_dup]


fout = open(out_fname,'wb')
for i in np.arange(0,len(sub_index_vec_dup)): # a gene in the network having pc
    [fout.write(str(Gene_Gene_Zscore_mat_dup[i][j]) + ' ') for j in np.arange(0,len(sub_index_vec_dup))]
    fout.write('\n')
fout.close()


