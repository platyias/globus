import sys
import os
#from scipy import odeint
import numpy as np

modelID = sys.argv[1] # query genome name '1035191.3.fas'
context_type = sys.argv[2] # 'pc' # or 'gc' 
#swissprotv = int(sys.argv[3]) # 1 (old) or 2 (newer) 

Dir1 = '/data/results/' + modelID + '/'

#if swissprotv == 1: 
#    phylo_corr_fname = Dir1 + modelID + '.' + context_type + '_sw1'
#    out_fname = Dir1 + modelID + '.' + context_type + '.Z_sw1'
#else:
phylo_corr_fname = Dir1 + modelID + '.' + context_type
out_fname = Dir1 + modelID + '.' + context_type + '.Z'

f = open(phylo_corr_fname,'r') # read the first line 
line = f.readline()
line = line.rstrip()
if context_type == 'gc':
    line_entry = line.split()
    num_genes = int(line_entry[-3]) + 1
else: # 'pc'
    gene_list = line.split()
    num_genes = len(gene_list)

Gene_Gene_context_mat = np.zeros([num_genes,num_genes]) 

if context_type == 'gc':
    with open(phylo_corr_fname,'r') as f:  
        line = f.readline()
        line = line.rstrip()
        gene_list = line.split()
        curr_gene = ''
        row=-1
        col=0
        for line in f:
            line = line.rstrip()
            line_entry = line.split()
            if curr_gene != line_entry[0]:
                row += 1
                col=0
                curr_gene = line_entry[0]
            Gene_Gene_context_mat[row][row+1+col] = float(line_entry[-1])
            Gene_Gene_context_mat[row+1+col][row] = float(line_entry[-1])
            col += 1
else: # 'pc'
    with open(phylo_corr_fname,'r') as f:  
        line = f.readline()
        line = line.rstrip()
        gene_list = line.split()
        row = 0
        for line in f:
            line = line.rstrip()
            line_entry = line.split()
            for i in np.arange(len(line_entry)):
                if line_entry[i] == 'nan':
                    Gene_Gene_context_mat[row][row+1+i] = np.nan #-20000
                    Gene_Gene_context_mat[row+1+i][row] = np.nan #-20000
                else:
                    Gene_Gene_context_mat[row][row+1+i] = float(line_entry[i])
                    Gene_Gene_context_mat[row+1+i][row] = float(line_entry[i])
            row += 1


# Gene_Gene_context_mat = symmetric num_genes x num_genes, Must ignore (i,i) entry
gene_context_mean_vec = np.zeros(num_genes)
gene_context_std_vec = np.zeros(num_genes)
for i in np.arange(0,num_genes): # a gene in the network having pc
    x_gene_corr = Gene_Gene_context_mat[i]
    x_gene_corr = x_gene_corr[ ~np.isnan(x_gene_corr) ]
    gene_context_mean_vec[i] = np.mean(x_gene_corr)
    gene_context_std_vec[i] = np.std(x_gene_corr)


# Gene_Gene_context_mat = symmetric num_genes x num_genes, **** -3 at (i,i) entry ****
Gene_Gene_Zscore_mat = np.zeros([num_genes,num_genes]) 

for i in np.arange(0,num_genes): # a gene in the network having pc
    Gene_Gene_Zscore_mat[i][i] = -3
        
    for j in np.arange(i+1,num_genes):
        if ~np.isnan(Gene_Gene_context_mat[i][j]) and ~np.isnan(gene_context_mean_vec[i]) and ~np.isnan(gene_context_std_vec[i]) and ~np.isnan(gene_context_mean_vec[j])  and ~np.isnan(gene_context_std_vec[j]):
            z1 = (Gene_Gene_context_mat[i][j] - gene_context_mean_vec[i])/gene_context_std_vec[i]
            z2 = (Gene_Gene_context_mat[i][j] - gene_context_mean_vec[j])/gene_context_std_vec[j]
            z_score = np.sqrt(z1*z1 + z2*z2) ##/2.0 GP Error in calculating Z-score
        else:
            z_score = -2
        Gene_Gene_Zscore_mat[i][j] = z_score
        Gene_Gene_Zscore_mat[j][i] = z_score


fout = open(out_fname,'wb')
for i in np.arange(0,num_genes): # a gene in the network having pc
    [fout.write(str(Gene_Gene_Zscore_mat[i][j]) + ' ') for j in np.arange(0,num_genes)]
    fout.write('\n')
fout.close()

