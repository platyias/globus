
# coding: utf-8

import sys
import numpy as np
import pandas as pd
import re

modelID_1 = sys.argv[1] #'1035191.3.fas'
#swissprotv = int(sys.argv[2]) # 1 (old) or 2 (newer) 
#global_netv = int(sys.argv[3]) # 1 (old) or 2 (newer) 

Dir1 = '/data/results/' + modelID_1 + '/'
Dir2 = '/data/database/'

bbh_out_fname = Dir1 + modelID_1 + '.blast.bbh.out'

#if global_netv == 2:
global_network_fname = Dir2 + 'GRGENEECnet_2014.txt'
#else:
#    global_network_fname = Dir2 + 'GRGENEECnet.org.txt'


#if swissprotv == 1:
#    SwProt_gene_fname = Dir1 + modelID_1 + '.gene_ec_link_only_sw1'
#    ort_out_fname = Dir1 + modelID_1 + '.ort_sw1'
#else:
SwProt_gene_fname = Dir1 + modelID_1 + '.gene_ec_link_only'
ort_out_fname = Dir1 + modelID_1 + '.ort'

if False:
    '''
    05-Get_Orthology.pl 
    /results/$1/$1.blast.bbh.out : write(best_gen_name + '\t' + query + '\t' +  str(query_BS_dict[best_gen_name]) + '\n')
    /OLD_CONTEXT/KEGG2Uniprot.txt 
    /OLD_CONTEXT/Swissprot2Uniprot.txt 
    /OLD_CONTEXT/Swissprot2EC.txt 
    /database/GRGENEECnet.txt 
    /results/$1/$1.gene
    $1 

    > /results/$1/$1.ort
    '''

GN_ec_index_dict = {} # global network
# key = EC number, val = EC number index
with open(global_network_fname,'r') as f: 
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        if GN_ec_index_dict.has_key(line_entry[1]):
            print 'error: same ec numbers: ' + line_entry[1]
        else:
            GN_ec_index_dict[line_entry[1]] = line_entry[0]

bbh_ortholog_dict = {}
# key = gene in genome x species hit name as a tuple, 
# val = the number of such gene-species hit appear in bbh_out_fname file. 
with open(bbh_out_fname,'r') as f: 
    for line in f:
        line = line.rstrip()
        line_entry = line.split() # 0=gene, 1=species hit, 2=bitscore (gene -> species alignment)
        if bbh_ortholog_dict.has_key( (line_entry[0], line_entry[1]) ):
            bbh_ortholog_dict[ (line_entry[0], line_entry[1]) ] += 1
        else:
            bbh_ortholog_dict[ (line_entry[0], line_entry[1]) ] = 0


gene_ecIndex_dict = {} # key = a gene index in genome, val = list of ec index in the global network
gene_in_order = [] # list of unique gene names. Gene whose ec # is in the global network
#SwProt_gene_fname = Dir1 + modelID_1 + '.gene_ec_link_only'
with open(SwProt_gene_fname,'r') as f: # 0=index, 1=gene, 2=# of ec annotated to the gene, 3~=ec numbers
    for line in f:
        line = line.rstrip()
        line_entry = line.split()
        gene_in_order.append(line_entry[1])
        for i in np.arange(int(line_entry[2])):
            if gene_ecIndex_dict.has_key(line_entry[0]):
                gene_ecIndex_dict[line_entry[0]].append(line_entry[3+i])
            else:
                gene_ecIndex_dict[line_entry[0]] = [line_entry[3+i]]


# key = UniProt (P21215), val = SwissProt (12AH_CLOS4)
UniProt_SwProt_dict = np.load(Dir2+'UniProt_SwProt_dict.npy').item()

# key = swissprot , val = Kegg gene 
SwProt_Kegg_dict = np.load(Dir2+'SwProt_kegg_dict.npy').item()

EC_keggGene_dict = {} 
# key = A tuple (EC number index in the global network x kegg gene), 
#val = count the number of such tuple, which doesn't mean anything. 
#This dictionary is to keep all the unique pairs of a EC number and a gene

SwProt_EC_dict = np.load(Dir2+'SwProt_EC_dict.npy').item()
for swprot in SwProt_EC_dict.keys():
    if SwProt_Kegg_dict.has_key(swprot):
        keggGene = SwProt_Kegg_dict[swprot]
        if len(SwProt_EC_dict[swprot]) == 1:
            ec = SwProt_EC_dict[swprot][0]
            if GN_ec_index_dict.has_key(ec):
                ec_index = GN_ec_index_dict[ec]
                if EC_keggGene_dict.has_key( (ec_index, keggGene) ):
                    EC_keggGene_dict[(ec_index, keggGene)] += 1
                else:
                    EC_keggGene_dict[(ec_index, keggGene)] = 0


fout = open(ort_out_fname,'wb')
#
# TODO: **** Make sure the order of EC numbers matches with the order in .gene file.
# i.e., the order of gene-EC homology score matches with gene-EC orthology  
#

#bbh_ortholog_dict = {}
# key = gene in genome x species hit name as a tuple, 
# val = the number of such gene-species hit appear in bbh_out_fname file. 

for i in np.arange(len(gene_in_order)): # gene_in_order=list of unique gene names. 
    # Gene whose ec # is in the global network
    ecIndices = gene_ecIndex_dict[str(i)] # the same ordered list of ECs for the gene
    fout.write(str(i) + '\t' + gene_in_order[i] + '\t' + str(len(ecIndices)) + '\t')
    [fout.write(ecInd+'\t') for ecInd in ecIndices]
    fout.write('///\t')
    
    
    for ec in ecIndices:
        kegg_genes_for_ec = [item[1] for item in EC_keggGene_dict.keys() if item[0] == ec] 
        # When the first key matches to the ec, the corresponding genes (item[1]).
        # genes for the ec from the EC_keggGene table.
        
        # If gene_in_order[i] gene has orthologs (species hits) to any of those genes (genes_for_ec)
        orth_exist = False
        for gene in kegg_genes_for_ec:
            if re.search('=', gene_in_order[i]):
                multi_d = re.search('=', gene_in_order[i])
                gene_d = gene_in_order[i][0:multi_d.start()]
            else:
                gene_d = gene_in_order[i]

            if bbh_ortholog_dict.has_key((gene_d, gene)):
                orth_exist = True
                break
                
        if orth_exist == True:
            fout.write('1\t')
        else:
             fout.write('0\t')
    fout.write('///\n')
    
    # yeast answers are repeated from .gene -> Do not write here.
    if False:
        if yeast_gene_ec_dict.has_key(gene_in_order[i]):
            yeast_ecIndices = yeast_gene_ec_dict[gene_in_order[i]]
            fout.write(len(yeast_ecIndices) + '\t')
            [fout.write(ecInd+'\t') for ecInd in yeast_ecIndices]
            fout.write('///\n')
        else:
            fout.write('0\t///\n')
        
fout.close()        

