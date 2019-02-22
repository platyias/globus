
# coding: utf-8

import sys
import os
import numpy as np

num_genomes = 71

modelID = sys.argv[1] #'1035191.3.fas' 

Dir1 = '/data/results/' + modelID + '/profile/'

out_fname = Dir1 + modelID + '.profile.txt'

parsed_file_name_list = []
for file in os.listdir(Dir1):
    if file.endswith(".parsed"):
        parsed_file_name_list.append(os.path.join(Dir1, file))

parsed_file_name_list.sort()

if len(parsed_file_name_list) != num_genomes:
    print 'error the number of genomes is ' +  str(len(parsed_file_name_list)) + '. Has to be ' + str(num_genomes)


phylo_corr_table = {} 
# key = query genome gene, val = a list(length = 71) of 0 & 1s for existance of orthologs in 71 genomes.

index_genome = 0
for parsedf in parsed_file_name_list:
    with open(parsedf,'r') as f: 
        for line in f:
            line = line.rstrip()
            line_entry = line.split()
            if phylo_corr_table.has_key(line_entry[0]):
                phylo_corr_table[line_entry[0]][index_genome] = line_entry[1]
                #phylo_corr_table[line_entry[0]].append(line_entry[1])
            else:
                #phylo_corr_table[line_entry[0]] = ['-1']*num_genomes
                phylo_corr_table[line_entry[0]] = ['0']*num_genomes
                phylo_corr_table[line_entry[0]][index_genome] = line_entry[1]
                #phylo_corr_table[line_entry[0]] = [line_entry[1]]
    index_genome += 1

fout = open(out_fname,'wb')
for key in phylo_corr_table.keys(): # query gene
    fout.write(key + ' ')
    [fout.write(cor + ' ') for cor in phylo_corr_table[key]]
    fout.write('\n')

