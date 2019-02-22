
# coding: utf-8

import sys #only needed to determine Python version number
#from Bio import SeqIO, Seq, SeqRecord
import numpy as np
import pandas as pd

modelID1 = sys.argv[1] #'1035191.3.fas'
modelID2 = sys.argv[2] #'a.aeolicus.pep'

Dir1 = '/data/results/' + modelID1 + '/profile/'
Dir2 = '/data/database/'

eval_cufoff = 0.001

#if len(sys.argv)>3:
#    eval_cufoff = sys.argv[4] #sys.argv[0] = xxx.py

# blastp_fname are from blast (C-BlastProfileGenomes.sh)
blastp_fname = Dir1 + modelID1 + '.TO.' + modelID2 + '.blastp'

out_fname = Dir1 + modelID1 + '.TO.' + modelID2 + '.parsed'

outfmt6_cols = 'qseqid sseqid evalue'.strip().split(' ')

df = pd.read_table(blastp_fname, header=None, names=outfmt6_cols)

fout = open(out_fname,'wb')
gp_df = df.groupby('qseqid')
keys = gp_df.groups.keys() # qseqid
keys.sort()
for query in keys: # each query
    q_gp_df = gp_df.get_group(query)
    
    # sort by e-val in descending order
    q_df = q_gp_df.sort_values(by='evalue', ascending=True, inplace=False)
    if q_df.iloc[0]['evalue']<=eval_cufoff:
        fout.write( query + '\t1\n')
    else:
        fout.write( query + '\t0\n')
        
fout.close()

