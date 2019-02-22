
# coding: utf-8

import sys #only needed to determine Python version number
import numpy as np
import pandas as pd

modelID_1 = sys.argv[1] # .fas
speID = sys.argv[2] # spe
# modelID_1 = '1035191.3.fas'
# speID = 'sce' # no organism called atc

Dir1 = '/data/results/' + modelID_1 + '/blast/'

Dir2 = '/data/database/'

blastp_fname_1 = Dir1 + speID + '.vs.' + modelID_1 + '.blastp'
blastp_fname_2 = Dir1 + modelID_1 + '.vs.' + speID + '.blastp'

out_filename = Dir1 + modelID_1 + '.vs.' + speID + '.bbhs'


if False:
    '''
    04-Get_BBHs.pl  
    /results/$b/blast/${b}.vs.${a}.blastp # b=genome, a=organism
    /results/$b/blast/${a}.vs.${b}.blastp 
    > /results/$b/blast/${b}.vs.${a}.bbhs

    05-Get_Orthology.pl 
    /results/$1/$1.blast.bbh.out 
    /OLD_CONTEXT/KEGG2Uniprot.txt 
    /OLD_CONTEXT/Swissprot2Uniprot.txt 
    /OLD_CONTEXT/Swissprot2EC.txt 
    /database/GRGENEECnet.txt 
    /results/$1/$1.gene
    $1 

    > /results/$1/$1.ort
    '''

#default_outfmt6_cols =  'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'.strip().split(' ')
#default_outfmt6_cols =  'qseqid sseqid qlen slen pident length mismatch gaps qstart qend sstart send evalue bitscore'.strip().split(' ')
default_outfmt6_cols =  'qseqid sseqid qlen slen pident mismatch gaps qstart qend sstart send evalue bitscore'.strip().split(' ')


# df is a pandas.DataFrame and blastp_fname won't have a header
df_gen2org = pd.read_table(blastp_fname_2, header=None, names=default_outfmt6_cols)
df_org2gen = pd.read_table(blastp_fname_1, header=None, names=default_outfmt6_cols)


# Kepp the best hits with the highest bitscores only
drop_index_df_gen2org = []
gp_df = df_gen2org.groupby('qseqid')
for query in gp_df.groups.keys():
    q_gp_df = gp_df.get_group(query) # sub-df of a query
    q_sorted_df = q_gp_df.sort_values(by='bitscore', ascending=False, inplace=False)
    [drop_index_df_gen2org.append(index) for index in q_sorted_df.index[1:]]

df_gen2org = df_gen2org.drop(drop_index_df_gen2org)


#plt.hist(df_gen2org['send'] - df_gen2org['sstart'])
#plt.hist(df_gen2org['qend'] - df_gen2org['qstart'])


### TODO ###
# Strict threshold (combination of the length of the alignment, e-value, and bitscore )
length_cr = (df_gen2org['send'] - df_gen2org['sstart'] + 1 >= 40) & (df_gen2org['qend'] - df_gen2org['qstart'] + 1 >= 40)
frac_cr = (df_gen2org['send'] - df_gen2org['sstart'] + 1>= 0.7*df_gen2org['slen'] ) & (df_gen2org['qend'] - df_gen2org['qstart'] + 1>= 0.7*df_gen2org['qlen'] )
bitsc_cr = df_gen2org['bitscore'] > 26.4 # from bitscore distribution.
eval_cr = df_gen2org['evalue'] < 2 # from evalue distribution.


df_gen2org_filtered = df_gen2org[(length_cr | frac_cr) & bitsc_cr & eval_cr]

df_gen2org_filtered = df_gen2org[df_gen2org['bitscore'] > 26.4]


query_hit_dict = {} # key = query, val = hit name
query_BS_dict = {} # key = query, val = bit score >26.4
# This is for df_gen2org only with the best bitscores ? or with stricter criteria?
for index in df_gen2org_filtered.index:
    query = df_gen2org_filtered.loc[index]['qseqid']
    query_hit_dict[query] = df_gen2org_filtered.loc[index]['sseqid']
    query_BS_dict[query] = df_gen2org_filtered.loc[index]['bitscore']


#Finding BBH 
fout = open(out_filename,'wb')
gp_df = df_org2gen.groupby('qseqid') # org to genome
for query in gp_df.groups.keys():
    q_gp_df = gp_df.get_group(query) # sub-df of a query
    q_sorted_df = q_gp_df.sort_values(by='bitscore', ascending=False, inplace=False)
    best_gen_name = q_sorted_df['sseqid'].iloc[0]
    # if Best hits match in both directions
    if query_hit_dict.has_key(best_gen_name) and query_hit_dict[best_gen_name] == query: 
        fout.write(best_gen_name + '\t' + query + '\t' + str(query_BS_dict[best_gen_name]) + '\n') 
        # gene name, species hit name, bitscore of genome to species alginment
        

