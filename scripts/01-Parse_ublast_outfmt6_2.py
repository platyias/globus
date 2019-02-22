import os
import sys #only needed to determine Python version number
import numpy as np
import pandas as pd
#from Bio import SeqIO, Seq, SeqRecord


modelID = sys.argv[1]
#swissprotv = int(sys.argv[2])  # 1 (old) or 2 (newer)

Dir1 = '/data/results/' + modelID + '/'
Dir2 = '/data/database/'

#if swissprotv == 1:
#    #swissprot.1
#    blastp_fname = Dir1 + modelID + '.VsSprot_5e-2_out6_seq.blastp'
#    MultipleECdb_fname = Dir2 + 'SPEmulti.merge'
#    out_filename = Dir1 + modelID + '_M_sw1.txt'
#    out_filename_hit = Dir1 + modelID + '_hit_M_sw1.txt'
#else:
#swissprot.2
blastp_fname = Dir1 + modelID + '.VsSprot2_5e-2_out6_seq.blastp'
MultipleECdb_fname = Dir2 + 'SPEmulti.merge2'
out_filename = Dir1 + modelID + '_M.txt'
out_filename_hit = Dir1 + modelID + '_hit_M.txt'


# Read MultipleECdb and build a 2D dictionary: protein name -> EC# -> s/e positions
MultipleECdb = {}
with open(MultipleECdb_fname,'r') as f:
        for line in f:
            line = line.rstrip()
            line_entry = line.split()
            # 0= protein name, 1=ec|ec|.... 2. one of the ecs 3. start position 4. end position
            startP = line_entry[3]
            endP = line_entry[4]
            if startP.isdigit() or (startP.startswith('-') and startP[1:].isdigit()):
                if endP.isdigit() or (endP.startswith('-') and endP[1:].isdigit()):
                    startP = int(startP)
                    endP = int(endP)
                    if MultipleECdb.has_key(line_entry[0]):
                        if MultipleECdb[line_entry[0]].has_key(line_entry[2]):
                            print 'error: same ec numbers within a protein:' + line_entry[0]
                        else:
                            MultipleECdb[line_entry[0]][line_entry[2]] = [startP, endP]
                    else:
                        MultipleECdb[line_entry[0]] = {line_entry[2]:[startP, endP]}
                else:
                    print line
            else:
                print line

# the outfmt 6 columns
default_outfmt6_cols = 'qseqid sseqid salltitles qlen slen pident mismatch gaps qstart qend sstart send evalue bitscore qseq sseq'.strip().split(' ')
#df.columns = default_outfmt6_cols

# df is a pandas.DataFrame and blastp_fname won't have a header
df = pd.read_table(blastp_fname, header=None, names=default_outfmt6_cols, sep='\s+')

# qseqid = Query Seq-id
# sseqid = Subject Seq-id
# salltitles = All Subject Title(s), separated by a '<>'
# qlen = Query sequence length
# slen = Subject sequence length
# pident = Percentage of identical matches
# length = Alignment length
# mismatch = Number of mismatches
# gaps = Number of gaps
# qstart = Start of alignment in query
# qend = End of alignment in query
# sstart = Start of alignment in subject
# send = End of alignment in subject
# evalue = Expect value
# bitscore = Bit score
# qseq = Aligned part of query sequence
# sseq = Aligned part of subject sequence

df['pident'] = df['pident'].apply(lambda x:round(x,1)) 
#df['salltitles']=df['sseqid'];


length_cr = (df['send'] - df['sstart'] + 1 >= 100) & (df['qend'] - df['qstart'] + 1 >= 100)
frac_cr = (df['send'] - df['sstart'] + 1>= 0.9*df['slen'] ) & (df['qend'] - df['qstart'] + 1>= 0.9*df['qlen'] )
df = df[length_cr | frac_cr]

df['salltitles'] = df['salltitles'].apply(lambda x:x.split()[0]) # 2.2.1.9|4.2.1.113|4.2.99.20

def find_shifted_posi_in_aligned(h_i, hit_seq, query_seq):
    cnt_r = 0 # residue
    posi = 0
    while(cnt_r<h_i):
        if hit_seq[posi] != '-':
            cnt_r += 1
        posi += 1
    q_posi = posi - query_seq[0:posi].count('-')
    return q_posi

init_data = {}
init_data['query'] = []
init_data['ec'] = []
init_data['start'] = []
init_data['end'] = []
init_data['pid'] = []
init_data['sseqid'] = []

gp_df = df.groupby('qseqid')
keys = gp_df.groups.keys()
keys.sort()
for query in keys:
    all_ec_for_query = []
    saved_hit_names = []
    
    q_gp_df = gp_df.get_group(query)
    #hit_gp_df = q_df.groupby('sseqid')
    #hit_keys = hit_gp_df.groups.keys() # hit sseqid grouped
    
    # sort by e-val in descending order
    q_df = q_gp_df.sort_values(by='evalue', ascending=True, inplace=False)
    
    # for each query, store the best e-value ec among all hit-hps
    for ind in np.arange(len(q_df)):
        
        # if hit and ec are not reported yet, read a line from the q_df
        #if ((q_df['sseqid'].iloc[ind] not in saved_hit_names)
        #    or (q_df['salltitles'].iloc[ind] not in all_ec_for_query)):
        if (q_df['salltitles'].iloc[ind] not in all_ec_for_query):
            if q_df['sseqid'].iloc[ind] in MultipleECdb.keys():
                ec_list = MultipleECdb[q_df['sseqid'].iloc[ind]] # 'RLMKL_PSYA2': {'2.1.1.173': [21, 386], '2.1.1.264': [455, 765]},
                for ec in ec_list:
                    if ec not in all_ec_for_query:

                        if (MultipleECdb[q_df['sseqid'].iloc[ind]][ec][0] <= q_df['sstart'].iloc[ind] and                            MultipleECdb[q_df['sseqid'].iloc[ind]][ec][1] >= q_df['sstart'].iloc[ind])                            or (MultipleECdb[q_df['sseqid'].iloc[ind]][ec][0] >= q_df['sstart'].iloc[ind] and                            MultipleECdb[q_df['sseqid'].iloc[ind]][ec][0] <= q_df['send'].iloc[ind]):

                                db_ec_len = MultipleECdb[q_df['sseqid'].iloc[ind]][ec][1] -                                        MultipleECdb[q_df['sseqid'].iloc[ind]][ec][0] + 1

                                if(q_df['sstart'].iloc[ind] < MultipleECdb[q_df['sseqid'].iloc[ind]][ec][0]):
                                    h_s = MultipleECdb[q_df['sseqid'].iloc[ind]][ec][0]
                                else:
                                    h_s = q_df['sstart'].iloc[ind]
                                    #print 'down'
                                if(q_df['send'].iloc[ind] < MultipleECdb[q_df['sseqid'].iloc[ind]][ec][1]):
                                    h_e = q_df['sstart'].iloc[ind]
                                else:
                                    h_e = MultipleECdb[q_df['sseqid'].iloc[ind]][ec][1]

                                # q_s and q_e are indices start from 0
                                # original position: q_[se] + q_df[[sq]start][ind]
                                q_s = find_shifted_posi_in_aligned(h_s - q_df['sstart'].iloc[ind], 
                                            q_df['sseq'].iloc[ind], q_df['qseq'].iloc[ind])
                                q_e = find_shifted_posi_in_aligned(h_e - q_df['sstart'].iloc[ind], 
                                            q_df['sseq'].iloc[ind], q_df['qseq'].iloc[ind])

                                if( (q_e-q_s+1>=100 and h_e-h_s>=100) or 
                                   (q_e-q_s+1>= 0.9*db_ec_len and h_e-h_s+1>= 0.9*db_ec_len) ):

                                    all_ec_for_query.append(ec)
                                    
                                    init_data['query'].append(q_df['qseqid'].iloc[ind])
                                    init_data['ec'].append(ec)
                                    init_data['start'].append(q_df['qstart'].iloc[ind])
                                    init_data['end'].append(q_df['qend'].iloc[ind])
                                    init_data['pid'].append(q_df['pident'].iloc[ind])
                                    init_data['sseqid'].append(q_df['sseqid'].iloc[ind])
            else: #not in the MultipleECdb
                ec_list = q_df['salltitles'].iloc[ind].split('|')
                for ec in ec_list:
                    if ec not in all_ec_for_query:
                        all_ec_for_query.append(ec)
                        
                        init_data['query'].append(q_df['qseqid'].iloc[ind])
                        init_data['ec'].append(ec)
                        init_data['start'].append(q_df['qstart'].iloc[ind])
                        init_data['end'].append(q_df['qend'].iloc[ind])
                        init_data['pid'].append(q_df['pident'].iloc[ind])
                        init_data['sseqid'].append(q_df['sseqid'].iloc[ind])
            
            if q_df['salltitles'].iloc[ind] not in all_ec_for_query:
                all_ec_for_query.append(q_df['salltitles'].iloc[ind])

df_dom = pd.DataFrame(init_data)
gp_q_df = df_dom.groupby('query')
q_names = gp_q_df.groups.keys() 

# Read BlastP results ( bacteria genome e.g. 1035191.3.fas to Swist-Prot enzyme sequence )
# for each gene, among all ec # hits, Keep only ones that are in the global network
fout = open(out_filename,'wb')
fout_hit = open(out_filename_hit,'wb')

for query in q_names:
    q_df = gp_q_df.get_group(query)
    # sort by start position
    q_df = q_df.sort_values(by='start', ascending=True, inplace=False)
    num_hits = len(q_df)
    lst_gp = [1]
    gp = 1
    for i in np.arange(num_hits-1):
        if q_df.iloc[i]['end'] < q_df.iloc[i+1]['start']:
            gp = gp + 1
        lst_gp.append(gp)
    
    q_df['group'] = lst_gp
    if gp > 1: # multiple domains in the query
        for ind in np.arange(num_hits):
            fout.write(
                q_df['query'].iloc[ind] + '=' + str(q_df['group'].iloc[ind]) + '\t' +
                q_df['ec'].iloc[ind] + '\t' +
                str('%s' %(q_df['start'].iloc[ind])) + '\t' +
                str('%s' %(q_df['end'].iloc[ind])) + '\t' +
                str('%s' %(q_df['pid'].iloc[ind])) + '\n') 
                #str('%g' %(q_df['evalue'].iloc[ind])) + '\t' +
                #q_df['sseqid'].iloc[ind] + '\t' +
                #q_df['salltitles'].iloc[ind] + '\n' )
            fout_hit.write(
                 q_df['query'].iloc[ind] + '=' + str(q_df['group'].iloc[ind]) + '\t' +
                 q_df['ec'].iloc[ind] + '\t' +
                 str('%s' %(q_df['start'].iloc[ind])) + '\t' +
                 str('%s' %(q_df['end'].iloc[ind])) + '\t' +
                 str('%s' %(q_df['pid'].iloc[ind])) + '\t' +
                 str('%s' %(q_df['sseqid'].iloc[ind])) + '\n')

    else:
        for ind in np.arange(num_hits):
            fout.write(
                q_df['query'].iloc[ind] + '\t' +
                q_df['ec'].iloc[ind] + '\t' +
                str('%s' %(q_df['start'].iloc[ind])) + '\t' +
                str('%s' %(q_df['end'].iloc[ind])) + '\t' +
                str('%s' %(q_df['pid'].iloc[ind])) + '\n') 
                #str('%g' %(q_df['evalue'].iloc[ind])) + '\t' +
                #q_df['sseqid'].iloc[ind] + '\t' +
                #q_df['salltitles'].iloc[ind] + '\n' )
            fout_hit.write(
                 q_df['query'].iloc[ind] + '\t' +
                 q_df['ec'].iloc[ind] + '\t' +
                 str('%s' %(q_df['start'].iloc[ind])) + '\t' +
                 str('%s' %(q_df['end'].iloc[ind])) + '\t' +
                 str('%s' %(q_df['pid'].iloc[ind])) + '\t' +
                 str('%s' %(q_df['sseqid'].iloc[ind])) + '\n')


fout.close()

