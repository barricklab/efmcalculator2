import itertools
from Bio import SeqIO
from Bio import Align
from Bio.Seq import Seq
import pandas as pd
import re
import time

_start_time = time.time()

filename = "l-6-10_plasmid_bba.fasta"
#filename = "test.fasta"
sequences = SeqIO.parse(filename, 'fasta')
data = {'Sequence': [], 'Size': [], 'Distance': [], 'Start-Pos':[], 'End-Pos':[]}
df= pd.DataFrame(data)
unique_seq = {'Sequence': []}
usdf = pd.DataFrame(unique_seq)
for record in sequences:
    seq_len = len(record.seq)
    seq = record.seq
    for i,letter in enumerate(seq):
        #print("\n next nucleotide:",seq[i:i+1])
        
        for j in range(17,6,-1):
        
            if len(seq[i:i+j]) > 6:
                sub_seq = seq[i:i+j]
                rem_seq = seq[i+1:]
                #print(sub_seq,':',rem_seq)
                found_overlap = False
                
                #rem_sub_seq = seq[k:k+j]
                count = seq.count_overlap(sub_seq)
                if count > 1:
                   distance  =  1
                   start_pos  = 0
                   rem_start  = 0
                   end_pos = 0
                   prv_end_pos = 0
                   for q in range(0,count):
                        start_pos = seq.find(sub_seq,rem_start)
                        end_pos = start_pos + len(sub_seq)
                        start_pos += 1
                        distance = start_pos - prv_end_pos  
                        
                        print("Sequence = {} : Size = {} : Distance = {} : Start-Pos = {} : End-Pos = {}".format(str(sub_seq), len(sub_seq),distance,start_pos,end_pos))
                        df.loc[str(sub_seq)+ '~' + str(start_pos)] =[str(sub_seq), len(sub_seq),distance,start_pos,end_pos] 
                        rem_start = start_pos + 1 
                        prv_end_pos = end_pos
                
                   

                    
print(df)    
# Converting to a CSV file because the sequence is too long
df.sort_values(['Size'],ascending=[True]).to_csv("match_new.csv", index=False)        
t_sec = round(time.time() - _start_time)
(t_min, t_sec) = divmod(t_sec,60)
(t_hour,t_min) = divmod(t_min,60) 
print('\nTime taken: {}hour:{}min:{}sec'.format(t_hour,t_min,t_sec))       
            

def tic():
    global _start_time 
    _start_time = time.time()

def tac():
    t_sec = round(time.time() - _start_time)
    (t_min, t_sec) = divmod(t_sec,60)
    (t_hour,t_min) = divmod(t_min,60) 
    print('Time passed: {}hour:{}min:{}sec'.format(t_hour,t_min,t_sec))
