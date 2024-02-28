from Bio import SeqIO
from Bio.Seq import Seq
import pandas as pd
import time
import numpy
import short_seq_scan as sss

#set debug to false if print and time not needed
_start_time = time.time()
_debug = False

#iteration objects
MIN_SHORT_SEQ_LEN = 6
MAX_SHORT_SEQ_LEN = 17

#debugging tool prints out the srs when found + time
def main (filename,output,debug):
    global _debug
    _debug = debug
    data = {'Sequence': [], 'Size': [], 'Distance': [], 'Start-Pos': [], 'End-Pos': []}
    df = pd.DataFrame(data)
    seq_parser(filename,output,df)
    #debug timer
    if _debug:        
        print(df)
        t_sec = round(time.time() - _start_time)
        t_min, t_sec = divmod(t_sec, 60)
        t_hour, t_min = divmod(t_min, 60)
        print('\nTime taken: {}hour:{}min:{}sec'.format(t_hour, t_min, t_sec))

 # records and parses through the fasta file 
def seq_parser(filename,output,df):
    sequences = SeqIO.parse(filename, 'fasta')
    for record in sequences:
            seq_len = len(record.seq)
            seq = record.seq
            build_sub_seq_from_seq(seq,df)
            #csv file sorted by ascending values
    df.sort_values(['Size'], ascending=[True]).to_csv(output, index=False)
  
#iterating through 6-17
def build_sub_seq_from_seq(seq,df):
    #visited sequences stores combinations when visited or found
    global visited_sequences
    visited_sequences = set()
    for i, letter in enumerate(seq):
        for j in range(MAX_SHORT_SEQ_LEN, MIN_SHORT_SEQ_LEN, -1):
            if len(seq[i:i + j]) > MIN_SHORT_SEQ_LEN:
                sub_seq = seq[i:i + j]
                rem_seq = seq[i + 1:]
                find_short_seq(seq,sub_seq,df)

    #circular short sequence look up
    k = MIN_SHORT_SEQ_LEN
    for j in (MIN_SHORT_SEQ_LEN, MAX_SHORT_SEQ_LEN,1):
        if (j > k):
            for l in range(len(seq) % k):
                #finding repeats at the seams via rolls
                roll_arr = numpy.roll([c for c in seq],l)
                roll_seq = ''.join(str(x) for x in roll_arr)
                sub_seq = Seq(roll_seq[0:k])
                find_short_seq(seq,sub_seq,df)

#stores the srs if count is >1 and not in visted sequences 
def find_short_seq(seq,sub_seq,df):
    count = seq.count_overlap(sub_seq)
    if count > 1 and str(sub_seq) not in visited_sequences:
        #comment out this portion if every combination possible is needed 
        df['Sequence'] = df['Sequence'].astype(str)
        df_filter = df[df['Sequence'].str.contains(str(sub_seq))]
        if df_filter.empty:
            for seq_attr in sss.scan_short_sequence(seq,sub_seq):
                #debugging 
                if _debug:
                    print("Sequence = {} : Size = {} : Distance = {} : Start-Pos = {} : End-Pos = {}".format(
                    seq_attr.sub_seq, seq_attr.length, seq_attr.distance, seq_attr.start_pos, seq_attr.end_pos))
                df.loc[seq_attr.sub_seq + '~' + str(seq_attr.start_pos)] = [seq_attr.sub_seq, seq_attr.length, seq_attr.distance,
                                                                      seq_attr.start_pos, seq_attr.end_pos]
                #adding to the csv file if not in vistited_sequences
                if count <= 1 or str(sub_seq) not in visited_sequences:
                    visited_sequences.add(str(sub_seq))

#cli call function
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Find Short Sequences from Fasta File')
    parser.add_argument('--input_file_name', metavar='path', required=True,
                        help='the path to fasta file')
    parser.add_argument('--output_csv_file_name', metavar='path', required=True,
                        help='path to out csv file name')
    parser.add_argument('--debug', metavar='path', required=False, action=argparse.BooleanOptionalAction,
                        help='--debug/--no-debug prints debug information')
    
    args = parser.parse_args()
    main(filename=args.input_file_name, output = args.output_csv_file_name, debug = args.debug)



