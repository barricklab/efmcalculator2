class SeqAttr:
    def __init__(self,sub_seq,distance, start_pos, end_pos):
        self.sub_seq = str(sub_seq)
        self.length = len(sub_seq)
        self.distance = distance
        self.start_pos = start_pos
        self.end_pos = end_pos

def scan_short_sequence(seq,sub_seq):
        
        distance = 1
        start_pos = 0
        rem_start = 0
        end_pos = 0
        prv_end_pos = 0
        
        start_pos = seq.find(sub_seq, rem_start)
        while (start_pos != -1):
            
            end_pos = start_pos + len(sub_seq)
            start_pos += 1
            distance = start_pos - prv_end_pos
            yield SeqAttr(sub_seq,distance,start_pos,end_pos)
            rem_start = start_pos + 1
            prv_end_pos = end_pos
            start_pos = seq.find(sub_seq, rem_start)
        
        
            
