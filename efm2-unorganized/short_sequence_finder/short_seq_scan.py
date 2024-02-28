class SeqAttr:
    def __init__(self, sub_seq, distance, start_pos, end_pos):
        self.sub_seq = str(sub_seq)
        self.length = len(sub_seq)
        self.distance = distance
        self.start_pos = start_pos
        self.end_pos = end_pos


def scan_short_sequence(seq, sub_seq, seq_len,isCircular, count):
    distance = 1
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_end_pos = 0
    prv_start_pos = 0

    start_pos = seq.find(sub_seq, rem_start)
    while (start_pos != -1):

        end_pos = start_pos + len(sub_seq)
        if end_pos >= seq_len:
            end_pos = end_pos - seq_len
        start_pos += 1
        distance = start_pos - prv_end_pos

        # fixes distance for 
        if isCircular == True:
            if distance > seq_len / 2:
                distance = (seq_len + prv_start_pos) - end_pos
            if start_pos < seq_len and not (count == 2 and end_pos <= 20):
                yield SeqAttr(sub_seq, distance, start_pos, end_pos)
        else:
            yield SeqAttr(sub_seq, distance, start_pos, end_pos)
        rem_start = start_pos + 1
        prv_end_pos = end_pos
        prv_start_pos = start_pos
        start_pos = seq.find(sub_seq, rem_start)



