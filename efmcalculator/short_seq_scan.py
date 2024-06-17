import logging

class SeqAttr:
    def __init__(self, sub_seq, distance, start_pos, end_pos, note):
        self.sub_seq = str(sub_seq)
        self.length = len(sub_seq)
        self.distance = distance
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.note = note

logger = logging.getLogger(__name__)
logging.getLogger(__name__).addHandler(logging.NullHandler())

def _build_seq_attr(sub_seq, seq_len, start_positions, isCircular, count):
    distance = 1
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_end_pos = 0
    prv_start_pos = 0


    for start_pos in start_positions:
        note = ""
        end_pos = start_pos + len(sub_seq)
        #if end_pos >= seq_len:
         #   end_pos = end_pos - seq_len
        start_pos += 1
        distance = start_pos - prv_end_pos

        # fixes distance for
        if distance == 1:
            note = "SSR"
        if isCircular == True:
            if distance > seq_len / 2:
                distance = (seq_len + prv_start_pos) - end_pos
            if start_pos < seq_len and end_pos > seq_len: # if repeat wraps around
                # fix end_pos
                end_pos = end_pos - seq_len
        # if overlapping
        if distance < 1:
            note = "skip for SSR"
            if count == 2:
                note = "skip"
        yield SeqAttr(sub_seq, distance, start_pos, end_pos, note)
        rem_start = start_pos + 1
        prv_end_pos = end_pos
        prv_start_pos = start_pos


def _find_repeat_positions(seq, sub_seq, seq_len, isCircular, count):
    '''testing function for sequence scanning'''
    start_pos = 0
    rem_start = 0
    end_pos = 0
    prv_start_pos = 0
    start_positions = []

    start_pos = seq.find(sub_seq, rem_start)
    while start_pos != -1:
        start_positions.append(start_pos)
        start_pos += 1
        rem_start = start_pos + 1

        if len(sub_seq) == 1:
            start_pos = seq.find(sub_seq, rem_start-1)
        else:
            start_pos = seq.find(sub_seq, rem_start)
    return start_positions


