# Jun01, 2022, ms
# python-blast-imager.py

import sys
from collections import defaultdict


def fmt6parser(fmt6):
    """
    reading file and storing info
    - get query name
    - store subject ids in a list
        sbject ids are in following bucket as keys
        this is to keep their order in the fmt6 file and call them in the order
    - store {sid: [[q_s, q_e, s_s, s_e, perc_identity]...]}
        i've read python will keep item order as dict is populated from
        some point of version, but i pretend i do not know that

    # blast+ fmt6 columns (assuming it is not customized because it can be)
    0: qseqid  ---> query id, qid
    1: sseqid  ---> subject id, sid
    2: pident  ---> percent identity, p_i
    3: length  ---> hit length
    4: mismatch  ---> mm count
    5: gapopen  ---> gap count
    6: qstart  ---> q_s
    7: qend  ---> q_e
    8: sstart  ---> s_s
    9: send  ---> s_e
    10: evalue
    11: bitscore
    """
    col_labels = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
        'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore'
        ]

    qid = ''
    sids = []  # subject ids in order
    hsp_bucket = defaultdict(list)
    hsp_count = 0

    with open(fmt6, mode='r') as FMT6:
        for line in FMT6:
            parts = line.rstrip().split('\t')
            # get what i want
            qid = parts[col_labels.index('qseqid')]
            sid = parts[col_labels.index('sseqid')]
            p_i = parts[col_labels.index('pidebt')]
            q_s = parts[col_labels.index('qstart')]
            q_e = parts[col_labels.index('qend')]
            s_s = parts[col_labels.index('sstart')]
            s_e = parts[col_labels.index('semd')]
            # populate buckets
            sids.append(sid)
            hsp_bucket[sid].append((q_s, q_e, s_s, s_e, p_i))
            hsp_count += 1

    return qid, sids, hsp_bucket, hsp_count


def mainStory():
    # get fmt6 file path
    fmt6 = sys.argv[1]
    query, sids, hsp_bucket, hsp_count = fmt6parser(fmt6)
    #





if __name__ == "__main__":
    mainStory()
