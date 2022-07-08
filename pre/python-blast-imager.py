# Jun01, 2022, ms
# python-blast-imager.py

import sys
from collections import defaultdict
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager


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
    # for query start min and query end max in hsp
    qs_min = 1e20
    qe_max = 0

    with open(fmt6, mode='r') as FMT6:
        for line in FMT6:
            parts = line.rstrip().split('\t')
            # get what i want
            qid = parts[col_labels.index('qseqid')]
            sid = parts[col_labels.index('sseqid')]
            p_i = float(parts[col_labels.index('pident')])
            q_s = int(parts[col_labels.index('qstart')])
            q_e = int(parts[col_labels.index('qend')])
            s_s = int(parts[col_labels.index('sstart')])
            s_e = int(parts[col_labels.index('send')])
            # populate buckets
            if sid not in hsp_bucket:
                # `not in sids` will do the same.  I followed original pl code.
                sids.append(sid)
            hsp_bucket[sid].append((q_s, q_e, s_s, s_e, p_i))  # all in str
            hsp_count += 1
            # min, max
            if q_s < qs_min:
                qs_min = q_s
            if q_e > qe_max:
                qe_max = q_e

    return qid, sids, hsp_bucket, hsp_count, qs_min, qe_max


def colormap(percentage, colors):
    if percentage >= 100:
        idx = 0
    else:
        idx = int((109 - percentage) / 10)

    return colors[idx]


def scale(x, qs_min, scale, left):
    scale = (x - qs_min) * scale + left
    return scale


def mainStory():
    # get fmt6 file path
    fmt6 = sys.argv[1]
    # parse hit info from fmt6
    query, sids, hsp_bucket, hsp_count, qs_min, qe_max = fmt6parser(fmt6)

    # graph parameter setup
    LEFT, BODY, RIGHT = (150, 550, 50)
    HEADER_H, FOOTER_H = (40, 20)
    LINE_W, HSP_SEP, SUBJ_SEP = (3, 14, 18)
    # drawing canvas size
    v_size = HEADER_H + FOOTER_H + (HSP_SEP * hsp_count)
    v_size += (SUBJ_SEP * len(hsp_bucket.keys()))
    if v_size < 100:  # floor
        v_size = 100
    if v_size > 4000:  # ceiling
        v_size = 4000
    h_size = sum((LEFT, BODY, RIGHT))
    SCALE = BODY / (qe_max - qs_min + 1)
    # colors
    colors = (
        (0, 0, 0),  # black
        (196, 0, 255), (0, 0, 255), (0, 255, 255), (0, 255, 0),
        (255, 255, 0), (255, 196, 0), (255, 0, 0), (128, 128, 128),
        (255, 255, 255)  # white
    )

    # creating ImageDraw object
    img = Image.new("RGB", (h_size, v_size), colors[-1])  # white is the last
    draw = ImageDraw.Draw(img)
    # fonts
    font = font_manager.FontProperties(family='monospace')
    font_b = font_manager.FontProperties(family='monospace', weight='bold')
    font_file = font_manager.findfont(font)
    font_file_b = font_manager.findfont(font_b)

    font_tiny = ImageFont.truetype(font_file, size=8)
    font_small = ImageFont.truetype(font_file, size=10)
    font_medium_b = ImageFont.truetype(font_file_b, size=12)

    # header part
    # query name
    draw.text(
        (5, HEADER_H - 8), query[:18], font=font_medium_b, fill=colors[0]
    )
    # query as a line
    draw.line([(LEFT, HEADER_H), (LEFT + BODY, HEADER_H)], fill=colors[0])
    # query hit start and end
    draw.text(
        (LEFT, HEADER_H - 20), str(qs_min),
        font=font_small, fill=colors[0]
    )
    draw.text(
        (LEFT + BODY, HEADER_H - 20), str(qe_max),
        font=font_small, fill=colors[0]
    )

    # percent identity color scale
    draw.text((670, 5), "% Identity", font=font_small, fill=colors[0])
    # color belt
    for perc in range(20, 101, 10):
        x = LEFT + BODY/2 + perc*2
        draw.rectangle([(x, 5), (x + 10, 15)], fill=colormap(perc, colors))
        draw.text((x, 17), str(perc), font=font_tiny, fill=colors[0])

    # alignments
    depth = defaultdict(int)
    v = 0
    for sid in sids:
        v += SUBJ_SEP
        draw.text(
            (10, HEADER_H + v + 9), sid[:18], font=font_small, fill=colors[0]
        )
        for hsp in hsp_bucket[sid]:
            v += HSP_SEP
            (q_s, q_e, s_s, s_e, p_i) = hsp
            if s_s < s_e:
                strand = '+'
            else:
                strand = '-'
            # map start end location in the drawing
            x1 = scale(q_s, qs_min, SCALE, LEFT)
            x2 = scale(q_e, qs_min, SCALE, LEFT)
            print(x1, x2)
            y = HEADER_H + v
            # for x in range(x1, x2 + 1):
            #     depth[x] += 1
            # draw hsp line, location on both end
            draw.rectangle(
                [(x1, y), (x2, y + LINE_W)], fill=colormap(p_i, colors)
            )
            draw.text(
                (x1 - 5*len(str(q_s)), y - 5), str(q_s),
                font=font_tiny, fill=colors[0]
            )
            draw.text(
                (x2 + 2, y - 5), str(q_e),
                font=font_tiny, fill=colors[0]
            )
            draw.text(
                (x1 - 5*len(str(s_s)), y + 2), str(s_s),
                font=font_tiny, fill=colors[0]
            )
            draw.text(
                (x2 + 2, y + 2), str(s_e) + ' ' + strand,
                font=font_tiny, fill=colors[0]
            )

    # alignments depth in the header


    img.save('py-b-imgr-test.png')



if __name__ == "__main__":
    mainStory()
