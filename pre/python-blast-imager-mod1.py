# Jun12, 2022, ms
# python-blast-imager-mod1.py
# modify the *-pl2py.py version for
# - better readability (move some global variables to different scope)
# - automatic naming of the output file

import os
import sys
from collections import defaultdict
from PIL import Image, ImageDraw, ImageFont
from matplotlib import font_manager


# === GLOBALS ===
# drawing parameter setup
LEFT, BODY, RIGHT = (150, 550, 50)
HEADER_H, FOOTER_H = (40, 20)
LINE_W, HSP_SEP, SUBJ_SEP = (3, 14, 18)

# color palette
COLORS = (
    (0, 0, 0),  # black
    (196, 0, 255), (0, 0, 255), (0, 255, 255), (0, 255, 0),
    (255, 255, 0), (255, 196, 0), (255, 0, 0), (128, 128, 128)
)
WHITE = (255, 255, 255)
BLACK = COLORS[0]

# fonts
# this part is to minimize platform difference with font_manager
# from matplotlib
font = font_manager.FontProperties(family='monospace')
font_b = font_manager.FontProperties(family='monospace', weight='bold')
font_file = font_manager.findfont(font)
font_file_b = font_manager.findfont(font_b)
# try to mimic gdTinyFont, gdSmallFont, gdMediumBoldFont in the perl GD
TINY = ImageFont.truetype(font_file, size=8)
SMALL = ImageFont.truetype(font_file, size=10)
MEDIUM_B = ImageFont.truetype(font_file_b, size=12)


# === functions ===
def generateOutFilePath(fmt6, outdir_switch=0):
    """
    generate output file path based on the input file path
    outdir_switch:
        0: current working dir (default)
        non zero: the same as fmt6
    """
    head, tail = os.path.split(fmt6)
    stem = os.path.splitext(tail)[0]
    out_file = stem + '-imager.png'
    if outdir_switch == 0:
        return out_file
    else:
        return os.path.join(head, out_file)


def getDrawSize(hsp_count, hsp_bucket):
    """
    set the drawing canvas size
    """
    v_size = HEADER_H + FOOTER_H + (HSP_SEP * hsp_count)
    v_size += (SUBJ_SEP * len(hsp_bucket.keys()))
    # floor, ceil the v_size
    if v_size < 100:
        v_size = 100
    if v_size > 4000:
        v_size = 4000
    h_size = sum((LEFT, BODY, RIGHT))

    return h_size, v_size


def fmt6parser(fmt6):
    """
    reading file and storing info
    - get query name
    - store subject ids in a list
        this is to keep their order in the fmt6 file and call them in the order
    - store {sid: [[q_s, q_e, s_s, s_e, p_i],,,]}
        i've read python will keep item order as dict is populated from
        some point of version, but i pretend i do not know that

    # blast+ fmt6 columns (assuming it is not customized)
    0: qseqid  ---> query id, qid
    1: sseqid  ---> subject id, sid
    2: pident  ---> percent identity, p_i
    3: length
    4: mismatch
    5: gapopen
    6: qstart  ---> query start, q_s
    7: qend  ---> query end, q_e
    8: sstart  ---> subject start, s_s
    9: send  ---> subject end, s_e
    10: evalue
    11: bitscore
    """
    col_labels = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch',
        'gapopen', 'qstart', 'qend', 'sstart', 'send',
        'evalue', 'bitscore'
    ]

    # info containers
    qid = ''
    sids = []  # subject ids in order in the file (from top to down)
    hsp_bucket = defaultdict(list)
    # key: sid
    # value: list of [q_s, q_e, s_s, s_e, p_i]
    hsp_count = 0
    # for query start min and query end max in hsp
    qs_min = 1e20  # use an extremely big number
    qe_max = 0

    with open(fmt6, mode='r') as FMT6:
        for line in FMT6:
            parts = line.rstrip().split('\t')
            # get what I want
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


def colormap(percentage):
    """
    map percent identity to pre-defined color palette
    """
    if percentage >= 100:
        idx = 0
    else:
        idx = int((109 - percentage) / 10)

    return COLORS[idx]


def scale(x, qs_min, scale, left):
    """
    map q_s or q_e (arg x) of a HSP to a corresponding x coordinate
    in the drawing
    """
    scale = (x - qs_min) * scale + left

    return scale


def mainStory():
    # parse the blast result
    # get fmt6 file path
    fmt6 = sys.argv[1]
    # parse hit info from fmt6
    query, sids, hsp_bucket, hsp_count, qs_min, qe_max = fmt6parser(fmt6)

    # let's start drawing!
    # determine drawing board size
    h_size, v_size = getDrawSize(hsp_count, hsp_bucket)
    # ratio for drawing width and max hit length
    scale_ratio = BODY / (qe_max - qs_min + 1)

    # creating ImageDraw object
    img = Image.new("RGB", (h_size, v_size), WHITE)
    draw = ImageDraw.Draw(img)

    # header part
    # query name
    draw.text(
        (5, HEADER_H - 8), query[:18], font=MEDIUM_B, fill=BLACK
    )
    # query line
    draw.line([(LEFT, HEADER_H), (LEFT + BODY, HEADER_H)], fill=BLACK)
    # query hit start and end
    draw.text(
        (LEFT, HEADER_H - 20), str(qs_min), font=SMALL, fill=BLACK
    )
    draw.text(
        (LEFT + BODY, HEADER_H - 20), str(qe_max), font=SMALL, fill=BLACK
    )

    # percent identity color scheme
    draw.text((645, 5), "% Identity", font=SMALL, fill=BLACK)
    # draw color palette
    for perc in range(20, 101, 10):
        x = LEFT + BODY/2 + perc*2
        draw.rectangle([(x, 5), (x + 10, 15)], fill=colormap(perc))
        draw.text((x, 17), str(perc), font=TINY, fill=BLACK)

    # alignments
    depth = defaultdict(int)
    v = 0
    for sid in sids:
        v += SUBJ_SEP
        draw.text((10, HEADER_H + v + 9), sid[:18], font=SMALL, fill=BLACK)
        for hsp in hsp_bucket[sid]:
            v += HSP_SEP
            (q_s, q_e, s_s, s_e, p_i) = hsp
            if s_s < s_e:
                strand = '+'
            else:
                strand = '-'
            # map start end location in the drawing
            x1 = scale(q_s, qs_min, scale_ratio, LEFT)
            x2 = scale(q_e, qs_min, scale_ratio, LEFT)
            # print(x1, x2)
            # print(int(x1), int(x2))
            # make x1 and x2 into int
            x1 = int(x1)
            x2 = int(x2)
            for x in range(x1, x2 + 1):
                depth[x] += 1

            y = HEADER_H + v
            # draw hsp line, location on both end
            draw.rectangle(
                [(x1, y + 1), (x2, y + 1 + LINE_W)], fill=colormap(p_i)
            )
            # query start
            draw.text(
                (x1 - 5*len(str(q_s)), y - 5), str(q_s), font=TINY, fill=BLACK
            )
            # query end
            draw.text(
                (x2 + 2, y - 5), str(q_e), font=TINY, fill=BLACK
            )
            # subject start
            draw.text(
                (x1 - 5*len(str(s_s)) + 1, y + 2), str(s_s),
                font=TINY, fill=BLACK
            )
            # subject end with strand
            draw.text(
                (x2 + 2, y + 2), str(s_e) + ' ' + strand,
                font=TINY, fill=BLACK
            )

    # alignments depth in the header
    maxDepth = max(depth.values())
    dScale = int(maxDepth/10) + 1
    draw.text(
        (LEFT + BODY + 2, HEADER_H + 2), str(dScale) + "/line",
        font=TINY, fill=BLACK
    )
    # piling dots to show which location hits more
    for i in depth.keys():
        level = int(depth[i]/dScale + 1)
        for j in range(level):
            # draw a line, but actually it is a dot
            draw.line(
                [(i, HEADER_H + (j * 2)), (i, HEADER_H + (j * 2))], fill=BLACK
            )

    # saving png file
    out_path = generateOutFilePath(fmt6)  # out to current dir

    print(out_path, "generated.")
    img.save(out_path)


if __name__ == "__main__":
    mainStory()
