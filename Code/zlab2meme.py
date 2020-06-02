"""
Conver motifs from the Zlab to MEME
"""

import sys

if len(sys.argv)<3:
    print("\nUsage: python zlab2meme.py <file-input> <motif> \n")
    sys.exit(1)
meme=sys.argv[1]
#output=sys.argv[2]
motif=sys.argv[2]
with open(meme) as en:
    no=0
    n=1
    m=0
    for line in en:
        if "letter-probability" in line:
            no=int(line.split(" ")[5])
            n=0
            m+=1
            print "\nMOTIF "+motif+"_ENCODE_ChIP-seq_"+str(m)+"\n"
            print line,
            continue
        if n<no:
            print line,
            n+=1
            continue
