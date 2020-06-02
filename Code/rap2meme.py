"""
Convert RAP formatted motifs to meme
"""
import sys

#pref=sys.argv[5]

if len(sys.argv)<4:
    print("\nUsage: python rap2meme.py <rap-input> <meme-out> <motifname> \n")
    sys.exit(1)
rap=sys.argv[1]
output=sys.argv[2]
#memehead=sys.argv[3]
motifname=sys.argv[3]

motif=("MOTIF %s %s\n\n"% (motifname,motifname))
with open(rap) as f1:
    with open(output, "w") as f3:
        Lines=f1.readlines()
        for i, line in enumerate(Lines):
            if line.startswith("A:"):
                As=line.split()
            if line.startswith("C:"):
                Cs=line.split()
            if line.startswith("G:"):
                Gs=line.split()
            if line.startswith("T:"):
                Ts=line.split()
        f3.write('''MEME version 4.4\n
ALPHABET= ACGT\n\nstrands: + -\n
Background letter frequencies (from uniform background):\n
A 0.25000 C 0.25000 G 0.25000 T 0.25000 \n\n''')
        f3.write(motif)
        header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(As)-1))) 
        f3.write(header)
        for i in range(1,len(As)):
            out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %(float(As[i]),float(Cs[i]),float(Gs[i]),float(Ts[i])))
            f3.write(out)

