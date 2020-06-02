"""
Convert POUR formatted motif to meme format
"""
import sys


As=[]
Cs=[]
Gs=[]
Ts=[]
motifname="motifname"
with open("../Motifs/Kheradpour-2013.meme", "w") as f3:
    f3.write('''MEME version 4.4\n
ALPHABET= ACGT\n\nstrands: + -\n
Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000''')
    with open("/home/kipkurui/Project/Motif_Assessment/Data/Motifs/Kheradpour-2013.txt") as f1:
        for line in f1:
            if line.startswith(">"):
                if len(As)>0:
                    motif=("\nMOTIF %s\n\n"% (motifname))
                    f3.write(motif)
                    header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(As)))) 
                    f3.write(header)
                    for i in range(0,len(As)):
                        out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %(float(As[i]),float(Cs[i]),float(Gs[i]),float(Ts[i])))
                        f3.write(out)
                else:
                    f3.write("\n")
                As=[]
                Cs=[]
                Gs=[]
                Ts=[]
                motifname=line.strip().split(">")[1]
            else:
                As+=[line.strip().split()[1]]
                Cs+=[line.strip().split()[2]]
                Gs+=[line.strip().split()[3]]
                Ts+=[line.strip().split()[4]]

