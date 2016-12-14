"""
These are some codes used in analsys in the notebooks that are 
are many and deserve their own modules
"""
from math import exp

def gomeroccupancyscore(pwm_dictionary, seq):
    """
    Takes as input a PWM dictionary, and a sequences and
    computes the gomer score
    """
    if "N" in seq:
        return 0
    else:
        pwm_length = len(pwm_dictionary["A"])
        gomer_occupancy = 1
        area_pwm_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(pwm_length - 1, 1, -1):
            prod_gomer = 1
            prod_gomer_rc = 1
            for j in range(pwm_length):
                if j <= i:
                    prod_gomer *= 0.25
                    prod_gomer_rc *= 0.25
                elif (j + i) > len(seq)-1:
                    prod_gomer *= 0.25
                    prod_gomer_rc *= 0.25
                else:
                    # print "got to else"
                    s = seq[j + i]
                    prod_gomer *= pwm_dictionary[s][j]
                    prod_gomer_rc *= area_pwm_rc[s][j]
            gomer_occupancy *= (1 - prod_gomer) * (1 - prod_gomer_rc)
        for i in range(len(seq) - 1):
            prod_gomer = 1
            prod_gomer_rc = 1
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq)-1:
                    prod_gomer *= 0.25
                    prod_gomer_rc *= 0.25
                else:
                    prod_gomer *= pwm_dictionary[seq[j + i]][j]
                    prod_gomer_rc *= area_pwm_rc[seq[j + i]][j]
            gomer_occupancy *= (1 - prod_gomer) * (1 - prod_gomer_rc)
        gomer_occupancy = 1 - gomer_occupancy
        return gomer_occupancy


def energyscore(pwm_dictionary, seq):
    """
    Score sequences using the beeml energy scoring approach.

    Borrowed greatly from the work of Zhao and Stormo

    P(Si)=1/(1+e^Ei-u)

    Ei=sumsum(Si(b,k)e(b,k))

    Previous approaches seem to be using the the minimum sum of the
    energy contribution of each of the bases of a specific region.

    This is currently showing some promise but further testing is
    needed to ensure that I have a robust algorithm.
    """
    if "N" in seq:
        return 0
    else:
        energy_list = []
        pwm_length = len(pwm_dictionary["A"])
        pwm_dictionary_rc = rc_pwm(pwm_dictionary, pwm_length)
        for i in range(len(seq) - 1):
            energy = 0
            energy_rc = 0
            for j in range(pwm_length - 1):
                if (j + i) >= len(seq):
                    energy += 0.25
                    energy_rc += 0.25
                else:
                    energy += pwm_dictionary[seq[j + i]][j]
                    energy_rc += pwm_dictionary_rc[seq[j + i]][j]

                energy_list.append(1 / (1 + (exp(energy))))
                energy_list.append(1 / (1 + (exp(energy_rc))))
        energy_score = min(energy_list)
        return energy_score

    
def get_from_cluster_key(cluster_key, top_10):
    """
    Given top ten motifs based on 
    """
    extract =[]
    for i in top_10:
        with open(cluster_key) as cluster:
            for line in cluster:
                if i in line:
                    lines = line.split()
                    if len(lines)> 1:
                        adds = lines[1].split(",")
                        extract.extend(adds)
                    else:
                        extract.extend(lines)
    return extract

def extract_pfm(pfm_in, pfm_out, extract):
    with open(pfm_out, "w") as pfm_o:
        stat = 0
        for mot in extract:
            with open(pfm_in) as pfm_i:
                for line in pfm_i:
                    if ">" in line:
                        if mot == line.strip().split(">")[1]:
                            pfm_o.write(line)
                            #print line,
                            stat = 1
                            continue
                    if ">" in line:
                        if mot in line:
                            continue
                        else:
                            stat = 0
                    if stat == 1:
                        pfm_o.write(line)

def extract_meme(meme_in, motif, meme_out, raw_dict):
    with open(meme_in) as f1:
        with open(meme_out, "a") as out_meme:
            lines = f1.readlines()
            for i, line in enumerate(lines):
                head = line.split()
                if motif in line and motif == head[1]:
                    k = i
                    out_meme.write("\n"+lines[i].strip()+"_"+raw_dict[motif]+"\n\n")
                    if "log-odds" in lines[i+1]:
                        odds = lines[k+2].split()
                        for j in range(2, (int(odds[5])+3)*2):
                            out_meme.write(lines[i+j]),
                    elif "letter-probability" in lines[i+2]:
                        out_meme.write(lines[i+2])
                        odds = lines[k+2].split()
                        for j in range(0, (int(odds[5]))):
                            out_meme.write(lines[i+3+j])
                    else:
                        pass


def get_dict(raw_file):
    from collections import OrderedDict
    raw_dict = OrderedDict()
    with open(raw_file) as raw_in:
        for line in raw_in:
            raw_dict[line.split()[0]] = line.split()[-1]
    return raw_dict


def extract_scored_meme(meme_in, out_meme, raw_dict):
    meme_head = """MEME version 4.4\n\nALPHABET= ACGT\n\nstrands: + -\n
Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000\n"""
    with open(out_meme, "w") as meme_out:
        meme_out.write(meme_head)
    for key in raw_dict.iterkeys():
        extract_meme(meme_in, key, out_meme, raw_dict)

def extract_meme_list(meme_in, out_meme, mot_list):
    """
    Use this to extract motifs in a given list from the main file. This differes
    from extract_scored since it does not add the scores of the motifs to the 
    meme header name
    """
    meme_head = """MEME version 4.4\n\nALPHABET= ACGT\n\nstrands: + -\n
Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000\n"""
    with open(out_meme, "w") as meme_out:
        meme_out.write(meme_head)
    for mot in mot_list:
        extract_meme_clean(meme_in, mot, out_meme)

def extract_meme_clean(meme_in, motif, meme_out):
    """
    Get the given motif from file
    """
    with open(meme_in) as f1:
        with open(meme_out, "a") as out_meme:
            lines = f1.readlines()
            for i, line in enumerate(lines):
                head = line.split()
                if motif in line and motif == head[1]:
                    k = i
                    out_meme.write("\n"+lines[i].strip()+"\n\n")
                    if "log-odds" in lines[i+1]:
                        odds = lines[k+2].split()
                        for j in range(2, (int(odds[5])+3)*2):
                            out_meme.write(lines[i+j]),
                    elif "letter-probability" in lines[i+2]:
                        out_meme.write(lines[i+2])
                        odds = lines[k+2].split()
                        for j in range(0, (int(odds[5]))):
                            out_meme.write(lines[i+3+j])
                    else:
                        pass

def meme2gimme(meme, gimme):
    with open(meme) as motif:
        with open(gimme, 'w') as gmotif:
            for line in motif:
                if line.startswith("MOTIF"):
                    if len(line.split(" ")) > 2:
                        gmotif.write(">"+line.split(" ")[1]+"\n")
                    else:
                        gmotif.write(">"+line.split(" ")[1])
                elif line.startswith('letter-probability'):
                    continue
                elif line.startswith('  '):
                    a = line.split()
                    if len(a)>0:
                        gmotif.write(a[0]+'\t'+a[1]+'\t'+a[2]+'\t'+a[3]+'\n')
                    else:
                        continue
                elif line.startswith("\n"):
                    continue
                else:
                    continue

def gimme2meme(gimme_in,meme_out):

    As=[]
    Cs=[]
    Gs=[]
    Ts=[]
    motifname="motifname"
    with open(meme_out, "w") as f3:
        f3.write("""MEME version 4.4\n\nALPHABET= ACGT\n\n
strands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n
        """)
        with open(gimme_in) as f1:
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
                    As+=[line.strip().split()[0]]
                    Cs+=[line.strip().split()[1]]
                    Gs+=[line.strip().split()[2]]
                    Ts+=[line.strip().split()[3]]
        motif=("\nMOTIF %s\n\n"% (motifname))
        f3.write(motif)
        header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(As)))) 
        f3.write(header)
        for i in range(0,len(As)):
            out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %(float(As[i]),float(Cs[i]),float(Gs[i]),float(Ts[i])))
            f3.write(out)



def gimmepfm2meme(pfm_in, meme_out):
    
    import sys
    As=[]
    Cs=[]
    Gs=[]
    Ts=[]
    in_file = pfm_in

    motifname="motifname"
    with open(meme_out, "w") as f3:
        f3.write("""MEME version 4.4\n\nALPHABET= ACGT\n\n
strands: + -\n\nBackground letter frequencies (from uniform background):\nA 0.25000 C 0.25000 G 0.25000 T 0.25000\n
        """)
        with open(in_file) as f1:
            for line in f1:
                if line.startswith(">"):
                    if len(As)>0:
                        motif=("\nMOTIF %s\n\n"% (motifname))
                        f3.write(motif)
                        header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(As)))) 
                        f3.write(header)
                        for i in range(0,len(As)):
                            totals = float(As[i])+float(Cs[i])+float(Gs[i])+float(Ts[i])
                            out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %
                                 (float(As[i])/totals,float(Cs[i])/totals,float(Gs[i])/totals,float(Ts[i])/totals))
                            f3.write(out)
                    else:
                        f3.write("\n")
                    As=[]
                    Cs=[]
                    Gs=[]
                    Ts=[]
                    motifname=line.strip().split(">")[1]
                else:
                    As+=[line.strip().split()[0]]
                    Cs+=[line.strip().split()[1]]
                    Gs+=[line.strip().split()[2]]
                    Ts+=[line.strip().split()[3]]
            motif=("\nMOTIF %s\n\n"% (motifname))
            f3.write(motif)
            header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(As)))) 
            f3.write(header)
            for i in range(0,len(As)):
                totals = float(As[i])+float(Cs[i])+float(Gs[i])+float(Ts[i])
                out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %
                     (float(As[i])/totals,float(Cs[i])/totals,float(Gs[i])/totals,float(Ts[i])/totals))
                f3.write(out)
