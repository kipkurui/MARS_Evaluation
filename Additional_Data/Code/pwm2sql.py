import collections

def getmotif(meme, motif="MOTIF"):
    """
    Extract a motif from meme file given a unique motif
    name and create dictionary for sequence scoring

    Default motif name is keyword MOTIF for single motif files. 
    """

    areapwm = collections.OrderedDict()
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    flag = 0
    check = 0
    ID=2459
    with open(meme, "r") as f1:
        for line in f1:
            if str(motif)in line and str(motif)==line.split()[0]:
                BASE_ID=line.split()[1]
                NAME=line.split()[1]
                flag = 1
                areapwm = collections.OrderedDict()
                areapwm["A"] = []
                areapwm["C"] = []
                areapwm["G"] = []
                areapwm["T"] = []
                check = 0
            if "letter-probability" in line and flag == 1:
                w = line.split(" ")[5]
                flag += 1
                continue
            if flag == 2 and int(check) < int(w):
                if line == "\n":
                    continue
                else:
                    words = line.split()
                    areapwm["A"].append(float(words[0]))
                    areapwm["C"].append(float(words[1]))
                    areapwm["G"].append(float(words[2]))
                    areapwm["T"].append(float(words[3]))
                    check += 1
            #if flag==2 and "URL" in line:
                #URL=line.split()[1]
                #flag=3
            if flag==2:
                if int(check)==int(w):
                    ID+=1
                    check = 0
                    flag = 1
                    for key in areapwm:
                        for i in range(len(areapwm[key])):
                            print "%i\t%s\t%i\t%f\n" % (ID,key,i+1,areapwm[key][i]),
            #if flag==2:
                #if int(check)==int(w):
                    #ID+=1
                    #check = 0
                    #flag = 1
                    #for key in areapwm:
                        #for i in range(len(areapwm[key])):
                            #print "%i\t%s\t%i\t%f\n" % (ID,key,i,areapwm[key][i]),



def geturl(meme):
    flag=0
    ID=2459
    with open(meme) as f1:
        for line in f1:
            if "MOTIF" in line:
                ID+=1
                BASE_ID=line.split()[1]
                NAME=line.split()[2]
                print "%i\t%s\t%s\t%i\n" % (ID,BASE_ID,NAME,ID),
                
meme="../Motifs/Clean_DB/HOCOMOCOv9.meme"
getmotif(meme)
#geturl(jas)
