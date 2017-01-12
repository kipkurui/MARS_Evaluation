"""
Fetch is python code that can be used to dierectly fetch and convert motifs
to PWM format. 

"""
import MySQLdb
import os
import sys
mydb= MySQLdb.connect(host='localhost',
		      user = 'root',
		      passwd ='#',
		      db = 'MAT_db')
cur=mydb.cursor()

def getmotif(ID):
    statement="SELECT * FROM MATRIX_DATA WHERE MATRIX_ID='%i'" % ID
    command = cur.execute(statement)
    A=[]
    C=[]
    G=[]
    T=[]
    row=cur.fetchall()
    for i in row:
        if i[2]=="A":
            A.append(i[4])
        elif i[2]=="C":
            C.append(i[4])
        elif i[2]=="G":
            G.append(i[4])
        else:
            T.append(i[4])
    statement="SELECT * FROM MATRIX WHERE ID='%i'" % ID
    command = cur.execute(statement)
    det=cur.fetchall()
    mid=det[0][1]+"."+str(det[0][4])
    name=det[0][2]

    motif=("\nMOTIF %s %s\n\n"% (mid,name))
    print motif,
    
    #l=A[0]+C[0]+G[0]+T[0] #use this to convert to PWM from PFM
    header=("letter-probability matrix: alength= 4 w= %s nsites= 20 E= 0\n" % (str(len(A)-1)))
    print header,
    for i in range(1,len(A)):
        out=("  %.6f\t  %.6f\t  %.6f\t  %.6f\t\n" %(A[i],float(C[i]),float(G[i]),float(T[i])))
        print out,

def forasseessmotifs(ID):
    statement="SELECT * FROM MATRIX_DATA WHERE MATRIX_ID='%i'" % ID
    command = cur.execute(statement)
    areapwm = {}
    areapwm["A"] = []
    areapwm["C"] = []
    areapwm["G"] = []
    areapwm["T"] = []
    row=cur.fetchall()
    for i in row:
        if i[1]=="A":
            areapwm["A"].append(i[3])
        elif i[1]=="C":
            areapwm["C"].append(i[3])
        elif i[1]=="G":
            areapwm["G"].append(i[3])
        else:
            areapwm["T"].append(i[3])
    statement="SELECT * FROM MATRIX WHERE ID='%i'" % ID
    command = cur.execute(statement)
    det=cur.fetchall()
    mid=det[0][2]+"."+str(det[0][3])
    name=det[0][4]
    return areapwm

def getchip(ID):
    statement="SELECT CHIP_ID FROM CHIP_SEQ WHERE TF_ID='%s'" % ID
    command = cur.execute(statement)
    tfid=cur.fetchall()

    statement="SELECT AT_100 FROM CHIP_DATA WHERE CHIP_ID='%i'" % tfid[0]
    command = cur.execute(statement)
    tfid=cur.fetchall()
    tflist=[]
    for i in tfid:
        tflist.append(i[0])
    #print tflist
    return tflist
    
   
    
#fetch all max motifs
#command = cur.execute("SELECT * FROM JASPAR_test.MATRIX WHERE lower(NAME) LIKE '%max%'")
#fetch all max motifs except those which are dimeric
tf="Egr1" #sys.argv[1]
print '''MEME version 4.4\n
ALPHABET= ACGT\n\nstrands: + -\n
Background letter frequencies (from uniform background):
A 0.25000 C 0.25000 G 0.25000 T 0.25000 \n''',

command = cur.execute("SELECT TF_ID FROM MATRIX WHERE lower(MOTIF_NAME) LIKE '%"+tf.lower()+"%' AND lower(MOTIF_NAME) NOT LIKE '%:%' ")
IDs=[]
a=cur.fetchall()
ID=a[0][0]
#tflist=getchip(ID) 
command = cur.execute("SELECT ID FROM MATRIX WHERE TF_ID LIKE '%s'" % ID)
a=cur.fetchall()
for i in a:
    IDs.append(i[0])
for i in IDs:
    #print i,
    #what=forasseessmotifs(i)
    getmotif(i)

