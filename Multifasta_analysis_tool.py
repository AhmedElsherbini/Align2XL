# -*- coding: utf-8 -*-
"""
@author: ahmed_elsherbini (drahmedsherbini@yahoo.com)
created in : 7-June 2020
This is a copy right for the author - do not distrbute
dependacies: see below
update: 25/11/2020
"""
#import

import os as os
import sys
from io import StringIO
from Bio import Phylo
from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import Phylo
from Bio.SeqUtils import GC
from Bio import AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from difflib import SequenceMatcher 
import re
from collections import Counter
import copy




print("Hi, let's start. Just answer questions PRECISELY!")

###################################################################################
#%% #in SPYDER this help to chunk the code!

f = input ("1-do you want to extract genes using its start/end postion in aligned_file.afa (more accurate) or fasta  file? y/n: ")
if (f =="y"):
        er = input ("what is the name of the file you want to extract from?")
        ah = input ("start position:")
        med = input ("end position:")
        with open ("filter_by_position_%s.fasta"%(er), "w") as f:
             for seq_record in SeqIO.parse(er, "fasta"):
                  f.write(">"+str(seq_record.id)+"\n")
                  f.write(str(seq_record.seq[int(ah)-1:int(med)+1])+ "\n")
        print("here you are the file : Extract_by_position.fasta")
###########################################################################################
lui = input("2-do you want to extract gene from fasta/multifasta file using upstream and downstream seq? y/n:")
if (lui == "y"):
       sherb = input("what is the name of your fasta/multifasta file:")
       fin = "exctracted_genes.fasta"
       out = open(fin,"w")
       print("make sure you give me from 10 to 15 bp length seq")
       x = input("what is sequnece  the 5'/upstream of the gene ?:")
       y = input("what is sequnece the 3'/downstream of the gene ?:")
       for seq_record in SeqIO.parse(sherb, "fasta"):
           ss = re.findall('%s(.+)%s'%(x,y), str(seq_record.seq))
           if ss != []:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(x)+str(ss[0])+str(y)+"\n")
       out.close()
       print("here you are the file:exctracted_genes.fasta")
#########################################################################################
tui = input("3-do you want to inlucde certain sequences from a multifasta file using a seq pattern (ex: genes with certain pattern)? y/n:")
if (tui == "y"):
       sherb = input("what is the name of your fasta/multifasta file:")
       homdemha = input("what is pattern you want to keep its genes ATTGCGTGTGTGT or KOPGTLSTTSG :")
       fin = "fasta_filtred_by_inculsion.fasta"
       out = open(fin,"w")
       for seq_record in SeqIO.parse(sherb, "fasta"):
           if (seq_record.seq.find(homdemha) != -1):
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq.ungap("-"))+"\n")
     
       print("here you are the file : fasta_filtred_by_seq_inclusion.fasta")
       out.close()


########################################################################################
tui = input("4-do you want to exlude sequences in a multifasta file using sequneces pattern (ex:N,X)? y/n:")
if (tui == "y"):
       sherb = input("what is the name of your fasta/multifasta file:")
       homdemha = input("what is sequnece pattern you want to exclude:")
       fin = "fasta_filtred_by_exclusion.fasta"
       out = open(fin,"w")
       for seq_record in SeqIO.parse(sherb, "fasta"):
           if (seq_record.seq.find(homdemha) != -1):
               continue
           else:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq.ungap("-"))+"\n")

       print("here you are the file : fasta_filtred_by_seq_Exclusion.fasta")
       out.close()
           
#########################################################################################

f = input("5-do you want to  print all > ID headers in your multifasta? y/n:")
if (f == "y"):
    sherb = input("what is the name of the multifasta file?")
    fin = "your_file_headers.txt"
    out = open(fin, "w")
    for seq_record in SeqIO.parse(sherb, "fasta"):
        out.write(">"+str(seq_record.id)+"\n")
    out.close()
    print("here you are the file:your_file_headers.txt")
###########################################################################################
f = input("6-do you want to include sequences using a pattern in their > ID headers (ex:-2019-)? y/n:")
if (f == "y"):
    sherb = input("what is the name of the fasta file?")
    homdemha = input("what is the pattern in seq header you want extract?")
    fin = "extracted_by_header.fasta"
    out = open(fin, "w")
    for seq_record in SeqIO.parse(sherb, "fasta"):
        if (seq_record.id.find(homdemha) != -1):
            out.write(">"+str(seq_record.id)+"\n")
            out.write(str(seq_record.seq)+"\n")
    out.close()
    print("here you are the file:extracted_by_header_inclusion.fasta")

########################################################################################
tui = input("7-do you want to exlude sequences in a multifasta file using pattern in their > ID header? y/n:")
if (tui == "y"):
       sherb = input("what is the name of your fasta/multifasta file:")
       homdemha = input("what is sequnece pattern you want to exclude:")
       fin = "fasta_filtred_by_exclusion.fasta"
       out = open(fin,"w")
       for seq_record in SeqIO.parse(sherb, "fasta"):
           if (seq_record.id.find(homdemha) != -1):
               continue
           else:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq.ungap("-"))+"\n")

       print("here you are the file : fasta_filtred_by_header_Exclusion.fasta")
       out.close()
           

################################################
a = input("8-do you want to translat DNA fasta file on its 1 frame? y/n:")
if (a == "y"):
    zeze = input("what is the name of DNA file?")
    with open ("translated_%s_file.fasta"%(zeze) , "w") as aa_fa:
        for dna_record in SeqIO.parse(zeze, "fasta"):
            aa_fa.write(">"+dna_record.id+ "\n")
            aa_fa.write(str(dna_record.seq.translate(to_stop=True))+"\n")
        aa_fa.close()
        print("here you are the file:translated_%s_file.fasta"%(zeze))

 
    
############################################################################################
#%%
#gc conent and At and number of unkown bases (extra work)
u = input("9-do you want to know GC content and N bases content of your Mmultifasta file? press y/n:")
if (u == "y"):
    file_path_out = input("what is the name of your file?")
    k = [("ID","GC content%")]
    w = [("ID", "unknown bases%")]
    for seq_record in SeqIO.parse(file_path_out,"fasta"):
        k.append((seq_record.id,GC(seq_record.seq)))
        w.append((seq_record.id,(((float(str(seq_record.seq).count("N" or "n"))/len(seq_record.seq)))*100)))
        #I care about unkown bases percentage as i care about GC%
   
    GFG = pd.ExcelWriter("GC_content.xlsx")
    n_bases = pd.ExcelWriter("N_bases.xlsx")
    df = pd.DataFrame(k)
    nf = pd.DataFrame(w)
    df.to_excel(GFG, index = False)
    nf.to_excel(n_bases, index = False)
    GFG.save()
    n_bases.save()
    print("here you are: GC_content.xlsx and N_bases.xlsx sheets")
    #I have provided the output in a list of tupules but as you can convert easily to dic 


################################################################################    
#%%
#to extract conserved _muation from protein
 
w = input("10-do you want to extract the longest conserved seq & the mutations inside your clustal_file.aln? press y/n:")
#C:/Users/ahmed/Downloads/merged_file.aln #kindly know that this code is not adapted to clustal files only
if (w == "y"):
    print("make sure you input file.aln does not have any outliers,indels and outgroups")
    zizo = input("what is the name of your file.aln:")
    aln = AlignIO.read(zizo, "clustal")
    #C:/Users/ahmed/Downloads/protein_alignment.aln
    #salah=int(input("which % of conservation you want to extract you seq (ex:100,90,..)?:"))
    #mz = int(float(salah/100)*int(len(aln)))
    mz = int(len(aln))
   
    A = ("A"*mz)
    G = ("G"*mz)
    C = ("C"*mz)
    T = ("T"*mz)
    R = ("R"*mz)
    N = ("N"*mz)
    D = ("D"*mz)
    Q = ("Q"*mz)
    E = ("E"*mz)
    H = ("H"*mz)
    I = ("I"*mz)
    L = ("L"*mz)
    K = ("K"*mz)
    M = ("M"*mz)
    F = ("F"*mz)
    P = ("P"*mz)
    O = ("O"*mz)
    S = ("S"*mz)
    U = ("U"*mz)
    W = ("W"*mz)
    Y = ("Y"*mz)
    V = ("V"*mz)
    astr = ("*"*mz) # stop codon translated into *
    XUN = ("X"*mz) #some papers use X as X unknown
    
    yr = [] 
    for xp in range(aln.get_alignment_length()): 
        if str(aln[:,xp]).upper() == A : 
            yr.append((aln[:,xp][0])) 
        elif str(aln[:,xp]).upper() == T :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == C :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == G :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == R :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == N :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == D :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == Q :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == E :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == H :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == I :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == L :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == M :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == K :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == F :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == P :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == O :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == S :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == U :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == W :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == Y :
            yr.append((aln[:,xp][0]))
        elif str(aln[:,xp]).upper() == V :
            yr.append((aln[:,xp][0]))
        
#kindly remeber that sometimes clustal files come in upper/low case letters so i did str().upper()
    yr = "".join(yr)
    s1 = str(yr) #the conserved bases from your your file
    s2 = str(aln[0].seq) #the first seq in the clustal file
    len1, len2 = len(s1), len(s2)
    ir, jr = 0, -1 #the best solution for longest common substring problem (https://rosettacode.org/wiki/Longest_Common_Substring#Python)
    for i1 in range(len1): #takes around 5 min for lenght 30 kb and number of 2000 sequences
        i2 = s2.find(s1[i1])
        while i2 >= 0:
            j1, j2 = i1, i2
            while j1 < len1 and j2 < len2 and s2[j2] == s1[j1]:
                if j1-i1 >= jr-ir:
                    ir, jr = i1, j1
                j1 += 1; j2 += 1
            i2 = s2.find(s1[i1], i2+1)
    rs = str(s1[ir:jr+1])

    Result = open("longest_conserved_%s_file.fasta"%(zizo), "w")
    Result.write('>Longest_conserved_seq_in_%s\n'%(zizo))
    Result.write(rs) #very useful for PCR primer designing or protein domain search
    Result.close()
    treka = ((len(yr)/len(aln[0])*100))
    #to qc my work,
    #len(yr) == str(aln.column_annotations).count("*")
    #if you get true, you are in the right track!
    print ("%f percent of conserved bases" %(treka))
    print("here you are: longest_conserved_seq_%s_file.fasta"%(zizo))
    
   ####################################################
    ##lets get unconserved basis
    hass =[]
    gogo = []
    for xe in range(aln.get_alignment_length()):
        if (str(aln[:,xe]).upper() != A) and (str(aln[:,xe]).upper() != T) and (str(aln[:,xe]).upper() != C) and (str(aln[:,xe]).upper() != G) and (str(aln[:,xe]).upper() != R) and (str(aln[:,xe]).upper() != N) and (str(aln[:,xe]).upper() != Q) and (str(aln[:,xe]).upper() != H) and (str(aln[:,xe]).upper() != E)  and (str(aln[:,xe]).upper() != I) and (str(aln[:,xe]).upper() != L) and (str(aln[:,xe]).upper() != K) and (str(aln[:,xe]).upper() != M) and (str(aln[:,xe]).upper() != F) and (str(aln[:,xe]).upper() != P) and (str(aln[:,xe]).upper() != O) and (str(aln[:,xe]).upper() != S) and (str(aln[:,xe]).upper() != U) and (str(aln[:,xe]).upper() != W) and (str(aln[:,xe]).upper() != Y) and (str(aln[:,xe]).upper() != D) and (str(aln[:,xe]).upper() != V) and (str(aln[:,xe]).upper() != astr) and (str(aln[:,xe]).upper() != XUN):
            cc = Counter(str(aln[:,xe]))
            cc = cc.most_common()
            hass.append((cc,int((xe)+1))) 
            gogo.append(cc) #this libe to get the pattern of the muatations
   
    yo = []
    ramy = []
    youssef = len(aln[:,xe])
    for elem in hass:
        tarb = float(100*(int(int(elem[0][0][1]))/youssef))
        ramy.append(100 -tarb)
        if len(elem[0]) == 2:
            yo.append("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]))
        if len(elem[0]) == 3:
            yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0])))
        if len(elem[0]) > 3:
            yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][3][0])))

                      
    fawzia = []
    for ghon in hass:
        fawzia.append(ghon[1])
   

    GRG = pd.ExcelWriter("mutations_%s_file.xlsx"%(zizo))
    df = pd.DataFrame({"Position_in_%s"%(zizo):fawzia,"{“WT”:frequency ,“mutation(S)”:frequency}":gogo ,"mutations(xpt if X or > 3 mutation in same pos)":yo,"prev%_of_mutations":ramy})
    df.to_excel(GRG, index = False)
    GRG.save()

    kp =[]
    for xx in hass:#here i want to extract the mutation position like (8,19,200..) from the list of tupules
        kp.append(xx[1])
    kp.pop(0) # to remove the first string in kp list #take care if you repeat this you will lose one element
    lp = list(range(0, len(kp))) #the index of list will represent the mutation number
    print("Done, you have %s mutations/UNconserved bases" %(len(kp)))
    plt.scatter(kp,lp ,color = 'red' , s =.1,lw=4)
    plt.xlabel('length of your genome , gene or protein')
    plt.ylabel('number of unconserved bases/mutations in %s'%(zizo))
    plt.title('distrubtion of mutations/unconserved bases')
    plt.savefig("mutations_map_%s_file.jpg"%(zizo)) #save your file!
    plt.close()
    
    plt.scatter(fawzia,ramy)
    for i, txt in enumerate(yo):
        plt.annotate(txt, (fawzia[i], ramy[i]))
    plt.ylabel('prevelance_of_mutations_%')
    plt.xlabel('length of your genome , gene or protein')
    plt.savefig("mutations_frequency_%s_file.jpg"%(zizo)) #save your file!
   

        
    print("here you are: mutations_file.xlsx,mutations_map and frequency plots ")
    print("#This the end of our journey :( , hope to see you again! ")
    print ("If you encounter any issues using me, feel free to contact Ahmed")
    #the idea here i want to know where my mutations or my unconervead bases positions in geneme
    #do not forget, if you have gap in your align mutation position will be shifted, remove gap from clusta
    ##########################################################################################


