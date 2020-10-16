# -*- coding: utf-8 -*-
"""
@author: ahmed_elsherbini (drahmedsherbini@yahoo.com)
created in : June 2020
This is a copy right for the author - do not distrbute
dependacies: see below
update: 16/10/2020
please cite my page if you used this script
"""
#import
print("Hi, #1_to work make sure you have python 3.xx , biopython , ClustalW, Muscle and MAFFT installed on your PC #2_this is an semi automated tool, just run it and answer questions PRECISELY!")
import os as os
import sys
from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio import AlignIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import MafftCommandline
from Bio import Phylo
from Bio.SeqUtils import GC
from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Blast import NCBIWWW
import matplotlib.pyplot as plt
from Bio.Cluster import pca
import pandas as pd
import numpy as np
from difflib import SequenceMatcher 
import re

###################################################################################
#%% #in SPYDER this help to chunk the code!

f = input ("1_do you want to extract genes using its start/end postion in aligned_file.afa (more accurate) or fasta  file? y/n: ")
if (f =="y"):
        er = input ("what is the name of the file you want to extract?")
        ah = input ("start position:")
        med = input ("end position +1:")
        with open ("filter_by_position.fasta", "w") as f:
             for seq_record in SeqIO.parse(er, "fasta"):
                  f.write(">"+str(seq_record.id)+"\n")
                  f.write(str(seq_record.seq[int(ah)-1:int(med)-1])+ "\n")
        print("here you are the file : filter_by_position.fasta")
########################################################################################
tui = input("2_do you want to exlude sequences in a multifasta file using sequneces pattern (ex:NNN,XX)? y/n:")
if (tui == "y"):
       sherb = input("what is the name of your fasta/multifasta file:")
       homdemha = input("what is pattern you want to exclude:")
       fin = "fasta_filtred_by_Exclusion.fasta"
       out = open(fin,"w")
       for seq_record in SeqIO.parse(sherb, "fasta"):
           if (seq_record.seq.find(homdemha) != -1):
               continue
           else:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq.ungap("-"))+"\n")
           out.close()
           print("here you are the file : fasta_filtred_by_Exclusion.fasta")
        
#########################################################################################
f = input("3-do you want to  print all headers in your multifasta? y/n:")
if (f == "y"):
    sherb = input("what is the name of the multifasta file?")
    fin = "your_file_headers.txt"
    out = open(fin, "w")
    for seq_record in SeqIO.parse(sherb, "fasta"):
        out.write(">"+str(seq_record.id)+"\n")
    out.close()
    print("here you are the file:your_file_headers.txt")
###########################################################################################
f = input("4-do you want to  extract sequences using a pattern in their headers/metadate (ex:-2019-)? y/n:")
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
    print("here you are the file:extracted_by_header.fasta")


##########################################################
lui = input("5-do you want to extract gene from fasta/multifasta file using upstream and downstream seq? y/n:")
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
################################################
a = input("6-do you want to translat multi/DNA fasta file on the 1 frame? y/n:")
if (a == "y"):
    zeze = input("what is the name of DNA file?")
    with open ("translated_file.fasta" , "w") as aa_fa:
        for dna_record in SeqIO.parse(zeze, "fasta"):
            aa_fa.write(">"+dna_record.id+ "\n")
            aa_fa.write(str(dna_record.seq.translate(to_stop=True))+"\n")
        aa_fa.close()
        print("here you are the file:translated_file.fasta")

#########################################################################
#bring your data

f = input("7_if u want NCBI efetch (press y), if you want to merge all files in your Dir (press m), to skip (press any key):")
#if you have alreay a merged file skip and press any key
if (f == "y"):
    records = []
    Entrez.email = input("Enter your email:")
    db = input("Enter which DB (nucleotide or protein) you want:")
    id = input("Enter your list (comma sep) of your sequence access no.:") 
    #NC_004718.3,NC_019843.3,NC_045512.2,NC_038312.1   #mers,sars,covid, rinho virus
    #MT520385,NC_045512,MT520395, MT520399 #strains of covid-19
    #very useful if you have genes or proteins acesss IDs
    with Entrez.efetch(db=db, id=id, rettype='fasta', retmode='text') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            print("The sequence ID is %s and its length is %d." %(seq_record.id,len(seq_record.seq)))
            rec = SeqRecord(seq_record.seq, id=seq_record.id)
            records.append(rec)
            SeqIO.write(records, "fetched_sequences.fasta", format='fasta')
            print("here you are: fetched_sequences.fasta")

    

#here be aware, merging fasta of large files takes much time and merge lines sometimes!
elif (f =="m"): 
    print("Be alarmed, I will merge all txt (fasta) files in your dir to make one fasta file")
    dir = input("where is the direcroy you want to merge all files in?:") # my advice is to use linux seqkit tool if you have genomes or big files
    #C:\Users\ahmed\Downloads\
    oh = open("merged_sequneces.fasta", 'w')
    for f in os.listdir(dir):
        for seq_record in SeqIO.parse(f, "fasta"):
            oh.write(">"+str(seq_record.id)+"\n")
            oh.write(str(seq_record.seq)+"\n")
       
    oh.close()
    print("here you are: merged_sequences.fasta")

 
    
############################################################################################
#%%
#gc conent and At and number of unkown bases (extra work)
u = input("8_do you want to know GC content of your DNA seq? press y/n:")
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
######################################################################################
#%%
#alignmnent

x = input("9-To align, if for muscle press m, for Mafft press f,press any key to skip:")

if (x == "m"):
    file_path_out = input("what is name of the fasta file you would like to align?")
    a = ("muscle_aligned.aln")
     #C:\Users\ahmed\Downloads\merged_file.afa
    m = ("tree.phy")
    #C:\Users\ahmed\Downloads\merged_file.phy
    muscle_cline = MuscleCommandline(input=file_path_out ,out = a , tree1 = m)
    print(muscle_cline) # great advice: you can run this output in you shell/cmd
    stdout, stderr = muscle_cline() #if you have alot/big files you will wait so much , back to great advice
    print("here you are : musvle_aligned.aln and tree.phy ")

elif (x == "f"):
    file_path_out = input("what is name of the fasta file you would like to align?")
    mafft_cline = MafftCommandline(input=file_path_out)
    print(mafft_cline)
    stdout, stderr = mafft_cline() #mafft is super fast 
    #C:\Users\ahmed\Downloads\mafft_aligned.aln
    with open("Mafft_aligned.aln", "w") as handle:
        handle.write(stdout)
        handle.close()
    print("here you are: Mafft_aligned.aln")
#################################################################################
#%%
#phylogentic tree

tr = input("10_do you want to draw a tree.dnd? y/n?")
if (tr == "y"):
    sh = input("To draw a tree, enter your/path , directory/ merged_file.dnd: ")
    #C:/Users/ahmed/Downloads/merged_file.dnd
    tree = Phylo.read(sh, 'newick')
    tree.ladderize()   # Flip branches so deeper clades are displayed at top
    Phylo.draw(tree) 

################################################################################    
#%%
#to extract the longest conserved sequneces  mutations from gene , genome or ptotein 

 
w = input("11_do you want to extract the longest conserved & the mutations from  your clustal_file.aln? press y/n:")
#C:/Users/ahmed/Downloads/merged_file.aln #kindly know that this code is not adapted to clustal files only
if (w == "y"):
    print("make sure you input file does not have any outliers, outgroups")
    xp = input("what is the name of your file.aln:")
    aln = AlignIO.read(xp, "clustal")
    #C:/Users/ahmed/Downloads/protein_alignment.aln
    mz = int(len(aln))
    # 7 , 9 , 8000 #the advatge here is the script is friendely with computional resoursces
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
    ir, jr = 0, -1 #the best solution for longest common subtring problem (https://rosettacode.org/wiki/Longest_Common_Substring#Python)
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

    Result = open("longest_conserved_seq.fasta", "w")
    Result.write('>Longest_conserved_seq\n')
    Result.write(rs) #very useful for PCR primer designing
    Result.close()
    treka = ((len(yr)/len(aln[0])*100))
    #to qc my work,
    #len(yr) == str(aln.column_annotations).count("*")
    #if you get true you are in the right track!
    print ("%f percent of conserved bases" %(treka))
    print("here you are: longest_conserved_seq.fasta ")
    
   ####################################################
    ##lets get unconserved basis
    yo =[(("mutation_pattern"),("position"))]
    for xe in range(aln.get_alignment_length()):
        if (str(aln[:,xe]).upper() != A) and (str(aln[:,xe]).upper() != T) and (str(aln[:,xe]).upper() != C) and (str(aln[:,xe]).upper() != G) and (str(aln[:,xe]).upper() != R) and (str(aln[:,xe]).upper() != N) and (str(aln[:,xe]).upper() != Q) and (str(aln[:,xe]).upper() != H) and (str(aln[:,xe]).upper() != E)  and (str(aln[:,xe]).upper() != I) and (str(aln[:,xe]).upper() != L) and (str(aln[:,xe]).upper() != K) and (str(aln[:,xe]).upper() != M) and (str(aln[:,xe]).upper() != F) and (str(aln[:,xe]).upper() != P) and (str(aln[:,xe]).upper() != O) and (str(aln[:,xe]).upper() != S) and (str(aln[:,xe]).upper() != U) and (str(aln[:,xe]).upper() != W) and (str(aln[:,xe]).upper() != Y) and (str(aln[:,xe]).upper() != D) and (str(aln[:,xe]).upper() != V) and (str(aln[:,xe]).upper() != astr) and (str(aln[:,xe]).upper() != XUN):
            yo.append((set(aln[:,xe]),int(xe)+1)) # this step to give you where really the uncoserved/mutations part is, +1 as python indexing start from 0

    
    GRG = pd.ExcelWriter("mutations_file.xlsx")
    df = pd.DataFrame(yo)
    df.to_excel(GRG, index = False)
    GRG.save()
    
    kp =[]
    for xx in yo:#here i want to extract the mutation position like (8,19,200..) from the list of tupules
        kp.append(xx[1])
    kp.pop(0) # to remove the first string in kp list #take care if you repeat this you will lose one element
    lp = list(range(0, len(kp))) #the index of list will represent the mutation number
    print("Done, you have %s mutations/UNconserved bases" %(len(kp)))
    plt.scatter(kp,lp ,color = 'red' , s =.1,lw=4)
    plt.xlabel('length of genome , gene or protein')
    plt.ylabel('number of unconserved bases/mutations')
    plt.title('distrubtion of mutations/unconserved bases')
    plt.savefig("mutations_graph.jpg") #save your file!
    print("here you are: mutations_file.xlsx,mutations_graph.jpg ")
    #the idea here i want to know where my mutations or my unconervead bases positions in geneme
    #do not forget, if you have gap in your align mutation position will be shifted, remove gap from clusta

    ##########################################################################################
