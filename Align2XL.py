#!/usr/bin/env python3

"""
@uthor: ahmed_elsherbini (drahmedsherbini@yahoo.com)
created in :07/06/2020
last_update:31/05/2021
This is a copy right for the author !
dependacies: see below !
"""
#import

import os as os
import sys
from io import StringIO
from Bio import Phylo
from Bio import Entrez
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import GC
from Bio import AlignIO
from Bio.Align import AlignInfo
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from difflib import SequenceMatcher 
from collections import Counter
import copy
import fnmatch
from Bio.Align.Applications import MafftCommandline
import regex
import re
####################################################################################
#some_instructions
print("*********************************************************************")
print("Align2XL.0.1.:2020-2021")
print("In-silico-PCR,multifasta manipulation and call mutations")
print("Here, my rules, (answer y for working on a single file) , (b for to batch analysis for all files with certain extension in this folder) and (press any other key to skip)")
print("For b/Batch mode, make sure all of your files have the same extension (example: .fasta)")
print("For (Non) aligned files,I accept only Fasta files. For aligned files,I accept Clustal,Fasta,Phylip and Stockholm")
print("This pipeline is crafted by Ahmed M. A. Elsherbini, drahmedsherbini@yahoo.com")
print("*********************************************************************")

###################################################################################

f = input("1-do you want efetch your sequences using their access numbers (Example: NC_004718.3,NC_019843.3) from ncbi? y/any other key to skip:")
#if you have alreay a merged file skip and press any key
if (f == "y"):
    records = []
    Entrez.email = "ahmed@yahoo.com"
    db = input("Enter which DB (nucleotide or protein) you want:")
    id_f = input("Enter your list (comma seperated) of your sequence access numbers(Example:NC_004718.3,NC_019843.3):") 
    #very useful if you have genes or proteins acesss IDs
    with Entrez.efetch(db=db, id=id_f, rettype='fasta', retmode='text') as handle:
        for seq_record in SeqIO.parse(handle, 'fasta'):
            print("The sequence ID is %s and its length is %d." %(seq_record.id,len(seq_record.seq)))
            rec = SeqRecord(seq_record.seq, id=seq_record.id)
            records.append(rec)
    SeqIO.write(records,'new_fetched_file.fasta' ,'fasta')
###################################################################################
#works only on LINUX OS
f = input("2-do you want to add your reference sequence on top of your sequence file(s)?y/b/any other key to skip:")
if (f == "y"):
      inpt_ref = input("What is the name of your reference sequence file?")
      inpt_seq = input("What is the name of your sequence file?")
      cmd1 = "cat %s %s > ref_added_to_%s.fa" %(inpt_ref,inpt_seq,inpt_seq) 
      os.system(cmd1)

elif (f == "b"): #Here the same reference file will be added to all *.xxx files.
      type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
      inpt_ref = input("What is the name of your refernce sequence file?")
      for inpt_seq in os.listdir():
        if fnmatch.fnmatch(inpt_seq,type_exten):
                cmd1 = "cat %s %s > ref_added_to_%s.fa" %(inpt_ref,inpt_seq,inpt_seq) 
                os.system(cmd1)
###################################################################################
f = input ("3-do you want to know the coordinates + extract a gene/protein in (an) aligned file(s) using few first and last few sequences ? y/b/any other key to skip: ")
if (f =="y"):
        inpt_f = input ("What is the name of the aligned file you want to extract from?")
        type_format = input("What is the format of your aligned file (Example:clustal,fasta,phylip,stockholm):")
        fi_seq= list(AlignIO.read(inpt_f,type_format))
        fi_seq = str(fi_seq[0].seq).upper() #let get the top reference sequence here to know where our coordinates are.
        start = input("What is sequence of the first few (like 18) bp/aa of 5'/N-term of the gene/protein ?:").upper()
        end =   input("What is sequence of the last few (like 18) bp/aa  of 3'/C-term  of the gene/protein ?:").upper()
        result_o = int(fi_seq.index(start)+1) #we will search for the index of first b/aa of the gene/protein using the starting few seq. 
        result_t = int(fi_seq.index(end)+len(end)) #we will get the index of the last b/aa of the gene/protein using the last few seq.
        z = input("What is the name of the gene/protein  you want to extract (Example:mgrb, spike, ORF1ab,..):")
        print("Your gene/protein starts at position %s  and end at position %s in %s" %(result_o,result_t,inpt_f))
        with open ("%s_from_%s.fasta"%(z,inpt_f), "w") as f:
             for seq_record in AlignIO.read(inpt_f, type_format):#lets loop on every seq in your file!
                  f.write(">"+str(seq_record.id)+"\n")
                  f.write(str(seq_record.seq[result_o-1:result_t])+ "\n") #the pro of alignment file is the coodinates are the same
        print("Here you are the file")

elif (f == "b"):
        start = input("What is sequence of the first few (like 18) bp/aa of 5'/N-term of the gene/protein ?:").upper()
        end =   input("What is sequence of the last few (like 18) bp/aa  of 3'/C-term  of the gene/protein ?:").upper()
        type_exten = input("What is the (*.extension) of your files (Example : *.fasta , *.afa ,*.phy ,*.stk, *.aln)?:")
        type_format = input("What is the format of your aligned files (Example:clustal,fasta,phylip,stockholm):")
        z = input("What is the name of the gene/protein you want to extract (Example:mgrb, spike, ORF1ab,..):")

        for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                 f = open("%s_from_%s_XFiles.faa" %(z,inpt_f),"w")
                 fi_seq= list(AlignIO.read(inpt_f,type_format))
                 fi_seq = str(fi_seq[0].seq).upper()
                 result_o , result_t = 0,1 #to avoid a broken loop if one we could not find the coodinates (example:mutations in flanking primers)
                 result_o = int(fi_seq.index(start)+1)
                 result_t = int(fi_seq.index(end)+len(end))
                 print("Your gene/protein starts at position %s  and end at position %s in %s" %(result_o,result_t,inpt_f))
                 if (result_o != 0) and (result_t != 1): 
                     for seq_record in AlignIO.read(inpt_f, type_format):
                                  f.write(">"+str(seq_record.id)+"_from_%s"%(inpt_f)+"\n") #this is important if you want to merge your files.
                                  f.write(str(seq_record.seq[result_o-1:result_t])+ "\n")
        f.close()
        ask = input ("Do you want to merge your output files into one file?: y/n")
        if (ask == "y"):
              cmd1 = "cat *_XFiles.faa > %s_from_all_files.fasta" %(z) #this is a strange solution. But,I found it the fastest way to merge output from large sequence files
              cmd2 = "rm -r *_XFiles.faa" # I am very sure that the user has no file with extension _XFiles.faa !!
              os.system(cmd1)
              os.system(cmd2) 
        print("Here you are the files.Please, change the name and extension of the your file")    

###################################################################################

f = input ("4-do you want to extract gene/protein using its manually co-ordinates (start/end) postion in (an) aligned file(s)? y/b/any other key to skip: ")
print("You can extract coordinates using any sequence viewer ex:(CLC seq viewer/Jalview/..)")
if (f =="y"):
        inpt_f = input ("What is the name of the file you want to extract from?")
        ah = input ("start position:")
        med = input ("end position:")
        type_format = input("What is the format of your aligned file (Example:clustal,fasta,phylip,stockholm):")
        z = input("What is the name of the gene/protein you want to extract (Example:mgrb, spike, ORF1ab,..):")
        with open ("%s_from_%s.fasta"%(z,inpt_f), "w") as f:
             for seq_record in AlignIO.read(inpt_f, type_format):
                  f.write(">"+str(seq_record.id)+"\n")
                  f.write(str(seq_record.seq[int(ah)-1:int(med)+1])+ "\n") 
        print("Here you are the file : %s_from_%s.fasta"%(z,inpt_f))

elif (f == "b"): #It is quite rare that user will loop on every file using the same coordinates
    ah = input ("start position:")
    med = input ("end position:")
    type_exten = input("What is the (*.extension) of your files (Example : *.fasta , *.afa ,*.phy ,*.stk, *.aln)?:")
    type_format = input("What is the format of your aligned files (Example:clustal,fasta,phylip,stockholm):")
    z = input("What is the name of the gene/protein  you want to extract (Example:mgrb, spike, ORF1ab,..):")
    for inpt_f in os.listdir():
        if fnmatch.fnmatch(inpt_f,type_exten):
            f = open("%s_from_%s_XFiles.afa" %(z,inpt_f),"w")
            for seq_record in AlignIO.read(inpt_f, type_format):
                 f.write(">"+str(seq_record.id)+"_from_%s"%(inpt_f)+"\n")
                 f.write(str(seq_record.seq[int(ah)-1:int(med)+1])+ "\n")
    f.close()
    ask = input ("Do you want to merge output files your files into one file?y/n:")
    if (ask == "y"):
          cmd1 = "cat *_XFiles.faa > %s_from_all_files.fna" %(z) 
          cmd2 = "rm -r *_XFiles.faa" 
          os.system(cmd1)
          os.system(cmd2)
    elif (ask =="n"):
        print("Kindly,remove XFiles from the name of your files")
        

    print("Here you are the file(s), Please, change the name and extension of the your files")


    
###########################################################################################
f = input("5-do you want to extract a gene/protein from (a) NON-aligned or aligned fasta/Contigs/multifasta file(s) using its first few and last few sequences? y/b/any other key to skip:")
print("This function can work even if want to allow  mismatch in your primers")
if (f == "y"):
       inpt_f = input("What is the name of your fasta/multifasta file:")
       x = input("What is sequence of the first few (like 18) bp/aa of 5' or N-term of the gene/protein ?:").upper()
       y = input("What is sequence of the last few (like 18) bp/aa  of 3' or C-term  of the gene/protein ?:").upper()
       z = input("Give a name for the gene/protein you want to extract (Example:mgrb, spike, ORF1ab,..):")
       out = open("%s_exctracted_genes_%s.fasta" %(z,inpt_f),"w") 
       mis_ma = int(input("what is the value of mismatch you want to allow in your primers (example:0,1,2,3)? "))
       out = open("%s_exctracted_genes_%s.fasta" %(z,inpt_f),"w")
       for seq_record in SeqIO.parse(inpt_f, "fasta"):
           ss = regex.findall('(%s(.+)%s){s<=%d}'%(x,y,mis_ma), str(seq_record.seq.ungap("-")).upper())
           if ss != []:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(ss[0][0])+"\n")
               
       out.close()
       print("Here you are the file")

elif (f == "b"):
     x = input("What is sequence of the first few (like 18) bp/aa of 5' or N-term of the gene/protein ?:").upper()
     y = input("What is sequence of the last few (like 18) bp/aa  of 3' or C-term  of the gene/protein ?:").upper()
     z = input("Give a name for the gene/protein you want to extract (Example:mgrb, spike, ORF1ab,..):")
     type_exten = input("What is the (*.extension) of your files (Example : *.fasta , *.fna , *.afa. *.fa)?:")
     mis_ma = int(input("what is the value of mismatch you want to allow in your primers (example:0,1,2,3)? ")) 

     for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                  out = open("extracted_%s_from_%s_XFiles.faa" %(z,inpt_f),"w")
                  for seq_record in SeqIO.parse(inpt_f, "fasta"):
                      ss = regex.findall('(%s(.+)%s){s<=%d}'%(x,y,mis_ma), str(seq_record.seq.ungap("-")).upper())
                      if ss != []:
                          out.write(">"+str(seq_record.id)+"\n")
                          out.write(str(ss[0][0])+"\n")
                          
                  out.close()
     ask = input ("Do you want to merge your output files into one file?: y/n:")
     if (ask == "y"):
          cmd1 = "cat *_XFiles.faa > %s_from_all_files.fna" %(z) 
          cmd2 = "rm -r *_XFiles.faa" 
          os.system(cmd1)
          os.system(cmd2) 
     elif (ask =="n"):
        print("Kindly,remove XFiles from the name of your files")
        
   
     print("Here you are the file.Please, change the name and extension of the your file")

############################################################################################################
f = input("6-do you want to do Step.4.But, on the reverse complement strand of your DNA file(s) ? y/b/any other key to skip:")
print("This function can work even if want to allow  mismatch in your primers")
if (f == "y"):
       inpt_f = input("What is the name of your fasta/multifasta file:")
       x = input("What is sequence of the first few (like 10) bp of 5' of the gene ?:").upper()
       y = input("What is sequence of the last few (like 10) bp of 3' of the gene ?:").upper()
       mis_ma = int(input("what is the value of mismatch you want to allow in your primers (example:0,1,2,3)? "))
       z = input("Give name of the gene you want to extract (Example:mgrb, spike, ORF1ab,..):")
       out = open("%s_from_reverse_%s.fasta" %(z,inpt_f),"w")
       
       for seq_record in SeqIO.parse(inpt_f, "fasta"):
           seq_record2 = seq_record.reverse_complement()
           ss = regex.findall('(%s(.+)%s){s<=%d}'%(x,y,mis_ma), str(seq_record2.seq.ungap("-")).upper())
           if ss != []:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(ss[0][0])+"\n")
       out.close()
       print("Here you are the file")

elif (f == "b"):
     x = input("What is sequence of the first few (like 10) bp of 5' of the gene ?:").upper()
     y = input("What is sequence of the last few (like 10 ) bp of 3' of the gene ?:").upper()
     z = input("Give a name for the gene you want to extract (Example:mgrb, spike, ORF1ab,..):")
     type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
     mis_ma = int(input("What is the value of mismatch you want to allow in your primers (example:0,1,2,3)? "))
     for inpt_f in os.listdir():
       if fnmatch.fnmatch(inpt_f,type_exten):
          out = open("exctracted_%s_from_%s_rev_XFiles.faa" %(z,inpt_f),"w")
          for seq_record in SeqIO.parse(inpt_f, "fasta"):
              seq_record2 = seq_record.reverse_complement()
              ss = regex.findall('(%s(.+)%s){s<=%d}'%(x,y,mis_ma), str(seq_record2.seq.ungap("-")).upper())
              if ss != []:
                  out.write(">"+str(seq_record.id)+"_from_%s"%(inpt_f)+"\n")
                  out.write(str(ss[0][0])+"\n")
          out.close()
     ask = input ("Do you want to merge your output files into one file? y/n:")
     if (ask == "y"):
        cmd1 = "cat *_XFiles.faa > %s_from_rev_all_files.fna" %(z)
        cmd2 = "rm -r *_XFiles.faa" # I am very sure that the user has no file with extension _XFiles.faa !!
        os.system(cmd1)
        os.system(cmd2)
     elif (ask =="n"):
        print("Kindly,remove XFiles from the name of your files")
        
      

     print("Here you are the file.Please, change the name and extension of the your files")

############################################################################################################
f = input("7-do you want to know GC content % and N bases content % of your DNA (a) OR X amino acids in multifasta/contigs file(s)? press y/b/any other key to skip:")

if (f == "y"):
    ch = input("Is it a nucleotide or protein file? nt/pro:")
    if (ch == "nt"):
        file_path_out = input("What is the name of your file?")
        k = [("ID","GC content%")]
        w = [("ID", "unknown bases%")]
        for seq_record in SeqIO.parse(file_path_out,"fasta"):
            k.append((seq_record.id,GC(seq_record.seq)))
            w.append((seq_record.id,(((float(str(seq_record.seq).count("N" or "n"))/len(seq_record.seq.ungap("-"))))*100))) #it is important to ungap here as the tool 
       
        GFG = pd.ExcelWriter("GC_content_%s.xlsx"%(file_path_out))
        n_bases = pd.ExcelWriter("N_bases_%s.xlsx"%(file_path_out))
        df = pd.DataFrame(k)
        nf = pd.DataFrame(w)
        df.to_excel(GFG, index = False)
        nf.to_excel(n_bases, index = False)
        GFG.save()
        n_bases.save()
        
    elif (ch == "pro"):
         file_path_out = input("What is the name of your file?")
         w = [("ID", "unknown amino acids %")]
         for seq_record in SeqIO.parse(file_path_out,"fasta"):
            w.append((seq_record.id,(((float(str(seq_record.seq).count("X" or "x"))/len(seq_record.seq.ungap("-"))))*100))) #it is important to ungap here as the tool 
         
         GFG = pd.ExcelWriter("X_AA_%s.xlsx"%(file_path_out))
         nf = pd.DataFrame(w)
         nf.to_excel(GFG, index = False)
         GFG.save()
        
elif (f == "b"): 
    type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
    ch = input("Are they nucleotide or protein files? nt/pro:")
    for file_path_out in os.listdir():
            if fnmatch.fnmatch(file_path_out,type_exten):
                 if (ch == "nt"):
                        k = [("ID","GC content%")]
                        w = [("ID", "unknown bases%")]
                        for seq_record in SeqIO.parse(file_path_out,"fasta"):
                            k.append((seq_record.id,GC(seq_record.seq)))
                            w.append((seq_record.id,(((float(str(seq_record.seq).count("N" or "n"))/len(seq_record.seq.ungap("-"))))*100))) #it is important to ungap here as the tool 
                       
                        GFG = pd.ExcelWriter("GC_content_%s.xlsx"%(file_path_out))
                        n_bases = pd.ExcelWriter("N_bases_%s.xlsx"%(file_path_out))
                        df = pd.DataFrame(k)
                        nf = pd.DataFrame(w)
                        df.to_excel(GFG, index = False)
                        nf.to_excel(n_bases, index = False)
                        GFG.save()
                        n_bases.save()  
                        
                 elif (ch == "pro"):
                     w = [("ID", "unknown amino acids %")]
                     for seq_record in SeqIO.parse(file_path_out,"fasta"):
                        w.append((seq_record.id,(((float(str(seq_record.seq).count("X" or "x"))/len(seq_record.seq.ungap("-"))))*100))) #it is important to ungap here as the tool 
                   
                     GFG = pd.ExcelWriter("N_bases_%s.xlsx"%(file_path_out))
                     nf = pd.DataFrame(w)
                     nf.to_excel(GFG, index = False)
                     GFG.save()
    print("Here you are the files")
#################################################################################################################
f = input("8-do you want to extract sequences from (a) multifasta file(s) by exclusion of certain percentage of ambiguous bases/aminoacids (ex:N,X)? y/b/any other key to skip:")

if (f == "y"):
       ch = input ("Is it a nucleotide or protein sequences? nt/pro:")
       inpt_f = input("What is the name of your fasta/multifasta file:")
       user_n = input("what is the cut-off percentage you want to exclude sequences above (1,2,10)?:")
       if (ch == "nt"):
           out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
           for seq_record in SeqIO.parse(inpt_f, "fasta"):
               per_n = float(((str(seq_record.seq).count("N" or "n"))/len(seq_record.seq.ungap("-")))*100)
               if (per_n >= float(user_n)):
                   continue
               else:
                   out.write(">"+str(seq_record.id)+"\n")
                   out.write(str(seq_record.seq)+"\n")
    
           out.close()
       elif (ch == "pro"):
           out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
           for seq_record in SeqIO.parse(inpt_f, "fasta"):
               per_n = float(((str(seq_record.seq).count("X" or "x"))/len(seq_record.seq.ungap("-")))*100)
               if (per_n >= float(user_n)):
                   continue
               else:
                   out.write(">"+str(seq_record.id)+"\n")
                   out.write(str(seq_record.seq)+"\n")
                   
elif (f == "b"):
    ch = input ("Is it a nucleotide or protein sequences? nt/pro:")
    type_exten = input("What is the (*.extension) of your files (Example : *.fasta , *.fna , *.afa. *.fa)?:")
    user_n = input("what is the cut-off percentage you want to exclude sequences above (1,2,10)?:")
    for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                 if (ch == "nt"):
                   out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
                   for seq_record in SeqIO.parse(inpt_f, "fasta"):
                       per_n = float(((str(seq_record.seq).count("N" or "n"))/len(seq_record.seq.ungap("-")))*100)
                       if (per_n >= float(user_n)):
                           continue
                       else:
                           out.write(">"+str(seq_record.id)+"\n")
                           out.write(str(seq_record.seq)+"\n")
            
                   print("here you are the file : fasta_filtred_by_seq_Exclusion.fasta")
                   out.close()
                 elif (ch == "pro"):
                   out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
                   for seq_record in SeqIO.parse(inpt_f, "fasta"):
                       per_n = float(((str(seq_record.seq).count("X" or "x"))/len(seq_record.seq.ungap("-")))*100)
                       if (per_n >= float(user_n)):
                           continue
                       else:
                           out.write(">"+str(seq_record.id)+"\n")
                           out.write(str(seq_record.seq)+"\n")
                  
    print("here you are the files")                           
    
    
#########################################################################################
f = input("9-do you want to extract sequences from (a) multifasta file(s) by inclusion of certain sequence pattern (ex: genes with certain mutations)? y/b/any other key to skip:")
if (f == "y"):
       inpt_f = input("What is the name of your fasta/multifasta file:")
       inpt_patt = input("What is pattern you want use it to extract (Example:ATTGCGTGTGTGT or ANNVLKOPGTLSTTSG):").upper()
       out = open("fasta_filtred_by_inculsion.fasta","w")
       for seq_record in SeqIO.parse(inpt_f, "fasta"):
           if (seq_record.seq.find(inpt_patt) != -1):
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq)+"\n")
     
       print("here you are the file : fasta_filtred_by_seq_inclusion.fasta")
       out.close()
       
elif (f == "b"):
   inpt_patt = input("What is pattern you want use it to extract (Example:ATTGCGTGTGTGT or ANNVLKOPGTLSTTSG):").upper()
   type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
   for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                   out = open("fasta_filtred_by_inculsion%s.fasta" %(inpt_f),"w")
                   for seq_record in SeqIO.parse(inpt_f, "fasta"):
                       if (seq_record.seq.find(inpt_patt) != -1):
                           out.write(">"+str(seq_record.id)+"\n")
                           out.write(str(seq_record.seq)+"\n")
                 
                   out.close()
   print("here you are the file : fasta_filtred_by_seq_inclusion.fasta")




########################################################################################
f = input("10-do you want to extract sequences from (a) multifasta file(s) by exclusion of certain sequence pattern (ex:N,X)? y/b/any other key to skip:")
if (f == "y"):
       inpt_f = input("What is the name of your fasta/multifasta file:")
       inpt_patt = input("What is sequence pattern you want to exclude from your file:").upper()
       out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
       for seq_record in SeqIO.parse(inpt_f, "fasta"):
           if (seq_record.seq.find(inpt_patt) != -1):
               continue
           else:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq)+"\n")

       print("here you are the file : fasta_filtred_by_seq_Exclusion.fasta")
       out.close()
       
elif (f == "b"):
    inpt_patt = input("What is sequence pattern you want to exclude from you files:").upper()
    type_exten = input("What is the (*.extension) of your files (Example : *.fasta , *.fna , *.afa. *.fa)?:")
    for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                   out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
                   for seq_record in SeqIO.parse(inpt_f, "fasta"):
                       if (seq_record.seq.find(inpt_patt) != -1):
                           continue
                       else:
                           out.write(">"+str(seq_record.id)+"\n")
                           out.write(str(seq_record.seq)+"\n")
            
                   out.close()
    print("here you are the file : fasta_filtred_by_seq_Exclusion.fasta")

           
#########################################################################################
f = input("11-do you want to  print all >_ID_headers in your multifasta file(s)? y/b/any other key to skip:")
if (f == "y"):
    inpt_f = input("What is the name of the multifasta file?")
    out = open("your_file_headers_%s.txt" %(inpt_f), "w")
    for seq_record in SeqIO.parse(inpt_f, "fasta"):
        out.write(">"+str(seq_record.id)+"\n")
    out.close()
    print("here you are the file:your_file_headers.txt")
    
elif (f == "b"):
     type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
     for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                out = open("your_file_headers_%s.txt" %(inpt_f), "w")
                for seq_record in SeqIO.parse(inpt_f, "fasta"):
                    out.write(">"+str(seq_record.id)+"\n")
                out.close()
     print("here you are the file:your_file_headers.txt")
                
###########################################################################################
f = input("12-do you want to extract sequences from your multifasta file(s) by inclusion of a certain pattern in their >_ID_header ex(2019)? y/b/any other key to skip:")
if (f == "y"):
    inpt_f = input("What is the name of the fasta file?")
    inpt_patt = input("What is the pattern in seq header you want extract? (Take-care:CASE senstive)")
    out = open("extracted_by_header_%s.fasta" %(inpt_f), "w")
    for seq_record in SeqIO.parse(inpt_f, "fasta"):
        if (seq_record.id.find(inpt_patt) != -1):
            out.write(">"+str(seq_record.id)+"\n")
            out.write(str(seq_record.seq)+"\n")
    out.close()
    print("Here you are the file:extracted_by_header_inclusion.fasta")

elif (f == "b"):
     inpt_patt = input("What is the pattern in seq header you want extract? (Take-care:CASE senstive)")
     type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
     for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                out = open("extracted_by_header_%s.fasta" %(inpt_f), "w")
                for seq_record in SeqIO.parse(inpt_f, "fasta"):
                    if (seq_record.id.find(inpt_patt) != -1):
                        out.write(">"+str(seq_record.id)+"\n")
                        out.write(str(seq_record.seq)+"\n")
                out.close()
     print("Here you are the file:extracted_by_header_inclusion.fasta")


########################################################################################
f = input("13-do you want to extract sequences from (a) multifasta file(s) by exclusion of a certain pattern in their >_ID_header ex(Italy)? y/b/any other key to skip:")
if (f == "y"):
       inpt_f = input("What is the name of your fasta/multifasta file:")
       inpt_patt = input("What is ID pattern you want to exclude:")
       out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
       for seq_record in SeqIO.parse(inpt_f, "fasta"):
           if (seq_record.id.find(inpt_patt) != -1):
               continue
           else:
               out.write(">"+str(seq_record.id)+"\n")
               out.write(str(seq_record.seq)+"\n")

       print("here you are the file : fasta_filtred_by_header_Exclusion.fasta")
       out.close()
       
elif (f == "b"):
     inpt_patt = input("What is ID pattern you want to exclude:")
     type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
     for inpt_f in os.listdir():
            if fnmatch.fnmatch(inpt_f,type_exten):
                   out = open("fasta_filtred_by_exclusion_%s.fasta" %(inpt_f),"w")
                   for seq_record in SeqIO.parse(inpt_f, "fasta"):
                       if (seq_record.id.find(inpt_patt) != -1):
                           continue
                       else:
                           out.write(">"+str(seq_record.id)+"\n")
                           out.write(str(seq_record.seq)+"\n")
            
                   out.close()
     print("here you are the files : fasta_filtred_by_header_Exclusion.fasta")
    

    
    #I have provided the output in a list of tupules but as you can convert easily to dic 
##################################################################################
f = input("14-do you want to translate (a) DNA multifasta file(s) on its (theirs) 1st reading frame ? y/b/any other key to skip :")
if (f == "y"):
    inp_file = input("What is the name of DNA file?")
    with open ("translated_%s_file.fasta"%(inp_file) , "w") as aa_fa:
        for dna_record in SeqIO.parse(inp_file, "fasta"):
            aa_fa.write(">"+dna_record.id+ "\n")
            aa_fa.write(str(dna_record.seq.ungap("-").translate(to_stop=True))+"\n")
        aa_fa.close()
        print("Here you are the file:translated_%s_file.fasta"%(inp_file))

elif (f =="b"): 
    type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
    for f in os.listdir():
        if fnmatch.fnmatch(f,type_exten):
            with open ("translated_%s_file.fasta"%(f) , "w") as aa_fa:
                for dna_record in SeqIO.parse(f, "fasta"):
                    aa_fa.write(">"+dna_record.id+ "\n")
                    aa_fa.write(str(dna_record.seq.ungap("-").translate(to_stop=True))+"\n")
                aa_fa.close()
    print("here you are the file:translated_files.fasta")


#################################################################################
f = input("15-do you want to align using MAFFT? y/b/any other key to skip:")

if (f == "y"): 
    a = input("What is name of your multifasta file?")
    print("running....")
    mafft_cline = MafftCommandline(input=a)
    print(mafft_cline)
    stdout, stderr = mafft_cline() #mafft is super fast 
    with open('mafft_aligned_%s.afa'%(a), "w") as handle:
        handle.write(stdout)
        handle.close()
    print("Done,check for your aligned.aln in your dir")
        
elif (f == "b"):
    type_exten = input("What is the (*.extension) of your files (Example: *.fasta , *.fna , *.afa. *.fa)?:")
    for inp_f in os.listdir():
        if fnmatch.fnmatch(inp_f,type_exten):
                        mafft_cline = MafftCommandline(input=inp_f)
                        print(mafft_cline)
                        print("running....")
                        stdout, stderr = mafft_cline() 
                        with open('mafft_aligned_%s.afa'%(inp_f), "w") as handle:
                            handle.write(stdout)
                            handle.close()
    print("Done,check for your aligned.aln in your dir")        
        
######################################################################################
f  = input("17-do you want to know if there are any INDEL(s) in your alignment file(s) y/b/any other key to skip:" )

if (f == "y"):
    print("Here, per align file, you will get 3 files -> 2 regarding insertion + 1 regarding deletions")
    inp_f = input("What is the name of your file:")
    type_format = input("What is the format of your aligned file (Example:clustal,fasta,phylip,stockholm):")
    aln = AlignIO.read(inp_f,type_format) 
    get_ins = []
    ins_p = []
    get_del = [] 
    substring = "-"
    #str(aln[0].seq) #first raw in alignment file
    #str(aln[:,0]) #first column in alignment
    
    q_ref = input("did you have/put your reference seq as the first/top  the alignment? y/n:")
    if q_ref == "y":
     
        for xp in range(aln.get_alignment_length()):
                col_len= len(str(aln[:,xp])) #len of the align column
                col_str = str(aln[:,xp]).upper()  #seq of the align column
                ins_len = col_str.count(substring) #count of -
                cc = list(Counter(str(aln[:,xp]))) #ordered count of columns with -

                if substring == str(aln[0].seq[xp]): 
                    freq_i =int(len(aln))-int(ins_len)
                    prev_i = 100 * ((freq_i/len(aln)))
                    get_ins.append((int(xp)+1,"%s_ins"%(cc[1]),freq_i,prev_i))
                    ins_p.append(int(xp)+1)
                    
                elif substring != str(aln[0].seq[xp]) and re.findall(substring, col_str) :
                    prev_d = 100 * (int(ins_len)/len(aln))
                    get_del.append((int(xp)+1,"%s_del"%(cc[0][0]),int(ins_len),prev_d))
                    
        GRG_3 = pd.ExcelWriter("insertion_%s.xlsx"%(inp_f))
        df_3 = pd.DataFrame(get_ins, columns=["position", "inserted_base/aa","frequency","prev_%"])
        df_3.to_excel(GRG_3, index = False)
        print("Using %s as a reference seq in %s file...."%(aln[0].id,inp_f))
        print("You have %s INSERTED bases/aa in %s!" %(int(len(df_3)),inp_f)) 
        GRG_3.save()
        plt.scatter(df_3["position"],df_3["prev_%"],color = 'green' , s =.1,lw=4)
        plt.ylabel('prevelance_of_insertion_%')
        plt.ylim(0,100)
        plt.xlabel('length of the protein')
        plt.xlim(0,len(aln[0].seq))
        plt.title('insertion_in_%s'%(inp_f))
        plt.savefig("prev_of_insertion_%s_file.jpg"%(inp_f)) 
        plt.close()
        
        GRG_4 = pd.ExcelWriter("deletion_%s.xlsx"%(inp_f))
        df_4 = pd.DataFrame(get_del, columns=["position", "deleted_base/aa","frequency","prev_%"])
        df_4.to_excel(GRG_4, index = False)
        print("You have %s DELETED bases/aa in %s !"%(int(len(df_4)),inp_f))
        GRG_4.save()
        plt.scatter(df_4["position"],df_4["prev_%"],color = 'blue' , s =.1,lw=4)
        plt.ylabel('prevelance_of_deletion_%')
        plt.ylim(0,100)
        plt.xlabel('length of the protein')
        plt.xlim(0,len(aln[0].seq))
        plt.title('Deletion_in_%s'%(inp_f))
        plt.savefig("prev_of_deletions_%s_file.jpg"%(inp_f)) 
        plt.close()
        
        seq_id = []
        for y in aln:
            for x in ins_p:
                x= int(x)
                if str(y.seq)[x-1] != substring:
                    seq_id.append(y.id)
        df_99 = pd.DataFrame({">seq_with_insertions":seq_id})
        df_99 = df_99.groupby([">seq_with_insertions"]).sum()
        df_99.to_excel("insertion_per_SEQ_%s.xlsx"%(inp_f))


    if q_ref == "n":
    
        for xp in range(aln.get_alignment_length()):
                col_len= len(str(aln[:,xp]))
                col_str = str(aln[:,xp]).upper() 
                ins_len = col_str.count(substring)
                cc = list(Counter(str(aln[:,xp])).most_common())
                
                if (int(ins_len) >= int(cc[0][1])) and re.findall(substring, col_str):
                    freq_i = len(aln)-int(ins_len)
                    prev_i = 100 * ((freq_i/len(aln)))
                    get_ins.append((int(xp)+1,"%s_ins"%(cc[1][0]),freq_i,prev_i))
                    ins_p.append(int(xp)+1)

                elif (int(ins_len) < int(cc[0][1])) and re.findall(substring, col_str):
                    prev_d = 100 * (int(ins_len)/len(aln))
                    get_del.append((int(xp)+1,"%s_del"%(cc[0][0]),int(ins_len),prev_d))

        GRG_3 = pd.ExcelWriter("insertion_%s.xlsx"%(inp_f))
        df_3 = pd.DataFrame(get_ins, columns=["position", "inserted_base/aa","frequency","prev_%"])
        df_3.to_excel(GRG_3, index = False)
        print("You have %s INSERTED bases/aa in %s!" %(int(len(df_3)),inp_f)) 
        GRG_3.save()
        plt.scatter(df_3["position"],df_3["prev_%"],color = 'green' , s =.1,lw=4)
        plt.ylabel('prevelance_of_insertion_%')
        plt.ylim(0,100)
        plt.xlabel('length of the protein')
        plt.xlim(0,len(aln[0].seq))
        plt.title('Insertions_in_%s'%(inp_f))
        plt.savefig("prev_of_insertion_%s_file.jpg"%(inp_f)) 
        plt.close()

        GRG_4 = pd.ExcelWriter("deletion_%s.xlsx"%(inp_f))
        df_4 = pd.DataFrame(get_del, columns=["position", "deleted_base/aa","frequency","prev_%"])
        df_4.to_excel(GRG_4, index = False)
        print("You have %s DELETED bases/aa in %s !"%(int(len(df_4)),inp_f))
        GRG_4.save()
        plt.scatter(df_4["position"],df_4["prev_%"],color = 'blue' , s =.1,lw=4)
        plt.ylabel('prevelance_of_deletion_%')
        plt.ylim(0,100)
        plt.xlabel('length of the protein')
        plt.xlim(0,len(aln[0].seq))
        plt.title('Deletion_in_%s'%(inp_f))
        plt.savefig("prev_of_deletions_%s_file.jpg"%(inp_f)) 
        plt.close()
        
        seq_id = []
        for y in aln:
            for x in ins_p:
                x= int(x)
                if str(y.seq)[x-1] != substring:
                    seq_id.append(y.id)
        df_99 = pd.DataFrame({">seq_with_insertions":seq_id})
        df_99 = df_99.groupby([">seq_with_insertions"]).sum()
        df_99.to_excel("insertion_per_SEQ_%s.xlsx"%(inp_f))

if (f == "b"):
    type_format = input("What is the format of your aligned files (Example:clustal, fasta ,phylip,stockholm):")
    type_exten = input("What is the (*.extension) of your files (Example: *.aln, *.afa ,*.fa,*.fasta,*.fna, *.phy,*.sth)?:")
    q_ref = input("did you have/put your reference seq as the first/top  the alignment? y/n:")
    
    for inp_f in os.listdir():
        if fnmatch.fnmatch(inp_f, type_exten):
                    aln = AlignIO.read(inp_f,type_format) 
                    get_ins = []
                    ins_p = [] 
                    get_del = [] 
                    substring = "-"
                    #str(aln[0].seq) #first raw in alignment file
                    #str(aln[:,0]) #first column in alignment
                    

                    if q_ref == "y":
                     
                        for xp in range(aln.get_alignment_length()):
                                col_len= len(str(aln[:,xp])) #len of the align column
                                col_str = str(aln[:,xp]).upper()  #seq of the align column
                                ins_len = col_str.count(substring) #count of -
                                cc = list(Counter(str(aln[:,xp]))) #ordered count of columns with -
                
                                if substring == str(aln[0].seq[xp]): 
                                    freq_i =int(len(aln))-int(ins_len)
                                    prev_i = 100 * ((freq_i/len(aln)))
                                    get_ins.append((int(xp)+1,"%s_ins"%(cc[1]),freq_i,prev_i))
                                    ins_p.append(int(xp)+1)
                                    
                                elif substring != str(aln[0].seq[xp]) and re.findall(substring, col_str) :
                                    prev_d = 100 * (int(ins_len)/len(aln))
                                    get_del.append((int(xp)+1,"%s_del"%(cc[0][0]),int(ins_len),prev_d))
                                    
                        GRG_3 = pd.ExcelWriter("insertion_%s.xlsx"%(inp_f))
                        df_3 = pd.DataFrame(get_ins, columns=["position", "inserted_base/aa","frequency","prev_%"])
                        df_3.to_excel(GRG_3, index = False)
                        print("Using %s as a reference seq in %s file...."%(aln[0].id,inp_f))
                        print("You have %s INSERTED bases/aa in %s!" %(int(len(df_3)),inp_f)) 
                        GRG_3.save()
                        plt.scatter(df_3["position"],df_3["prev_%"],color = 'green' , s =.1,lw=4)
                        plt.ylabel('prevelance_of_insertion_%')
                        plt.ylim(0,100)
                        plt.xlabel('length of the protein')
                        plt.xlim(0,len(aln[0].seq))
                        plt.title('insertion_in_%s'%(inp_f))
                        plt.savefig("prev_of_insertion_%s_file.jpg"%(inp_f)) 
                        plt.close()
                        
                        GRG_4 = pd.ExcelWriter("deletion_%s.xlsx"%(inp_f))
                        df_4 = pd.DataFrame(get_del, columns=["position", "deleted_base/aa","frequency","prev_%"])
                        df_4.to_excel(GRG_4, index = False)
                        print("You have %s DELETED bases/aa in %s !"%(int(len(df_4)),inp_f))
                        GRG_4.save()
                        plt.scatter(df_4["position"],df_4["prev_%"],color = 'blue' , s =.1,lw=4)
                        plt.ylabel('prevelance_of_deletion_%')
                        plt.ylim(0,100)
                        plt.xlabel('length of the protein')
                        plt.xlim(0,len(aln[0].seq))
                        plt.title('Deletion_in_%s'%(inp_f))
                        plt.savefig("prev_of_deletions_%s_file.jpg"%(inp_f)) 
                        plt.close()
                        
                        seq_id = []
                        for y in aln:
                            for x in ins_p:
                                x= int(x)
                                if str(y.seq)[x-1] != substring:
                                    seq_id.append(y.id)
                        df_99 = pd.DataFrame({">seq_with_insertions":seq_id})
                        df_99 = df_99.groupby([">seq_with_insertions"]).sum()
                        df_99.to_excel("insertion_per_SEQ_%s.xlsx"%(inp_f))
                
                
                    if q_ref == "n":
                    
                        for xp in range(aln.get_alignment_length()):
                                col_len= len(str(aln[:,xp]))
                                col_str = str(aln[:,xp]).upper() 
                                ins_len = col_str.count(substring)
                                cc = list(Counter(str(aln[:,xp])).most_common())
                                
                                if (int(ins_len) >= int(cc[0][1])) and re.findall(substring, col_str):
                                    freq_i = len(aln)-int(ins_len)
                                    prev_i = 100 * ((freq_i/len(aln)))
                                    get_ins.append((int(xp)+1,"%s_ins"%(cc[1][0]),freq_i,prev_i))
                                    ins_p.append(int(xp)+1)

                                elif (int(ins_len) < int(cc[0][1])) and re.findall(substring, col_str):
                                    prev_d = 100 * (int(ins_len)/len(aln))
                                    get_del.append((int(xp)+1,"%s_del"%(cc[0][0]),int(ins_len),prev_d))
                
                        GRG_3 = pd.ExcelWriter("insertion_%s.xlsx"%(inp_f))
                        df_3 = pd.DataFrame(get_ins, columns=["position", "inserted_base/aa","frequency","prev_%"])
                        df_3.to_excel(GRG_3, index = False)
                        print("You have %s INSERTED bases/aa in %s!" %(int(len(df_3)),inp_f)) 
                        GRG_3.save()
                        plt.scatter(df_3["position"],df_3["prev_%"],color = 'green' , s =.1,lw=4)
                        plt.ylabel('prevelance_of_insertion_%')
                        plt.ylim(0,100)
                        plt.xlabel('length of the protein')
                        plt.xlim(0,len(aln[0].seq))
                        plt.title('Insertions_in_%s'%(inp_f))
                        plt.savefig("prev_of_insertion_%s_file.jpg"%(inp_f)) 
                        plt.close()
                
                        GRG_4 = pd.ExcelWriter("deletion_%s.xlsx"%(inp_f))
                        df_4 = pd.DataFrame(get_del, columns=["position", "deleted_base/aa","frequency","prev_%"])
                        df_4.to_excel(GRG_4, index = False)
                        print("You have %s DELETED bases/aa in %s !"%(int(len(df_4)),inp_f))
                        GRG_4.save()
                        plt.scatter(df_4["position"],df_4["prev_%"],color = 'blue' , s =.1,lw=4)
                        plt.ylabel('prevelance_of_deletion_%')
                        plt.ylim(0,100)
                        plt.xlabel('length of the protein')
                        plt.xlim(0,len(aln[0].seq))
                        plt.title('Deletion_in_%s'%(inp_f))
                        plt.savefig("prev_of_deletions_%s_file.jpg"%(inp_f)) 
                        plt.close()
                        for y in aln:
                            for x in ins_p:
                                x= int(x)
                                if str(y.seq)[x-1] != substring:
                                    seq_id.append(y.id)
                        df_99 = pd.DataFrame({">seq_with_insertions":seq_id})
                        df_99 = df_99.groupby([">seq_with_insertions"]).sum()
                        df_99.to_excel("insertion_per_SEQ_%s.xlsx"%(inp_f))

#####################################################################################
f =  input("18-do you want to remove any insertion in your data so position prediction in next steps will be precise? y/any other key to skip:")

if f == "y":
    print("Take care this step accepts only FASTA file with no dulpicated ID headers")
    q_inst = input("Did you install CIAlign tool before? y/n:") #This is a new tool for removing insertions from alignment files
    if q_inst == "n":
        os.system("pip3 install cialign")
    inpt_f  = input("What is the name of your input file?")
    rem_ins = "CIAlign --infile %s --remove_insertions  --insertion_min_size 1" %(inpt_f)
    os.system(rem_ins)
    print("Find CIAlign_cleaned.fasta in your folder and chnage the name of the file if you want")
    print("For CIAalign, cite (Tumescheit, Firth & Brown, 2020) (https://doi.org/10.1101/2020.09.14.291484)")
    #no batch mode here as the CIAlign export the same output name for any input file 
################################################################################    
#%%
 
f = input("19-do you want to extract the possible longest conserved sequence & call variants in each position inside your aligned file(s)? press y/b/any other key to skip:")

if (f == "y"):
    print("Make sure your aligened file does not have any insertions.")
    print("if you have a reference, make sure it is the first sequence on top of your alignment.If you do not have ref, we can still do the job ")
    inp_f = input("What is the name of your file:")
    type_format = input("What is the format of your aligned file (Example:clustal,fasta,phylip,stockholm):")
    aln = AlignIO.read(inp_f,type_format)
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
        
    yr = "".join(yr) 
    s1 = str(yr).upper() 
    s2 = str(aln[0].seq).upper() #the first seq in the align. file #think that I compare collected starts*****, with the first seq.
    len1, len2 = len(s1), len(s2)
    ir, jr = 0, -1 #the best solution for longest common substring problem (https://rosettacode.org/wiki/Longest_Common_Substring#Python)
    for i1 in range(len1): #takes around 5 min for lenght 30 kb and number of 2000 sequences
        i2 = s2.find(s1[i1]) #with my own adaptions
        while i2 >= 0:
            j1, j2 = i1, i2
            while j1 < len1 and j2 < len2 and s2[j2] == s1[j1]:
                if j1-i1 >= jr-ir:
                    ir, jr = i1, j1
                j1 += 1; j2 += 1
            i2 = s2.find(s1[i1], i2+1)
    rs = str(s1[ir:jr+1])
    
    
    Result = open("longest_conserved_%s_file.fasta"%(inp_f), "w")
    Result.write('>the_possible_longest_conserved_seq_in_%s\n'%(inp_f))
    Result.write(rs) #very useful for protein domain search
    Result.close()
    cons_p = ((len(yr)/len(aln[0])*100))
    #to qc my work,
    #len(yr) == str(aln.column_annotations).count("*")
    #if you get true, you are in the right track!
    print("here you are: longest_conserved_seq_%s_file.fasta"%(inp_f))
    
   ####################################################
    ##lets get unconserved basis
    q_ref = input("did you have your reference seq as the first on top of the alignment? y/n:")
   
    if (q_ref == "y"):
        
            hass =[]
            mutat_patt = []
            refre = []
            for xe in range(aln.get_alignment_length()):
                if (str(aln[:,xe]).upper() != A) and (str(aln[:,xe]).upper() != T) and (str(aln[:,xe]).upper() != C) and (str(aln[:,xe]).upper() != G) and (str(aln[:,xe]).upper() != R) and (str(aln[:,xe]).upper() != N) and (str(aln[:,xe]).upper() != Q) and (str(aln[:,xe]).upper() != H) and (str(aln[:,xe]).upper() != E)  and (str(aln[:,xe]).upper() != I) and (str(aln[:,xe]).upper() != L) and (str(aln[:,xe]).upper() != K) and (str(aln[:,xe]).upper() != M) and (str(aln[:,xe]).upper() != F) and (str(aln[:,xe]).upper() != P) and (str(aln[:,xe]).upper() != O) and (str(aln[:,xe]).upper() != S) and (str(aln[:,xe]).upper() != U) and (str(aln[:,xe]).upper() != W) and (str(aln[:,xe]).upper() != Y) and (str(aln[:,xe]).upper() != D) and (str(aln[:,xe]).upper() != V) and (str(aln[:,xe]).upper() != astr) and (str(aln[:,xe]).upper() != XUN):
                    cc = Counter(str(aln[:,xe].upper()))
                    cc = list(Counter(cc).items())
                    hass.append((cc,int((xe)+1))) 
                    mutat_patt.append(cc) #this line to get the pattern of the muatations
                    refre.append(str(aln[0].seq[xe].upper()))
           
            yo = []
            prev_list = []
            total_len = len(aln[:,xe])
            for elem in hass:
                tarb = float(100*(int(int(elem[0][0][1]))/total_len))
                prev_list.append(100 -tarb)
                if len(elem[0]) == 2:
                    yo.append("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]))
                if len(elem[0]) == 3:
                    yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0])))
                if len(elem[0]) > 3:
                    yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][3][0])))
        
                              
            pos_list = []
            for ghon in hass:
                pos_list.append(ghon[1])
            
           
        
            GRG = pd.ExcelWriter("mutations_%s_file.xlsx"%(inp_f))
            df = pd.DataFrame({"Position_in_%s"%(inp_f):pos_list,"{“WT”:frequency ,“mutation(S)”:frequency}":mutat_patt ,"reference":refre,"Top mutations":yo,"prev%_of_mutations":prev_list})
            df.to_excel(GRG, index = False)
            GRG.save()
        
            kp =[]
            for xx in hass:#here i want to extract the mutation position like (8,19,200..) from the list of tupules
                kp.append(xx[1])
            kp.pop(0) # to remove the first string in kp list 
            lp = list(range(0, len(kp))) #the index of list will represent the mutation number
            print("Using %s as a reference seq...."%(aln[0].id))
            print("Done, you have %s mutations" %(len(kp)))
            plt.scatter(kp,lp ,color = 'red' , s =.1,lw=4)
            plt.xlabel('length of the protein')
            plt.xlim(0,len(aln[0].seq))
            plt.ylabel('number of mutations')
            plt.title('Mutations in your query seq (red dot = 1 mutation)')
            plt.savefig("mutations_map_%s_file.jpg"%(inp_f)) #save your file!
            plt.close()
            
            plt.scatter(pos_list,prev_list)
            for i, txt in enumerate(yo):
                plt.annotate(txt, (pos_list[i], prev_list[i]))
            plt.ylabel('prevelance_of_mutations_%')
            plt.xlabel('length of the protein')
            plt.xlim(0,len(aln[0].seq))
            plt.savefig("mutations_per_position_%s_file.jpg"%(inp_f)) #save your file!
            plt.close()
            print("here you are: mutations_file.xlsx,mutations_map and frequency plots")
           
          
            ##########################################################################################
    elif (q_ref =="n"):
        print("Then, I will consider the most common (base/aa) as the WT")
        hass =[]
        mutat_patt = []
        for xe in range(aln.get_alignment_length()):
            if (str(aln[:,xe]).upper() != A) and (str(aln[:,xe]).upper() != T) and (str(aln[:,xe]).upper() != C) and (str(aln[:,xe]).upper() != G) and (str(aln[:,xe]).upper() != R) and (str(aln[:,xe]).upper() != N) and (str(aln[:,xe]).upper() != Q) and (str(aln[:,xe]).upper() != H) and (str(aln[:,xe]).upper() != E)  and (str(aln[:,xe]).upper() != I) and (str(aln[:,xe]).upper() != L) and (str(aln[:,xe]).upper() != K) and (str(aln[:,xe]).upper() != M) and (str(aln[:,xe]).upper() != F) and (str(aln[:,xe]).upper() != P) and (str(aln[:,xe]).upper() != O) and (str(aln[:,xe]).upper() != S) and (str(aln[:,xe]).upper() != U) and (str(aln[:,xe]).upper() != W) and (str(aln[:,xe]).upper() != Y) and (str(aln[:,xe]).upper() != D) and (str(aln[:,xe]).upper() != V) and (str(aln[:,xe]).upper() != astr) and (str(aln[:,xe]).upper() != XUN):
                cc = Counter(str(aln[:,xe].upper()))
                cc = cc.most_common() #here lets make the common base/aa as the WT
                hass.append((cc,int((xe)+1))) 
                mutat_patt.append(cc) #this line to get the pattern of the muatations
       
        yo = []
        prev_list = []
        total_len = len(aln[:,xe])
        for elem in hass:
            tarb = float(100*(int(int(elem[0][0][1]))/total_len))
            prev_list.append(100 -tarb)
            if len(elem[0]) == 2:
                yo.append("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]))
            if len(elem[0]) == 3:
                yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0])))
            if len(elem[0]) > 3:
                yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][3][0])))
    
                          
        pos_list = []
        for ghon in hass:
            pos_list.append(ghon[1])
       
    
        GRG = pd.ExcelWriter("mutations_%s_file.xlsx"%(inp_f))
        df = pd.DataFrame({"Position_in_%s"%(inp_f):pos_list,"{“WT”:frequency ,“mutation(S)”:frequency}":mutat_patt ,"Top mutations":yo,"prev%_of_mutations":prev_list})
        df.to_excel(GRG, index = False)
        GRG.save()
    
        kp =[]
        for xx in hass:#here i want to extract the mutation position like (8,19,200..) from the list of tupules
            kp.append(xx[1])
        kp.pop(0) # to remove the first string in kp list #take care if you repeat this you will lose one element
        lp = list(range(0, len(kp))) #the index of list will represent the mutation number
        print("Done, you have %s mutations" %(len(kp)))
        plt.scatter(kp,lp ,color = 'red' , s =.1,lw=4)
        plt.xlabel('length of the protein')
        plt.xlim(0,len(aln[0].seq))
        plt.ylabel('Number of mutations')
        plt.title('Mutations in your query seq (red dot = 1 mutation)')
        plt.savefig("mutations_map_%s_file.jpg"%(inp_f)) 
        plt.close()
        
        plt.scatter(pos_list,prev_list)
        for i, txt in enumerate(yo):
            plt.annotate(txt, (pos_list[i], prev_list[i]))
        plt.ylabel('prevelance_of_mutations_%')
        plt.xlabel('length of protein')
        plt.xlim(0,len(aln[0].seq))
        plt.savefig("mutations_per_position_%s_file.jpg"%(inp_f))
        plt.close()
        print("here you are: mutations_file.xlsx,mutations_map and frequency plots ")

elif (f =="b"): 
    type_exten = input("What is the (*.extension) of your files (Example : *.aln, *.afa ,*.fa,*.fasta,*.fna, *.phy,*.sth)?:")
    type_format = input("What is the format of your aligned files (Example:clustal, fasta ,phylip,stockholm):")
    q_ref = input("did you have/put your reference sequence on the top of alignment? y/n:")
    for inp_f in os.listdir():
        if fnmatch.fnmatch(inp_f, type_exten):
                aln = AlignIO.read(inp_f,type_format) 
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
                    
            #kindly remeber that sometimes align files come in upper/low case letters so i did str().upper()
                yr = "".join(yr)
                s1 = str(yr).upper() 
                s2 = str(aln[0].seq).upper() #the first seq in the align file
                len1, len2 = len(s1), len(s2)
                ir, jr = 0, -1 #the best solution for longest common substring problem (https://rosettacode.org/wiki/Longest_Common_Substring#Python)
                for i1 in range(len1): #takes around 5 min for lenght 30 kb and number of 2000 sequences
                    i2 = s2.find(s1[i1]) #with my own adaptions
                    while i2 >= 0:
                        j1, j2 = i1, i2
                        while j1 < len1 and j2 < len2 and s2[j2] == s1[j1]:
                            if j1-i1 >= jr-ir:
                                ir, jr = i1, j1
                            j1 += 1; j2 += 1
                        i2 = s2.find(s1[i1], i2+1)
                rs = str(s1[ir:jr+1])
                
               
                
            
                Result = open("longest_conserved_%s_file.fasta"%(inp_f), "w")
                Result.write('>Longest_conserved_seq_in_%s\n'%(inp_f))
                Result.write(rs) #very useful for PCR primer designing or protein domain search
                Result.close()
                #len(yr) == str(aln.column_annotations).count("*") #qc step
                #if you get true, you are in the right track!
                print("Here you are: longest_conserved_seq_%s_file.fasta"%(inp_f))
                
               ####################################################
                
               
                if (q_ref == "y"):
                    
                        hass =[]
                        mutat_patt = []
                        refre = []
                        for xe in range(aln.get_alignment_length()):
                            if (str(aln[:,xe]).upper() != A) and (str(aln[:,xe]).upper() != T) and (str(aln[:,xe]).upper() != C) and (str(aln[:,xe]).upper() != G) and (str(aln[:,xe]).upper() != R) and (str(aln[:,xe]).upper() != N) and (str(aln[:,xe]).upper() != Q) and (str(aln[:,xe]).upper() != H) and (str(aln[:,xe]).upper() != E)  and (str(aln[:,xe]).upper() != I) and (str(aln[:,xe]).upper() != L) and (str(aln[:,xe]).upper() != K) and (str(aln[:,xe]).upper() != M) and (str(aln[:,xe]).upper() != F) and (str(aln[:,xe]).upper() != P) and (str(aln[:,xe]).upper() != O) and (str(aln[:,xe]).upper() != S) and (str(aln[:,xe]).upper() != U) and (str(aln[:,xe]).upper() != W) and (str(aln[:,xe]).upper() != Y) and (str(aln[:,xe]).upper() != D) and (str(aln[:,xe]).upper() != V) and (str(aln[:,xe]).upper() != astr) and (str(aln[:,xe]).upper() != XUN):
                                cc = Counter(str(aln[:,xe].upper()))
                                cc = list(Counter(cc).items())
                                hass.append((cc,int((xe)+1))) 
                                mutat_patt.append(cc) #this line to get the pattern of the muatations
                                refre.append(str(aln[0].seq[xe].upper()))
                       
                        yo = []
                        prev_list = []
                        total_len = len(aln[:,xe])
                        for elem in hass:
                            tarb = float(100*(int(int(elem[0][0][1]))/total_len))
                            prev_list.append(100 -tarb)
                            if len(elem[0]) == 2:
                                yo.append("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]))
                            if len(elem[0]) == 3:
                                yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0])))
                            if len(elem[0]) > 3:
                                yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][3][0])))
                    
                                          
                        pos_list = []
                        for ghon in hass:
                            pos_list.append(ghon[1])
                        
                       
                    
                        GRG = pd.ExcelWriter("mutations_%s_file.xlsx"%(inp_f))
                        df = pd.DataFrame({"Position_in_%s"%(inp_f):pos_list,"{“WT”:frequency ,“mutation(S)”:frequency}":mutat_patt ,"reference":refre,"Top mutations":yo,"prev%_of_mutations":prev_list})
                        df.to_excel(GRG, index = False)
                        GRG.save()
                    
                        kp =[]
                        for xx in hass:#here i want to extract the mutation position like (8,19,200..) from the list of tupules
                            kp.append(xx[1])
                        kp.pop(0) # to remove the first string in kp list #take care if you repeat this you will lose one element
                        lp = list(range(0, len(kp))) #the index of list will represent the mutation number
                        print("Using %s as a reference seq in_%s file...."%(aln[0].id,inp_f))
                        print("Done, you have %s mutations" %(len(kp)))
                        plt.scatter(kp,lp ,color = 'red' , s =.1,lw=4)
                        plt.xlabel('length of the protein')
                        plt.xlim(0,len(aln[0].seq))
                        plt.ylabel('number of mutations in %s')
                        plt.title('Mutations in your query (red dot = 1 mutation)')
                        plt.savefig("mutations_map_%s_file.jpg"%(inp_f)) 
                        plt.close()
                        
                        plt.scatter(pos_list,prev_list)
                        for i, txt in enumerate(yo):
                            plt.annotate(txt, (pos_list[i], prev_list[i]))
                        plt.ylabel('prevelance_of_mutations_%')
                        plt.xlabel('length of the protein')
                        plt.xlim(0,len(aln[0].seq))
                        plt.savefig("mutations_per_position_%s_file.jpg"%(inp_f)) 
                        plt.close()
                        print("here you are: mutations plots")
                       
                    
               
                        ##########################################################################################
                elif (q_ref =="n"):
                    hass =[]
                    mutat_patt = []
                    for xe in range(aln.get_alignment_length()):
                        if (str(aln[:,xe]).upper() != A) and (str(aln[:,xe]).upper() != T) and (str(aln[:,xe]).upper() != C) and (str(aln[:,xe]).upper() != G) and (str(aln[:,xe]).upper() != R) and (str(aln[:,xe]).upper() != N) and (str(aln[:,xe]).upper() != Q) and (str(aln[:,xe]).upper() != H) and (str(aln[:,xe]).upper() != E)  and (str(aln[:,xe]).upper() != I) and (str(aln[:,xe]).upper() != L) and (str(aln[:,xe]).upper() != K) and (str(aln[:,xe]).upper() != M) and (str(aln[:,xe]).upper() != F) and (str(aln[:,xe]).upper() != P) and (str(aln[:,xe]).upper() != O) and (str(aln[:,xe]).upper() != S) and (str(aln[:,xe]).upper() != U) and (str(aln[:,xe]).upper() != W) and (str(aln[:,xe]).upper() != Y) and (str(aln[:,xe]).upper() != D) and (str(aln[:,xe]).upper() != V) and (str(aln[:,xe]).upper() != astr) and (str(aln[:,xe]).upper() != XUN):
                            cc = Counter(str(aln[:,xe].upper()))
                            cc = cc.most_common() #here lets make the common base/aa as the WT
                            hass.append((cc,int((xe)+1))) 
                            mutat_patt.append(cc) #this step to get the pattern of the muatations
                   
                    yo = []
                    prev_list = []
                    total_len = len(aln[:,xe])
                    for elem in hass:
                        tarb = float(100*(int(int(elem[0][0][1]))/total_len))
                        prev_list.append(100 -tarb)
                        if len(elem[0]) == 2:
                            yo.append("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]))
                        if len(elem[0]) == 3:
                            yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0])))
                        if len(elem[0]) > 3:
                            yo.append(("%s%s%s"%(elem[0][0][0],elem[1],elem[0][1][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][2][0]),"%s%s%s"%(elem[0][0][0],elem[1],elem[0][3][0])))
                
                                      
                    pos_list = []
                    for ghon in hass:
                        pos_list.append(ghon[1])
                   
                
                    GRG = pd.ExcelWriter("mutations_%s_file.xlsx"%(inp_f))
                    df = pd.DataFrame({"Position_in_%s"%(inp_f):pos_list,"{“WT”:frequency ,“mutation(S)”:frequency}":mutat_patt ,"Top mutations":yo,"prev%_of_mutations":prev_list})
                    df.to_excel(GRG, index = False)
                    GRG.save()
                
                    kp =[]
                    for xx in hass:#here i want to extract the mutation position like (8,19,200..) from the list of tupules
                        kp.append(xx[1])
                    kp.pop(0) # to remove the first string in kp list #take care if you repeat this you will lose one element
                    lp = list(range(0, len(kp))) #the index of list will represent the mutation number
                    print("Done, you have %s mutations" %(len(kp)))
                    plt.scatter(kp,lp ,color = 'red' , s =.1,lw=4)
                    plt.xlabel('length of the protein')
                    plt.xlim(0,len(aln[0].seq))
                    plt.ylabel('number of mutations')
                    plt.title('Mutations in you query seq (red dot = 1 mutation)')
                    plt.savefig("Mutations_%s_file.jpg"%(inp_f)) 
                    plt.close()
                    prev_list
                    plt.scatter(pos_list,prev_list)
                    for i, txt in enumerate(yo):
                        plt.annotate(txt, (pos_list[i], prev_list[i]))
                    plt.ylabel('prevelance_of_mutations_%')
                    plt.xlabel('length of the protein')
                    plt.xlim(0,len(aln[0].seq))
                    plt.savefig("mutations_per_position_%s_file.jpg"%(inp_f)) 
                    plt.close()
                    print("here you are: mutations plots")

########################################################################################
f  = input("20-do you want to extract mutations per sequence and the frequnet mutations combinations? y/b/any key to skip:" )
if (f == "y"):
    inp_f = input("What is the name of your file:")
    type_format = input("What is the format of your aligned file(Example:clustal,fasta,phylip,stockholm):")
    aln = AlignIO.read(inp_f,type_format)
    f_2  = input("Do you have the reference sequence on top of your alignmnet? y/n:" )
    if (f_2 == "y"):
        zfo = []
        zxg = []
        for y in list(aln):
          for x in range(len(aln[0].seq)):
              if (str((aln[0].seq)[x].upper()) != str((y.seq)[x].upper())):
                  zxg.append(y.id)
                  zfo.append("%s%d%s," %(str((aln[0].seq)[x]),int(x+1),str((y.seq)[x])))
        df_99 = pd.DataFrame({"Seq_ID":zxg,"mutation":zfo})
        df_100 = df_99.groupby(['Seq_ID']).sum()
        df_101 = df_100.mutation.value_counts()
        df_100.to_excel("mutations_per_seq_%s.xlsx"%(inp_f))
        df_101.to_excel("mutations_combination_freq_%s.xlsx"%(inp_f))
        print("here you are: mutations_file.xlsx")

    
    elif (f_2 == "n"):
        print("Then, I will consider the most common (base/aa) as the WT")
        zfo = []
        zxg = []
        for y in list(aln):
          for x in range(len(aln[0].seq)):
              cc = Counter(str(aln[:,x].upper()))
              cc = cc.most_common() #here lets make the common base/aa as the WT
              if (str(cc[0][0].upper()) != str((y.seq)[x].upper())):
                  zxg.append(y.id)
                  zfo.append("%s%d%s," %(str(cc[0][0]),int(x+1),str((y.seq)[x])))
        df_99 = pd.DataFrame({"Seq_ID":zxg,"mutation":zfo})
        df_100 = df_99.groupby(['Seq_ID']).sum()
        df_101 = df_100.mutation.value_counts()
        df_100.to_excel("mutations_per_seq_%s.xlsx"%(inp_f))
        df_101.to_excel("mutations_combination_freq_%s.xlsx"%(inp_f))
        print("here you are: mutations_files.xlsx")

      

elif (f =="b"): 
    type_format = input("What is the format of your aligned files (Example:clustal, fasta ,phylip,stockholm):")
    type_exten = input("What is the (*.extension) of your files (Example : *.aln, *.afa ,*.fa,*.fasta,*.fna, *.phy,*.sth)?:")
    f_2  = input("do you have the reference sequence on top of your alignmnet? y/n:" )
    for inp_f in os.listdir():
           if fnmatch.fnmatch(inp_f, type_exten):
                if (f_2 == "y"):
                    zfo = []
                    zxg = []
                    aln = AlignIO.read(inp_f,type_format)
                    for y in list(aln):
                      for x in range(len(aln[0].seq)):
                          if (str((aln[0].seq)[x].upper()) != str((y.seq)[x].upper())):
                              zxg.append(y.id)
                              zfo.append("%s%d%s," %(str((aln[0].seq)[x]),int(x+1),str((y.seq)[x])))
                    df_99 = pd.DataFrame({"Seq_ID":zxg,"mutation":zfo})
                    df_100 = df_99.groupby(['Seq_ID']).sum()
                    df_101 = df_100.mutation.value_counts()
                    df_100.to_excel("mutations_per_seq_%s.xlsx"%(inp_f))
                    df_101.to_excel("mutations_combination_freq_%s.xlsx"%(inp_f))
                    print("here you are: mutations plots")
            
                
                elif (f_2 == "n"): #no reference
                    print("Then, I will consider the most common (base/aa) as the WT")
                    zfo = []
                    zxg = []
                    aln = AlignIO.read(inp_f,type_format)
                    for y in list(aln):
                      for x in range(len(aln[0].seq)):
                          cc = Counter(str(aln[:,x].upper()))
                          cc = cc.most_common() 
                          if (str(cc[0][0].upper()) != str((y.seq)[x].upper())): 
                              zxg.append(y.id)
                              zfo.append("%s%d%s," %(str(cc[0][0]),int(x+1),str((y.seq)[x])))
                    df_99 = pd.DataFrame({"Seq_ID":zxg,"mutation":zfo})
                    df_100 = df_99.groupby(['Seq_ID']).sum()
                    df_101 = df_100.mutation.value_counts()
                    df_100.to_excel("mutations_per_seq_%s.xlsx"%(inp_f))
                    df_101.to_excel("mutations_combination_freq_%s.xlsx"%(inp_f))
                    print("here you are: mutataion plots")
    

print("Thank you for this journey, \ ͡ ͜ʖ ͡/ Ahmed M. A. Elsherbini")
    

        
           
