![align2xl](https://user-images.githubusercontent.com/49863685/110328892-3580b880-801c-11eb-83f7-8b646556fa18.jpg)

Since 30/12/2020: the script is now updating contact me to get the new release

This tool is used to filter and manipulate  multifasta files (DNA or protein). Therefore, you can extract/ exlude sequneces from the file. Also, you can convert  to get the longest conserved seq among your culstal files and you will get you a list and plots of the variants.

## Installation

Make sure that you have python 3.7 or > and all the dependacies (Biopython, matplotlib, pandas , numpy , difflib, re , collections )

The tool was written to be used on Linux command line or any python IDE !

## Usage
For your own comfort, make sure you have the files in the same directory as the tool!
Type in your command line, then just answer questions !


```python
python3 Align2XL.py

```

## This tool can be used if..
1- you want to extract a gene using its start/end postion in aligned_file.afa (more accurate) or multi/fasta file

2- you want to inlucde certain sequences in your multifasta file with a certain pattern (ex: gene with certain mutation)

3- you want to exlude sequences in a multifasta file using a pattern in the sequneces (ex:NNN,XX)

4- you want to  print all  ID headers in your multifasta

5- you want to include/exlcude sequences using a pattern in their headers/metadate (ex:> E.coli-2019-)

6- you want extract a gene from fasta/multifasta file using it's upstream and downstream sequence

7- you want to translate DNA multifasta file or ( a batch of files at onece) on it's (theirs) 1st frame

8- you want to know GC content and N bases content of your multifasta seq file

10- (Align2xlsx) you want to extract the longest conserved sequnece and you want to call variants between your aligned genomes , genes or proteins from a clustal_file.aln and the output is just a fasta for your conserved seq and the xlsx file of your variants with a graph which maps your mutations are!


Extra: you can Align2xlsx for a batch of clustal files at once

## Contributing
Pull requests are very welcome. 


For major changes, please open an issue (or Email me) first to discuss what you would like to change.

The tool is still under development !

Contact me directly on email : drahmedsherbini@yahoo.com
## License
This tool is a copy right of the author and  part of a project related to Nile University 

Please, cite my page if this tool was useful for your work, until the paper is out!
