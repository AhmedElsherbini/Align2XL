# Multifasta_analysis_tool

This tool is used to filter and manipulate of multifasta files. Also you can extract genes from the files !

## Installation

Make sure that you have python 3.7 or above and all the dependacies (Biopython, matplotlib, pandas , numpy , difflib, re )

The tool was written to be used on Linux command line or any python IDE !

## Usage
For own comfort, make sure you have the files in the same directory as the tool!

```python
python3 Multifasta_analysis_tool.py

```
then just answer questions !

## This tool can be used if..
1- you want extract genes using its start/end postion in aligned_file.afa (more accurate) or fasta file

2- you want to exlude sequences in a multifasta file using sequneces pattern (ex:NNN,XX)

3- you want to  print all  ID headers in your multifasta

4- you want extract gene from fasta/multifasta file using it's upstream and downstream seq

5- you want to translat DNA fasta file on the its 1 frame

6- you want to efetch data from NCBI or you want to merge your fasta files 

7- you want to know GC content and N bases content of your DNA seq

8- you want to align muscle or mafft

9- you want to draw a tree.dnd using python matplotlib

10-you want to extract the longest conserved and call mutations between your genomes , gene proteins clustal_file.aln and the output is just xlsx file 

## Contributing
Pull requests are very welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.
the tool is still under development !
## License
Please, cite my page if this tool was useful for your work
