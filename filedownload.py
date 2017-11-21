import os
import requests

'''
This file contains simple functions that retrieve PDB files inside a given directory,
extract the PDB 4-letter ID to download the FASTA file or the PDB file
'''

def getPDB(filename):
    """
    Extract the 4-letter ID of a PDB file
    input : a PDB file
    ouput : the 4 letter ID
    """
    filename = filename.replace(".pdb", "")
    return(filename)

def parseDirectory(indir,outdir):
    """
    Parses a directory to retrieve PDB files inside
    And uses the downloadFASTA function to retrieve the FASTA sequence
    input1 : the path of the directory containing the PDB files
    input2 : the path of the directory containing the FASTA sequences
    """
    for files in os.listdir(indir):
        if files.endswith(".pdb") :
            downloadFASTA(getPDB(files),outdir)
          
def downloadFASTA(filename,outdir):
    """
    Retrieves the FASTA sequence of a protein subunit from the PDB website
    input1 : a PDB 4-letter ID
    input2 : the path of the directory containing the FASTA sequences
    ouptut : creates a file containing the FASTA sequence
    """
    filename = filename.upper()
    print("Retrieving FASTA sequence of"+filename)
    url = "https://www.rcsb.org/pdb/download/downloadFastaFiles.do?structureIdList="+filename+"&compressionType=uncompressed"
    r = requests.get(url)
    with open(os.path.join(outdir,filename+".fasta.txt"), "wb") as code:
        code.write(r.content)
    print('Finished downloading FASTA')

def downloadPDB(PDBid, outdir):
    """
    Download the PDB file corresponding to an ID
    input1 : a list of 4 letter PDB ID
    input2 : the path of the directory containing the PDB files of the interologs
    output : a PDB file
    """
    for pid in PDBid:
        print('Downloading PDB file'+pid)
        url = "https://files.rcsb.org/download/"+pid+".pdb"
        r = requests.get(url)
        with open(os.path.join(outdir,pid+".pdb"),"wb") as code:
            code.write(r.content)
        print('Finished downloading PDB file')
        

def FASTAfromPDB(PDBfile, indir, outdir):
    """
    Create a FASTA file from a PDB file reading the amino acids in the fourth column
    input1 : the name of the PDB file
    input2 : the directory containing the PDB file
    input3 : the directory where the FASTA file will be stored
    ouput : a FASTA file
	@originalauthor : github.com/jameslyons (source code modified)
    """
    
    input_file = open(os.path.join(indir,PDBfile),"r")
    output_file_path = os.path.join(outdir,PDBfile+".fasta")
    try:
        os.mknod(output_file_path)
    except:
        print("It seems that this file already exists")
    output_file = open(output_file_path,"r+")
    output_file.write(">"+PDBfile+"\n")

    letters = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H',
               'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
               'TYR':'Y','VAL':'V', 'BGLU' : 'E'}

    prev = '-1'
        
    for line in input_file:
        toks = line.split()
        if len(toks)<1: continue
        if toks[0] != 'ATOM': continue
        if len(toks[4]) == 1:
            if (toks[5] != prev) & (toks[3] in letters.keys()):
                output_file.write('%c' % letters[toks[3]])
                prev = toks[5]
        else:
             if (toks[4] != prev) & (toks[3] in letters.keys()):
                output_file.write('%c' % letters[toks[3]])
                prev = toks[4]           

    output_file.write('\n')
    input_file.close()
    output_file.close()
    return(output_file_path)

