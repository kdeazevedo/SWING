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
    with open(outdir+filename+".fasta.txt", "wb") as code:
        code.write(r.content)
    print('Finished downloading FASTA')

def downloadPDB(PDBid, outdir):
    """Download the PDB file corresponding to an ID
    input1 : a list of 4 letter PDB ID
    input2 : the path of the directory containing the PDB files of the interologs
    output : a PDB file
    """
    for id in PDBid:
        print('Downloading PDB file'+id)
        url = "https://files.rcsb.org/download/"+id+".pdb"
        r = requests.get(url)
        with open(outdir+id+".pdb","wb") as code:
            code.write(r.content)
        print('Finished downloading PDB file')