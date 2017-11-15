import argparse
import os.path
import pathlib
import logging
import subprocess
import re
import dockerasmus.pdb as pdb
import models.residue
import models.protein
from models.complex import Complex
from functions import angles_generator, angles_random
from askInterEvol import runAlign
import filedownload as fd

FORMAT = '%(asctime)s - %(name)s - %(levelname)8s : %(message)s'
logging.basicConfig(format=FORMAT,level=logging.DEBUG)
logger = logging.getLogger('main')

parser = argparse.ArgumentParser()
parser.add_argument('-rec',help='Recepter\'s file path')
parser.add_argument('-lig',help='Ligand\'s file path')
parser.add_argument('-n',help='Number of sampling')
parser.add_argument('-fasta_dir',default='FASTA',help='The directory where the FASTA sequences will be stored')
parser.add_argument('-inter_dir',default='INTER',help='The directory where the interologs will be stored')

logger.DEBUG('Start parsing arguments')
args = parser.parse_args()


PDB = r'(\w+).pdb'
# Check if files exist
assert os.path.isfile(args.rec), "Recepter file not found"
assert os.path.isfile(args.lig), "Ligand file not found"
LIG = re.match(PDB,pathlib.PurePath(args.lig).name).group(1)
REC = re.match(PDB,pathlib.PurePath(args.rec).name).group(1)
rec = pdb.Protein.from_pdb_file(args.rec)
lig = pdb.Protein.from_pdb_file(args.lig)
rec.name = REC
lig.name = LIG
cpx = Complex(rec,lig)

# Create if fasta or interolog directory doesn't exist
FASTADIR = args.fasta_dir
INTERDIR = args.inter_dir
pathlib.Path(FASTADIR).mkdir(parents=True, exist_ok=True) 
pathlib.Path(INTERDIR).mkdir(parents=True, exist_ok=True) 




################################
###  Download from InterEvol ###
################################
ligf = os.path.join(FASTADIR,LIG+'.fasta')
recf = os.path.join(FASTADIR,REC+'.fasta')
lig.write_fasta(ligf)
rec.write_fasta(recf)
assert False
dico = runAlign(recf,ligf,INTERDIR)
PDBid = [k for k in dico.keys()]
#fd.downloadPDB(PDBid,INTERDIR) #obsolete as file downloaded directly from InterEvol




################################
###         Sampling         ###
################################

subprocess.call(['cp',args.rec,'Proteins'])
#counter = 0
#ntotal = 0
#while counter<500:
#    l =angles_random()
#    A = cpx.rotations(l[0],l[1],l[2],l[3],l[4])
#    m = cpx.min_ca(A)
#    ntotal += 1
#    if m <6 and m>4:
#        counter += 1
#print(ntotal)

for l in angles_generator(args.n):
    A = cpx.rotations(l[0],l[1],l[2],l[3],l[4])
    m = cpx.min_ca(A)
    if m <10:
        counter += 1
print(counter)
assert False

lig.write_atoms('Proteins/1JP3_B_rota12.pdb')
subprocess.call(["python3", 'Minimizer/runMini.py', '-rec','Proteins/1JP3_A.pdb', '-lig','Proteins/1JP3_B_rota12.pdb'])


# Multiple templates
positions = download(rec,lig)

for pos in positions:
    for t in translation():
        for r in rotation():
            lig.trans_rotate(trans=t,rotate=r)
            lig.write_atoms('Proteins/ligand.pdb')
            subprocess.call(["python3", 'Minimizer/runMini.py', '-rec',rec.path(), '-lig',lig.path()])
