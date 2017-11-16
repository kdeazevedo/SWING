import argparse
import os.path
import pathlib
import logging
import subprocess
import re
import numpy as np
import dockerasmus.pdb as pdb
import models.residue
import models.protein
from models.complex import Complex
from functions import angles_generator, angles_random, vec_to_dist
from askInterEvol import runAlign
import filedownload as fd
from alignInterolog import runProfit

FORMAT = '%(asctime)s - %(name)s - %(levelname)8s : %(message)s'
logging.basicConfig(format=FORMAT,level=logging.DEBUG)
logger = logging.getLogger('main')

parser = argparse.ArgumentParser()
parser.add_argument('-rec',help='Recepter\'s file path')
parser.add_argument('-lig',help='Ligand\'s file path')
parser.add_argument('-n',help='Number of sampling')
parser.add_argument('-fasta_dir',default='FASTA',help='The directory where the FASTA sequences will be stored')
parser.add_argument('-inter_dir',default='INTER',help='The directory where the interologs will be stored')

logger.debug('Start parsing arguments')
args = parser.parse_args()


PDB = r'(\w+).pdb'
# Check if files exist
assert os.path.isfile(args.rec), "Recepter file not found"
assert os.path.isfile(args.lig), "Ligand file not found"
LIG = re.match(PDB,pathlib.PurePath(args.lig).name).group(1)
REC = re.match(PDB,pathlib.PurePath(args.rec).name).group(1)
subprocess.call(['cp',args.rec,'Proteins'])
subprocess.call(['cp',args.lig,'Proteins'])
rec = pdb.Protein.from_pdb_file(args.rec)
lig = pdb.Protein.from_pdb_file(args.lig)
rec.name = REC
rec.path = os.path.join('Proteins',REC+'.pdb')
lig.name = LIG
lig.path = os.path.join('Proteins',LIG+'.pdb')

# Create if fasta or interolog directory doesn't exist
FASTADIR = args.fasta_dir
INTERDIR = args.inter_dir
pathlib.Path(FASTADIR).mkdir(parents=True, exist_ok=True) 
pathlib.Path(INTERDIR).mkdir(parents=True, exist_ok=True) 

n_samples = int(args.n)
n_digit = len(args.n)


################################
###  Download from InterEvol ###
################################
ligf = os.path.join(FASTADIR,LIG+'.fasta')
recf = os.path.join(FASTADIR,REC+'.fasta')
lig.write_fasta(ligf)
rec.write_fasta(recf)
dico = runAlign(recf,ligf,INTERDIR)
PDBid = [k for k in dico.keys()]
print(dico)


################################
###  Alignment with Profit   ###
################################


for key in dico.keys():
	#print(key)
    #print(dico[key])
    liste = dico[key]
    #print(liste[3],liste[2])
    deg = min(int(liste[0][:-1]),int(liste[1][:-1]))
    runProfit(lig.path, rec.path, os.path.join(INTERDIR,key+".pdb"), liste[3], liste[2])




lig_aligned = pdb.Protein.from_pdb_file('Proteins/'+LIG+'_aligned.pdb')
lig_aligned.name = LIG
cpx = Complex(rec,lig_aligned)

################################
###         Sampling         ###
################################


for idx,l in enumerate(angles_generator(n_samples,deg=deg)):
    A = cpx.rotations(l[0],l[1],l[2],l[3],l[4])
    D = cpx.ca_dist(A)
    i,j = np.unravel_index(D.argmin(), D.shape)
    m = D[i,j]
    if m <5:
        print(m,'--')
        A = A+vec_to_dist(cpx.rec.get_ca()[i],A[cpx.lig.get_ca_ind()][j],5)
        print(np.min(cpx.ca_dist(A)))
    cpx.lig.write_atoms('Proteins/B'+'{:05d}'.format(idx)+'.pdb',A)
    

for k in range(n_samples):
    subprocess.call(["python3", 'Minimizer/runMini.py', '-rec','Proteins/'+REC+'.pdb', '-lig','Proteins/B'+'{:05d}'.format(k)+'.pdb'])

