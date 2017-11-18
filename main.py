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
parser.add_argument('-rec',required=True,help='Recepter\'s file path')
parser.add_argument('-lig',required=True,help='Ligand\'s file path')
parser.add_argument('-n',help='Number of sampling')
parser.add_argument('-o',default='out',help='The directory where program\'s output sotres')

logger.debug('Start parsing arguments')
args = parser.parse_args()


PDB = r'(\w+).pdb'
# Check if files exist
assert os.path.isfile(args.rec), "Recepter file not found"
assert os.path.isfile(args.lig), "Ligand file not found"

# Create output folders
pathlib.Path(args.o).mkdir(parents=True, exist_ok=True)
fld_lst = [('FASTA','Fasta'),('INTER','Inter'),('CPX','Complex'),('PRO','Proteins'),('PMIN','pdb_mini'),('GOUT','global_out')]
FLD = {'OUT':os.path.abspath(args.o)}
for d,n in fld_lst:
    t = os.path.join(FLD['OUT'],n)
    pathlib.Path(t).mkdir(parents=True,exist_ok=True)
    FLD[d] = t

LIG = re.match(PDB,pathlib.PurePath(args.lig).name).group(1)
REC = re.match(PDB,pathlib.PurePath(args.rec).name).group(1)
subprocess.call(['cp',args.rec,FLD['PRO']])
subprocess.call(['cp',args.lig,FLD['PRO']])
rec = pdb.Protein.from_pdb_file(args.rec)
lig = pdb.Protein.from_pdb_file(args.lig)
rec.name = REC
rec.path = os.path.join(FLD['PRO'],REC+'.pdb')
lig.name = LIG
lig.path = os.path.join(FLD['PRO'],LIG+'.pdb')

n_samples = int(args.n)
n_digit = len(args.n)


################################
###  Download from InterEvol ###
################################
ligf = os.path.join(FLD['FASTA'],LIG+'.fasta')
recf = os.path.join(FLD['FASTA'],REC+'.fasta')
lig.write_fasta(ligf)
rec.write_fasta(recf)
dico = runAlign(recf,ligf,FLD['INTER'])


################################
###  Alignment with Profit   ###
################################


for key in dico.keys():
	#print(key)
    #print(dico[key])
    liste = dico[key]
    #print(liste[3],liste[2])
    deg = min(int(liste[0][:-1]),int(liste[1][:-1]))
    runProfit(lig.path, rec.path, os.path.join(FLD['INTER'],key+".pdb"), liste[3], liste[2])




lig_aligned = pdb.Protein.from_pdb_file(os.path.join(FLD['PRO'],LIG+'_aligned.pdb'))
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
        A = A+vec_to_dist(cpx.rec.get_ca()[i],A[cpx.lig.get_ca_ind()][j],5)
    cpx.lig.write_atoms(os.path.join(FLD['PRO'],'B{:05d}.pdb'.format(idx)),A)
    

for k in range(n_samples):
    subprocess.call(["python3", 'Minimizer/runMini.py', '-rec',cpx.rec.path, '-lig',os.path.join(FLD['PRO'],'B{:05d}.pdb'.format(k)),'-o',FLD['OUT']])
    out_file = os.path.join(FLD['OUT'],'pdb_mini/{}_B{:05d}_min1.pdb'.format(cpx.rec.name,k))
    cpx_out_file = os.path.join(FLD['CPX'],'{}_B{:05d}.pdb'.format(cpx.rec.name,k))
    subprocess.call("sed '1d' {} > tmp.txt; cat {} tmp.txt > {}; rm tmp.txt".format(out_file,cpx.rec.path,cpx_out_file),shell=True)
