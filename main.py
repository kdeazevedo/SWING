import argparse
import os.path
import logging
import dockerasmus.pdb as pdb
import models.protein
from models.complex import Complex
import subprocess
from functions import angles_generator, angles_random

logger = logging.getLogger('spam_application')

parser = argparse.ArgumentParser()
parser.add_argument('-rec',help='Recepter\'s file path')
parser.add_argument('-lig',help='Ligand\'s file path')
args = parser.parse_args()

# Check if files exist
assert os.path.isfile(args.rec), "Recepter file not found"
assert os.path.isfile(args.lig), "Ligand file not found"

subprocess.call(['cp',args.rec,'Proteins'])
rec = pdb.Protein.from_pdb_file(args.rec)
lig = pdb.Protein.from_pdb_file(args.lig)

cpx = Complex(rec,lig)
counter = 0
ntotal = 0
while counter<500:
    l =angles_random()
    A = cpx.rotations(l[0],l[1],l[2],l[3],l[4])
    m = cpx.min_ca(A)
    ntotal += 1
    if m <6 and m>4:
        counter += 1
print(ntotal)

assert False
for l in angles_generator(1000):
    A = cpx.rotations(l[0],l[1],l[2],l[3],l[4])
    m = cpx.min_ca(A)
    if m <10:
        counter += 1
print(counter)

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
