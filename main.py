import argparse
import os.path
import logging
import dockerasmus.pdb as pdb
import models.protein
import subprocess

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
