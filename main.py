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
from alignInterolog import runPymolAlignment
import runMini as mini
from argParser import parser
import json
from datetime import datetime

# Define arguments from argParser.py
args = parser.parse_args()


PDB = r'(\w+).pdb'
# Check if files exist
assert os.path.isfile(args.rec), "Recepter file not found"
assert os.path.isfile(args.lig), "Ligand file not found"
# Create output folders:
# out
#  - Fasta        : store proteins' fasta files
#  - Inter        : store Interologs dowloaded from InterEvol, their alignements and config files
#  - Complex      : store complex's pdb files after sampling (recepter + rotated ligand)
#  - Proteins     : store all proteins' pdb files (original recepter, original ligand, rotated ligand, ... )
#  - pdb_mini     : store minimiser's pdb output
#  - global_out   : store minimiser's global output
pathlib.Path(args.o).mkdir(parents=True, exist_ok=True)
fld_lst = [('FASTA','Fasta'),('INTER','Inter'),('CPX','Complex'),('PRO','Proteins'),
    ('PMIN','pdb_mini'),('GOUT','global_out')]
FLD = {'OUT':os.path.abspath(args.o)}
for d,n in fld_lst:
    t = os.path.join(FLD['OUT'],n)
    pathlib.Path(t).mkdir(parents=True,exist_ok=True)
    FLD[d] = t

################################
###          Loggers         ###
################################

# A filter that returns child name of logger
class LastPartFilter(logging.Filter):
    def filter(self, record):
        record.name_last = record.name.rsplit('.', 1)[-1]
        return True
FORMAT = '%(asctime)s - %(name_last)10s - %(levelname)8s - %(message)s'
logger = logging.getLogger('main')
logger.setLevel(logging.DEBUG)
# Define file handler, log file is stored in OUTPUT, level = INFO
fh = logging.FileHandler(os.path.join(FLD['OUT'],'log_{:%Y%m%d_%H%M%S}.txt'.format(datetime.now())))
fh.setLevel(logging.INFO)
# Define stream handler, level = DEBUG
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter(FORMAT)
fh.setFormatter(formatter)
ch.setFormatter(formatter)
logger.addHandler(fh)
logger.addHandler(ch)
fh.addFilter(LastPartFilter())
ch.addFilter(LastPartFilter())

# Extract recepter's and ligand's information
# Create recepter object
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


################################
###  Download from InterEvol ###
################################

# Download interologs if command is run or download
if args.cmd == 'run' or args.cmd == 'download':
    ligf = os.path.join(FLD['FASTA'],LIG+'.fasta')
    recf = os.path.join(FLD['FASTA'],REC+'.fasta')
    lig.write_fasta(ligf)
    rec.write_fasta(recf)
    logger.debug("Start download. Please wait...")
    dico = runAlign(recf,ligf,FLD['INTER'])
    # Write result in json format with prefix "Inter"
    with open(os.path.join(interDirectory,'Inter_{}.conf'.format(REC)),'w') as f:
        json.dump(dico,f,indent=2)
        logger.info("Write result into {}".format(f.name))


################################
###   Alignment with Pymol   ###
################################

# Start alignment with pymol if command is run or align
if args.cmd == 'run' or args.cmd == 'align':
    # Reconstruct dico from config file
    if args.cmd == 'align':
        with open(args.config,'r') as f:
            dico = json.load(f)
    logger.debug("Start of alignment with Pymol")
    # Run alignment on each interolog
    # Then, add aligned ligand's pdb file path into dico
    for key in dico.keys():
        liste = dico[key]
        res = runPymolAlignment(lig.path,rec.path,os.path.join(FLD['INTER'],key+".pdb"),liste['idn_l'],liste['idn_r'])
        dico[key]['lig_aligned']=res
    # Write dico into config file in json format with prefix "Samples"
    with open(os.path.join(FLD['INTER'],'Samples_{}.conf'.format(REC)),'w') as f:
        json.dump(dico,f,indent=2)
        logger.info("Write result into {}".format(f.name))
    logger.debug("End of alignment with Pymol")
    
################################
###         Sampling         ###
################################

# Run sampling if command is run or samples
if args.cmd == 'run' or args.cmd == 'samples':
    # Reconstruct dico from config file
    if args.cmd == 'samples':
        with open(args.config,'r') as f:
            dico = json.load(f)
    # Run sampling for each interolog
    for key in dico.keys():
        liste = dico[key]
        # Degree = min of interolog's recepter and ligand's identity
        deg = min(int(liste['idn_r'][:-1]),int(liste['idn_l'][:-1]))
        lig_aligned = pdb.Protein.from_pdb_file(liste['lig_aligned'])
        lig_aligned.name = LIG
        # Create Complex object for recepter and aligned ligand
        cpx = Complex(rec,lig_aligned)
        logger.debug('Start sampling of template {}'.format(key))
        # Rotate ligand for each generated angles
        for idx,l in enumerate(angles_generator(args.n,deg=deg)):
            logger.info('Ligand {:>6}\'s rotatation No. {:06d}({})'.format(cpx.lig.name,idx,l))
            A = cpx.rotations(l[0],l[1],l[2],l[3],l[4])
            # Move lignad if the minimum distance between two carbon alpha is less than 5
            D = cpx.ca_dist(A)
            i,j = np.unravel_index(D.argmin(), D.shape)
            m = D[i,j]
            if m <5:
                A = A+vec_to_dist(cpx.rec.get_ca()[i],A[cpx.lig.get_ca_ind()][j],5)
            # Write rotated ligand file in pdb 
            if args.minimizer:
                cpx.lig.write_atoms(os.path.join(FLD['PRO'],'{:06d}.pdb'.format(idx)),A)
            else:
                cpx.lig.write_atoms(os.path.join(FLD['PRO'],'{}_{}_{:06d}.pdb'.format(cpx.lig.name,key,idx)),A)
        logger.debug("End sampling of template {}".format(key))
        # Run minimizer for each rotation
        if args.minimizer: 
            logger.debug('Start of minimizer')
            for k in range(args.n):
                conf = "{}_{:06d}".format(key,k)
                lig_pro = os.path.join(FLD['PRO'],'{:06d}.pdb'.format(k))
                lig_pro_o = os.path.join(FLD['PRO'],'{}_{}.pdb'.format(cpx.lig.name,conf))
                successed = mini.run(cpx.rec.path,lig_pro,FLD['OUT'],conf=conf)
                # Rename rotated ligand pdb file
                subprocess.call(["mv",lig_pro,lig_pro_o])
                # Rename ligand file after minimizer and create complex pdb file
                if successed:
                    out_f_1 = os.path.join(FLD['OUT'],'pdb_mini/{}_min.pdb'.format(conf))
                    out_f = os.path.join(FLD['OUT'],'pdb_mini/{}_{}_min.pdb'.format(cpx.lig.name,conf))
                    subprocess.call(["mv",out_f_1,out_f])
                    cpx_out_file = os.path.join(FLD['CPX'],'{}_cpx_{}.pdb'.format(cpx.rec.name,conf))
                    subprocess.call("sed '1d' {} > tmp.txt; cat {} tmp.txt > {}; rm tmp.txt".format(out_f,cpx.rec.path,cpx_out_file),shell=True)
            logger.debug('End of minimizer')
logger.debug('End of program')
