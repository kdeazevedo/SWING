#!/usr/bin/env python

import glob, shutil, re, tarfile, subprocess
import sys, os, time
import logging

logger = logging.getLogger('main.minimizer')
WORKDIR = ''
FILEDIR = ''


def writeMini(inlines, receptor, ligand):
    """ Purpose: writes the script for the minimization from the template (wd/templates/runMini.sh) and args: 
             modifies prot names and paths
        Input: template content for the minimizer script stored in inlines, recpetor's name, ligand's name
        Output: runMini.sh 
    """
    if ligand == receptor:
        sep = ""
        ligand = (receptor, "_1")
        ligand = sep.join(ligand)
    outfile = open("%s/minimizer/run_mini/runMini.sh"%(WORKDIR), "w+")
    for line in inlines:
        if line[0:10] == "set SEPDIR":
            outfile.write("set SEPDIR = {}".format(os.path.join(FILEDIR,'Proteins')))
        elif line[0:13] == "foreach PROT ":
            outfile.write("foreach PROT (\'\\ls $SEPDIR/%s\')\n"%(receptor))
        elif line[0:13] == "foreach PROTT":
            outfile.write("foreach PROTT (\'\\ls $SEPDIR/%s\')\n"%(ligand))
        #elif line[0:10] == "set SEPDIR":
        #    outfile.write("set SEPDIR = %s/Proteins"%WORKDIR)
        else:
            outfile.write(line)
    outfile.close()
    os.system("chmod +x %s/minimizer/run_mini/runMini.sh"%(WORKDIR))

def writeBuilder(inlines, receptor, ligand):
    """ Purpose: writes the script for the PDB reconstruction from the template (wd/templates/rctrPDB.sh) and args: 
             modifies prot names and paths
        Input: template content for the reconstruction script stored in inlines, recpetor's name, ligand's name
        Output: rctrPDB.sh 
    """
    logger.debug(receptor)
    if ligand == receptor:
        sep = ""
        ligand = (receptor, "_1")
        ligand = sep.join(ligand)
    outfile = open("%s/minimizer/run_builder/rctrPDB.sh"%(WORKDIR), "w+")
    for line in inlines:
        if line[0:11] == "set PROTDIR":
            outfile.write("set PROTDIR = {}".format(os.path.join(FILEDIR,'Proteins')))
        elif line[0:13] == "foreach PROT ":
            outfile.write("foreach PROT (%s)\n"%(receptor))
        elif line[0:13] == "foreach PROTT":
            outfile.write("foreach PROTT (%s)\n"%(ligand))
        #elif line[0:11] == "set PROTDIR":
        #    outfile.write("set PROTDIR = %s/Proteins"%WORKDIR)
        else:
            outfile.write(line)
    outfile.close()
    os.system("chmod +x %s/minimizer/run_builder/rctrPDB.sh"%(WORKDIR))


############################ ARGS ################################


def run(rec,lig,filedir,conf=1):
    logger.debug(rec)
    global WORKDIR
    global FILEDIR
    FILEDIR = filedir
    WORKDIR = os.path.join(os.path.abspath(os.path.join(__file__,os.pardir)),'Minimizer')
    logger.debug('mini {}'.format(WORKDIR))
    
    ############################ END ARGS ################################
    
    #### preparation of executables
    os.system("chmod +x %s/minimizer/progs_MAXDo/simulmain.out"%(WORKDIR))
    os.system("chmod +x %s/minimizer/progs_MAXDo/Getarea2.out"%(WORKDIR))
    os.system("chmod +x %s/minimizer/progs_builder/Interface.out"%(WORKDIR))
    
    
    #### template file for minimzation
    infileMini = open("%s/templates/runMini.sh"%(WORKDIR))
    inlinesM = infileMini.readlines()
    infileMini.close()
    
    #### template file for reconstruction
    infileBuild = open("%s/templates/rctrPDB.sh"%(WORKDIR))
    inlinesB = infileBuild.readlines()
    infileBuild.close()
    
    #### Creates script for launching minimization and pdb reconstruction
    writeMini(inlinesM, rec, lig)
    writeBuilder(inlinesB, rec, lig)
    
    #### Executes minimization
    os.chdir("%s/minimizer/run_mini/"%(WORKDIR))
    logger.debug("ready to start minimization1 \n============================\nworking directory: {}".format(os.getcwd()))
    os.system("./runMini.sh")
    
    #### Unfortunately, needs to modify format of output for rctrPDB.sh 
    globfile = open("%s/minimizer/run_mini/global.dat"%(WORKDIR), "r")
    infile = globfile.readlines()
    globfile.close()
    
    newglob = open("%s/minimizer/run_builder/global.dat"%(WORKDIR), "w+")
    for line in infile:
        globalelement = line.split()
        linetowrite = ('{0[0]:>5s} {0[1]:>3s} {0[2]:>13s} {0[3]:>12s} {0[4]:>12s} {0[5]:>12s} {0[6]:>12s} {0[7]:>12s} {0[8]:>12s} {0[9]:>12s} {0[10]:>12s}\n'.format(globalelement))
        newglob.write(linetowrite)
    newglob.close()
    
    #### Executes pdb reconstruction
    os.chdir("%s/minimizer/run_builder/"%(WORKDIR))
    logger.debug("working directory: {}".format(os.getcwd()))
    os.system("./rctrPDB.sh")
    
    #### Move newly created pdb file to pdb conf directory as well as the global file containing the energy information
    recid = os.path.basename(rec).split(".")[0]
    logger.debug(recid)
    ligid = os.path.basename(lig).split(".")[0]
    logger.debug(ligid)
    rootrecid = os.path.basename(rec)[:6]
    logger.debug(rootrecid)
    rootligid = os.path.basename(lig)[:6]
    logger.debug(rootligid)
    
    # PDB
    pdbfile = "%s/minimizer/run_builder/pdb_files/%s-%s.min1.pdb"%(WORKDIR, rootrecid, rootligid)   
    pdbpath = "{}/pdb_mini/{}_{}_min_{:05d}.pdb".format(FILEDIR, recid, ligid, conf)
    # global file
    inpath = "%s/minimizer/run_mini/global.dat"%(WORKDIR)
    outpath = "{}/global_out/global_{}_{}_min_{:05d}.dat".format(FILEDIR, recid, ligid, conf)
    try:
        assert os.path.isfile(pdbfile)
        assert os.path.isfile(inpath)
        subprocess.call(["mv", pdbfile, pdbpath])
        subprocess.call(["mv", inpath, outpath])
        return True
    except:
        logger.warn('Failed on ligand No. {}'.format(ligid))
        return False
    
    
