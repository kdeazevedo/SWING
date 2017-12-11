'''
Created on 2 dec. 2017

@author: marine

Purpose : use clusco to cluster pdb files from interEvolAlign and sampling before and after minimization. Require six argument :
	- pdbdir : the directory of pdb to cluster
Output : a list of clustered pdbs
Use (example) : python cluster.py -pdbdir test
'''

import os
import sys
import subprocess
import logging

logger = logging.getLogger('main.cluster')

def createList (pdb_directory, step):
    '''
    create a pdb list from the directory pdb's content
    input :
        - pdb_directory : the name of the repository which contains the pdb to cluster
        - step : indicate if the clustering is about ligand initial position (init) or sampled ligand conformations (sample)
    output : one
    '''    
    if step == "init" :
        cmd="ls "+os.path.join(pdb_directory,"*_aligned.pdb")+" > pdb_list"
    elif step == "sample" :
        cmd="ls "+os.path.join(pdb_directory,"*.pdb")+" > pdb_list"
    else :
        logger.warn("Invalid step argument for pdb list creation before clustering")
    subprocess.call(cmd, shell=True)
    
    return None
    

def runClusco(pdbListName) :
    '''
    cluster pdbs from a pdb list using clusco hierachical option
    input :
        - pdbListName : the name of the pdb_directory : the name of the pdb which contains the pdb to cluster
    output : a list of pdbs
    '''    
    # Test that the pdb list is not empty
    clusteredPDB=[]
    with open(pdbListName,'r') as f:
        count = sum(1 for line in f)
    if count == 0:
        logger.warn("pdb list to cluster is empty")
    # Cluster only if more than one pdb
    elif count > 2:
        cmdClusco="clusco -l "+pdbListName+" -s rmsd -o cluster_matrix.log -d 0.02 3"
        subprocess.call(cmdClusco, shell=True)
        print ("")
    
        integers=["0","1","2","3","4","5","6","7","8","9"]
        with open (pdbListName+".clustering0", "r") as cluster :
            for line in cluster :
                if line[3:4] in integers :
                    top=line.split()
                    clusteredPDB.append(top[2])

    return clusteredPDB
    


'''    
if __name__ == '__main__':

#Try cluster
    try :
        pdbdir = sys.argv[sys.argv.index("-pdbdir")+1]
    except :
        print("ERROR : specified pdb directory does not exist\n")

                
createList(pdbdir, step="init")
pdblist=runClusco(pdbListName="pdb_list")
print ("")

dico={
  "1iy9_CB": {
    "idn_r": "17%",
    "lig_aligned": "/home/kazevedo/Documents/M2/MeetU/2017-2018_Equipe10/out/Proteins/1FFT_B_sep_1xme_AB_aligned.pdb",
    "chn_l": "B",
    "idn_l": "18%",
    "chn_r": "A"
  },
  "3hb3_AB": {
    "idn_r": "19%",
    "lig_aligned": "/home/kazevedo/Documents/M2/MeetU/2017-2018_Equipe10/out/Proteins/1FFT_B_sep_3hb3_AB_aligned.pdb",
    "chn_l": "B",
    "idn_l": "38%",
    "chn_r": "A"
  },
  "2gsm_AB": {
    "idn_r": "16%",
    "lig_aligned": "/home/kazevedo/Documents/M2/MeetU/2017-2018_Equipe10/out/Proteins/1FFT_B_sep_2gsm_AB_aligned.pdb",
    "chn_l": "B",
    "idn_l": "38%",
    "chn_r": "A"
  },
  "1v54_AB": {
    "idn_r": "16%",
    "lig_aligned": "/home/kazevedo/Documents/M2/MeetU/2017-2018_Equipe10/out/Proteins/1FFT_B_sep_1v54_AB_aligned.pdb",
    "chn_l": "B",
    "idn_l": "39%",
    "chn_r": "A"
  }
}

for item in dico.keys() :
    print (item)
    if not any(item in s for s in pdblist):
        del dico[item]
print ("\n")
for item in dico.keys() :
    print (item)

'''
