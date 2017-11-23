'''
Created on 13 nov. 2017

@author: marine

Purpose : use Profit to align two complexes from PDB files. Require six argument :
    - dir : the pdb file directory
    - interolog : the PDB directory of the interolog complex
    - ligand : the PDB file of the ligand
    - receptor : the PDB file of the receptor
    - fasta : the fasta directory
Output : a pdb file of the ligand in its initial position
Use (example) : python alignInterolog.py -ligand test_data/1AY7_l_sep.pdb -receptor test_data/1AY7_r.pdb -interolog test_data/2za4.pdb -chainLigand B -chainReceptor A
'''

import os
import sys
import subprocess
import askInterEvol
import filedownload
   
    
def runProfit (ligand, receptor, interolog, list_interEvol):
    '''
    use Profit to align and fit the pdb structures
    input :
        - ligand : the ligand of the complex to solve
        - receptor : the static receptor of the complex to solve
        - interolog : the interolog complex structure 
        - chainLigand : the name of the corresponding ligand in the interolog complex
        - chainReceptor : the name of the corresponding receptor in the interolog complex
    output : a pdb file of the initial position of the complex to solve
    '''    
    
    #move to the pdb file directory
    #os.chdir(dir) useless with the absolute path
    

    #name the future aligned pdb files
    intermediate = interolog.replace(".pdb","_aligned.pdb")
    output = ligand.replace(".pdb", "_aligned.pdb")
    
    #write the script used by profit
    with open("profit_script", "w") as script :    
    
        script.write("REFERENCE "+receptor+"\nMOBILE "+interolog+"\nALIGN "+chainReceptor+"*:A*\nITERATE\nFIT\nWRITE "+intermediate+"\nREFERENCE "+intermediate+"\nMOBILE "+ligand+"\nALIGN B*:"+chainLigand+"*\nITERATE\nFIT\nWRITE "+output)
    
    script.close()
    
    #Configure Profit variables
    proFit=subprocess.check_output("find ~ -type d -name .\* -prune -or -name ProFitV3.1 -print", shell=True)[:-1]
    env = os.environ.copy()
    env["HELPDIR"] = proFit
    env["DATADIR"] = proFit
	
    #launch profit
    #subprocess.call("profit < profit_script",shell=True)
    script = open("profit_script","r")
    p=subprocess.Popen(["profit"],stdin=subprocess.PIPE,stdout=sys.stdout,stderr=sys.stderr,env=env)
    p.communicate(input=script.read()+"QUIT")
    #p.communicate(input=script.encode(encoding="UTF-8"))
    
    return output

    

    
   
def runPymolAlignment(ligand, receptor, interolog, chainLigand, chainReceptor):
    '''
    use pymol to align the ligand on the interolog ligand
    input :
        - ligand : the ligand of the complex to solve
        - receptor : the static receptor of the complex to solve
        - interolog : the interolog complex structure 
        - chainLigand : the name of the corresponding ligand in the interolog complex
        - chainReceptor : the name of the corresponding receptor in the interolog complex
    output : a pdb file of the initial position of the complex to solve
    '''   
    
    #move to the pdb file directory
    #PDBDirectory = os.path.dirname(ligand)
    #os.chdir(PDBDirectory)
    
	#name the alignment actors
    output = ligand.replace(".pdb", "_aligned.pdb")
    ligandName = os.path.basename(ligand).replace(".pdb","")
    receptorName = os.path.basename(receptor).replace(".pdb","")
    interologName = os.path.basename(interolog).replace(".pdb","")
    interologRec = interologName.replace(".pdb","_rec")
    interologLig = interologName.replace(".pdb","_lig")
    
    
    #move to the pdb file directory
    #os.chdir(dir)
    
    #write the script used by profit
    with open("pymol_script.pml", "w") as script :    
    
        script.write("load "+interolog+", "+interologName+"\nload "+ligand+", "+ligandName+"\nload "+receptor+", "+receptorName+"\nselect "+interologRec+", "+interologName+" and chain "+chainReceptor+"\nselect "+interologLig+", "+interologName+" and chain "+chainLigand+"\ncealign "+receptorName+", "+interologRec+"\ncealign "+interologLig+", "+ligandName+"\nselect tosave, "+ligandName+"\nsave "+output+", tosave\nquit")
    
    script.close()

	#launch pymol_script
    subprocess.call("pymol -c pymol_script.pml", shell=True)
	
    return output





    
if __name__ == '__main__':
    '''
    obsolete
    try :
        dir = sys.argv[sys.argv.index("-dir")+1]
    except :
        print("ERROR : specified directory does not exist\n")
    '''
    try :
        lig = sys.argv[sys.argv.index("-ligand")+1]
    except :
        print("ERROR : specified ligand file does not exist\n")
        
    try:
        rec = sys.argv[sys.argv.index("-receptor")+1]
    except :
        print("ERROR : specified receptor file does not exist\n")
        
    try :
        inter = sys.argv[sys.argv.index("-interolog")+1]
    except :
        print("ERROR : specified interolog directory does not exist\n")
        
    try :
        chainRec = sys.argv[sys.argv.index("-chainReceptor")+1]
    except :
        print("ERROR : unspecified interolog receptor chain\n")
        
    try :
        chainLig = sys.argv[sys.argv.index("-chainLigand")+1]
    except :
        print("ERROR : unspecified interolog ligand chain\n")    
        
    runPymolAlignment(lig, rec, inter, chainLig, chainRec)
    
'''
Try with interEvol

    try :
        fasta = sys.argv[sys.argv.index("-fasta")+1]
    except :
        print("ERROR : specified fasta directory does not exist\n")

    ligf = filedownload.FASTAfromPDB(lig, dir, fasta)
    recf = filedownload.FASTAfromPDB(rec, dir, fasta)
    
    dico = askInterEvol.runAlign(recf, ligf, inter)
    for key in dico.keys():
        print(key)
        print(dico[key])
        liste = dico[key]
        print(liste[3],liste[2])
        runProfit(dir, lig, rec, os.path.join(inter,key+".pdb"), liste[3], liste[2])     
        
'''

    
        
                
