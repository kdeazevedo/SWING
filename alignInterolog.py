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
Use (example) : python alignInterolog.py -dir pdb_input -ligand 1AY7_l_sep.pdb -receptor 1AY7_r.pdb -interolog absolutepath -fasta absolutepath
'''
import os
import sys
import subprocess
import askInterEvol
import filedownload
   
    
def runProfit (dir, ligand, receptor, interolog, chainLigand, chainReceptor):
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
    proFit=subprocess.check_output("find ~ -type d -name .\* -prune -or -name ProFitV3.1 -print", shell=True)
    env = os.environ.copy()
    env["HELPDIR"] = proFit
    env["DATADIR"] = proFit
	#print(env)
	

    #launch profit
    #subprocess.call("profit < profit_script",shell=True)
    script = open("profit_script","r")
    p=subprocess.Popen(["profit"],stdin=subprocess.PIPE,stdout=sys.stdout,stderr=sys.stderr,env=env)
    p.communicate(input=script.read()+"QUIT")

    

    
if __name__ == '__main__':
    
    try :
        dir = sys.argv[sys.argv.index("-dir")+1]
    except :
        print("ERROR : specified directory does not exist\n")
    
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
        
