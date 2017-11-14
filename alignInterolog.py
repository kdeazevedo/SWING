'''
Created on 13 nov. 2017

@author: marine

Purpose : use Profit to align two complexes from PDB files. Require six argument :
    - dir : the pdb file directory
    - interolog : the PDB file of the interolog complex
    - ligand : the PDB file of the ligand
    - receptor : the PDB file of the receptor
    - chainReceptor : the name of the receptor chain in the interolog complex
    - chainLigand : the name of the ligand chain in the interolog complex
Output : a pdb file of the ligand in its initial position
Use (example) : python alignInterolog.py -dir pdb_input -ligand 1AY7_l_sep.pdb -receptor 1AY7_r.pdb -interolog 2za4.pdb -chainLigand B -chainReceptor A
'''
import os
import sys
import subprocess
   
    
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
    os.chdir(dir)
    
    #name the future aligned pdb files
    intermediate = interolog.replace(".pdb","_aligned.pdb")
    output = ligand.replace(".pdb", "_aligned.pdb")
    
    #write the script used by profit
    with open("profit_script", "w") as script :    
    
        script.write("REFERENCE "+receptor+"\nMOBILE "+interolog+"\nALIGN "+chainReceptor+"*:A*\nITERATE\nFIT\nWRITE "+intermediate+"\nREFERENCE "+intermediate+"\nMOBILE "+ligand+"\nALIGN B*:"+chainLigand+"\nITERATE\nFIT\nWRITE "+output)
    
    script.close()
    
    #launch profit
    subprocess.call("profit", stdin="profit_script")

    
    
    
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
        print("ERROR : specified interolog file does not exist\n")
        
    try :
        chainLig = sys.argv[sys.argv.index("-chainLigand")+1]
    except :
        print("ERROR : unspecified ligand chain\n")
        
    try :
        chainRec = sys.argv[sys.argv.index("-chainReceptor")+1]
    except :
        print("ERROR : unspecified receptor chain\n")
    
    runProfit(dir, lig, rec, inter, chainLig, chainRec)
        
        
    
    
    
