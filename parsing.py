from ForceField3 import chargePDB, epsilon_vdw_PDB
from RunScore import compEner
from path import path
import os

dcharge = chargePDB()
dvdw, depsilon = epsilon_vdw_PDB()

def parserPDB(infile):
	"""
	Parse a PDB file to extract the chains, the residues, the atoms, and 3D coordinates of the atoms
	Input : a PDB file
	Output : a dictionnary containing the variables of interest
	"""
	fichier = open(infile,"r")
	lines = fichier.readlines()

	d_proteine= {}
	d_proteine["nchaine"]=[]

	for i in lines:
		if i[0:4]=="ATOM": #To make sure it is the line of an atom
			
			nchaine=i[21] #Chain number
			#Making sure that the chain doesn't already exist to avoid overriding
			if not nchaine in d_proteine["nchaine"]:
				d_proteine["nchaine"].append(nchaine)
				d_proteine[nchaine]={}
				d_proteine[nchaine]["position"]=[]
			
			#For each chain, keep the residue
			position=i[22:26].strip()
			if not position in d_proteine[nchaine]["position"]:
				d_proteine[nchaine]["position"].append(position)
				d_proteine[nchaine][position]={}
				d_proteine[nchaine][position]["residu"]=i[17:20].strip()
				d_proteine[nchaine][position]["atome"]=[]
			
			#For each residue, the atom
			atome=i[12:16].strip()
			if not atome in d_proteine[nchaine][position]["atome"]:
				d_proteine[nchaine][position]["atome"].append(atome)
				d_proteine[nchaine][position][atome]={}
			
			#For each atom, the 3D Coord
			d_proteine[nchaine][position][atome]["x"]=float(i[30:38])
			d_proteine[nchaine][position][atome]["y"]=float(i[38:46])
			d_proteine[nchaine][position][atome]["z"]=float(i[46:54])
			d_proteine[nchaine][position][atome]["ID"]=int(i[6:11])
	fichier.close()	
	return d_proteine

def assignParams(dPDB, dcharge, dvdw, depsilon):
    '''
    Assign to a PDB dictionnary the charge, vdw, and epsilon
    Input : the PDB dictionnary, the dictionnary of the charge, vdw, and epsilon of the atom
    '''
    first = True
    
    for chain in dPDB["nchaine"]:
		for resi in dPDB[chain]["position"] :
        
			# means this is the Nter residue
			if first and (dPDB[chain][resi].has_key("H1") or dPDB[chain][resi].has_key("H2") or dPDB[chain][resi].has_key("H3")) :

					for atomi in dPDB[chain][resi]["atome"] :
						if atomi == "H1" :
							dPDB[chain][resi]["H1"]["charge"] = 0.1984
							dPDB[chain][resi]["H1"]["vdw"] = 0.6
							dPDB[chain][resi]["H1"]["epsilon"] = 0.0157


						elif atomi == "H2" :
							dPDB[chain][resi]["H2"]["charge"] = 0.1984
							dPDB[chain][resi]["H2"]["vdw"] = 0.6
							dPDB[chain][resi]["H2"]["epsilon"] = 0.0157

						elif atomi == "H3" :
							dPDB[chain][resi]["H3"]["charge"] = 0.1984   
							dPDB[chain][resi]["H3"]["vdw"] = 0.6
							dPDB[chain][resi]["H3"]["epsilon"] = 0.0157

						elif atomi == "N" :
							dPDB[chain][resi]["N"]["charge"] = 0.1592
							dPDB[chain][resi]["N"]["vdw"] = 1.875
							dPDB[chain][resi]["N"]["epsilon"] =  0.17

						elif atomi == "CA" :
							dPDB[chain][resi]["CA"]["charge"] = 0.0221
							dPDB[chain][resi]["CA"]["vdw"] = 1.9080
							dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

						elif atomi == "HA" :
							dPDB[chain][resi]["HA"]["charge"] = 0.116
							dPDB[chain][resi]["HA"]["vdw"] = 1.1
							dPDB[chain][resi]["HA"]["epsilon"] = 0.0157

						else:
							dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["residu"]][atomi]
							dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["residu"]][atomi]
							dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["residu"]][atomi]
                        
                    
					first = False
                
			# means this is the Cter residue   
			elif first == False and (dPDB[chain][resi].has_key("OXT")):

				for atomi in dPDB[chain][resi]["atome"] :

					if atomi == "CA" :
						dPDB[chain][resi]["CA"]["charge"] = -0.2493
						dPDB[chain][resi]["CA"]["vdw"] = 1.9080
						dPDB[chain][resi]["CA"]["epsilon"] = 0.1094

					elif atomi == "C" :
						dPDB[chain][resi]["C"]["charge"] = 0.7231
						dPDB[chain][resi]["C"]["vdw"] = 1.9080
						dPDB[chain][resi]["C"]["epsilon"] = 0.0860
                    
					elif atomi == "O" :
						dPDB[chain][resi]["O"]["charge"] = -0.7855
						dPDB[chain][resi]["O"]["vdw"] = 1.6612
						dPDB[chain][resi]["O"]["epsilon"] = 0.2100

					elif atomi == "OXT" :
						dPDB[chain][resi]["OXT"]["charge"] = -0.7855
						dPDB[chain][resi]["OXT"]["vdw"] = 1.6612
						dPDB[chain][resi]["OXT"]["epsilon"] = 0.2100
                    
					else:
						dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["residu"]][atomi]
						dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["residu"]][atomi]
						dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["residu"]][atomi]

        
			else :
				# for all the other residues
				for atomi in dPDB[chain][resi]["atome"] :

					# particular case of Histidine
					if dPDB[chain][resi]["residu"] == "HIS" and dPDB[chain][resi].has_key("HD1") :
						# means this is a HID
						dPDB[chain][resi]["residu"] = "HID"
                    
					elif dPDB[chain][resi]["residu"] == "HIS" and dPDB[chain][resi].has_key("HE2") :
						# means this is a HIE
						dPDB[chain][resi]["residu"] = "HIE"

					# general case                    
					dPDB[chain][resi][atomi]["charge"] = dcharge[dPDB[chain][resi]["residu"]][atomi]
					dPDB[chain][resi][atomi]["vdw"] =  dvdw[dPDB[chain][resi]["residu"]][atomi]
					dPDB[chain][resi][atomi]["epsilon"] =  depsilon[dPDB[chain][resi]["residu"]][atomi]

					


def preparePDB(infile) :
    """
    assigning all the atomic params (e.g. charges, vdw and epsilon) for all the atoms
    of the infile PDB
    Input: filename of the pdb
    Output: dico dPDB with 3D coords and charges, vdw and epsilon params for each atom
    """
    
    dPDB = parserPDB(infile)
    assignParams(dPDB, dcharge, dvdw, depsilon)

    return dPDB


def parseDirectory(indir, dico_Rec):
	"""
	Parse the directory contaning all the conformations
	Input : name of the directory, and dictionnary of the receptor
	Ouput : a list of all file numbers in the directory and a list of all scores computed from the files
	"""
	n=0
	listFiles=[]
	listScores=[]
	for files in os.listdir(indir):
		if files.startswith("1BRS") : #To avoid hidden files
			pathfile=os.path.join(indir,files)
			n+=1
			togo=949-n
			print "Dealing with file "+files
			print "Number of files remaning : "+str(togo)
			dico_confH = preparePDB(pathfile)
			a, b, c, d, e, f, g=files.split("_")
			f=int(f) #f correspond to the number of the file
			listFiles.append(f)
			listScores.append(compEner(dico_Rec, dico_confH, 'B', 'D'))
			dico_confH = {}
	
	return listFiles, listScores
