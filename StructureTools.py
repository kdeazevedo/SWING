import os
import string
import math
import numpy
from path import path
from itertools import islice

def distancePoints((x1,y1,z1),(x2,y2,z2)):
    """
    Computes the distance between the two sets of coordinates
    input: 2 tuples with the corresponding coordinates 
    output: distance
    """
    x = (x1-x2)
    y = (y1-y2)
    z = (z1-z2)
    return math.sqrt(x*x+y*y+z*z)


def centerMassOfResidue(dPDB):
    """
    Calculates the center of mass of each residue contained in dPDB
    Input : dPDB a dictionary
    Output : dPDB with XCM, YCM and ZCM the coordinates of the center of mass of each residue
    """

    reslist = dPDB["position"]
    
        
    for res in reslist :        
        x = y = z = 0.0
        
        # looping over the current residue atoms
        for atom in dPDB[res]["atome"] :
            x += dPDB[res][atom]["x"]
            y += dPDB[res][atom]["y"]
            z += dPDB[res][atom]["z"]
        
        #Calculation of the center of mass of the current residue    
        Xcm = float(x)/len(dPDB[res]["atome"]) 
        Ycm = float(y)/len(dPDB[res]["atome"])
        Zcm = float(z)/len(dPDB[res]["atome"])
        dPDB[res]["XCM"] = Xcm
        dPDB[res]["YCM"] = Ycm
        dPDB[res]["ZCM"] = Zcm
        

def computeDist_dico(d_res1, d_res2, mode) :
    """
    Compute the distance between to residues
    Input : d_res1, d_res2 are dictionary corresponding to residue 1 and residue 2 respectively, mode atom or center
    Output : the distance between d-res1 and d_res2
    """

	#Compute the distance between atoms of a couple of residues
    if mode == "atom" :
        minval = 1000000
        for atom1 in d_res1["atome"] :
            coord1 = [d_res1[atom1]["x"], d_res1[atom1]["y"], d_res1[atom1]["z"]]
            for atom2 in d_res2["atome"] :
                coord2 = [d_res2[atom2]["x"], d_res2[atom2]["y"], d_res2[atom2]["z"]]
                dist = distancePoints((coord1[0], coord1[1], coord1[2]),(coord2[0],coord2[1], coord2[2]))
                if minval > dist :
                    minval = dist

	#Computes the distance between the CM of the 2 given residues
    elif mode == "center" : 
        dPDBtmp = {}
        dPDBtmp["position"] = ["res1", "res2"]
        dPDBtmp["res1"] = d_res1
        dPDBtmp["res2"] = d_res2

        centerMassOfResidue(dPDBtmp)
        minval = distancePoints((dPDBtmp["res1"]["XCM"],dPDBtmp["res1"]["YCM"],dPDBtmp["res1"]["ZCM"]),(dPDBtmp["res2"]["XCM"],dPDBtmp["res2"]["YCM"],dPDBtmp["res2"]["ZCM"]))
        
    return minval


def scorelist(listFiles, listScores, filename):
	"""
	Create in a file a ranking of the scores from best (minimum value) to worst (maximum value) with the correspond file number in front
    Input: the list of the files and the list of the scores, as well as the name wished for the file
    Output: the ranking in a file
    """
	dir_path = os.path.dirname(os.path.realpath(__file__))
	pathname=dir_path+"/scoring_Cornell"
	
	try:
		if not os.path.exists(pathname):
			os.makedirs(pathname)
	except:
		print "ERROR : You don't seem to have the rights to write in this directory"
		sys.exit()
		
	filepath = os.path.join(pathname, filename)
	fileid = open(filepath, 'w+')
	
	xarray= numpy.array(listFiles)
	yarray= numpy.array(listScores)
	
	data= numpy.array([xarray,yarray])
	data= data.T
	data= data[data[:, 1].argsort()]
	
	numpy.savetxt(fileid,data,fmt=['%d','%d'])
	fileid.close()
	
	return data[0][0]


def writePDB(dPDBrec, dPDBlig, prediction) :
	"""
	Write the PDB file of the predicted complex
	Input: dictionnaries of the receptor and solution ligand, as well as the name of the PDB file
	"""
	
	pred = open(prediction, "w")
	
	for chain in dPDBrec["nchaine"]:
		for res in dPDBrec[chain]["position"] :
			for atom in dPDBrec[chain][res]["atome"] :
				pred.write("ATOM  %5d  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X \n"%(dPDBrec[chain][res][atom]["ID"], atom, dPDBrec[chain][res]["residu"],chain, res,dPDBrec[chain][res][atom]["x"], dPDBrec[chain][res][atom]["y"],dPDBrec[chain][res][atom]["z"] ))
	
	for chain in dPDBlig["nchaine"]:
		for res in dPDBlig[chain]["position"] :
			for atom in dPDBlig[chain][res]["atome"] :
				pred.write("ATOM  %5d  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  1.00 X X\n"%(dPDBlig[chain][res][atom]["ID"], atom, dPDBlig[chain][res]["residu"],chain, res,dPDBlig[chain][res][atom]["x"], dPDBlig[chain][res][atom]["y"],dPDBlig[chain][res][atom]["z"] ))
                    
	pred.close()

def writeInterfacePDB(dPDBcompl, prediction) :
	"""
	Add the Bfactor in the PDB file of the predicted complex
    Input: dictionnary containing the Bfactor as well as the name of the new file
    """
    
	pred = open(prediction, "w")
	
	for chain in dPDBcompl["nchaine"]:
		for res in dPDBcompl[chain]["position"] :
			for atom in dPDBcompl[chain][res]["atome"] :
				pred.write("ATOM  %5d  %-4s%3s %s%4s    %8.3f%8.3f%8.3f  1.00  %d X X\n"%(dPDBcompl[chain][res][atom]["ID"], atom, dPDBcompl[chain][res]["residu"],chain, res,dPDBcompl[chain][res][atom]["x"], dPDBcompl[chain][res][atom]["y"],dPDBcompl[chain][res][atom]["z"], dPDBcompl[chain][res]["bfactor"] ))
	
	pred.close()

def extract100Bests(indir):
	"""
	Extract the 100 bests scores following the first scoring method to apply the second scoring method to them
	Input: the directory containing the ranking of the scores
	Output: a list containing the number of the best files, a list containing the name of the best files, a list containing the best scores
	"""
	infile = open("scoring_Cornell/Scoring1.txt")
	head = list(islice(infile, 100))
	
	listBestHits=[]
	listBestFiles=[]
	listBestScores=[]
	for i in head:
		num, score = i.split(" ")
		num=str(num)
		number=int(num)
		score=score.strip('\n')
		score=float(score)
		listBestHits.append(indir+"/1BRS_A_1BRS_B_allatom_"+num+"_DP.pdb")
		listBestScores.append(score)
		listBestFiles.append(number)
	
	infile.close()
	return listBestHits, listBestScores, listBestFiles
