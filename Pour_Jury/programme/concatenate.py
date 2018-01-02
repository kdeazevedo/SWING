import pathlib
import sys
import subprocess
import os

def concatenate(indir, recfile):
	for files in os.listdir(indir):
		ligfile=os.path.join(indir,files)
		os.system("cat "+recfile+" >> "+ligfile)

if __name__ == '__main__':
	
	try :
		indir = sys.argv[sys.argv.index("-dir")+1]
	except :
		print("ERROR : specified directory does not exist\n")
        
	try:
		rec = sys.argv[sys.argv.index("-rec")+1]
	except :
		print("ERROR : specified receptor file does not exist\n")
	
	concatenate(indir, rec)
