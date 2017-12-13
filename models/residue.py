#Modified dockerasmus.pdb.residue.py
#Original Author: Martin Larralde (althonos)
#Modified by Hua-Ting Yao ([anthony]htyao)
import dockerasmus.pdb as pdb

# Overwrite dockerasmus.pdb.Residue name function
@property
def name(self):
    if self._name is not None:
        return self._name

pdb.Residue.name = name
