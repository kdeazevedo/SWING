import dockerasmus.pdb as pdb

# Overwrite dockerasmus.pdb.Residue name function
@property
def name(self):
    if self._name is not None:
        return self._name

pdb.Residue.name = name
