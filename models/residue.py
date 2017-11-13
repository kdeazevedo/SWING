import dockerasmus.pdb as pdb

@property
def name(self):
    if self._name is not None:
        return self._name

pdb.Residue.name = name
