import dockerasmus.pdb as pdb

def iter_chain_atoms(self):
    """
    Yield every atom and its chain's id in ``self``. Similar to self.iteratoms()
    """
    for chain in self.itervalues():
        for residue in chain.itervalues():
            for atom in sorted(residue.itervalues(), key=lambda a: a.id):
                yield atom, chain.id

pdb.Protein.iter_chain_atoms = iter_chain_atoms

def write_atoms(self,path):  
    """
    A naive protein's atoms' records printer.

    Keyword arguments:
    path -- the output file's path
    """
    atom_template = '{:7}{:4d}  {:4}{:3} {:1} {:3d}      {:5.3f}  {:5.3f}  {:5.3f}'
    out = open(path,'w')
    for atom, c in self.iter_chain_atoms():
        print(atom_template.format("ATOM",atom.id,atom.name,atom.residue.name,c,atom.residue.id,atom.x,atom.y,atom.z), file=out)
    print("END", file=out)
    out.close()


def trans_rotate(self,trans=None,rotate=None):
    """
    The function takes translateion and rotation's parameters, retruns atoms' positions translated and rotated matrix.
    """
    pass


pdb.Protein.write_atoms = write_atoms
pdb.Protein.trans_rotate = trans_rotate
