import dockerasmus.pdb as pdb

def write_atoms(self,path):  
    """
    A naive protein's atoms' records printer.

    Keyword arguments:
    path -- the output file's path
    """
    # TODO : atom format uncompleted. Atom's residue should be writen. 
    atom_template = '{:7}{:4d}  {:4}{:3} {:1} {:3d}      {:5.3f}  {:5.3f}  {:5.3f}'
    out = open(path,'w')
    for atom in self.iteratoms():
        print(atom_template.format("ATOM",atom.id,atom.name,"ARG","B",123,atom.x,atom.y,atom.z), file=out)
    print("END", file=out)
    out.close()

def trans_rotate(self,trans=None,rotate=None):
    """
    The function takes translateion and rotation's parameters, retruns atoms' positions translated and rotated matrix.
    """
    pass


pdb.Protein.write_atoms = write_atoms
pdb.Protein.trans_rotate = trans_rotate
