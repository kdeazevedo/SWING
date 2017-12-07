# 2017-2018_Equipe10, Sampling
Template based protein docking written in Python3. The program is developed by Kévin De Azevedo, Marine Duhamel, and Hua-Ting Yao (Paris 11)
## Abstract
SWING stands for Sampling With INteroloG and aims at producing ligand conformations for complex structure identification. Our program use the [InterAlign](http://biodev.cea.fr/interevol/) website in order to find interologs of the complex whose the structure is to identify. Interologs are used to determine a starting point around which the program will sample.

## Installation
Python dependent packages can be found in `requirements.txt`.
Use `pip3` to install.
```
pip3 install -r requirements.txt
```
### The third party software
`Minimiser`, `Pymol` (https://pymol.org/2/) and `ClusCo` (https://bitbucket.org/mjamroz/clusco/overview) are used in this program.

Please copy `minimiser` folder (from [here](https://github.com/meetU-MasterStudents/2017-2018_partage/tree/master/Codes/Minimizer)) into `Minimiser` folder.

Then, install `Pymol` following instructions from its official page ([here](http://www.pymol.org/install))

Finally, download `ClusCo` from [here](https://bitbucket.org/mjamroz/clusco/downloads/) and follow instructions from [here](https://bitbucket.org/mjamroz/clusco/overview) for installation. Beware, you must also install tclap and cmake.

For ubuntu users :
```
sudo apt-get install libtclap-dev
sudo apt-get install cmake
```
For mac users :
```
brew install tclap
brew install cmake
```

## Quick start:
To run the program from A to Z :
```bash
python3 main.py run -rec test/1AY7_r.pdb -lig test/1AY7_l_sep.pdb -n 10 --minimizer
```
However, for a given pair of receptor and ligand, request from [InterEvol](http://biodev.cea.fr/interevol/) and pymol alignement could only be executed once for ever. Hence, three parts of program can be run separately with proper positional arguments.

The command above is equivelent to the followings :
```bash
python3 main.py download -rec test/1AY7_r.pdb -lig test/1AY7_1_sep.pdb
python3 main.py align -rec test/1AY7_r.pdb -lig test/1AY7_l_sep.pdb -c out/Inter/Inter.conf
python3 main.py samples -rec test/1AY7_r.pdb -lig test/1AY7_l_sep.pdb -n 10 -c out/Inter/Samples.conf --minimizer
```

For more arguments' detail:
```bash
python3 main.py -h
```
