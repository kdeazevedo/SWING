# 2017-2018_Equipe10, Sampling
Template based protein docking written in Python3. The program is developed by Kévin De Azevedo, Marine Duhamel, and Hua-Ting Yao (Paris 11)
## Abstract

## Installation
Python dependent packages can be found in `requirements.txt`.
Use `pip3` to install.
```
pip3 install -r requirements.txt
```
### The third party software
`Minimiser` and [ProFit](http://www.bioinf.org.uk/software/profit/index.html) are used in this program.

Please copy `minimiser` folder(from [here](https://github.com/meetU-MasterStudents/2017-2018_partage/tree/master/Codes/Minimizer)) into `Minimiser` folder.

Then, download `ProFit` from its official page ([here](http://www.bioinf.org.uk/software/swreg.html))
## Quick start:
To run the program from A to Z :
```bash
python3 main.py run -rec test/1AY7_r.pdb -lig test/1AY7_l_sep.pdb -n 10 --minimizer
```
However, for a given pair of receptor and ligand, request from [InterEvol](http://biodev.cea.fr/interevol/) and ProFit alignement could only be execuated once for ever. Hence, three parts of program can be run separately with proper positional arguments.

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
