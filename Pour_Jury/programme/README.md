#For the examiners
We chose to integrate our sampling method ([SWING](https://github.com/meetU-MasterStudents/2017-2018_Equipe10)) with the scoring methods of our choice, eiter [MeetDockOne](https://github.com/meetU-MasterStudents/2017-2018_Equipe1) or [DeNovo](https://github.com/meetU-MasterStudents/2017-2018_Equipe11) through a pipeline made with [SnakeMake](http://snakemake.readthedocs.io/en/latest/), a workflow management system used in life science. The configuration file (config.json) given with the pipeline let the user choose which scoring method they want to use with SWING, along with the possibility to change several parameters.
Details on how to install SnakeMake can be found on the official website. However, since MeetDockOne requires the use of a [conda](https://conda.io) environment, we highly recommend to install it via conda with the command line :

```
conda install -c bioconda snakemake
```

##Prerequisites

These prerequisites only concern the pipeline and do not take into account those specifically needed by SWING and the scoring programs. Please, be sure to first meet these prerequisites, available on the README.md in the respective Git repositories, before moving on.

###Conda environment

Because the MeetDockOne method make use of conda environment, our pipeline does too. The environment file required is not the same as the one supplied by the MeetDockOne team : it is available in this directory under the name environment.yml, and contains everything needed in order to run a complete docking with SWING and either scoring method.

To create the environment, please follow these instructions from the MeetDockOne readme :

#### 1- create the environment from the environment.yml file

`conda env create -f ./path/to/environment.yml`

#### 2- activate the environment

`source activate meetu`

#### 3- go to pipeline folder
`cd ./path/to/programme/`

###Recommended directory arborescence

One of the points of a SnakeMake pipeline is the fact that no specific directory arborescence is required : the program will
run in its entierety wherever the inputs and scripts are so long as their path are specified in the configuration file. However, we recommend you use the customary arborescence of a SnakeMake pipeline, which is the one already provided in this directory (name of files in input directory are examples) :

```
programme/
├── DeNovo/
│   └── ...
├── MeetDockOne/
│   └── ...
├── outputs
│   └── sampling
│   └── scoring
│── SWING/
│   └── ...
├── test/
│   ├── 1AY7_r.pdb
│   └── 1AY7_l_sep.pdb
│   └── Inter_1AY7_r.conf
│   └── Samples_1AY7_r.conf
│   └── ...
│── blast.conf
│── concatenate.py
│── config.json
│── environment.yml
│── parameters.conf
│── Snakefile
```

The `outputs` directory is created durin the deployment of the pipeline. For more informations on the `blast.conf` and `parameters.conf` files, please see the wiki of the SWING repository [here](https://github.com/meetU-MasterStudents/2017-2018_Equipe10/wiki/Parameters-of-InterEvol).

##Running the pipeline

After setting up the environment in the directory where the Snakefile is, the user just needs to write the command `snakemake` in a terminal to deploy the pipeline. With everything provided in the directory, the user is able to run the full docking program on the 1AY7 complex. If they want to re-deploy it again with the same complex, they should write `touch` instead. However, should the user wish to test it on other complex, they will have to modify the config.json file given with the directory. This file contains the data given in to the program in order for it to run. Here is an explanation of each term :
```
{
    "swingdir" : the directory of the SWING program
    "inputdir" : directory where all the input files are stored (must contain ligand and receptor files at least)
    "scoringdir" : the directory of the scoring program,
    "whichscoring" : which method of scoring you wish to use. Valid arguments : "MeetDockOne" or "DeNovo"
    "mode" : by which end the sampling the program will start. Valid arguments : "all", "align", "samples"
    "recfile" : the name of the receptor file
    "ligfile" : the name of the ligand file
    "interfile" : the name of the Inter.conf generated at the end of the download part
    "samplefile" : the name of the Samples.conf generated at the end of the alignment part
    "nconf" : number of confirmation given
    "recChain" : the chain of the receptor for the MeetDockOne scoring
    "ligChain" : the chain of the ligand for the MeetDockOne scoring
    "depth" : "naccess" or "msms" (see MeetDockOne scoring method)
    "ph" : pH (see MeetDockOne scoring method)
    "dist" : threshold for interface determination (see MeetDockOne scoring method)
}
```

To run the program with another complex, the receptor and ligand files **have** to be in the input directory. The mode argument basically tells the program which parts should be ommitted :
* "samples" omits the alignment and download part if you already have the corresponding Samples.conf post-alignment
* "align" omits the download part if you have already downloaded the interologs and have the corresponding Inter.conf
* "all" runs the entierety of the program
For more informations on what the different parts of the sampling entails, please read the README in the SWING git repository.
