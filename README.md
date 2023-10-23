# PROTACable
![Project Logo](./assets/logo/project-logo.png)
![Project Scheme](./assets/logo/project-scheme.png)

PROTACable is an end-to-end in-silico design toolkit for novel PROTACs

## Table of Contents

- [Preliminary](#preliminary)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Stage 1: Docking POI and Elaborating Variations](#stage-1-docking-poi-and-elaborating-variations)
- [Stage 2: POI-E3 Docking and Pose Filtering](#stage-2-poi-e3-docking-and-pose-filtering)
- [Stage 3: Linker Ligation and Pose Filtering](#stage-3-linker-ligation-and-pose-filtering)
- [Stage 4: SE(3) Transformer Network Score Prediction](#stage-4-se-3-transformer-network-score-prediction)
- [Stage 5: Ternary Complex Minimization](#stage-5-ternary-complex-minimization)
- [References](#references)
- [License](#license)

## Preliminary

## Requirements

Different stages of PROTACable pipeline have variable set of requirements.
### Conda environment:
* See [Installation](#installation) below.
### Linux OS:
* This toolkit collection had been tested on Centos 8 OS.
* SLURM scheduler is needed for Stage 2 (protein-protein docking) due to the way the code was implemented.
### Software:
* Stage 1: GNINA v1.0 (optional).
* Stage 2: ProPOSE v2022.
* Stage 5: Maestro Schrödinger (Preferred to have).
### Hardware:
* CUDA >= 11.0 (optional)


## Installation

To get started with PROTACable, follow these instructions:

```
git clone https://github.com/giaguaro/PROTACable.git 
cd PROTACable
export PROTACable=${pwd}
```
> **Note:** You need to have the environment variable defined at all times during running any part of the code.

### Conda environment:

```
conda env create -f PROTACable.yml
```
> [Follow NVIDIA's instructions](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/#axzz4TWipdwX1) to install the latest version of CUDA (>= 11.0 is required). Make sure `nvcc` is in your PATH.

### Software:

- **GNINA**: Refer to GNINA official repository - https://github.com/gnina/gnina
  
Once GNINA is installed place it in the ```PROTACable_stage_1/utilities``` directory.
 
- **ProPOSE**: Refer to ProPOSE official website for academic licensing - https://mm.nrc-cnrc.gc.ca/ccbpub/propose_main.php

Once ProPOSE software and license is obtained place them in the ```PROTACable_stage_2/ProPOSE``` directory.

*Preferred*

- **Maestro**: Refer to Maestro official website for download and tokens - https://www.schrodinger.com/products/maestro

Only Maestro's PrepWizard will be used and will be used only in Stage 5. This tool will serve the purpose of minimization of the selected ternary complexes modelled throughout the pipeline. 



## Usage

Due to the complexity of the pipeline, each stage was kept operating separately such that there is a finer control over the output before the subsequent stage is activated. The majority of the stages in the pipeline do more than one task. To get familiar with the tasks refer to our paper: [**doi link here**].

First, you need to define the PROTACable environment:

```
export PROTACable=/path/to/PROTACable/directory
```

 Then you need to validate installation by running:

```
sh $PROTACable/test_installation.sh
```

### Stage 1: Docking POI and Elaborating Variations

This stage is broken into two parts: docking and creating R-group variations for linker ligation.

+ You can optionally dock your POI with the cognate ligand (if it's not already complexed with the ligand pose):

```
sh $PROTACable/PROTACable_stage_1/dock_poi.sh -i <input PDB file> -l <input ligand in SDF format> -w <number of workers>
```

+ It is recommended to generate variations of the R-group for linker ligation using the docked PDB ligand (generated from above) or by extracting it from an existing PDB file of the POI-Lig complex. You can generated variations programmatically using the PDB file of the ligand and an exit vector. The exit vector here refers to the index which corresponds to the atom you want the linker to be attached to. 

```
python $PROTACable/PROTACable_stage_1/make_variations.py <ligand in PDB format> <exit vector> <output prefix>
```
> **Tip:** Once you decide on the atom you think is apt for linker ligation (usually an atom that is in maximally solvent exposed moiety), you may take note of the atom number (atom label). You can then find the corresponding index by opening the PDB file in a text viewer and locating the atom with the sought atom number. The exit vector will be exactly that atom's serial number. Refer to [PDB format guide](https://www.biostat.jhsph.edu/~iruczins/teaching/260.655/links/pdbformat.pdf) for further clarification.

This will generate two variations with carboxyl ```C(=O)O``` ${outout_prefix}_carboxyl.pdb and amide ```C(=O)N``` ${outout_prefix}_amide.pdb files as well the original unmodified ligand pose ${outout_prefix}_original.pdb.


```
```
from project_title import AwesomeFeature

# Initialize
feature = AwesomeFeature()

# Use the feature
result = feature.do_something_cool()
print(result)
\```

(Note: Remove the extra backticks around the code blocks above.)

## Features

### Feature 1

![Feature 1 Screenshot](https://path-to-your-image/feature1.png)

A brief description of this feature.

### Feature 2

Sub-headings can be used to describe sub-features or components:

#### Component 1
Details about this component.

#### Component 2
Details about this component.

## References

1. [Reference Title](https://www.example.com) - A short description of this reference.
2. Author Name, "Paper Title", Journal Name, Year. [Link to paper](https://www.example-paper.com)

## License

This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.
