# PyAPM: Python Automated Polymer Modeler for LAMMPS

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

PyAPM is an open-source Python tool for automated generation of amorphous polymer models from SMILES strings. It supports both **all-atom (AA)** and **coarse-grained (CG)** modeling, enabling multiscale simulations of polymers with LAMMPS.

**Key capabilities**:
- Full automation from monomer SMILES to LAMMPS input files
- Support for 8 polymer architectures: oligomers, homopolymers, random/alternating/block/sequential copolymers, branched polymers, and crosslinked networks
- Automated crosslinking reaction setup using AutoMapper
- Coarse-graining via Iterative Boltzmann Inversion (IBI) method
- Built-in DREIDING force field with optimized LJ parameters
- Modular Python design for easy extension

## Workflow Overview

PyAPM consists of six core modules:

| Module | Description |
|--------|-------------|
| `ChainBuilder` | Converts SMILES of repeating units to polymer SMILES |
| `MoleculeBuilder` | Builds single‑chain structure from polymer SMILES |
| `AtmosphereBuilder` | Packs multiple chains into amorphous box (PACKMOL + Moltemplate) |
| `CrosslinkBuilder` | Handles crosslinking reactions (AutoMapper + LAMMPS fix bond/react) |
| `IBIBuilder` | Maps AA trajectories to CG models via IBI |
| `CGBuilder` | Generates CG polymer systems with prescribed bond/angle distributions |

## Installation

### Requirements
- Python 3.8 or higher
- LAMMPS (2 Aug 2023 – Update 1) with REACTION, KSPACE, MOLECULE packages
- PACKMOL, Moltemplate, AutoMapper (see below)
- Python packages: numpy, rdkit-pypi, openbabel, MDAnalysis, scipy, pandas, natsort

### Step 1: Install Anaconda and Python packages
```bash
conda install numpy
conda install -c conda-forge rdkit openbabel mdanalysis pandas natsort
```
### Step 2: Install external tools
Moltemplate:
`git clone https://github.com/jewettaij/moltemplate` and follow its README.
```bash
pip install moltemplate
```

AutoMapper:
An improved version is included in the AutoMapper/ folder of PyAPM. (Original: https://github.com/m-bone/AutoMapper)

PACKMOL:
Download from https://m3g.github.io/packmol/, compile with make, and add to PATH.

Such as:
```bash
export PACKMOL_EXEC=/home/python/packmol-20.14.4-docs1/packmol
```
Or define in '~/.bashrc'

### Step 3: Install LAMMPS

Download LAMMPS 2 Aug 2023 – Update 1 from https://www.lammps.org/download.html.
Compile with `-D PKG_REACTION=yes -D PKG_KSPACE=yes -D PKG_MOLECULE=yes`.

Ensure all executables (packmol, moltemplate.sh, cleanup_moltemplate.sh, lmp_mpi or lmp) are in your PATH or set appropriate environment variables.

### Step 4: Install PyAPM

```bash
git clone https://github.com/Tianyi_Wang/PyAPM.git
cd PyAPM
pip install -e .
```

## Quick Start

Here is a minimal example to build a polyethylene homopolymer (10 repeat units) and generate LAMMPS input files:

```python
from PyAPM import build_homopolymer
from pkg_resources import resource_filename
import os
import shutil

if __name__ == "__main__":
  os.mkdir('test1')
  dreiding_path = resource_filename('PyAPM', 'dreiding.lt')
  shutil.copy(dreiding_path,'test1')
  os.chdir('test1')
  smi = ["CC(*)(C*)C(=O)OC12CC3CC(CC(C3)C1)C2"]
  n_moltype = 1
  n_mol = [10]
  length = [25]
  box=[50,50,50]
  build_homopolymer(smi, n_moltype, n_mol, box, length)
```

```bash
source /home/python/anaconda3/etc/profile.d/conda.sh
conda activate pyapm_env
export PACKMOL_EXEC=/home/python/packmol-20.14.4-docs1/packmol
python test_homopolymer.py
```


After execution, the following files are created:

    system.data – LAMMPS data file

    system.in.init and system.in.settings – LAMMPS input scripts

    Polymer.pdb – PDB file of the packed system

For crosslinked systems, please refer to the examples/ folder.

## Documentation

A full API documentation is under construction. For now, please refer to the docstrings in the source code and the `examples/` folder.

## Citation
If you use PyAPM in your work, please cite:

- **Software**: Tianyi Wang. PyAPM: Python Automated Polymer Modeler for LAMMPS. Version 1.0.0. GitHub: https://github.com/TYWangEcust/PyAPM.git

- **Paper** (to appear)

## Licence
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

