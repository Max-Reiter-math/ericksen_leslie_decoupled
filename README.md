# ericksen_leslie_decoupled
[//]: # "show preview in vs code by ctrl+shift+v"
FEniCS implementation of a decoupled numerical FEM scheme for the Ericksen-Leslie Equations equipped with Dirichlet Energy:

In this package we consider two types of ericksen-leslie models. One is a rather simplified one (submodel = "simple") while the second one is of rather general nature (submodel = "general"). Both are equipped with the Dirichlet energy.

We offer a numerical schemes which decouples almost all (only velocity and pressure are coupled) quantities and rely on a fixpoint iteration. The goal of this project is to make the implementation of an academically used algorithm transparent and accessible. 
Mathematical background about the algorithm and its application will later be explained in an according paper.

## Table of contents
* [Installation](#installation)
* [Usage](#usage)
* [File structure](#file-structure)
* [Citation](#citation)

## Installation

Clone the repository
```sh
git clone https://github.com/Max-Reiter-math/ericksen_leslie_decoupled.git
```

### Dependencies
FEniCS

Matplotlib

numpy

## Usage

sim.py allows the execution of a single simulation via command line inputs. For details about the necessary and optional arguments, simply type 
```
python sim.py -h
```
example usage:
```
python sim.py -m decoupled_fp_solver -e annihilation_2 -s "simple" -d 3 -dt 0.01
```

## File structure:
[//]: # "type tree /F into windows terminal"

```bash
│   CITATION.cff
│   example_config.json
│   lib.py
│   README.md
│   run_experiment.py
│   sim.py
│   __init__.py
│
├───accessories
│   │   communicate.py
│   │   logger.py
│   │   my_plots.py
│   │   postprocess.py
│   │   __init__.py
│
├───experiments
│   │   annihilation_2.py
│   │   smooth_2d.py
│   │   velocity_driven_flow_2.py
│   │   __init__.py
│ 
├───model_library
│   │   decoupled_fp_solver.py
│   │   __init__.py
│   │
│   ├───basic_models
│   │   │   assist_funcs.py
│   │   │   basemodel.py
│   │   │   basemodel_fp.py
│   │   │   basemodel_linear_decoupled.py
│   │   │   basemodel_linear_fp_decoupled.py
│   │   │   __init__.py
```

## Citation
For citation we refer to the CITATION.cff file in this directory.