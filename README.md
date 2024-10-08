![MolParse](https://github.com/mwinokan/MolParse/blob/master/graphics/molparse-01.png?raw=true)

![GitHub Tag](https://img.shields.io/github/v/tag/mwinokan/molparse?include_prereleases&label=PyPI&link=https%3A%2F%2Fpypi.org%2Fproject%2Fmolparse%2F)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mwinokan/molparse/python-publish.yml?label=publish&link=https%3A%2F%2Fgithub.com%2Fmwinokan%2FMolParse%2Factions%2Fworkflows%2Fpython-publish.yml)
![GitHub Actions Workflow Status](https://img.shields.io/github/actions/workflow/status/mwinokan/molparse/black.yaml?label=lint&link=https%3A%2F%2Fgithub.com%2Fmwinokan%2FMolParse%2Factions%2Fworkflows%2Fblack.yaml)
[![Documentation Status](https://readthedocs.org/projects/hippo-db/badge/?version=latest)](https://molparse.winokan.com/en/latest/?badge=latest)
![GitHub last commit](https://img.shields.io/github/last-commit/mwinokan/molparse)
![GitHub Issues or Pull Requests](https://img.shields.io/github/issues/mwinokan/molparse)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

A python package for parsing, modifying, and analysis of molecular structure files.

## Installation

Easiest way to install is from [PyPI](https://pypi.org/project/MolParse/):

`pip install molparse`

## Usage

### Python Module

MolParse is primarily a python module which can be used interactively, or within (batch) scripts:

Use `pydoc` to see help on the `molparse` module, or its methods & classes. E.g. from a shell:

`pydoc molparse`

`pydoc molparse.System`

### Binaries and command-line programs

#### moltree

In addition to the python module, an interactive command-line interface is available with the binary `moltree`. Pass a
PDB or GRO file as follows:

`moltree <FILE>`

![moltree](https://github.com/mwinokan/MolParse/blob/master/graphics/moltree.png?raw=true)

Use the mouse to interact with buttons and CTR-C to exit.

#### molxvg

Gromacs produces data files in XVG format by default, these can be parsed using the `molparse.xvg.parseXVG` method from
within a python environment, alternatively a binary exists to access its basic functionality from the command line. Run
the following to open an interactive plotly graph of an xvg:

`molxvg [FILE.xvg] -s`

![moltree](https://github.com/mwinokan/MolParse/blob/master/graphics/molxvg.png?raw=true)

Other options can be found by running `molxvg --help`.

## Installation from source

### Requirements

* [ASE](#https://wiki.fysik.dtu.dk/ase/index.html)
* [MPyTools](#https://github.com/mwinokan/MPyTools)

### Installing ASE

* `pip install --upgrade --user ase`
* `export PATH=$PATH:~/.local/bin` to your `.bash_profile`
* `export PYTHONPATH=$PYTHONPATH:~/.local/lib/python3.X/site-packages` to your `.bash_profile` Where X is your python
  version.

### MPyTools

* `git clone https://github.com/mwinokan/MPyTools.git`
* Add `export MPYTOOLS=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MPYTOOLS` to your `.bash_profile`

### MolParse

* `git clone https://github.com/mwinokan/MolParse.git`
* Add `export MOLPARSE=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MOLPARSE` to your `.bash_profile`
