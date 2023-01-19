![MolParse](https://github.com/mwinokan/MolParse/blob/master/graphics/molparse-01.png?raw=true)

A python package for parsing, modifying, and analysis of molecular structure files. 

## Installation

### Requirements

*   [ASE](#https://wiki.fysik.dtu.dk/ase/index.html)
*   [MPyTools](#https://github.com/mwinokan/MPyTools)

### Installing ASE

*   `pip install --upgrade --user ase`
*   `export PATH=$PATH:~/.local/bin` to your `.bash_profile`
*   `export PYTHONPATH=$PYTHONPATH:~/.local/lib/python3.X/site-packages` to your `.bash_profile` Where X is your python version.

### MPyTools

* `git clone https://github.com/mwinokan/MPyTools.git`
* Add `export MWPYTPATH=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MWPYTPATH` to your `.bash_profile`

### MolParse

* `git clone https://github.com/mwinokan/MolParse.git`
* Add `export MWMPPATH=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MWMPPATH` to your `.bash_profile`

## Usage

try importing the module and using the builtin help method:

`import molparse`
`help(molparse)`
