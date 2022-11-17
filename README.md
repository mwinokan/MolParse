# AseMolPlot

Extension to ASE's plotting allowing for easy plotting and animations of trajectories with ray-tracing support through PoV-Ray. Intended for use on EUREKA.

## Installation

### Requirements

*   [ASE](#https://wiki.fysik.dtu.dk/ase/index.html)
*   [MPyTools](#https://github.com/mwinokan/MPyTools)

### ASE on EUREKA

*   `pip install --upgrade --user ase`
*   `export PATH=$PATH:~/.local/bin` to your `.bash_profile`
*   `export PYTHONPATH=$PYTHONPATH:~/.local/lib/python3.7/site-packages` to your `.bash_profile`

### AseMolPlot

* `git clone https://github.com/mwinokan/AseMolPlot.git`
* Add `export MWAMPPATH=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MWAMPPATH` to your `.bash_profile`

## Usage

try importing the module and using the builtin help method:

`import asemolplot as amp`
`help(amp)`
