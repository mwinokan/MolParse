# AseMolPlot

Extension to ASE's plotting allowing for easy plotting and animations of trajectories with ray-tracing support. Intended for use on EUREKA.

Initialised at /users/mw00368/py/AseMolPlot on mw00368@login7.swmgmt.eureka

</newgit.sh>

## Installation

* `git clone https://github.com/mwinokan/AseMolPlot.git`
* Add `export MWAMPPATH=/path/to/directory` to your `.bash_profile`

### POVRAY on EUREKA

**!!! Please do this on a debug node !!!**

* `cd`
* `git clone https://github.com/POV-Ray/povray.git`
* `cd povray`
* `source $MWAMPPATH/load_pov.sh`
* `cd unix`
* `./prebuild.sh`
* `cd ../`
* `./configure --prefix=$HOME COMPILED_BY="Your Name your@email.com"`
* `make clean`
* `make install -j16`

## Usage

`import asemolplot as amp`