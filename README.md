# AseMolPlot

Extension to ASE's plotting allowing for easy plotting and animations of trajectories with ray-tracing support. Intended for use on EUREKA.

Initialised at /users/mw00368/py/AseMolPlot on mw00368@login7.swmgmt.eureka

</newgit.sh>

## Installation

* `git clone https://github.com/mwinokan/AseMolPlot.git`
* Add `export MWAMPPATH=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MWAMPPATH` to your `.bash_profile`

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

The following examples have imported AMP as: `import asemolplot as amp`. See *test.py* for an example script.

Current functionality that can be used:

### Use PoV-Ray to render a PNG from an atomic image

`amp.makePovImage(filename,atoms,**parameters)`

*   `filename` String which will name the temporary and output files - i.e. filename.pov, filename.png
*   `atoms` An ASE Atoms object containing the molecule.
*   `**parameters` Parameters. Available options:
    -   `canvas_width`
    -   `radii`
    -   `rotation`
    -   `celllinewidth`

### Style templates for images

*amp.styles* contains several style templates that can be passed to makePovImage. For example: `amp.makePovImage(filename,atoms,**amp.styles.standard)`