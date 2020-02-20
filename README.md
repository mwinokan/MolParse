# AseMolPlot

Extension to ASE's plotting allowing for easy plotting and animations of trajectories with ray-tracing support. Intended for use on EUREKA.

## Installation

### Requirements

*   [ASE](#https://wiki.fysik.dtu.dk/ase/index.html)
*   [PoV-Ray](#https://github.com/POV-Ray/povray)
*   [MPyTools](#https://github.com/mwinokan/MPyTools)

### ASE on EUREKA

*   `pip install --upgrade --user ase`
*   `export PATH=$PATH:~/.local/bin` to your `.bash_profile`
*   `export PYTHONPATH=$PYTHONPATH:~/.local/lib/python3.7/site-packages` to your `.bash_profile`

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

### AseMolPlot

* `git clone https://github.com/mwinokan/AseMolPlot.git`
* Add `export MWAMPPATH=/path/to/directory` to your `.bash_profile`
* Add `export PYTHONPATH=$PYTHONPATH:$MWAMPPATH` to your `.bash_profile`

## Usage

The following examples have imported AMP as: `import asemolplot as amp`. See *test.py* for an example script.

*Use PoV-Ray to render a PNG from an atomic image*

`amp.makePovImage(filename,atoms,**style)`

*   `filename` String which will name the temporary and output files - i.e. filename.pov, filename.png
*   `atoms` An ASE Atoms object containing the molecule.
*   `**style` Parameters. Available options:
    -   `<int> canvas_width`
    -   `<float> radii`
    -   `<string> rotation`
    -   `<float> celllinewidth`
    -   `<bool> drawCell`

*Use PoV-Ray to render a PNG for every image in a trajectory*

`amp.makePovImages(filename,subdirectory="pov",interval=1,**style)`

*Use PoV-Ray to create an animation from a trajectory*

`amp.makePovAnimation(filename,subdirectory="pov",interval=1,**style)`

### Style templates for images

*amp.styles* contains several style templates that can be passed to makePovImage. For example: `amp.makePovImage(filename,atoms,**amp.styles.standard)`

*Custom styles*

Copy a style and modify one part:

```custom_style = amp.styles.standard
custom_style['celllinewidth'] = 0.07
```

