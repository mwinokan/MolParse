# AseMolPlot

Extension to ASE's plotting allowing for easy plotting and animations of trajectories with ray-tracing support through PoV-Ray. Intended for use on EUREKA.

*Ciprofloxacin with ASE's PNG Renderer:*

![ASE PNG Example](https://github.com/mwinokan/AseMolPlot/blob/master/amp.png "Standard ASE PNG Renderer")

*Ciprofloxacin with ASE & PoV-Ray*

![ASE POV Example](https://github.com/mwinokan/AseMolPlot/blob/master/pov.png "ASE & PoV-Ray Render")

## Installation

### Requirements

*   [ASE](#https://wiki.fysik.dtu.dk/ase/index.html)
*   [PoV-Ray](#https://github.com/POV-Ray/povray)
*   [MPyTools](#https://github.com/mwinokan/MPyTools)
*   `makePovAnimation` requires: *imageio* `conda install imageio`

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

`amp.makePovImage(filename,image,verbosity=1,rmPovFiles=True,bonds=False,bondradius=1.1,forceLoad=False,**style)`

*   `filename` String which will name the temporary and output files - i.e. filename.pov, filename.png
*   `image` An ASE Atoms object containing the molecule.
*   `verbosity` Verbosity level. Every nested function call will pass `verbosity=verbosity-1`.
*   `rmPovFiles` Remove the PoV-Ray temp and input files?
*   `bonds` Boolean to draw bonds in the image.
*   `bondRadius` Floating radius within which to search for neighbouring (bonded) atoms.
*   `forceLoad` Force the loading of the PoV-Ray dependencies?
*   `**style` A dictionary of parameters. See `amp.styles` and [ase.io.pov](#https://wiki.fysik.dtu.dk/ase/_modules/ase/io/pov.html#write_pov). Available options:
    -   `<int> canvas_width` Required.
    -   `<float> radii` Relative scaling of rendered atoms.
    -   `<string> rotation` Rotation sting: i.e. "-90x,65z"
    -   `<bool> drawCell` Draw the unit cell?
    -   `<float> celllinewidth` Width of the unit cell edges. 

*Use PoV-Ray to render a PNG for every image in a trajectory*

`amp.makePovImages(filename,subdirectory="pov",interval=1,verbosity=1,rmPovFiles=True,bonds=False,bondradius=1.1,filenamePadding=4,forceLoad=False,**style)`

*   `filename` File from which to read the system's trajectory.
*   `subdirectory` Subdirectory within which to save the generated images.
*   `interval` Create an image for every `n` timesteps.
*   `verbosity` Verbosity level. Every nested function call will pass `verbosity=verbosity-1`.
*   `rmPovFiles` Remove the PoV-Ray temp and input files?
*   `bonds` Draw bonds in the image?
*   `bondRadius` Floating radius within which to search for neighbouring (bonded) atoms.
*   `filenamePadding` Pad the filenames of the outputs up to `n` characters. i.e. with `n=4`: "0001.png"
*   `forceLoad` Force the loading of the PoV-Ray dependencies?
*   `**style` A dictionary of parameters. See `amp.styles`. All parameters from makePovImage() are included.

*Use PoV-Ray to create an animation from a trajectory file*

`amp.makePovAnimation(filename,subdirectory="pov",interval=1,gifstyle=styles.gif_standard,verbosity=1,forceLoad=False,**plotstyle)`

*   `filename` File from which to read the system's trajectory.
*   `subdirectory` Subdirectory within which to save the generated images.
*   `interval` Create an image for every `n` timesteps.
*   `gifstyle` A dictionary of parameters for making the gif. This is passed to `imageio.mimsave()`.
    -   `<string> background` Currently only "white" is supported. No specification is transparent.
    -   `<int> fps` Frames per second of animated gif.
    -   `<0 or 1> loop` Loop the gif?
*   `**style` A dictionary of parameters. See `amp.styles`. All parameters from makePovImage() with the addition of:
    -   `canvas_height` Set the height of the output images (this will be achieved by cropping.)
    -   `crop_xshift` Crop from `n` pixels in from the top left.
    -   `crop_yshift` Crop from `n` pixels in from the top left.

*Wrappers with verbosity control for ase.io.read and ase.io.write*

*amp* contains functions `read()` and `write()` which call the relevant ASE functions, except that they include a verbosity control. Useful for example when user output is desired when loading/writing large structure files.

`amp.write(filename,image,verbosity=1,**parameters)`
`amp.read(filename,index=None,verbosity=1,**parameters)`

### Style templates for images

*amp.styles* contains several style templates that can be passed to makePovImage. For example: `amp.makePovImage(filename,atoms,**amp.styles.standard)`

*Custom styles*

Copy a style and modify one part:

```custom_style = amp.styles.standard.copy()
custom_style['celllinewidth'] = 0.07
```

### Common Errors

Any PoV-Ray plotting function (`makePov*()`):

*   `OSError: Povray command povray after.ini 2> pov.log failed with error code 32512`
    -   Cause: The `povray` executable was called without the correct environment modules. This often happens if you call a `MakePov___()` function (loading the correct modules), then load/unload some modules for another package, and then call an AMP POV-Ray function again.
    -   Fix: call `amp.loadPov()` before making images again. Alternatively add `forceLoad=True` to `MakePov___()` function calls.

PoV-Ray animations (`makePovAnimation()`):

*   `ValueError: operands could not be broadcast together with shapes (211,300,4) (184,300,4) `
    -   Cause: `imagio.mimsave()` does not like processing images which aren't of the same size.
    -   Fix: generate the images using a style containing both a `canvas_width` and `canvas_height` so that the PoV-Ray images are of a uniform size. Canvas height is not usually set.