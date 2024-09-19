isPovLoaded = False


def loadPov(verbosity=1, purge=True, anaconda=False):
    import module  # https://github.com/mwinokan/MPyTools

    # load the correct modules for ASE/PoV-Ray
    if purge:
        module.module("purge")
    module.module("--expert", "load", "Boost/1.63.0-intel-2017a-Python-2.7.13")
    module.module("--expert", "load", "zlib/1.2.8-intel-2016a")
    module.module("--expert", "load", "libpng/1.6.24-intel-2016a")
    module.module("--expert", "load", "libjpeg-turbo/1.5.0-intel-2016a")
    module.module("--expert", "load", "LibTIFF/4.0.6-intel-2016a")
    if anaconda:
        module.module("--expert", "load", "anaconda3/2019.03")

    if verbosity > 0:
        mout.out(
            "PoV-Ray dependencies loaded.",
            printScript=True,
        )

    global isPovLoaded
    isPovLoaded = True


def makePovImage(
    filename,
    image,
    verbosity=1,
    rmPovFiles=True,
    bonds=False,
    bondradius=1.1,
    forceLoad=False,
    printScript=False,
    **style,
):
    from ase import io
    import mcol
    import mout

    if not isPovLoaded or forceLoad:
        loadPov(verbosity=verbosity - 1)

    if verbosity > 0:
        mout.out(
            "processing " + mcol.file + filename + ".pov" + mcol.clear + " ... ",
            printScript=printScript,
            end="",
        )  # user output

    if not style["drawCell"]:
        image.set_cell([0, 0, 0])
    del style["drawCell"]

    if "canvas_height" in style:
        del style["canvas_height"]

    if bonds:
        from ase.io.pov import get_bondpairs, set_high_bondorder_pairs

        bondpairs = get_bondpairs(image, radius=bondradius)
        if len(bondpairs) > 5000:
            mout.warningOut(
                "Too many bondpairs (" + str(len(bondpairs)) + "), not drawing bonds!",
                end=" ",
            )
        else:
            style["bondatoms"] = bondpairs

    io.write(
        filename + ".pov", image, run_povray=True, camera_type="perspective", **style
    )

    if rmPovFiles:
        import os

        os.system("rm " + filename + ".ini")
        os.system("rm " + filename + ".pov")

    if verbosity > 0:
        mout.out("Done.")  # user output


def makePovImages(
    filename,
    subdirectory="pov",
    interval=1,
    verbosity=1,
    rmPovFiles=True,
    bonds=False,
    bondradius=1.1,
    filenamePadding=4,
    printScript=False,
    forceLoad=False,
    index=":",
    **style,
):
    from ase import io
    import mcol
    import mout

    if not isPovLoaded or forceLoad:
        loadPov(verbosity=verbosity - 1)

    import os
    import math

    os.system("mkdir -p " + subdirectory)
    # os.system("rm -v "+subdirectory+"/* ")
    os.system("rm " + subdirectory + "/* 2> /dev/null")

    if index != ":":
        image = io.read(filename, index=index)
        makePovImage(
            subdirectory + "/" + str(index).zfill(filenamePadding),
            image,
            verbosity=verbosity - 1,
            bonds=bonds,
            bondradius=bondradius,
            printScript=printScript,
            rmPovFiles=False,
            forceLoad=forceLoad,
            **style,
        )
    else:
        traj = io.read(filename, index=index)

        global num_traj_images
        global num_frames
        num_traj_images = len(traj)
        num_frames = math.ceil(num_traj_images / interval)
        if verbosity > 0:
            mout.varOut("Trajectory", num_traj_images, unit="images")
            mout.varOut("Animation", num_frames, unit="frames")

        for n, image in enumerate(traj):
            if n % interval != 0 and n != 100:
                continue

            if verbosity == 1:
                mout.progress(
                    n + 1,
                    num_traj_images,
                    prepend="Creating images",
                    printScript=printScript,
                )

            makePovImage(
                subdirectory + "/" + str(n).zfill(filenamePadding),
                image,
                verbosity=verbosity - 1,
                bonds=bonds,
                bondradius=bondradius,
                printScript=printScript,
                rmPovFiles=False,
                forceLoad=forceLoad,
                **style,
            )

    if rmPovFiles:
        os.system("rm " + subdirectory + "/*.ini")
        os.system("rm " + subdirectory + "/*.pov")


# Using ImageMagick (artefacting!):
def makePovAnimationIM(filename, subdirectory="pov", interval=1, verbosity=1, **style):
    import mcol
    import mout

    if not isPovLoaded or forceLoad:
        loadPov(verbosity=verbosity - 1)

    import os
    import module  # https://github.com/mwinokan/MPyTools

    makePovImages(
        filename,
        subdirectory=subdirectory,
        interval=interval,
        verbosity=verbosity - 1,
        **style,
    )

    module.module("--expert", "load", "ImageMagick/7.0.3-1-intel-2016a")
    if verbosity > 0:
        mout.out("ImageMagick loaded.", printScript=True)

    if verbosity > 0:
        mout.out(
            "creating " + mcol.file + "animation.gif" + mcol.clear + " ... ",
            printScript=True,
            end="",
        )  # user output
    os.system(
        "convert -delay 10 "
        + subdirectory
        + "/*.png -fill white -opaque none -loop 1 "
        + subdirectory
        + ".gif"
    )
    if verbosity > 0:
        mout.out("Done.")  # user output


# Using imageio
# https://stackoverflow.com/questions/753190/programmatically-generate-video-or-animated-gif-in-python
def makePovAnimation(
    filename,
    subdirectory="pov",
    interval=1,
    gifstyle=None,
    verbosity=1,
    printScript=False,
    forceLoad=False,
    useExisting=False,
    dryRun=False,
    **plotstyle,
):
    import mcol
    import mout

    if not isPovLoaded or forceLoad:
        loadPov(verbosity=verbosity - 1)

    from . import styles

    if gifstyle is None:
        gifstyle = styles.gif_standard

    if plotstyle == {}:
        plotstyle = styles.standard
    # print(plotstyle)

    global num_traj_images
    global num_frames

    import os
    import imageio

    cropping = False
    shifting = False

    # Set canvas sizes:
    canv_w = plotstyle["canvas_width"]
    if "canvas_height" in plotstyle:
        canv_h = plotstyle["canvas_height"]
        del plotstyle["canvas_height"]
    else:
        mout.errorOut("No canvas_height specified in plotstyle!", fatal=True)

    # Check if cropping:
    if "crop_w" in plotstyle and "crop_w" in plotstyle:
        cropping = True
        crop_w = plotstyle["crop_w"]
        crop_h = plotstyle["crop_h"]
    if "crop_w" in plotstyle:
        del plotstyle["crop_w"]
    if "crop_h" in plotstyle:
        del plotstyle["crop_h"]

    # Check if crop offset:
    if "crop_x" in plotstyle:
        shifting = True
        crop_x = plotstyle["crop_x"]
        crop_y = plotstyle["crop_y"]
    if "crop_x" in plotstyle:
        del plotstyle["crop_x"]
    if "crop_y" in plotstyle:
        del plotstyle["crop_y"]

    mout.varOut("cropping", cropping)
    if cropping:
        mout.varOut("crop_w", crop_w)
        mout.varOut("crop_h", crop_h)
    mout.varOut("shifting", shifting)

    # Generate the PNG's
    if not useExisting:
        if verbosity > 0:
            mout.out(
                "generating "
                + mcol.file
                + subdirectory
                + "/*.png"
                + mcol.clear
                + " ... ",
                printScript=printScript,
                end="",
            )  # user output
        if verbosity > 1:
            mout.out(" ")

        if not dryRun:
            # Generate all the images
            makePovImages(
                filename,
                subdirectory=subdirectory,
                interval=interval,
                printScript=printScript,
                verbosity=verbosity - 1,
                **plotstyle,
            )
        else:
            # Generate just the first image
            makePovImages(
                filename,
                subdirectory=subdirectory,
                interval=interval,
                printScript=printScript,
                verbosity=verbosity - 1,
                index=0,
                **plotstyle,
            )

        if verbosity == 1:
            mout.out("Done.")

    # Load ImageMagick
    # if cropping or backwhite:
    import module  # https://github.com/mwinokan/MPyTools

    ret = module.module("--expert", "load", "ImageMagick/7.0.3-1-intel-2016a")
    if ret == 0 and verbosity > 0:
        mout.out("ImageMagick loaded.", printScript=printScript)

    # Combine the images
    if verbosity > 0:
        mout.out(
            "loading " + mcol.file + subdirectory + "/*.png" + mcol.clear + " ... ",
            printScript=printScript,
            end="",
        )  # user output

    images = []

    # loop over all files in the subdirectory:
    for file in sorted(os.listdir(subdirectory)):

        # get the relative path to the file:
        filename = subdirectory + "/" + file

        # check if the file is a PNG:
        if file.endswith(".png"):

            tempname = subdirectory + "/" + "temp.png"

            # run different IM commands depending on cropping and shifting:
            if not cropping and not shifting:
                # os.system("convert "+filename+
                #           " -background white -extent "+
                #           str(canv_w)+"x"+
                #           str(canv_h)+" "+
                #           filename)
                os.system("convert " + filename + " -flatten " + tempname)
            elif cropping and not shifting:
                os.system(
                    "convert "
                    + filename
                    + " -crop "
                    + str(crop_w)
                    + "x"
                    + str(crop_h)
                    + " -background white -extent "
                    + str(crop_w)
                    + "x"
                    + str(crop_h)
                    + " "
                    + tempname
                )
                print(
                    "convert "
                    + filename
                    + " -crop "
                    + str(crop_w)
                    + "x"
                    + str(crop_h)
                    + " -background white -extent "
                    + str(crop_w)
                    + "x"
                    + str(crop_h)
                    + " "
                    + tempname
                )
                os.system("ls NEB")
            elif shifting and not cropping:
                os.system(
                    "convert "
                    + filename
                    + " -crop +"
                    + str(crop_x)
                    + "+"
                    + str(crop_y)
                    + " -background white -extent "
                    + str(canv_w)
                    + "x"
                    + str(canv_h)
                    + " "
                    + tempname
                )
            else:
                os.system(
                    "convert "
                    + filename
                    + " -crop "
                    + str(crop_w)
                    + "x"
                    + str(crop_h)
                    + "+"
                    + str(crop_x)
                    + "+"
                    + str(crop_y)
                    + " -background white -extent "
                    + str(crop_w)
                    + "x"
                    + str(crop_h)
                    + " "
                    + tempname
                )

            os.system("mv " + tempname + " " + filename)

            # Read in the image and append to the image array
            image = imageio.imread(filename)
            images.append(image)

            if verbosity > 0 and not dryRun:
                mout.progress(
                    len(images),
                    num_frames,
                    prepend="Cropping & loading images",
                    printScript=printScript,
                )

    if (verbosity > 0) and not useExisting:
        mout.out("Done.")  # user output

    if verbosity > 0:
        mout.out(
            "creating " + mcol.file + subdirectory + ".gif" + mcol.clear + " ... ",
            printScript=printScript,
            end="",
        )  # user output

    # Generate the animated GIF:
    imageio.mimsave(subdirectory + ".gif", images, **gifstyle)

    if verbosity > 0:
        mout.out("Done.")  # user output


def crop(filename, width=500, height=500, xshift=0, yshift=0, verbosity=1):
    import mcol
    import mout
    import os

    import module  # https://github.com/mwinokan/MPyTools

    module.module("--expert", "load", "ImageMagick/7.0.3-1-intel-2016a")
    global isPovLoaded
    isPovLoaded = False
    if verbosity > 0:
        mout.out("ImageMagick loaded.", printScript=printScript)

    os.system(
        "convert "
        + filename
        + " -crop "
        + str(width)
        + "x"
        + str(height)
        + "+"
        + str(xshift)
        + "+"
        + str(yshift)
        + " "
        + filename
    )
