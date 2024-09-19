def import_checks():
    import importlib

    # tight dependencies
    specs = ["mout", "mcol", "mplot", "mwin"]

    for spec in specs:
        if not importlib.util.find_spec(spec):
            print(f"Fatal Error: Missing '{spec}' package")
            print(
                "Clone the latest version of MPyTools and ensure its in your PYTHONPATH"
            )
            print("https://github.com/mwinokan/MPyTools")
            exit(1)

    # python version
    import sys

    if sys.version_info < (3, 10):
        import mout

        mout.warningOut("Not using python 3.10 or newer! This may lead to errors.")

    # loose dependencies
    dependecies = {
        "ase": [
            "various molecular manipulations and I/O",
            "https://wiki.fysik.dtu.dk/ase/",
        ],
        "numpy": ["countless numerical methods"],
        "scipy": ["curve fitting and analysis "],
        "curses": ["moltree CLI program"],
        "json": ["I/O to the CJSON format"],
        "matplotlib": ["Plotting graphs"],
        "plotly": ["Interactive graphs"],
        # 'imageio': ["creating gif animations (deprecated)"],
    }

    for key in dependecies:
        if not importlib.util.find_spec(key):
            import mout

            mout.warningOut(f"Missing {key} package! This may lead to errors.")
            mout.out(f"MolParse uses {key} for: {dependecies[key][0]}")
            if len(dependecies[key]) > 1:
                mout.out(f"get it from: {dependecies[key][1]}")
