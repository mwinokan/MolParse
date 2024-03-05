#!/usr/bin/env python3

import molparse as mp

sys = mp.parse("/Users/mw00368/OneDrive/Leverhulme/Projects/ArcherUmbrellaHelicase/gro/2CR_dry.gro")

sys.name = "2CR_dry"

sys.tree()

# mp.tree(sys)
