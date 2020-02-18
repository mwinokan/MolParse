from ase import io

import mcol              # https://github.com/mwinokan/MPyTools
import mout              # https://github.com/mwinokan/MPyTools
import asemolplot as amp # https://github.com/mwinokan/AseMolPlot

atoms = io.read('traj.traj', index=-10)

amp.makePovImage('test',atoms,**amp.styles.standard)