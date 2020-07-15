from pymol.cgo import *
from pymol import cmd
import glob
import os

lista = ["e200k_coulomb_mean.dx"]
#        "E200K_clust-cutoff02_2.dx", \
#        "E200K_clust-cutoff02_3.dx"]

lista = glob.glob("wt*.dx")

for idx, name in enumerate(lista):
    cmd.load(name)
    basename = os.path.splitext(name)[0]
    cmd.isosurface(basename+'_sp', basename, +0.24000)
    cmd.color('blue', basename+'_sp')
    cmd.isosurface(basename+'_sm', basename, -0.16000)
    cmd.color('red', basename+'_sm')

