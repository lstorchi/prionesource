from pymol.cgo import *
from pymol import cmd


cmd.load('test_apbs.dx')
cmd.isosurface('sp', "test_apbs", +2.0)
cmd.color('blue', 'sp')
cmd.isosurface('sm', "test_apbs", -2.0)
cmd.color('red', 'sm')

