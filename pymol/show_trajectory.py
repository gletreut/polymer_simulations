############################################################################
## imports
############################################################################
from pymol.cgo import *
from pymol import cmd, stored
import os,sys,shutil
import yaml

# directory of current .py
filepath=os.path.realpath(sys.argv[0])
origin=os.path.dirname(filepath)
print ("{:<20s}{:<s}".format("origin",origin))

# import pymol code
#py = os.path.join(pymoldir,"py")
sys.path.append(origin)
from pymol_functions import *
from boundingbox import *

############################################################################
## global parameters
############################################################################
prefix='state'
quit=False

############################################################################
## main script
############################################################################
if (len(sys.argv) != 2):
    print ("syntax: <.xyz directory>")
    sys.exit()

# trajectory directory
trajdir=sys.argv[1]
if not (os.path.isdir(trajdir)):
    raise ValueError("Directory {} does not exist!".format(trajfile))
print ("{:<20s}{:<s}".format("traj dir",trajdir))

# load states in trajectory
load_xyz_states(trajdir, prefix=prefix)

# color first 3 monomers
color_123()

# draw box
#if ( ('Lx' in params) and ('Ly' in params) and ('Lz' in params) ):
#    print "making box!"
#    Lx=float(params['Lx'][0])
#    Ly=float(params['Ly'][0])
#    Lz=float(params['Lz'][0])
#    drawBox(-Lx/2.,-Ly/2.,-Lz/2.,Lx/2.,Ly/2.,Lz/2., r=0.5,g=0.5,b=0.5)

if (quit):
    cmd.quit()
