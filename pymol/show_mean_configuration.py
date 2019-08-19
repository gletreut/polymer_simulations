############################################################################
## imports
############################################################################
#from __future__ import print_function
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
## functions
############################################################################
def check_file(myfile):
    if not os.path.isfile(myfile):
        sys.exit("File doesn't exist: {:s}".format(myfile))
    return

############################################################################
## global parameters
############################################################################
# general
prefix='state'
quit=True
nc=80       # index of the centromer
zoom=20
view=None
#view="-0.467148870,0.704388320,-0.534410000,0.881349385,0.322722554,-0.345051497,-0.070578903,-0.632196784,-0.771575153,0.000000000,0.000000000,-164.520812988,-2.193372726,-3.725057602,-2.012924194,129.709457397,199.332168579,20.000000000" # scale 1

# for movie making
makemovie=False
fname_moviestate="state_avg_fluctuations.gif"
nmin=0
nmax=1000
fps=20

# for state movie making
makestatemovie=True
state_sel=1

############################################################################
## main script
############################################################################
if (len(sys.argv) < 2):
    print ("syntax: <file.xyz> [fluctuations.dat]")
    sys.exit()

# trajectory file
trajfile=sys.argv[1]
check_file(trajfile)
print ("{:<20s}{:<s}".format("trajfile",trajfile))
trajdir = os.path.dirname(trajfile)

# fluctuations file
try:
    flucfile=sys.argv[2]
    check_file(flucfile)
    print ("{:<20s}{:<s}".format("trajfile",flucfile))
except IndexError:
    flucfile=None

# load states in trajectory
load_single_xyz(trajfile, file_surface_radius=flucfile)

# set view
if not view is None:
    set_view_states(view)

# set zoom
if not zoom is None:
    cmd.center()
    cmd.zoom("all", zoom)

# color first 3 monomers
#color_123()

# emphasize the centromere region
#emphasize_centromer(nc)

# draw box
#if ( ('Lx' in params) and ('Ly' in params) and ('Lz' in params) ):
#    print "making box!"
#    Lx=float(params['Lx'][0])
#    Ly=float(params['Ly'][0])
#    Lz=float(params['Lz'][0])
#    drawBox(-Lx/2.,-Ly/2.,-Lz/2.,Lx/2.,Ly/2.,Lz/2., r=0.5,g=0.5,b=0.5)

if (makemovie):
    make_movie(npts=1001,nmin=nmin,nmax=nmax, imgpx_w=800, imgpx_h=600, framerate=fps)

if (makestatemovie):
    make_state_movie(state=state_sel,fileout=os.path.join(trajdir,fname_moviestate),nframes=150,imgpx_w=800, imgpx_h=600)

if (quit):
    cmd.quit()
