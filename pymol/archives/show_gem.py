from pymol.cgo import *
from pymol import cmd, stored
import os,sys,shutil
import yaml

############################################################################
## imports
############################################################################
# directory of current .py
pymolpath=os.path.realpath(sys.argv[0])
pymoldir=os.path.dirname(pymolpath)
print "{:<20s}{:<s}".format("pymol dir",pymoldir)

# import pymol code
py = os.path.join(pymoldir,"py")
sys.path.append(py)
from pymol_functions import *
from functions import *
from boundingbox import *

############################################################################
## input
############################################################################
if (len(sys.argv) != 3):
    print "syntax: <param.yml> <.xyz file>"
    sys.exit()

# param file
paramfile = sys.argv[1]
if not (os.path.isfile(paramfile)):
    raise ValueError("File {} does not exist!".format(paramfile))
if not (os.path.splitext(paramfile)[1] == '.yml'):
    raise ValueError("Wrong parameter file: {}".format(paramfile))

# trajectory
trajfile=sys.argv[2]
if not (os.path.isfile(trajfile)):
    raise ValueError("File {} does not exist!".format(trajfile))

trajfile=os.path.realpath(trajfile)
simudir=os.path.dirname(trajfile)

dest = os.path.basename(paramfile)
dest = os.path.join(simudir,dest)
shutil.copyfile(paramfile,dest)
paramfile = dest

fin = open(paramfile,"r")
params = yaml.load(fin)
fin.close()

# data file
datafile=None
if ('datafile' in params):
    datafile = params['datafile']
    datafile = os.path.join(simudir,datafile)
    if not (os.path.isfile(datafile)):
        print ("File {} does not exist! Continuing without.".format(datafile))
        datafile = None

############################################################################
##  MAIN
############################################################################
# check extension
if not (os.path.splitext(trajfile)[1] in params['exts']):
    print "File {} is not a trajectory file!".format(trajfile)
    sys.exit(1)
print "{:<20s}{:<s}".format("trajfile",trajfile)

# draw box
#if ( ('Lx' in params) and ('Ly' in params) and ('Lz' in params) ):
#    print "making box!"
#    Lx=float(params['Lx'][0])
#    Ly=float(params['Ly'][0])
#    Lz=float(params['Lz'][0])
#    drawBox(-Lx/2.,-Ly/2.,-Lz/2.,Lx/2.,Ly/2.,Lz/2., r=0.5,g=0.5,b=0.5)

# load polymer
polymer_colloids_standard([trajfile], **params['polymer_colloids_standard'])

# load bonds
if (datafile != None):
    bonds=data_lammps_get_couplings(datafile)
    color_pairs(bonds,**params['color_pairs'])

make_movie_bool=bool(params['make_movie_bool'])
if (make_movie_bool):
    make_movie(**params['make_movie'])
make_state_movie_bool=bool(params['make_state_movie_bool'])

if (make_state_movie_bool):
    make_state_movie(**params['make_state_movie'])

quit=bool(params['quit'])
if (quit):
    cmd.quit()
