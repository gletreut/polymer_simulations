#!/usr/bin/env python3
#################### imports ####################
# standard
import sys,os,glob
import numpy as np
import yaml
import argparse

# custom
origin=os.path.dirname(os.path.realpath(sys.argv[0]))
sys.path.append(origin)
from utils import write_map, read_traj_xyz

#################### global params ####################

# yaml formats
npfloat_representer = lambda dumper,value: dumper.represent_float(float(value))
nparray_representer = lambda dumper,value: dumper.represent_list(value.tolist())
float_representer = lambda dumper,value: dumper.represent_scalar(u'tag:yaml.org,2002:float', "{:<.8e}".format(value))
#unicode_representer = lambda dumper,value: dumper.represent_unicode(value.encode('utf-8'))
yaml.add_representer(float,float_representer)
yaml.add_representer(np.float64,npfloat_representer)
yaml.add_representer(np.ndarray,nparray_representer)
#yaml.add_representer(unicode,unicode_representer)

# parameters
params={}
params['delta_t'] = 0.01
params['istart'] = None
params['iend'] = None

mydict = {}
mydict['threshold'] = 1.0
mydict['dump_start'] = None
mydict['dump_end'] = None
mydict['ffac'] = 'gaussian'
params['compute_contacts_args']=mydict

#################### functions ####################
def compute_contacts(traj, outputdir='.', dump_start=None, dump_end=None, threshold=1.0, filename="cmat_xyz.dat", ffac='gaussian'):
    """
    Compute the contact probability matrix from a given ensemble of states.
    """

    # initialize place holders for average
    nstate = len(traj)
    state = traj[0]
    N = len(state)
    print ("Size of matrix N = {:d}".format(N))
    count = 0
    C = np.zeros((N,N), dtype=np.float64)

    # initialize form factor function
    # function takes a numpy array of distances as input
    if ffac == 'gaussian':
        print("Gaussian form factor")
        ffunc = lambda D: np.exp(-1.5*D**2/threshold**2)
    else:
        print("Theta form factor")
        ffunc = lambda D: (D <= threshold).astype(np.float64)

    # iterate over configuration
    for i in range(nstate)[dump_start:dump_end]:
        XYZ = traj[i]
        X = XYZ[:,0]
        Y = XYZ[:,1]
        Z = XYZ[:,2]
        count += 1
        for n in range(N):
            xn = X[n]
            yn = Y[n]
            zn = Z[n]
            norm = np.sqrt((X-xn)**2 + (Y-yn)**2 + (Z-zn)**2)
            C[n] += ffunc(norm)

        # end loop on monomers
    # end loop on states
    if count == 0:
        raise ValueError("No state selected for the contact map.")
    C /= count

    # write
    fileout = os.path.join(outputdir, filename)
    write_map(C, fileout)
#    with open(fileout,'w') as fout:
#        np.savetxt(fout, C)
    print ("Contact matrix written to: {}".format(fileout))

    return

#################### main ####################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Analysis to analyze a trajectory of xyz files")
    parser.add_argument('trajdir', type=str, help='Paths to directory containing .xyz files.')
    #parser.add_argument('--paramfile', '-f',  type=file, required=True, help='Path to param file.')
    parser.add_argument('--outputdir', '-d',  type=str, required=False, default=os.path.join('.', 'analysis_xyz'), help='Path to output directory.')
    parser.add_argument('--compute_contacts',  action='store_true', help='Compute the contact probability map.')

    # INITIALIZATION
    # load arguments
    namespace = parser.parse_args(sys.argv[1:])

    # trajdir
    if os.path.isdir(namespace.trajdir):
        print ("{:<20s}{:<s}".format("trajdir", namespace.trajdir))

        # trajectory
        trajfiles = glob.glob(os.path.join(namespace.trajdir,"state*.xyz"))
        trajfiles.sort()
        traj = []
        for f in trajfiles:
            if not os.path.isfile(f):
                print("File does not exist: {:s}".format(f))
                continue
            print ("loading file {}...".format(f))
            state = np.loadtxt(f, skiprows=2, ndmin=2)
            indices = state[:,0]
            XYZ = state[:,1:]
            if ( len(np.unique(np.diff(indices))) != 1 ):
                sys.exit("Indices are not consecutive in the file.")
            idx = np.argsort(indices)
            XYZ = XYZ[idx]
            traj.append(XYZ)

    elif os.path.isfile(namespace.trajdir):
        print ("{:<20s}{:<s}".format("trajfile", namespace.trajdir))
        traj, iters = read_traj_xyz(namespace.trajdir)
        traj = np.array(traj)
        traj = traj[:,:,1:]
        print(traj.shape)
    else:
        sys.exit("Directory does not exist!")
    nstates = len(traj)
    print ("nstates = {:d}".format(nstates))

    if (nstates == 0):
        sys.exit("Empty set of states!")

    # output directory
    outputdir = namespace.outputdir
    if not os.path.isdir(outputdir):
        os.makedirs(outputdir)
    print ("{:<20s}{:<s}".format("outputdir", outputdir))

    # parameter file
#    if namespace.paramfile is None:
#        paramfile="analysis_default.yml"
#        params = default_parameters()
#        with open(paramfile,'w') as fout:
#            yaml.dump(params,fout)
#    else:
#        paramfile = namespace.paramfile.name
#        params = yaml.load(namespace.paramfile)
#
#    dest = os.path.join(outputdir, os.path.basename(paramfile))
#    if (os.path.realpath(dest) != os.path.realpath(paramfile)):
#        shutil.copy(paramfile,dest)
#    paramfile = dest

    # OPERATIONS
    ## trajectories (of concentration + production rate)
    if namespace.compute_contacts:
        compute_contacts(traj, outputdir=outputdir, **params['compute_contacts_args'])
