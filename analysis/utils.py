import numpy as np
import sys
import re

def print_stats(vector):
    N = len(vector)
    mu = np.mean(vector)
    std = np.std(vector)
    cv = std/mu

    print ("{N:<18s}{mu:<18s}{std:<18s}{cv:<18s}".format(N='N',mu='mean',std='std',cv='cv'))
    print ("{N:<18d}{mu:<18.6f}{std:<18.6f}{cv:<18.1f}".format(N=N,mu=mu,std=std,cv=cv*100))
    return

def read_traj_xyz(filepath,istart=None,iend=None):
    """
    Load trajectory in xyz format.
    """

    iters = []
    traj = []
    #pattern="Timestep:(\d+)$"
    pattern="Timestep:(\d+).*$"
    fin = open(filepath,'r')
    state = []

    if (istart is None):
        istart = 0
    if (iend is None):
        iend = 99.9e99

    it_current=-1
    while True:
        line=fin.readline()
        if ((line == "")):
            print("Reached the end of file")
            break

        m = re.search(pattern,line)
        try:
            it_current = np.uint(m.group(1))
        except AttributeError:
            pass

        if (it_current < istart):
           continue

        if (it_current > iend):
            break

        tab = line.replace("\n","").split()
        try:
            si, x, y, z = tab
            si = np.uint(si)
            x = np.float_(x)
            y = np.float_(y)
            z = np.float_(z)
            coord = [si,x,y,z]
            state.append(coord)
        except ValueError as e:
            natoms = len(state)
            if (natoms > 0):
                iters.append(it_current)
                traj.append(state)
                #print(np.array(state))
                print("iter: {:d}. Configuration with {:d} atoms imported".format(it_current,natoms))
                state = []
    # end loop
    # add atoms if a state is loaded
    natoms = len(state)
    if (natoms > 0):
        traj.append(state)
        #print(np.array(state))
        print("iter: {:d}. Configuration with {:d} atoms imported".format(iters[-1],natoms))
        state = []


    print ("niters = {:d}    nconfs = {:d}".format(len(iters),len(traj)))
    return traj, iters






