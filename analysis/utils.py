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
    Load a trajectory in the xyz format written by `MDWorld::dump_xyz`:
    frames of a count line, an `Atoms. Timestep:<iter>` line, then
    `index x y z` lines. Returns the list of frames and their iterations.
    """

    iters = []
    traj = []
    pattern = re.compile(r"Timestep:\s*(\d+)")
    state = []

    if (istart is None):
        istart = 0
    if (iend is None):
        iend = float('inf')

    def flush():
        # store the frame being read, if any
        if state:
            iters.append(it_current)
            traj.append(list(state))
            print("iter: {:d}. Configuration with {:d} atoms imported".format(it_current,len(state)))
            state.clear()

    it_current = -1
    with open(filepath,'r') as fin:
        for line in fin:
            m = pattern.search(line)
            if m:
                it_current = int(m.group(1))

            if (it_current < istart):
                continue

            if (it_current > iend):
                break

            tab = line.split()
            try:
                si, x, y, z = tab
                state.append([int(si), float(x), float(y), float(z)])
            except ValueError:
                # header line: the previous frame is complete
                flush()
    # store the last frame
    flush()

    print ("niters = {:d}    nconfs = {:d}".format(len(iters),len(traj)))
    return traj, iters

def read_map(tfile,xy=False):
    # read matrices
    DATA=np.loadtxt(tfile)
    NN,D=DATA.shape

    if (D != 3):
        print ("wrong format: D={:d}!".format(D))
        sys.exit()

    N = np.sqrt(float(NN))
    if (int(N) != N):
        print ("Length is not a square root! sqrt(NN) = {:.4g}".format(N))
    N = int(N)

#    print "{:<20s}{:<d}".format("N",N)
#    print "{:<20s}{:<d}".format("D",D)

    MAT=np.reshape(DATA[:,2],(N,N))
    if not (xy):
        return MAT
    else:
        X=np.reshape(DATA[:,0],(N,N))
        Y=np.reshape(DATA[:,1],(N,N))
        return X,Y,MAT

def write_map(MAT,tfile,xy=True):
    if (len(MAT.shape) != 2):
        raise ValueError("wrong format: len(shape)={:d}!".format(len(MAT.shape)))

    N,N=MAT.shape

    fout = open(tfile,"w")
    for i in range(N):
        for j in range(N):
            val = MAT[i][j]
            if (xy):
                fout.write("{:<20d}{:<20d}{:<20.8e}\n".format(i,j,val))
            else:
                fout.write("{:<20.8e}\n".format(val))
    fout.close()
    return

