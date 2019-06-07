#################### imports ####################
# standard
import sys,io,os,glob
import numpy as np
import yaml
import argparse
import shutil
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
import matplotlib.ticker
from copy import deepcopy
from scipy.stats import chi2
from scipy.stats import gamma


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
yaml.add_representer(np.float_,npfloat_representer)
yaml.add_representer(np.ndarray,nparray_representer)
#yaml.add_representer(unicode,unicode_representer)

# matplotlib controls
plt.rcParams['svg.fonttype'] = 'none'  # to embed fonts in output ('path' is to convert as text as paths)
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams['axes.linewidth']=0.5

# parameters
params={}
params['delta_t'] = 0.01
#params['istart'] = 700*1e5
params['istart'] = None
#params['iend'] = 800*1e5
params['iend'] = None
mydict = {}
#mydict['iter_sel'] = (np.arange(10)+770)*1e5
mydict['iter_sel'] = (np.arange(2)+772)*1e5
mydict['delta_t'] = params['delta_t']
mydict['aratio'] = 5
mydict['xlo'] = -20.
mydict['xhi'] = 20.
mydict['binw'] = 1.
mydict['plot_npts'] =  1000
params['plot_density_profiles_args']=mydict

mydict = {}
mydict['threshold'] = 1.0
mydict['dump_start'] = None
mydict['dump_end'] = None
mydict['ffac'] = 'gaussian'
params['compute_contacts_args']=mydict

#################### functions ####################
def compute_linear_fit(X,Y):
    """
    Compute the linear fit Y=a + b. X
    Returns:
      * a
      * b
      * pvalue
    """

    npts = len(X)
    ndeg = npts-2
    sigmas = np.ones(npts)

    S = np.sum(1./sigmas**2)
    Sx = np.sum(X/sigmas**2)
    Sy = np.sum(Y/sigmas**2)
    Sxx = np.sum(X**2/sigmas**2)
    Sxy = np.sum(X*Y/sigmas**2)
    D = S*Sxx - Sx**2
    a = (Sxx*Sy - Sx*Sxy)/D
    b = (S*Sxy - Sx*Sy)/D
    xm=Sx/S
    ym=Sy/S

    lsq = np.sum((Y - a - b*X)**2/sigmas**2)

    pval = 1. - chi2(ndeg).cdf(lsq)
    r = np.sum((X-xm)*(Y-ym)/sigmas**2)/np.sqrt(np.sum((X-xm)**2/sigmas**2)*np.sum((Y-ym)**2/sigmas**2))

    # alternative
    #pval = 1. - gamma(0.5*npts).cdf(0.5*lsq)

    return a,b,pval,lsq,r

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
    C = np.zeros((N,N), dtype=np.float_)

    # initialize form factor function
    # function takes a numpy array of distances as input
    if ffac == 'gaussian':
        print "Gaussian form factor"
        ffunc = lambda D: np.exp(-1.5*D**2/threshold**2)
    else:
        print "Theta form factor"
        ffunc = lambda D: np.float_(np.uint(D <= threshold))

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
    C /= nstate

    # write
    fileout = os.path.join(outputdir, filename)
    write_map(C, fileout)
#    with open(fileout,'w') as fout:
#        np.savetxt(fout, C)
    print ("Contact matrix written to: {}".format(fileout))

    return

def plot_density_profiles(traj, iters, iter_sel, outputdir='.', xlo=-20., xhi=20., binw=1., delta_t=1., plot_npts=100, unit_length=None,aratio=4./3):
    """
    Plot density profile along cell length for configurations at input iterations:
    """

    # make list of configurations
    iter_sel = np.array(iter_sel, dtype=np.uint)
    nsel = len(iter_sel)
    confs = []
    for i in range(nsel):
        it = iter_sel[i]
        idx = iters.index(it)
        conf = traj[idx]
        confs.append(conf)

    # create color list
    cmap = matplotlib.cm.Set1
    norm = lambda n: n%(cmap.N)

    # compute distributions
    nbins = (xhi-xlo)/binw
    edges = np.arange(nbins+1,dtype=np.float_)*binw + xlo
    dists_free = []
    dists_tot = []
    for i in range(nsel):
        conf = confs[i]
        states = np.array([coord[0] for coord in conf])
        X = np.array([coord[1] for coord in conf])
        idx = (states == 0)
        Xfree = X[idx]

        hist, ed= np.histogram(X, bins=edges, density=False)
        dists_tot.append(hist)
        hist, ed= np.histogram(Xfree, bins=edges, density=False)
        dists_free.append(hist)

    # make figure
    ## figure init
    axh = 3.
    axw = axh*aratio
    nrows=2
    fig = plt.figure(num='none', facecolor='w', figsize=(axw,axh*nrows))
    gs = mgs.GridSpec(nrows,1)

    ## axes init
    axes=[]
    ax = fig.add_subplot(gs[0,0])
    axes.append(ax)
    ax = fig.add_subplot(gs[1,0],sharex=axes[0])
    axes.append(ax)
    naxes = len(axes)

    ## plots
    X = 0.5*(edges[:-1] + edges[1:])

    ### plot total concentration
    ax = axes[0]
    for i in range(nsel):
        Y = dists_tot[i]
        npart = np.uint(np.sum(Y))
        ax.plot(X, Y, '-', color=cmap(norm(i)), lw=0.5, label="iter {:d} / n = {:d}".format(iter_sel[i], npart))

    ax.set_title("total concentration", fontsize='large')
    ax.set_ylabel("pdf",fontsize="medium",labelpad=10,rotation='vertical')
    ax.legend(loc='best', frameon=False, fontsize='medium')

    ### plot free concentration
    ax = axes[1]
    for i in range(nsel):
        Y = dists_free[i]
        npart = np.uint(np.sum(Y))
        ax.plot(X, Y, '-', color=cmap(norm(i)), lw=0.5, label="iter {:d} / n = {:d}".format(iter_sel[i], npart))

    ax.set_title("free concentration", fontsize='large')
    ax.set_ylabel("pdf",fontsize="medium",labelpad=10,rotation='vertical')
    ax.legend(loc='best', frameon=False, fontsize='medium')

    ### additional tweaks
    for i in range(naxes):
        ax = axes[i]
        ax.set_ylim(0.,None)

        ax.tick_params(axis='both', labelbottom='off', labelsize='medium', length=4)
        #ax.set_xlim(0.,None)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['bottom'].set_smart_bounds(True)
        if (i  == naxes-1 ):
            ax.tick_params(labelbottom='on')

    # write figure
    gs.tight_layout(fig, w_pad=1.0)
    filename = 'density_profiles'
    exts=['.pdf', '.svg', '.png']
    for ext in exts:
        fileout = os.path.join(outputdir,filename+ext)
        fig.savefig(fileout, bbox_inches='tight', pad_inches=0)
        print ("Fileout: {:<s}".format(fileout))
    plt.close('all')
    return

#################### main ####################
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Analysis to analyze a trajectory of xyz files")
    parser.add_argument('trajdir', type=str, help='Paths to directory containing .xyz files.')
    #parser.add_argument('--paramfile', '-f',  type=file, required=True, help='Path to param file.')
    parser.add_argument('--outputdir', '-d',  type=str, required=False, default=os.path.join('.', 'analysis_xyz'), help='Path to output directory.')
    parser.add_argument('--compute_contacts',  action='store_true', help='Show density profiles.')

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
            with open(f,'r') as fin:
                state = np.loadtxt(f,skiprows=2)
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
        print traj.shape
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
    """
    if namespace.density_profiles:
        plot_density_profiles(traj, iters, outputdir=outputdir, **params['plot_density_profiles_args'])

    if namespace.showlineardependence:
        plot_linear_dependence(traj, outputdir=outputdir, **params['plot_linear_dependence_args'])
    """
