#! /bin/python3

import os, sys
from pathlib import Path
import argparse
import yaml
import numpy as np

##############################################################################
# GLOBAL
##############################################################################

##############################################################################
# MAIN
##############################################################################
########## I/O ##########
if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="Helper script to generate configuration file for polymers.")
    parser.add_argument('config', type=str, help='Path to configuration file to use.')
    parser.add_argument('--outputdir', '-d',  type=str, required=False, \
        default=Path('.'), help='Path to output directory.')
    parser.add_argument('--length',  type=int, default=0, help='Number of monomer per chain.')
    parser.add_argument('--density',  type=float, default=0, help='Density of chains.')
    parser.add_argument('--b',  type=float, default=None, help='Parameter b for a Gaussian chain potential.')
    parser.add_argument('--r0_ke',  type=float, nargs=2, default=None, help='Parameters r0 and ke for an harmonic chain potential.')
    parser.add_argument('--rc_ke_eps',  type=float, nargs=3, default=None, help='Parameters rc and ke for a FENE chain potential.')
    parser.add_argument('--lp',  type=float, default=None, help='Persistence length for a Kratky-Porod model.')

    # load arguments
    namespace = parser.parse_args(sys.argv[1:])

    # load configuration file
    fpath = Path(namespace.config)
    if not fpath.is_file():
      raise ValueError(f"File doesn't exist: {str(fpath)}")

    with open(fpath,'r') as fin:
      config = yaml.load(fin, Loader=yaml.loader.FullLoader)

    fpath_out = fpath.parent / (fpath.stem + '_helper.yaml')
    fout = open(fpath_out, 'w')

########## calculations ##########
    # get box size
    lx, ly, lz = config['MDWorld']['box']
    lx = float(lx)
    ly = float(ly)
    area = lx*ly
    print(f"lx = {lx}", f"ly = {ly}", f"area = {area}")

    # force fields
    fout.write('# forcefields\n')
    N = namespace.length
    print(f"N = {N}")
    if ( N == 0 ):
      raise ValueError("Must provide a non-zero length for the chains.")
    noffset = 0
    nchain = 0
    phi = 0.
    dphi = float(N)*np.pi*0.25 / area     # area fraction of one chain

    print(f"Target density = {namespace.density}")
    while phi < namespace.density:
      if namespace.b:  # Gaussian chain
        mydir = {}
        mydir['PolymerGaussian'] = {'offset': noffset, \
                                    'N': N, \
                                    'b': float(namespace.b) \
                                   }
        yaml.dump(mydir, fout)

      if namespace.r0_ke:  # Harmonic chain
        r0, ke = namespace.r0_ke

        mydir = {}
        mydir['PolymerHarmonic'] = {'offset': noffset, \
                                    'N': N, \
                                    'ke': float(ke), \
                                    'r0': float(r0) \
                                   }
        yaml.dump(mydir, fout)

      if namespace.rc_ke_eps:  # FENE chain
        rc, ke, eps = namespace.rc_ke_eps

        mydir = {}
        mydir['PolymerFENE'] = {'offset': noffset, \
                                    'N': N, \
                                    'ke': float(ke), \
                                    'rc': float(rc), \
                                    'sigma': 1., \
                                    'epsilon': eps \
                                   }
        yaml.dump(mydir, fout)

      if namespace.lp:  # Kratky-Porod potential
        mydir = {}
        mydir['PolymerKratkyPorod'] = {'offset': noffset, \
                                    'N': N, \
                                    'lp': float(namespace.lp) \
                                   }

        yaml.dump(mydir, fout)
      # increment
      nchain += 1
      noffset += N
      phi += dphi

    print("nchains = {:d}".format(nchain))
    print("nmonomers = {:d}".format(N*nchain))



