#!/usr/bin/env python3

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

########## calculations ##########
    # get box size
    lx, ly, lz = config['MDWorld']['box']
    lx = float(lx)
    ly = float(ly)
    area = lx*ly
    print(f"lx = {lx}", f"ly = {ly}", f"area = {area}")

    # number of chains
    N = namespace.length
    print(f"N = {N}")
    if ( N == 0 ):
      raise ValueError("Must provide a non-zero length for the chains.")
    print(f"Target density = {namespace.density}")
    if ( namespace.density <= 0 ):
      raise ValueError("Must provide a positive density of chains.")
    dphi = float(N)*np.pi*0.25 / area     # area fraction of one chain
    nchain = int(np.ceil(namespace.density / dphi))
    print("nchains = {:d}".format(nchain))
    print("nmonomers = {:d}".format(N*nchain))

    # force fields, one entry per chain and type. The simulation reads repeated
    # keys in `forcefields` as separate force fields, so the output repeats
    # keys and is written as text (PyYAML cannot dump repeated keys).
    def chain_forcefields(offset):
      forcefields = {}
      if namespace.b:  # Gaussian chain
        forcefields['PolymerGaussian'] = {'offset': offset, 'N': N, 'b': float(namespace.b)}
      if namespace.r0_ke:  # Harmonic chain
        r0, ke = namespace.r0_ke
        forcefields['PolymerHarmonic'] = {'offset': offset, 'N': N, 'ke': float(ke), 'r0': float(r0)}
      if namespace.rc_ke_eps:  # FENE chain
        rc, ke, eps = namespace.rc_ke_eps
        forcefields['PolymerFENE'] = {'offset': offset, 'N': N, 'ke': float(ke), 'rc': float(rc),
                                      'sigma': 1., 'epsilon': float(eps)}
      if namespace.lp:  # Kratky-Porod potential
        forcefields['PolymerKratkyPorod'] = {'offset': offset, 'N': N, 'lp': float(namespace.lp)}
      return forcefields

    lines = ['forcefields:']
    for c in range(nchain):
      for key, params in chain_forcefields(c*N).items():
        lines.append(f"  {key}: " + yaml.dump(params, default_flow_style=True, sort_keys=False).strip())

    with open(fpath_out, 'w') as fout:
      fout.write('\n'.join(lines) + '\n')
    print(f"Force fields written to: {fpath_out}")
