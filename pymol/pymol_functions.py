from pymol.cgo import *
from pymol import cmd, stored, util
import subprocess as sb
import os,sys,glob
import numpy as np

def load_xyz_states_list(files, monomer_color='tv_blue', bg_color='white', monomer_size=1.0, sphere_transparency=0.0, showspheres=True):
    """
    This function loads .xyz states in the input files and build a trajectory.
    """
    # initialization of pymol
    cmd.reinitialize()
    cmd.set("max_threads",1)
    nfiles = len(files)
    print ("objects to load: %d" %nfiles)

    # list names
    names = [os.path.splitext(os.path.basename(f))[0] for f in files]

    # load states
    for i in range(nfiles):
        n = names[i]
        f = files[i]
        cmd.load(f)
        cmd.color(monomer_color,n)

    # join states
    cmd.join_states('trajectory', " ".join(names), mode=0)
    cmd.delete(" ".join(names))

    # display parameters
    cmd.hide("everything","all")
    cmd.alter("all","vdw={:.1f}".format(monomer_size / 2.0))
    if (showspheres):
        cmd.show("spheres","all")
    cmd.bg_color(bg_color)
    cmd.set("orthoscopic",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("depth_cue",0)
    cmd.set("sphere_transparency",sphere_transparency)
    cmd.zoom("all",complete=0)
    return

cmd.extend("load_xyz_states_list",load_xyz_states_list)

def load_xyz_states(trajdir, prefix="state"):
    """
    This function loads .xyz states in the trajdir file and build a trajectory.
    The trajectory is ordered by alphabetical order.
    """
    # make list of state files
    lst = glob.glob(trajdir+prefix+"*.xyz")
    lst.sort()

    return load_xyz_states_list(lst)

cmd.extend("load_xyz_states",load_xyz_states)

def polymer_colloids_standard(path=["pol.xyz"], polymer_color='tv_blue', monomer_size=1.0, sphere_transparency=0.0, showspheres=True):
    ## LOADING
    cmd.reinitialize()
    cmd.set("max_threads",1)
    obj=[os.path.splitext(p)[0].split('/')[-1] for p in path]
    ntypes=len(obj)
    print ("objects to load: %d" %ntypes)

    # COLORS DEFINITION
    pol_color="black" # polymer is black
    col=[]
    for j in range(ntypes-1):
        hsv= (240-j*240/(ntypes-1),1,1)
        # convert to rgb
        rgb=hsv_to_rgb(hsv)
        # define color
        col_name="mycolor%d" %j
        cmd.set_color(col_name,rgb)
        col.append(col_name)

    # MAKE THE COLORING
    sel=obj[0]
    p=path[0]
    #sel_col="black"
    sel_col=polymer_color
    cmd.load(p)
    cmd.color(sel_col,sel)

    for j in range(1,ntypes):
        sel=obj[j]
        p=path[j]
        #sel_col=col[j-1]
        sel_col="tv_red"
        cmd.load(p)
        cmd.color(sel_col,sel)

    ## DISPLAY FORMATTING
    cmd.hide("everything","all")
    cmd.alter("all","vdw={:.1f}".format(monomer_size / 2.0))
    if (showspheres):
        cmd.show("spheres","all")
    cmd.bg_color("white")
    cmd.set("orthoscopic",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("depth_cue",0)
    cmd.set("sphere_transparency",sphere_transparency)
    cmd.zoom("all",complete=0)


cmd.extend("polymer_colloids_standard",polymer_colloids_standard)

