#from __future__ import print_function
from pymol.cgo import *
from pymol import cmd, stored, util
import subprocess as sb
import os,sys,glob
import numpy as np

def draw_sphere(x,y,z,radius,rgb=[0.,1.,0.], name="sphere", transparency=0.0):
    """
    Draw a sphere
    """
    cmd_lst = [
            COLOR, rgb[0], rgb[1], rgb[2],
            SPHERE, x, y, z, radius
            ]
    cmd.load_cgo(cmd_lst, name)
    cmd.show("cgo",name)
    cmd.set("cgo_transparency", transparency, name)
    return


def load_xyz_trajectory(trajfile, free_color='tv_blue', bound_color='tv_red', bg_color='white', monomer_size=1.0, target_radius=1.0, sphere_transparency=0.0, target_transparency=0.25):
    """
    This function loads a trajectory.xyz input file.
    """
    # initialization of pymol
    cmd.reinitialize()
    cmd.set("max_threads",1)
    print ("file to load: %s" %trajfile)

    # list names
    name = os.path.splitext(os.path.basename(trajfile))[0]

    # load
    cmd.load(trajfile)

    # color the bound atoms in red and the free ones in blue
    cmd.select("free", "name 0")
    cmd.select("bound", "name 1")
    cmd.color(free_color,"free")
    cmd.color(bound_color,"bound")

    # draw a sphere for the target
    name_target = "target"
    xtarget = 0.
    ytarget = 0.
    ztarget = 0.
    ctarget = [0.,1.,0.]
    draw_sphere(xtarget, ytarget, ztarget, target_radius, rgb=ctarget, name=name_target, transparency=target_transparency)

    # display parameters
    cmd.hide("everything","all")
    cmd.alter("all","vdw={:.1f}".format(monomer_size / 2.0))
    cmd.show("spheres","all")
    cmd.show("cgo",name_target)
    cmd.bg_color(bg_color)
    cmd.set("orthoscopic",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("depth_cue",0)
    cmd.set("sphere_transparency",sphere_transparency)
    cmd.zoom("all",complete=0)
    return

cmd.extend("load_xyz_trajectory",load_xyz_trajectory)


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
        #cmd.load(f, "lala")
        cmd.load(f)
        cmd.color(monomer_color,n)

    # join states
    if (nfiles > 1):
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

def load_single_xyz(file_xyz, file_surface_radius=None, monomer_color='tv_blue', bg_color='white', surface_color='gray', monomer_size=1.0, sphere_transparency=0.0, surface_transparency=0.5):
    """
    This function loads .xyz states in the input files and build a trajectory.
    """
    # initialization of pymol
    cmd.reinitialize()
    cmd.set("max_threads",1)

    # list names
    name = os.path.splitext(os.path.basename(file_xyz))[0]

    # load state
    cmd.load(file_xyz)
    cmd.color(monomer_color,name)

    # radius
    cmd.alter("{:s}".format(name),"vdw={:.1f}".format(monomer_size / 2.0))

    # case with fluctuation information
    if not file_surface_radius is None: # try to load fluctuation surface
        # load the radii
        with open(file_surface_radius,'r') as fin:
            radii = np.loadtxt(fin)

        if radii.shape[1] != 2:
            print(radii)
            print(radii.shape)
            raise ValueError("Wrong shape for the radii loaded!")
        natoms = len(radii)

        # re-load the main positions in a new object
        name_surface = "fluctuations"
        cmd.load(file_xyz, name_surface)
        cmd.color(surface_color,name_surface)

        # alter the radii
        for i in range(natoms):
            resi = int(radii[i,0])+1
            r = radii[i,1]
            print("resi {:d} with r = {:1f}".format(resi,r))
            cmd.alter("{:s} and resi {:d}".format(name_surface, resi),"vdw={:.1f}".format(r))


    # display parameters
    cmd.hide("everything","all")
    cmd.show("spheres",name)
    cmd.show("surface",name_surface)
    #cmd.show("spheres",name_surface)
    cmd.set("solvent_radius", 1.4*(2*monomer_size))       # if too small then surface is incomplete
    cmd.bg_color(bg_color)
    cmd.set("orthoscopic",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("depth_cue",0)
    cmd.set("sphere_transparency",sphere_transparency)  # for spheres
    cmd.set("transparency",surface_transparency)        # for surface
    cmd.set("surface_quality", 1)                       # for surface quality
    cmd.zoom("all",complete=1)


    return

cmd.extend("load_single_xyz",load_xyz_states_list)

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

def make_movie(npts=100,rdir='.',nmin=None,nmax=None, imgpx_w=800, imgpx_h=600, framerate=5):
    ## PARAMETERS
    dirpath=os.path.join(rdir,"output_mov")
    filename=os.path.join(dirpath,"mov")

    if ( not os.path.exists(dirpath) ):
        os.makedirs(dirpath)

    #cmd.viewport(320,240)
    cmd.set("ray_trace_mode", 0)

    ## STATES PROCESSING
    if (nmin==None): nmin=0
    if (nmax==None):nmax=cmd.count_states()-1
    print("number of states=%d" %nmax)
    nt=npts-1
    delta=(nmax-nmin)/float(nt)

    cmd.zoom("all", state=0, complete=1) # reajust view to fit all states (0 means all states)
    for i in range(npts):
        k=int(nmax-i*delta)
        print("state {:04d}".format(k))
        cmd.set("state",k+1)
        #cmd.zoom("type1", state=k+1, complete=1) # reajust view to fit the current state
#        cmd.zoom("all", state=k+1, complete=1) # reajust view to fit the current state
    	#cmd.ray(1200,900)
        cmd.png(filename+"%04d.png" %(npts-i), width=imgpx_w, height=imgpx_h, ray=1)
    bname=filename+"%04d.png"
    # for FFMPEG, see:
    # https://askubuntu.com/questions/610903/how-can-i-create-a-video-file-from-a-set-of-jpg-images
    # https://askubuntu.com/questions/745732/converting-png-files-to-a-movie
    # https://stackoverflow.com/questions/25869206/ffmpeg-png-to-mp4-black-screen
    cmdtab=[]
    cmdtab.append('ffmpeg')
    cmdtab.append('-y')
    cmdtab.append("-r  {:d}".format(framerate))
    cmdtab.append("-i {:s}".format(bname))
    cmdtab.append('-c:v')
    cmdtab.append('libx264')
    cmdtab.append('-pix_fmt')
    cmdtab.append('yuv420p')
    cmdtab.append('{:s}.mp4'.format(filename))
    command=' '.join(cmdtab)
    sb.call(command,shell=True)
    command="rm -f {:s}*.png".format(filename)
    sb.call(command,shell=True)

cmd.extend("make_movie",make_movie)

def make_state_movie(state=None,fileout='state.gif',nframes=150,imgpx_w=800, imgpx_h=600):
    """make a movie showing a view of the input state"""
    ## PARAMETERS
    if (state == None): state=cmd.count_states()
    rdir = os.path.dirname(os.path.realpath(fileout))
    dirstates_path=os.path.join(rdir,"output_mov_states")
    pngname=os.path.join(dirstates_path,"state")
    #filename=os.path.splitext(fileout)[0]

    if ( not os.path.exists(dirstates_path) ):
        os.makedirs(dirstates_path)

    ## BUILD MOVIE
    cmd.set("state",state)
    #cmd.zoom("all", state=state, complete=1) # reajust view to fit the current state
    cmd.mset("%d x%d" %(state, nframes))
    util.mroll(1,nframes,1)
    cmd.viewport(imgpx_w,imgpx_h)
    cmd.set("ray_trace_frames",1)
    cmd.set("cache_frames",0)
    cmd.mclear()
    cmd.mpng(pngname)

    # MAKE THE GIF AND DELETE PNGs
    #command="convert -delay %.1f -loop 0 %s*.png %s.gif" %(30./100,pngname,filename)
    command="convert -delay %.1f -loop 0 %s*.png %s" %(30./100,pngname,fileout)
    sb.call(command,shell=True)
    #command="rm -f %s*.png" %filename
    command="rm -rf %s" %dirstates_path
    sb.call(command,shell=True)

cmd.extend("make_state_movie",make_state_movie)

def set_view_states(view):
    nstates=cmd.count_states()

    for i in range(nstates):
        cmd.set("state",i+1)
        cmd.set_view(view)
    return

cmd.extend("set_view_states",set_view_states)

def color_123():
    """
    Color firsr 3 monomers
    """

    cmd.select("at0", "resi 1")
    cmd.select("at1", "resi 2")
    cmd.select("at2", "resi 3")
    cmd.deselect()
    cmd.color("red","at0")
    cmd.color("green","at1")
    cmd.color("yellow","at2")

cmd.extend("color_123",color_123)

def emphasize_centromer(nc, radius=1.0, color="green"):
    """
    Empasize the centromer
    """

    cmd.select("centromer", "resi {:d}".format(nc))
    cmd.deselect()
    cmd.color(color,"centromer")
    cmd.alter("centromer","vdw={:.1f}".format(radius))

cmd.extend("emphasize_centromer",emphasize_centromer)

