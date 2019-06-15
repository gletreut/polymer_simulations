from pymol.cgo import *
from pymol import cmd, stored, util
import subprocess as sb
import os,sys
from functions import get_constrains
import numpy as np

def polymer_colloids_standard(path=["pol.xyz"], polymer_color='tv_blue', monomer_size=1.0, sphere_transparency=0.0, showspheres=True):
    ## LOADING
    cmd.reinitialize()
    cmd.set("max_threads",1)
    obj=[os.path.splitext(p)[0].split('/')[-1] for p in path]
    ntypes=len(obj)
    print "objects to load: %d" %ntypes

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

def polymer_colloids(mydict,path="all_atoms.xyz"):
    """
    typesdict is a dictionary with LAMMPS type assigned to a given
    color index.
    e.g. mydict['types']={'1':1,'2':2,'3':2}
    """
    ## INPUT
    if ('types' in mydict):
        typesdict=mydict['types']
        ntypes=max(typesdict.values()) # number of color types
    else:
        typesdict={}
        ntypes=0

    if ('radius' in mydict):
        radiusdict=mydict['radius']
    else:
        radiusdict={}

    ## LOADING
    cmd.reinitialize()
    for p in path.split():
        cmd.load(p)

    ## COLORS DEFINITION
    base_color="black" # the first color is black
    cmd.color(base_color,'all')
    #col=[pol_color]
    col=[]
    for j in range(ntypes):
        hsv= (240-j*240/(ntypes-1),1,1)
        # convert to rgb
        rgb=hsv_to_rgb(hsv)
        # define color
        col_name="mycolor%d" %j
        cmd.set_color(col_name,rgb)
        col.append(col_name)

    ## MAKE THE COLORING
    for j in typesdict.keys():
        cind=typesdict[j]-1
        type_color=col[cind]
        name="type%d" %j
        name_cmd="name %d" %j
        cmd.select(name, name_cmd)
        cmd.color(type_color,name)

    cmd.deselect()

    ## RADIUS
    cmd.alter("all","vdw=0.5")
    for j in radiusdict.keys():
        name="type%d" %j
        radius=float(radiusdict[j])
        cmd.alter(name,"vdw=%.2f" %radius)


    ## DISPLAY FORMATTING
    cmd.hide("everything","all")
    cmd.show("spheres","all")
    cmd.bg_color("white")
    cmd.set("orthoscopic",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("depth_cue",0)
    cmd.zoom("all",complete=0)

"""
def polymer_colloids(ntypes,path="all_atoms.xyz"):
    ## LOADING
    cmd.reinitialize()
    #cmd.set("max_threads",1)
    #path="all_atoms.xyz"
    obj=[os.path.splitext(path)[0]]
    typ_ind=[i+2 for i in range(2*ntypes)]
    cmd.load(path)

    # COLORS DEFINITION
    pol_color="black" # polymer is black
    col=[]
    for j in range(ntypes):
        hsv= (240-j*240/(ntypes-1),1,1)
        # convert to rgb
        rgb=hsv_to_rgb(hsv)
        # define color
        col_name="mycolor%d" %j
        cmd.set_color(col_name,rgb)
        col.append(col_name)

    # MAKE THE COLORING
    spol="pol"
    name="name %d" %(1)
    cmd.select(spol, name)
    cmd.color(pol_color,spol)

    for j in range(ntypes):
        ind1=typ_ind[j]
        ind2=typ_ind[j+ntypes]
        scol="col%d" %(j+1)
        name="name %d+%d" %(ind1,ind2)
        col_color=col[j]
        cmd.select(scol,name)
        cmd.color(col_color,scol)
    cmd.deselect()

    ## DISPLAY FORMATTING
    cmd.hide("everything","all")
    cmd.alter("all","vdw=0.5")
    cmd.show("spheres","all")
    cmd.bg_color("white")
    cmd.set("orthoscopic",1)
    cmd.set("ray_trace_mode",0)
    cmd.set("ray_trace_fog",0)
    cmd.set("depth_cue",0)
    cmd.zoom("all",complete=0)
#"""

cmd.extend("polymer_colloids",polymer_colloids)

def color_pairs(bonds, line_width=1.0):
    """
    Given a list of bonds:
        bonds[i]=[at1,at2,num bond type]
    draw dashed lines in-between those atoms
    """

    forbid=[1]      # do not consider chain bonds
    nbonds=len(bonds)
    if (nbonds == 0):
        print "list of bonds is empty!"
        return

    # sort by increasing bond types
    idx=np.argsort(bonds[:,2])
    bonds=bonds[idx]

    # determine number values taken by the bonds
    btypes=np.unique(bonds[:,2])
    ntypes=len(btypes)
    groups={}
    for bt in btypes:
        if not (bt in forbid):
            groups[bt]=[]

    ## COLORS DEFINITION
    if (len(groups) == 0):
        return
    col={}
    keys = np.sort(groups.keys())
    ng = len(keys)
    for j,k in enumerate(keys):
#        hsv= (240-j*240/(len(groups)-1),1,1)
        hsv= (240-j*240/(ng-1),1,1)
        # convert to rgb
        rgb=hsv_to_rgb(hsv)
        # define color
        col_name="mycolor%d" %(k)
        cmd.set_color(col_name,rgb)
        col[k] = (col_name)

    ## MAKE THE COLORING
    cmd.unbond("all","all")
    for n in range(nbonds):
        at1=bonds[n,0]
        at2=bonds[n,1]
        bt = bonds[n,2]
        if not (bt in groups.keys()):
            continue
        groups[bt].append(at1)
        groups[bt].append(at2)
        kind=np.ravel(np.where(btypes==bonds[n,2]))[0]  # index in kval: 0 <= ind <= nkval-1
        cmd.bond("id {:d}".format(at1), "id {:d}".format(at2))
        #mycol=col[kind]
        #mycmd="id {:d} id {:d}".format(at1,at2)
        #print "id {:d} id {:d} col={:s}".format(at1,at2,mycol)
        #cmd.color(mycol,mycmd)

    for bt in groups:
        group = np.unique(groups[bt])
        myid = "+".join(["{:d}".format(id) for id in group])
        mygroup = "group{:d}".format(bt)
        mycol=col[bt]
        cmd.select(mygroup,"id {}".format(myid))
        cmd.set_bond("line_color",mycol,mygroup)
        #cmd.set_bond("line_width",line_width,mygroup)
        #cmd.set_bond("line_width",5,mygroup)
        cmd.show("lines",mygroup)

    cmd.deselect()
    return

cmd.extend("color_pairs",color_pairs)

def color_pairs_v2(bonds):
    """
    color the pairs given in the bonds vector:
    bonds[i]=[at1,at2,num bond type]
    """

    nbonds=len(bonds)
    if (nbonds == 0):
        print "list of bonds is empty!"
        return

    # sort by increasing bond types
    idx=np.argsort(bonds[:,2])
    bonds=bonds[idx]

    # determine number values taken by the bonds
    btypes=np.unique(bonds[:,2])
    ntypes=len(btypes)

    ## COLORS DEFINITION
    col=[]
    for j in range(ntypes):
        hsv= (240-j*240/(ntypes-1),1,1)
        # convert to rgb
        rgb=hsv_to_rgb(hsv)
        # define color
        col_name="mycolor%d" %(j)
        cmd.set_color(col_name,rgb)
        col.append(col_name)

    ## MAKE THE COLORING
    for n in range(nbonds):
        at1=bonds[n,0]
        at2=bonds[n,1]
        kind=np.ravel(np.where(btypes==bonds[n,2]))[0]  # index in kval: 0 <= ind <= nkval-1
        mycol=col[kind]
        mycmd="id {:d} id {:d}".format(at1,at2)
        #print "id {:d} id {:d} col={:s}".format(at1,at2,mycol)
        cmd.color(mycol,mycmd)

    cmd.deselect()

    return

cmd.extend("color_pairs_v2",color_pairs_v2)

def color_pairs_v1(params,refdir='.'):
    """
    color the pairs given in params
    if params contains a 'constrains' keyword
    """

    mydict={}
    # check that there is an input contact matrix
    if not ( ('cons_pairs' in params) and ('cons_values' in params) ):
        return

    # read pairs of constraints
    pairs=[]
    idx=[]
    filein=os.path.join(refdir,params['cons_pairs'][0])
    fin=open(filein)
    while True:
        line=fin.readline()
        if (line == ""):
            break
        else:
            try:
                tab=line.split()
                i=int(tab[0])
                j=int(tab[1])
                ind=int(tab[2])
                pairs.append(np.array([i,j]))
                idx.append(ind)
            except ValueError as e:
                print e
                pass
    fin.close()
    pairs=np.array(pairs)
    idx=np.array(idx)

    # read k values
    kval=[]
    filein=os.path.join(refdir,params['cons_values'][0])
    fin=open(filein)
    while True:
        line=fin.readline()
        if (line == ""):
            break
        else:
            tab=line.split()
            try:
                ind=int(tab[0])
                val=float(tab[1])
                kval.append(val)
            except ValueError as e:
                print e
                pass
    fin.close()
    kval=np.array(kval)

    npairs=len(pairs)
    nkval=len(kval)

    if (npairs == 0):
        return

    ## COLORS DEFINITION
    col=[]
    ntypes=nkval+1    # +1 to avoid blue as first color
    for j in range(ntypes):
        hsv= (240-j*240/(ntypes-1),1,1)
        # convert to rgb
        rgb=hsv_to_rgb(hsv)
        # define color
        col_name="mycolor%d" %(j)
        cmd.set_color(col_name,rgb)
        col.append(col_name)

    ## MAKE THE COLORING
    for j in range(npairs):
        n1=pairs[j][0]+1
        n2=pairs[j][1]+1
        ind=idx[j]+1  # index in kval: 0 <= ind <= nkval-1
        mycol=col[ind]
        mycmd="id %d" %(n1)
        cmd.color(mycol,mycmd)
        mycmd="id %d" %(n2)
        cmd.color(mycol,mycmd)

    cmd.deselect()

    return

cmd.extend("color_pairs_v1",color_pairs_v1)

def hsv_to_rgb(hsv):

        h = float(hsv[0])
        s = float(hsv[1])
        v = float(hsv[2])

        if( s == 0 ) :
                #achromatic (grey)
                r = g = b = v

        else:
                # sector 0 to 5
                h = h/60.
                i = int(h)
                f = h - i                       # factorial part of h
                #print h,i,f
                p = v * ( 1 - s )
                q = v * ( 1 - s * f )
                t = v * ( 1 - s * ( 1 - f ) )

                if i == 0:
                        (r,g,b) = (v,t,p)
                elif i == 1:
                        (r,g,b) = (q,v,p)
                elif i == 2:
                        (r,g,b) = (p,v,t)
                elif i == 3:
                        (r,g,b) = (p,q,v)
                elif i == 4:
                        (r,g,b) = (t,p,v)
                elif i == 5:
                        (r,g,b) = (v,p,q)
                else:
                        (r,g,b) = (v,v,v)
                        print "error, i not equal 1-5"

        return [r,g,b]

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
    print "number of states=%d" %nmax
    nt=npts-1
    delta=(nmax-nmin)/float(nt)

    cmd.zoom("all", state=0, complete=1) # reajust view to fit all states (0 means all states)
    for i in range(npts):
        k=int(nmax-i*delta)
	print "state %04d" %k
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

def make_state_movie(state=None,rdir='.',nframes=150,imgpx_w=800, imgpx_h=600):
    """make a movie showing a view of the input state"""
    ## PARAMETERS
    if (state == None): state=cmd.count_states()
    dirpath=os.path.join(rdir,"output_mov_state")
    filename=os.path.join(dirpath,"state")

    if ( not os.path.exists(dirpath) ):
        os.makedirs(dirpath)

    ## BUILD MOVIE
    cmd.set("state",state)
    cmd.zoom("all", state=state, complete=1) # reajust view to fit the current state
    cmd.mset("%d x%d" %(state, nframes))
    util.mroll(1,nframes,1)
    cmd.viewport(imgpx_w,imgpx_h)
    cmd.set("ray_trace_frames",1)
    cmd.set("cache_frames",0)
    cmd.mclear()
    cmd.mpng(filename)

    # MAKE THE GIF AND DELETE PNGs
    command="convert -delay %.1f -loop 0 %s*.png %s.gif" %(30./100,filename,filename)
    sb.call(command,shell=True)
    command="rm -f %s*.png" %filename
    sb.call(command,shell=True)

cmd.extend("make_state_movie",make_state_movie)
