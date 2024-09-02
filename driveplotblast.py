import numpy as np
import argparse
from argparse import RawTextHelpFormatter
import plotblast as pb

def main():

    parser = argparse.ArgumentParser(formatter_class=RawTextHelpFormatter)
    parser.add_argument("--dstart",type=int,default=0,help="start dump")
    parser.add_argument("--dend",type=int,default=0,help="end dump")
    parser.add_argument("--quant",type=str,default='dens',help='plot quantity')
    parser.add_argument("--run",type=str,default='test',help='run name')
    parser.add_argument("--log",action='store_true',default=False,help='plot logarithmic map')
    parser.add_argument("--save",action='store_true',default=False,help='save map')
    parser.add_argument("--minmax",type=float,nargs=2,default=None,help='min/max quantity')
    parser.add_argument("--axis",type=int,default=0,help='projection axis for 3D')
    parser.add_argument("--surf",action='store_true',default=False,help='average instead of midslice for 3D')
    parser.add_argument("--circ",type=float,default=-1.0,help='draw circle of given radius')
    parser.add_argument("--cmap",type=str,default='viridis',help='color map')
    parser.add_argument("--contour",type=int,default=-1,help='number of contours between minimum and maximum')
    parser.add_argument("--grid",action='store_true',default=False,help='plot grid lines')
    parser.add_argument("--title",action='store_true',default=False,help='show time stamp as title')
    parser.add_argument("--label",type=str,default='None',help='show label in upper right corner')
    parser.add_argument("--vector",type=str,default=None,help='plot vectors [vel/Bcc]')
    parser.add_argument("--samples",type=int,default=16,help='number of sample points for vector field')
    parser.add_argument("--cvmap",type=str,default='black',help='color for vectors/streamlines')
    args  = parser.parse_args()

    path = '../data/'
    #path   = '/21dayscratch/scr/f/h/fheitsch/expgrid/'
    #path   = '/nas/longleaf/home/fheitsch/athena++/staging/expgrid/'
    prob   = 'blast'
    run    = args.run
    field  = args.quant
    out    = 1
    dstart = args.dstart
    dend   = args.dend
    suffix = 'athdf'
    intype = 'hdf'

    stream   = args.vector
    level    = None
    average  = None
    colormap = args.cmap
    if args.minmax is not None:
        vmin     = args.minmax[0]
        vmax     = args.minmax[1]
    else:
        vmin     = None
        vmax     = None
    logc     = args.log
    midplane = not args.surf
    surface  = args.surf
    center   = False
    units    = False
    axis     = args.axis
    

    for d in range(dstart,dend+1):
        pb.map(path,prob,run,out,d,
               quantity=field,vector=args.vector,level=level,
               axis=axis,midplane=midplane,surface=surface,save=args.save,units=units,
               colormap=colormap,vmin=vmin,vmax=vmax,logc=logc,center=center,circ=args.circ,
               contour=args.contour,grid=args.grid,title=args.title,label=args.label,
               samples=args.samples,cvmap=args.cvmap)

main()
