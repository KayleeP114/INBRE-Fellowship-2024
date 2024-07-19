#!/usr/bin/env python3

#
# Script authored by FM Ytreberg <fmytreberg@gmail.com>
# Last modified Jun 2020
#
# Requires: mdtraj
#


#
# IMPORTS
#

import argparse
import itertools
import os
import re
import subprocess
import sys
import mdtraj as md
import numpy as np


#
# PARSE COMMAND LINE
#

parser = argparse.ArgumentParser(description = 'This script will analyze MD\
        trajectories using MDTraj.')

parser.add_argument('-i',
    metavar = 'TRJFILE',
    required = True,
    type = str,
    help = 'input trajectory file to analyze (.h5, .dcd, .pdb, .xtc, .inpcrd)')
parser.add_argument('-t',
    metavar = 'TOPFILE',
    type = str,
    default = None,
    help = 'input topology file, if needed for trajectory (.pdb, .cif, .psf,\
            .prmtop)')
parser.add_argument('-stride',
    metavar = 'STRIDE',
    type = int,
    default = None,
    help = 'load/analyze every stride-th frame')
parser.add_argument('-single',
    metavar = 'SINGLE',
    type = int,
    default = None,
    help = 'load/analyze a single frame (0 based)')
parser.add_argument('-bfac',
    metavar = 'BFAC',
    type = str,
    default = None,
    help = 'replace b-factor column with information from provided file;\
            only works with single frame')
parser.add_argument('-frames',
    metavar = 'FRAMES',
    type = str,
    help = 'frames to read in 0 based, e.g., "-frames 0:10" will analyze frames\
            1-10, "-frames -1" will analyze last frame, "-frames 100:110:2"\
            will analyze frames 100,102,...,110')
parser.add_argument('-select',
    metavar = 'SELECT',
    type = str,
    help = 'selection criteria for input, analysis and output, e.g.,\
            "all", "backbone", "protein and name CA", "chainid 0 to 5",\
            "chainid 0 2", "chainid 2 and not type H", "residue 18 to 93\
            and not resn ZNB", "name CA and residue 12 15 22"')
parser.add_argument('-super',
    metavar = 'SUPER_SELE',
    type = str,
    help = 'selection for superposition of trajectory (must be within -select)')
parser.add_argument('-ref',
    action = 'store_true',
    help = 'use topology from -t as reference structure for superposition, etc,\
            otherwise will use first frame of trajectory')
parser.add_argument('-superfr',
    metavar = 'SUPERFR',
    type = int,
    default = 0,
    help = '(default "0") snapshot to use as a reference for superposition\
            (e.g., "0", "2"); ignored if -ref is used')
parser.add_argument('-anchor',
    metavar = 'ANCHOR',
    type = str,
    default = None,
    help = 'selection to anchor for image_molecules; must be a single resid\
            or chainid, both 0 based (e.g., "chainid 0", "resid 3")')
parser.add_argument('-info',
    action = 'store_true',
    help = 'output info about trajectory and topology')
parser.add_argument('-hbond',
    metavar = ('SELE1','SELE2', 'FRAC'),
    nargs = 3,
    type = str,
    default = None,
    help = 'find hydrogen bonds between selections (e.g., -hbonds "chainid 0 1"\
            "chain 2"); for a trajectory will return a list with the bonds\
            and the fraction of frames in which they occur if > specified' )
parser.add_argument('-saltbr',
    metavar = ('SELE1','SELE2','FRAC'),
    nargs = 3,
    type = str,
    default = None,
    help = 'find salt bridges between selections (e.g., -saltbr "chainid 0 1"\
            "chain 2"); for a trajectory will return a list with the bonds\
            and the fraction of frames in which they occur if > specified' )
parser.add_argument('-dssp',
    action = 'store_true',
    help = 'output dssp info about trajectory: H=helix, E=beta, C=coil')
parser.add_argument('-rg',
    action = 'store_true',
    help = 'output radius of gyration of -select for trajectory')
parser.add_argument('-rmsd',
    metavar = 'RMSD_SELE',
    type = str,
    default = None,
    help = 'output RMSD of selection relative to -ref or -superfr')
parser.add_argument('-rmsf',
    metavar = 'RMSF_SELE',
    type = str,
    default = None,
    help = 'output RMSF of selection relative to -ref or -superfr')
parser.add_argument('-pi',
    action = 'store_true',
    help = 'output minimum distance to periodic image; uses gmx mindist -pi;\
            only gives meaningful results if you use -select "protein" and\
            choose -super and/or -anchor as needed to make structures whole')
parser.add_argument('-dist',
    metavar = ('SELE1', 'SELE2'),
    nargs = 2,
    type = str,
    help = 'output average, minimum and maximum distances between selections,\
            e.g., "chainid 0" "chainid 1"')
parser.add_argument('-contacts',
    metavar = ('SELE1', 'SELE2', 'CUTOFF', 'TYPE'),
    nargs = 4,
    type = str,
    help = 'output contacts between between selections, e.g., "chainid 0 and\
            name CA" "chainid 1 and name CA" "0.5" "num"; a contact is defined\
            as any atom in selections that is within specified cutoff in nm;\
            type can be "num" for number of contacts or "res" for a list of\
            residues')
parser.add_argument('-o',
    metavar = 'TRJOUT',
    type = str,
    help = 'output trajectory file (will not superpose unless using -super)')
parser.add_argument('-split',
    metavar = 'BASE:FMT',
    type = str,
    help = 'overwrites -o, splits into PDB frames, e.g., -split frame:3 will\
            split with names frame001.pdb, frame002.pdb, etc, -split out:1\
            will split into out1.pdb, out2.pdb, etc')

# Check args
if len(sys.argv) == 1:  # Print help if no arguments are given
    parser.print_help()
    sys.exit(1)
args = parser.parse_args()


#
# FUNCTIONS
#

def selFrames(traj, top):
    '''
    Allows use to select a specific range of frames from the trajectory.
    '''

    index = args.frames.split(':')
    if len(index) == 1:
        traj = traj[int(index[0])]
    elif len(index) == 2:
        traj = traj[int(index[0]):int(index[1])]
    else:
        traj = traj[int(index[0]):int(index[1]):int(index[2])]
    top = traj.topology
    return traj, top


def getHbonds(traj, top):
    '''
    Determines hydrogen bonds between two selections.
    '''

    grp1 = top.select(args.hbond[0])
    grp2 = top.select(args.hbond[1])
    frac_cut = float(args.hbond[2])
    
    '''
    print('# SELE1--SELE2 (donor --> acceptor)')
    hbonds = md.baker_hubbard(traj, freq=freq)
    for hbond in hbonds:
        don = hbond[0]
        acc = hbond[2]
        if don in grp1 and acc in grp2:
            print('{:s} --> {:s}'.format(str(top.atom(don)),
                    str(top.atom(acc))))
        elif acc in grp1 and don in grp2:
            print('{:s} <-- {:s}'.format(str(top.atom(acc)),
                    str(top.atom(don))))
    '''

    hbonds = md.wernet_nilsson(traj)
    hb_arr = []
    for frame in range(0, traj.n_frames):
        hbonds_slice = hbonds[frame]
        for hbond in hbonds_slice:
            don = hbond[0]
            acc = hbond[2]
            if don in grp1 and acc in grp2:
                hb_arr.append('{:d}:{:s} --> {:s}'.format(frame,
                        str(top.atom(don)), str(top.atom(acc))))
            elif acc in grp1 and don in grp2:
                hb_arr.append('{:d}:{:s} <-- {:s}'.format(frame,
                        str(top.atom(acc)), str(top.atom(don))))
    
    hb_arr = [x.split(':')[1] for x in hb_arr]
    hb, ct = np.unique(hb_arr, return_counts=True)
    hb_arr = []
    for h, c in zip(hb, ct):
        frac = c / traj.n_frames
        if frac > frac_cut:  # will not count bonds that are < frac_cut
            hb_arr.append('{:s}  {:.2f}'.format(h, frac))
    hb_sort = sorted(hb_arr, key=lambda x: x.split()[-1], reverse=True)
    
    print('# SELE1--SELE2  frac (donor --> acceptor)')
    for h in hb_sort:
        print(h)


def getSaltBridge(traj, top):
    '''
    Determines salt bridges by looking at distances between acidic and basic
    residues.
    From VMD tool:
    atomselect macro acidic "resname ASP GLU"
    atomselect macro basic "resname ARG HIS LYS HSP"
    '''

    cutoff = 0.4 # 4.0 Angstroms from https://doi.org/10.1093/nar/gkab375
    asp = ['OD1', 'OD2']
    glu = ['OE1', 'OE2']
    arg = ['NH1', 'NH2']
    his = ['ND1', 'NE2']
    lys = ['NZ']
    grp1 = top.select(args.saltbr[0])
    grp2 = top.select(args.saltbr[1])
    frac_cut = float(args.saltbr[2])
    
    # Generate indices for a single acidic O and basic N atoms for both
    # selections
    ac1 = []
    bs1 = []
    for a in grp1:
        atom_prop = str(top.atom(a))
        atom = str(atom_prop.split('-')[1])
        if 'ASP' in atom_prop and any(x in atom for x in asp):
            ac1.append(a)
        elif 'GLU' in atom_prop and any(x in atom for x in glu):
            ac1.append(a)
        elif 'ARG' in atom_prop and any(x in atom for x in arg):
            bs1.append(a)
        elif 'HIS' in atom_prop and any(x in atom for x in his):
            bs1.append(a)
        elif 'LYS' in atom_prop and any(x in atom for x in lys):
            bs1.append(a)
    ac2 = []
    bs2 = []
    for a in grp2:
        atom_prop = str(top.atom(a))
        atom = str(atom_prop.split('-')[1])
        if 'ASP' in atom_prop and any(x in atom for x in asp):
            ac2.append(a)
        elif 'GLU' in atom_prop and any(x in atom for x in glu):
            ac2.append(a)
        elif 'ARG' in atom_prop and any(x in atom for x in arg):
            bs2.append(a)
        elif 'HIS' in atom_prop and any(x in atom for x in his):
            bs2.append(a)
        elif 'LYS' in atom_prop and any(x in atom for x in lys):
            bs2.append(a)

    # Generate distances between acidic and basic atoms
    saltbr = []
    pairs = list(itertools.product(ac1, bs2)) # acidic grp1 - basic grp2
    if len(pairs) > 0:
        distmat = md.compute_distances(traj, pairs)
        row = 0
        for frame in distmat:
            row = row + 1
            for pair, dist in zip(pairs, frame):
                res1 = str(top.atom(pair[0])).split('-')[0]
                res2 = str(top.atom(pair[1])).split('-')[0]
                if res1 is not res2 and dist < cutoff:
                    saltbr.append('{:d}:{:s} --> {:s}'.format(row,
                            res1, res2))
            saltbr = list(set(saltbr))
    pairs = list(itertools.product(bs1, ac2))  # basic grp1 - acidic grp2
    if len(pairs) > 0:
        row = 0
        distmat = md.compute_distances(traj, pairs)
        for frame in distmat:
            row = row + 1
            for pair, dist in zip(pairs, frame):
                res1 = str(top.atom(pair[0])).split('-')[0]
                res2 = str(top.atom(pair[1])).split('-')[0]
                if res1 is not res2 and dist < cutoff:
                    saltbr.append('{:d}:{:s} <-- {:s}'.format(row,
                            res1, res2))
            saltbr = list(set(saltbr))

    # Print out salt bridge pairs
    saltbr = [x.split(':')[1] for x in saltbr]
    sb, ct = np.unique(saltbr, return_counts=True)
    saltbr = []
    for s, c in zip(sb, ct):
        frac = c / traj.n_frames
        if frac > frac_cut:  # will not count bonds that are < frac_cut
            saltbr.append('{:s}  {:.2f}'.format(s, frac))
    sb_sort = sorted(saltbr, key=lambda x: x.split()[-1], reverse=True)
    
    print('# SELE1--SELE2 (acidic --> basic  frac)')
    for s in sb_sort:
        print(s)


def getDSSP(traj, top):
    '''
    Calculates and prints out secondary structure information C=coil, H=helix,
    E=beta.
    '''

    dssp = md.compute_dssp(traj)
    for frame in range(0, traj.n_frames):
        dssp_str = ''
        begin = 0
        end = 0
        for chain in top.chains:
            end += chain.n_residues
            dssp_str += ''.join(str(x) for x in dssp[frame,begin:end])
            dssp_str += ' '
            begin = end
        print(dssp_str)


def getRg(traj):
    '''
    Calculates and prints out radius of gyration of the trajectory.
    '''

    rgs = md.compute_rg(traj)
    frame = 0
    print('# Frame     Rg(nm)')
    for rg in rgs:
        frame += 1
        print('{:6d}  {:10.4f} '.format(frame, rg))


def getRMSD(traj, top, ref):
    '''
    Calculates and prints out the RMSD of the trajectory.
    '''