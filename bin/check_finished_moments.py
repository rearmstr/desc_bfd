#! /usr/bin/env python

import glob
import os
import re
import argparse

parser = argparse.ArgumentParser(description='Check the number of finished moments')
parser.add_argument('--meas_dir',default='s16a_wide', help='rerun of meas files')
parser.add_argument('--moment_dir', default='s16a_wide_moments', help='rerun of moment files')
parser.add_argument('--root', default='/tigress/HSC/HSC/rerun/rearmstr/', help='root directory')
parser.add_argument('--verbose', default=False, type=bool, help='print more output')

args = parser.parse_args()

tracts = glob.glob('%s/%s/deepCoadd-results/HSC-I/*' % (args.root, args.meas_dir))
moment_dir = '%s/%s/deepCoadd-results/HSC-I' % (args.root, args.moment_dir)
patches=[]
for d in tracts:
    if args.verbose:
        print('Checking tract: ', d)
    m,tract=os.path.split(d)
    meas_files=glob.glob(m+'/'+tract+'/*/meas-*fits*')

    meas_patches=[]
    for file in meas_files:
        patch=re.search('\/(\d,\d)\/',file).groups()[0]
        meas_patches.append(patch) 
    meas_patches.sort()
    moment_files=glob.glob(moment_dir+'/'+tract+'/*/moment-*fits*')
    if args.verbose:
        print(' Found %d patches' % len(meas_patches))
        print(' Found %d moment files' % len(moment_files))

    if len(moment_files) == len(meas_patches):
        continue

    
    moment_patches=[]
    for file in moment_files:
        patch=re.search('\/(\d,\d)\/',file).groups()[0]
        moment_patches.append(patch)
    moment_patches.sort()

    missing=[]
    for p in meas_patches:
        if p not in moment_patches:
            missing.append(p)
    print(tract,len(meas_files),len(moment_files),missing)
