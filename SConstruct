# -*- python -*-
from lsst.sconsUtils import scripts, env
env.Append(CCFLAGS=['-DUSE_EIGEN', '-Wno-reorder', '-Wno-delete-non-virtual-dtor', '-DBALLTREE', '-fopenmp', 
    '-Wno-sign-compare',
    ])
scripts.BasicSConstruct("desc_bfd")

