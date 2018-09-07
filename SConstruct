# -*- python -*-
from lsst.sconsUtils import scripts, env
scripts.BasicSConstruct("desc_bfd")
env.Append(CCFLAGS = ['-DFFLOAT',# To make bfd code happy
                      '-fopenmp','-g',
                      '-Wno-reorder','-Wno-comment', '-Wno-sign-compare','-Wno-parentheses',# To supress warnings from bfd code
                      '-Wno-unused-but-set-variable',
                      '-Wno-delete-non-virtual-dtor', #KWeight warning
                      '-Wno-unused-variable', # annoying afw warnings
                      '-Wno-unused-local-typedefs',# boost warnings
                      '-Wno-misleading-indentation'
                      ])
