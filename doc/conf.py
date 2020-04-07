"""Sphinx configuration file for an LSST stack package.

This configuration only affects single-package Sphinx documenation builds.
"""

from documenteer.sphinxconfig.stackconf import build_package_configs
import lsst.desc.bfd


_g = globals()
_g.update(build_package_configs(
    project_name='desc_bfd',
    version=lsst.desc.bfd.version.__version__))
