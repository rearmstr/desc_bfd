#
# LSST Data Management System
# Copyright 2008-2017 LSST/AURA.
#
# This product includes software developed by the
# LSST Project (http://www.lsst.org/).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

# The Plugin classes here are accessed via registries, not direct imports.
#__all__ = ("BfdKMomentConfig")

#from .bfdkmoment import BfdKMomentControl, BfdKMoment

#from lsst.pex.config import makeConfigClass
#import lsst.meas.base


#BfdKMomentConfig = makeConfigClass(BfdKMomentControl)


#class BfdKMomentSingleFrameConfig(lsst.meas.base.SingleFramePluginConfig, BfdKMomentConfig):
#
#    def setDefaults(self):
#        lsst.meas.base.SingleFramePluginConfig.setDefaults(self)
#        BfdKMomentConfig.setDefaults(self)
#
#
#@lsst.meas.base.register("BfdKMoment", apCorrList=[])
#class BfdKMomentSingleFramePlugin(lsst.meas.base.SingleFramePlugin):
#    """Single-frame measurement interface for BfdKMomentAlgorithm.
#    This class simply provides __init__ and measure methods that matched the SingleFramePlugin signatures
#    and delegate to the BfdKMomentAlgorithm's methods.
#    """
#    ConfigClass = BfdKMomentSingleFrameConfig
#
#    @staticmethod
#    def getExecutionOrder():
#        return 3.0
#
#    def __init__(self, config, name, schema, metadata):
#        lsst.meas.base.SingleFramePlugin.__init__(self, config, name, schema, metadata)
#        self.algorithm = BfdKMoment(name, config.makeControl(), schema)
#
#    def measure(self, measRecord, exposure):
#        self.algorithm.measure(measRecord, exposure)
#
#    def fail(self, measRecord, error=None):
#        self.algorithm.fail(measRecord, error.cpp if error is not None else None)

#
#@lsst.meas.base.register("BfdKMoment", apCorrList=[])
#class BfdKMomentForcedPlugin(lsst.meas.base.ForcedPlugin):
#    """Forced measurement interface for BfdKMomentAlgorithm
#    This class simply provides __init__ and measure methods that matched the ForcedPlugin signatures
#    and delegate to BfdKMomentAlgorithm implementations.
#    The BfdKMoment algorithm currently cannot be run in forced mode when the measurement WCS is different
#    from the reference WCS (as is the case in CCD forced photometry).  This is a temporary limitation
#    that will be addressed on DM-5405.
#    BfdKMoment forced measurement when the measurement image is the same as the reference image should be
#    almost -- but not quite -- identical to unforced measurement.  The primary difference is that
#    the final fit region from the reference measurement will be used for the initial fit in forced mode
#    as well as the exp, dev, and combined exp+dev fits
#    """
#    ConfigClass = BfdKMomentConfig
#
#    @staticmethod
#    def getExecutionOrder():
#        return 3.0
#
#    def __init__(self, config, name, schemaMapper, metadata):
#        lsst.meas.base.ForcedPlugin.__init__(self, config, name, schemaMapper, metadata)
#        self.algorithm = BfdKMoment(name, config.makeControl(), schemaMapper)
#
#    def measure(self, measRecord, exposure, refRecord, refWcs):
#        if refWcs != exposure.getWcs():
#            raise lsst.meas.base.FatalAlgorithmError(
#                "BfdKMoment forced measurement currently requires the measurement image to have the same"
#                " Wcs as the reference catalog (this is a temporary limitation)."
#            )
#        self.algorithm.measure(measRecord, exposure)
#
#    def fail(self, measRecord, error=None):
#        self.algorithm.fail(measRecord, error.cpp if error is not None else None)
#
