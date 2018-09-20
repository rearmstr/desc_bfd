#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2013 LSST Corporation.
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

import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.table
from lsst.meas.base.baseMeasurement import (BaseMeasurementPluginConfig, BaseMeasurementPlugin,
                                            BaseMeasurementConfig, BaseMeasurementTask)

__all__ = ("BaseMeasureConfig", "BaseMeasureTask")

class BaseMeasureConfig(BaseMeasurementConfig):
    maxObjects = lsst.pex.config.Field(
        dtype=int,
        default=None,
        optional=True,
        doc="If not None, clip the catalog and process only this many objects (for fast-debug purposes)"
    )
    firstObject = lsst.pex.config.Field(
        dtype=int,
        default=0,
        optional=True,
        doc="If maxObjects not None, it will start on this entry in the catalog (for fast-debug purposes)"
    )

class BaseMeasureTask(BaseMeasurementTask):
    """An intermediate base class for top-level BFD moment tasks.  This is mostly copied from meas_multifit.
    Subclasses (for different kinds of input data) must implement three methods:
    readInputs(), prepCatalog(), makeLikelihood(), and writeOutputs().
    """
    ConfigClass = BaseMeasureConfig

    def __init__(self,  schema=None, **kwargs):
        """Initialize the measurement task, including the modelfits catalog schema,
        the model, prior, and calib objects, and the fitter subtask.
        """

        BaseMeasurementTask.__init__(self, **kwargs)
        if schema is None:
            self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()

    def run(self, dataRef):
        """Main driver
        """
        if self.initialCheck(dataRef) is False:
            return None

        self.log.info("Reading inputs")
        inputs = self.readInputs(dataRef)
        self.log.info("Preparing catalog")
        if self.config.maxObjects is not None:
            first = self.config.firstObject
            last = min(self.config.firstObject + self.config.maxObjects,len(inputs.sources))
            inputs.sources = inputs.sources[first:last]

        outCat = self.prepCatalog(inputs)

        self.log.info("Measuring catalog")
        self.runMeasure(inputs, outCat, dataRef) 
        self.log.info("Writing outputs")
        self.writeOutputs(dataRef, outCat)
        # I don't think I need to return anything here
        return None#lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

    def prep(self):
        pass

    def readInputs(self, dataRef):
        """Read task inputs using the butler.

        The form of the return value is subclass-dependent; it will simply be passed unmodified
        to prepCatalog() and makeLikelihood()
        """
        raise NotImplementedError()

    def initialCheck(self, dataRef):
        return True

    def prepCatalog(self, inputs):
        """Prepare and return the output catalog
        """
        raise NotImplementedError()

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        raise NotImplementedError()

    def runMeasure(self, inputs, outRecord, dataRef):
        """
        """
        raise NotImplementedError()
