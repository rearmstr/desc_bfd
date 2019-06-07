#!/usr/bin/env python
#
# LSST Data Management System
# Copyright 2008-2015 LSST Corporation.
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
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.    See the
# GNU General Public License for more details.
#
# You should have received a copy of the LSST License Statement and
# the GNU General Public License along with this program.  If not,
# see <http://www.lsstcorp.org/LegalNotices/>.
#

import lsst.pipe.base as pipeBase
import lsst.afw.table
from lsst.meas.base.pluginRegistry import PluginRegistry
#from lsst.meas.base.baseMeasurement import (BaseMeasurementPluginConfig, BaseMeasurementPlugin,
#                              BaseMeasurementConfig, BaseMeasurementTask)
from lsst.meas.base.noiseReplacer import NoiseReplacer, DummyNoiseReplacer
from lsst.meas.base.sfm import SingleFramePlugin
from .baseMeasure import BaseMeasureConfig, BaseMeasureTask
import lsst.pex.config as pexConfig
__all__ = ("MeasureImageConfig", "MeasureImageTask")


class MeasureImageConfig(BaseMeasureConfig):
    """!
    Config class for single frame measurement driver task.
    """

    plugins = SingleFramePlugin.registry.makeField(
        multi=True,
        default=["base_PixelFlags",
                 "base_SdssCentroid",
                 "base_NaiveCentroid",
                 "base_SdssShape",
                 "base_GaussianFlux",
                 "base_PsfFlux",
                 "base_CircularApertureFlux",
                 "base_SkyCoord",
                 "base_Variance",
                 "base_Blendedness",
                 "base_LocalBackground",
                 ],
        doc="Plugins to be run and their configuration"
    )
    algorithms = property(lambda self: self.plugins, doc="backwards-compatibility alias for plugins")
    prefix = pexConfig.Field(dtype=str, optional=True, default=None, doc="prefix for all measurement fields")
    nGrow = pexConfig.Field(dtype=int, optional=True, default=-1, doc="grow footprints")
    outputDir = pexConfig.Field(dtype=str, optional=True, default='.', doc="where to put output Files")
    outputName = pexConfig.Field(dtype=str, optional=True, default='MOMENTS', doc="prefix for data name")
    undeblended = SingleFramePlugin.registry.makeField(
        multi=True,
        default=[],
        doc="Plugins to run on undeblended image"
    )

## @addtogroup LSST_task_documentation
## @{
## @page MeasureImageTask
## @ref MeasureImageTask_ "MeasureImageTask"
## @copybrief MeasureImageTask
## @}


class MeasureImageTask(BaseMeasureTask):
    """!
    """

    ConfigClass = MeasureImageConfig

    NOISE_SEED_MULTIPLIER = "NOISE_SEED_MULTIPLIER"
    NOISE_SOURCE = "NOISE_SOURCE"
    NOISE_OFFSET = "NOISE_OFFSET"
    NOISE_EXPOSURE_ID = "NOISE_EXPOSURE_ID"

    def __init__(self, **kwds):
        """!
        Initialize the task. Set up the execution order of the plugins and initialize
        the plugins, giving each plugin an opportunity to add its measurement fields to
        the output schema and to record information in the task metadata.
        @param[in,out] schema      lsst.afw.table.Schema, to be initialized to include the
                                   measurement fields from the plugins already
        @param[in,out] algMetadata lsst.daf.base.PropertyList used to record information about
                                   each algorithm.  An empty PropertyList will be created if None.
        @param[in]     **kwds      Keyword arguments forwarded to lsst.pipe.base.Task.__init__
        """
        self.schema = lsst.afw.table.SourceTable.makeMinimalSchema()
        BaseMeasureTask.__init__(self, **kwds)
        self.config.slots.setupSchema(self.schema)	
        self.initializePlugins(schema=self.schema)


    @pipeBase.timeMethod
    def runMeasure(self, inputs, outCat, dataRef):

        noiseImage=None
        exposureId = None#, exposureId=None, beginOrder=None, endOrder=None
        beginOrder = None
        endOrder = None

        footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                      for measRecord in inputs.sources}

        # noiseReplacer is used to fill the footprints with noise and save heavy footprints
        # of the source pixels so that they can be restored one at a time for measurement.
        # After the NoiseReplacer is constructed, all pixels in the exposure.getMaskedImage()
        # which belong to objects in measCat will be replaced with noise
        if True:#self.config.doReplaceWithNoise:
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, inputs.exposure, footprints,
                                          noiseImage=noiseImage, log=self.log, exposureId=exposureId)
            algMetadata = inputs.sources.getMetadata()
            if algMetadata is not None:
                algMetadata.addInt(self.NOISE_SEED_MULTIPLIER, self.config.noiseReplacer.noiseSeedMultiplier)
                algMetadata.addString(self.NOISE_SOURCE, self.config.noiseReplacer.noiseSource)
                algMetadata.addDouble(self.NOISE_OFFSET, self.config.noiseReplacer.noiseOffset)
                if exposureId is not None:
                    algMetadata.addLong(self.NOISE_EXPOSURE_ID, exposureId)
        else:
            noiseReplacer = DummyNoiseReplacer()

        self.runPlugins(noiseReplacer, outCat, inputs.sources, inputs.exposure, beginOrder, endOrder)


    def readInputs(self, dataRef):

        """Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
        of the pixel variance and covariance of the matrix as determined by the Psf
        """
        exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True)
        return lsst.pipe.base.Struct(
            sources = dataRef.get(self.dataPrefix + "src", immediate=True),
            exposure = exposure
            )

    def prepCatalog(self, inputs):
        """Prepare and return the output catalog
        """
        outCat =  lsst.afw.table.SourceCatalog(self.schema)
        srcCat = inputs.sources

        exposurePsf = inputs.exposure.getPsf()
        exposureCalib = inputs.exposure.getCalib()

        # SchemaMapper will transfer ID, Coord, Footprint (we'll overwrite the Coord with that from RefCat)
        mapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
        mapper.addMinimalSchema(self.schema)

        for srcRecord in srcCat:

            outRecord = outCat.addNew()

            # Start by setting some miscellaneous calculated fields
            outRecord.assign(srcRecord, mapper)
            outRecord.setCoord(srcRecord.getCoord())
            outRecord.set('id', srcRecord.getId())

            # Next we determine the pixel region we want to fit.
            if self.config.nGrow > 0:
                outRecord.setFootprint(self.setupFitRegion(inputs.exposure, srcRecord, self.config.nGrow))

        return outCat

    def writeOutputs(self, dataRef, outCat):
        dataRef.put(outCat, "moment", flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)

    def selection(self, source, ref):
        return True

    def preMeasure(self, source, ref, exp):
        return True

    def preMeasureImage(self, source, ref, exp):
        return True

    def runPlugins(self, noiseReplacer, outCat, measCat, exposure, beginOrder=None, endOrder=None):
        """Function which calls the defined measument plugins on an exposure
        Parameters
        ----------
        noiseReplacer : lsst.meas.base.NoiseReplacer
            noiseReplacer to fill sources not being measured with noise.
        measCat : lsst.afw.table.SourceCatalog
            SourceCatalog to be filled with outputs. Must contain all the SourceRecords to be measured (with
            Footprints attached), and have a schema that is a superset of self.schema.
        exposure : lsst.afw.image.ExposureF
            Exposure contaning the pixel data to be measured and the associated PSF, WCS, etc.
        beginOrder : float
            beginning execution order (inclusive): measurements with executionOrder < beginOrder are not
            executed. None for no limit.
        endOrder : float
            ending execution order (exclusive): measurements with executionOrder >= endOrder are not
            executed. None for no limit.
        """
        # First, create a catalog of all parentless sources
        # Loop through all the parent sources, first processing the children, then the parent

         # run premeasure on full image
        measParentCat = measCat.getChildren(0)

        nMeasCat = len(measCat)
        nMeasParentCat = len(measParentCat)
        self.log.info("Measuring %d source%s (%d parent%s, %d child%s) ",
                      nMeasCat, ("" if nMeasCat == 1 else "s"),
                      nMeasParentCat, ("" if nMeasParentCat == 1 else "s"),
                      nMeasCat - nMeasParentCat, ("" if nMeasCat - nMeasParentCat == 1 else "ren"))

        for i, (source, ref) in enumerate(zip(outCat, measCat)):
            self.preMeasureImage(source, ref, exposure)

        for i, (source, ref) in enumerate(zip(outCat, measCat)):
            if not self.selection(source, ref):
                continue
            inserted = False
            try:
                noiseReplacer.insertSource(ref.getId())
                inserted = True
                
                self.preMeasure(source, ref, exposure)
                self.callMeasure(source, exposure, beginOrder=beginOrder, endOrder=endOrder)
            except Exception as e:
                if inserted:
                    noiseReplacer.removeSource(ref.getId())
                continue

            if inserted:
                noiseReplacer.removeSource(ref.getId())
            


        noiseReplacer.end()

        


# #!/usr/bin/env python
# #
# # LSST Data Management System
# # Copyright 2008-2013 LSST Corporation.
# #
# # This product includes software developed by the
# # LSST Project (http://www.lsst.org/).
# #
# # This program is free software: you can redistribute it and/or modify
# # it under the terms of the GNU General Public License as published by
# # the Free Software Foundation, either version 3 of the License, or
# # (at your option) any later version.
# #
# # This program is distributed in the hope that it will be useful,
# # but WITHOUT ANY WARRANTY; without even the implied warranty of
# # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# # GNU General Public License for more details.
# #
# # You should have received a copy of the LSST License Statement and
# # the GNU General Public License along with this program.  If not,
# # see <http://www.lsstcorp.org/LegalNotices/>.
# #

# import numpy
# import os
# import lsst.pex.config
# import lsst.pipe.base
# import lsst.afw.table
# from contextlib import contextmanager
# import lsst.meas.algorithms.algorithmsLib
# from lsst.meas.algorithms.algorithmRegistry import *
# from lsst.meas.algorithms.replaceWithNoise import *

# from .baseMeasure import BaseMeasureConfig, BaseMeasureTask

# __all__ = ("MeasureImageConfig", "MeasureImageTask")


# class MeasureImageConfig(BaseMeasureConfig):
#     algorithms = AlgorithmRegistry.all.makeField(
#         multi=True,
#         default=[],
#         doc="Algorithms that will be run.")
#     doReplaceWithNoise = pexConfig.Field(dtype=bool, default=True, optional=False,
#                                          doc='When measuring, replace other detected footprints with noise?')

#     replaceWithNoise = pexConfig.ConfigurableField(
#         target = ReplaceWithNoiseTask,
#         doc = ("Task for replacing other sources by noise when measuring sources; run when " +
#                "'doReplaceWithNoise' is set."),)

#     prefix = pexConfig.Field(dtype=str, optional=True, default=None, doc="prefix for all measurement fields")
#     nGrow = pexConfig.Field(dtype=int, optional=True, default=-1, doc="grow footprints")
#     outputDir = pexConfig.Field(dtype=str, optional=True, default='.', doc="where to put output Files")
#     outputName = pexConfig.Field(dtype=str, optional=True, default='MOMENTS', doc="prefix for data name")

#     def makeMeasureSources(self, schema, metadata=None):
#         """ Convenience method to make a MeasureSources instance and
#         fill it with the configured algorithms.

#         This is defined in the Config class to support use in unit tests without needing
#         to construct a Task object.
#         """
#         builder = algorithmsLib.MeasureSourcesBuilder(
#             self.prefix if self.prefix is not None else "",
#             False
#             )
#         builder.addAlgorithms(self.algorithms.apply())
#         return builder.build(schema, metadata)


# class MeasureImageTask(BaseMeasureTask):
#     """Copied from meas_multifit

#     Like ProcessImageTask, MeasureImageTask is intended to be used as a base
#     class with CCD and coadd derived-class specializations.  It is run
#     after process[Ccd|Coadd|Eimage].py, and generates a single output
#     catalog with the mapper name 'modelfits'.
#     """

#     ConfigClass = MeasureImageConfig
#     dataPrefix = ""

#     def __init__(self, **kwds):
#         BaseMeasureTask.__init__(self, **kwds)
#         self.measurer = self.config.makeMeasureSources(self.schema, None)
#         if self.config.doReplaceWithNoise:
#             self.makeSubtask('replaceWithNoise')

#     def readInputs(self, dataRef):

#         """Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
#         of the pixel variance and covariance of the matrix as determined by the Psf
#         """
#         exposure = dataRef.get(self.dataPrefix + "calexp", immediate=True)
#         return lsst.pipe.base.Struct(
#             sources = dataRef.get(self.dataPrefix + "src", immediate=True),
#             exposure = exposure
#             )

#     def prepCatalog(self, inputs):
#         """Prepare and return the output catalog
#         """
#         outCat =  lsst.afw.table.SourceCatalog(self.schema)
#         srcCat = inputs.sources

#         exposurePsf = inputs.exposure.getPsf()
#         exposureCalib = inputs.exposure.getCalib()

#         # SchemaMapper will transfer ID, Coord, Footprint (we'll overwrite the Coord with that from RefCat)
#         mapper = lsst.afw.table.SchemaMapper(lsst.afw.table.SourceTable.makeMinimalSchema())
#         mapper.addMinimalSchema(self.schema)

#         for srcRecord in srcCat:

#             outRecord = outCat.addNew()

#             # Start by setting some miscellaneous calculated fields
#             outRecord.assign(srcRecord, mapper)
#             outRecord.setCoord(srcRecord.getCoord())
#             outRecord.set('id',srcRecord.getId())

#             # Next we determine the pixel region we want to fit.
#             if self.config.nGrow > 0:
#                 outRecord.setFootprint(self.setupFitRegion(inputs.exposure, srcRecord, self.config.nGrow))

#         return outCat

#     def writeOutputs(self, dataRef, outCat):
#         dataRef.put(outCat, "moment", flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS)

#     def selection(self, source, ref):
#         return True

#     def preMeasure(self, source, ref, exp):
#         return True

#     def preMeasureImage(self, source, ref, exp):
#         return True

#     def runMeasure(self, inputs, outCat, dataRef):
#         """Need to change this to use contexts
#         """

#         # run premeasure on full image
#         for i, (source, ref) in enumerate(zip(outCat,inputs.sources)):
#                 self.preMeasureImage(source, ref, inputs.exposure)
#         noiseout = self.config.doReplaceWithNoise
#         if noiseout:
#             self.replaceWithNoise.begin(inputs.exposure, inputs.sources)
#         for i, (source, ref) in enumerate(zip(outCat,inputs.sources)):
#             with self.noiseContext(inputs, i, ref, noiseout) as inputs:
#                 if not self.selection(source, ref):
#                     continue

#                 if not self.preMeasure(source, ref, inputs.exposure):
#                     continue

#                 try:
#                     self.measurer.applyWithPeak(source, inputs.exposure, True, 0., float("inf"))
#                 except MemoryError:
#                     raise  # always let MemoryError propagate up, as continuing just causes
#                            # more problems
#                 except Exception as err:
#                     print ('error')
#                     #self.log.warn("Error measuring source %s at %s,%s: %s"
#                     #              % (source.getId(), ref.getX(), ref.getY(), err))
#                     self.log.warn("Error measuring source %s: %s"
#                                   % (source.getId(), err))

#         if noiseout:
#             # Put the exposure back the way it was
#             self.replaceWithNoise.end(inputs.exposure, inputs.sources)

#     @contextmanager
#     def noiseContext(self, inputs, i, ref, noiseout):
#         """Context manager that applies and removes gain
#         """
#         if noiseout:
#             self.replaceWithNoise.insertSource(inputs.exposure, i)

#         try:
#             yield inputs
#         finally:

#             if noiseout:
#                 self.replaceWithNoise.removeSource(inputs.exposure, inputs.sources, ref)


#     def setupFitRegion(self, exposure, source, nGrow, growIsotropic=True):
#             """Given a SourceRecord (with Footprint) and the Exposure it was detected on,
#             return a new Footprint containing the pixels that should be used in a model
#             fit of the given source.
#             """

#             fp = lsst.afw.detection.growFootprint(source.getFootprint(), nGrow, growIsotropic)
#             fp.clipTo(exposure.getBBox(lsst.afw.image.PARENT))
#             return fp

