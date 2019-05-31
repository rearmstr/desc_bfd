
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
# see <http://www.lsstcorp.org/LegalNotices/>.__
#

import lsst.pipe.base
import lsst.pex.config
import numpy as np
import os


import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom
from lsst.meas.base.noiseReplacer import NoiseReplacer, DummyNoiseReplacer

from .measureImage import MeasureImageTask

__all__ = ("MeasureMCShearConfig", "MeasureMCShearTask")

class MeasureMCShearConfig(MeasureImageTask.ConfigClass):
	"""Config for ProcessCoadd"""
	coaddName = lsst.pex.config.Field(
		doc = "coadd name: typically one of deep or goodSeeing",
		dtype = str,
		default = "deep",
	)
	oldHSC = lsst.pex.config.Field(
		doc = "allow ability to run on old HSC data",
		dtype = bool,
		default = False,
	)
	overWrite = lsst.pex.config.Field(
		doc = "Overwrite existing file?",
		dtype = bool,
		default = True,
	)
	maxArea = lsst.pex.config.Field(
		doc = "Maximum number of pixels in footprint to use",
		dtype = int,
		default = 10000,
	)
	multiplicity = lsst.pex.config.Field(
		doc = "Generate this number of sources for each object",
		dtype = int,
		default = 1,
	)
	outdir = lsst.pex.config.Field(
		doc = "Generate this number of sources for each object",
		dtype = str,
		default = './',
	)

class MeasureMCShearTask(MeasureImageTask):
	"""Specialization of MeasureImageTask for running on calexps
	"""

	_DefaultName = "MeasureMCShear"
	ConfigClass = MeasureMCShearConfig

	def __init__(self, **kwargs):
		MeasureImageTask.__init__(self, **kwargs)
		self.dataPrefix = self.config.coaddName + "Coadd_"
		


	@classmethod
	def _makeArgumentParser(cls):
		parser = lsst.pipe.base.ArgumentParser(name=cls._DefaultName)
		parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
							   ContainerClass=CoaddDataIdContainer)
		return parser

	def initialCheck(self, dataRef):
		"""Check to see if we should overwrite existing file
		If file already exists this will return false and exit the processing
		"""
		if self.config.overWrite is False:
			 if dataRef.datasetExists(self.dataPrefix+"moment"):
					self.log.info('Moment already exists %s. skipping' % dataRef.dataId)
					return False

	def readInputs(self, dataRef):
		"""Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
		of the pixel variance and covariance of the matrix as determined by the Psf
		"""
		name = self.config.coaddName + "Coadd_calexp"
		if self.config.oldHSC:
			name += "_hsc"
		exposure = dataRef.get(name, immediate=True)
		return lsst.pipe.base.Struct(
			sources=dataRef.get(self.dataPrefix + "meas", immediate=True),
			exposure=exposure
			)

	def selection(self, source, ref):
		childName = 'deblend_nchild'
		satName = 'flags_pixel_saturated_center'
		if self.config.oldHSC is False:
			childName = 'deblend_nChild'
			satName = 'base_PixelFlags_flag_saturated'
		if ref.getParent()==0 and ref.get(childName)>0:
			return False
		if ref.getFootprint().getArea() > self.config.maxArea:
			return False
		if ref.getFootprint().getArea() == 0:
			return False
		if ref.get(satName):
			return False

		return True

	def prepCatalog(self, inputs):
		"""Prepare and return the output catalog
		"""

		outCat =  lsst.afw.table.SourceCatalog(self.schema)
		#mcCat =  lsst.afw.table.SourceCatalog(self.schema)
		metadata = dafBase.PropertyList()
		outCat.getTable().setMetadata(metadata)
		srcCat = inputs.sources

		exposurePsf = inputs.exposure.getPsf()
		exposureCalib = inputs.exposure.getCalib()

		# SchemaMapper will transfer ID, Coord, Footprint
		mapper = lsst.afw.table.SchemaMapper(srcCat.schema, self.schema)

		# create output catalog
		iRecords = []
		for i, srcRecord in enumerate(srcCat):
			if self.selection(srcRecord, srcRecord) is False:
				continue
			for j in range(self.config.multiplicity):
				outRecord = outCat.addNew()
				iRecords.append(i)

				outRecord.assign(srcRecord, mapper)
				outRecord.setCoord(srcRecord.getCoord())
				outRecord.setId(srcRecord.getId())

		return outCat,iRecords
	
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

		outCat,indices = self.prepCatalog(inputs)

		self.log.info("Measuring catalog")
		self.runMeasure(inputs, outCat, dataRef, indices) 
		self.log.info("Writing outputs")
		self.writeOutputs(dataRef, outCat)
		# I don't think I need to return anything here
		return None#lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)


	def runMeasure(self, inputs, outCat, dataRef, indices):

		noiseImage=None
		exposureId = None
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

		self.runPlugins(noiseReplacer, outCat, indices, inputs.sources, inputs.exposure, beginOrder, endOrder)
	
	def runPlugins(self, noiseReplacer, outCat, indices, measCat, exposure, beginOrder=None, endOrder=None):
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
		self.log.info("Measuring %d source%s (%d parent%s, %d child%s) %d output",
					  nMeasCat, ("" if nMeasCat == 1 else "s"),
					  nMeasParentCat, ("" if nMeasParentCat == 1 else "s"),
					  nMeasCat - nMeasParentCat, ("" if nMeasCat - nMeasParentCat == 1 else "ren"),
					  len(outCat))
		
		for i, (source, index) in enumerate(zip(outCat, indices)):
			self.log.debug("Premeasure for object %d", i)
			self.preMeasureImage(source, measCat[index], exposure)

		for i, (source, index) in enumerate(zip(outCat, indices)):
			ref = measCat[index]
			self.log.debug("Measure object %d from ref %d", i, index)
			if not self.selection(source, ref):
				self.log.debug('Failed selection')
				continue
			inserted = False
			try:
				noiseReplacer.insertSource(ref.getId())
				inserted = True
				
				self.preMeasure(source, ref, exposure)
				self.callMeasure(source, exposure, beginOrder=beginOrder, endOrder=endOrder)
			except Exception as e:
				self.warn('Error %s', e)
				if inserted:
					noiseReplacer.removeSource(ref.getId())
				continue

			if inserted:
				noiseReplacer.removeSource(ref.getId())
			


		noiseReplacer.end()

	def preMeasureImage(self, source, ref, exp, setFlags=True):
		return True


	def writeOutputs(self, dataRef, outCat):
		#dataRef.put(outCat, self.dataPrefix+"moment", flags=afwTable.SOURCE_IO_NO_HEAVY_FOOTPRINTS)
		mask = np.isfinite(outCat['g1'])
		outCat = outCat.subset(mask)
		if os.path.exists(self.config.outdir) is False:
			os.makedirs(self.config.outdir)
		outCat.writeFits('%s/mc_sample_%s_%s_%s.fits'%(self.config.outdir, dataRef.dataId['tract'],
													dataRef.dataId['patch'],
													dataRef.dataId['filter']
			),flags=afwTable.SOURCE_IO_NO_HEAVY_FOOTPRINTS)

	def getPreviousTaskClass(self):
		return MeasureMCShearTask

	def _getConfigName(self):
		"""Return the name of the config dataset
		"""
		return "%s_measureCoadd_config" % (self.config.coaddName,)

	def _getMetadataName(self):
		"""Return the name of the metadata dataset
		"""
		return "%s_measureCoadd_metadata" % (self.config.coaddName,)





