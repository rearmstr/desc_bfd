
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
import numpy


import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.detection as afwDet
import lsst.afw.geom as afwGeom

def getMaskFromFootprint(mask, footprint, bitmask):
    """Function to get which of pixels have a specific mask in footprint
    """
    bbox = footprint.getBBox()
    bbox.clip(mask.getBBox(afwImage.PARENT))
    fp = afwImage.Mask(bbox)
    subMask = mask.Factory(mask, bbox, afwImage.PARENT)
    #fp.getFootprint().span.setMask(
    afwDet.setMaskFromFootprint(fp, footprint, 0x1)
    
    return ((subMask.getArray()&bitmask==0) & (fp.getArray())) > 0
    #return (subMask.getArray()==bitmask & fp.getArray()) > 0


from .measureImage import MeasureImageTask

try:
    from lsst.meas.mosaic import applyMosaicResults
except ImportError:
    applyMosaicResults = None

__all__ = ("MeasureCoaddConfig", "MeasureCoaddTask")

class MeasureCoaddConfig(MeasureImageTask.ConfigClass):
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
    noiseGrow = lsst.pex.config.Field(
        doc = "how large should the footprint be grown to measure the noise",
        dtype = int,
        default = 5,
    )
    calculateVariance = lsst.pex.config.Field(
        doc = "Should the variance be calculated for each individual object?",
        dtype = bool,
        default = True,
    )
    useCorrelatedNoise = lsst.pex.config.Field(
        doc = "Should we use an estimate of the correlated noise, where do we get this?",
        dtype = bool,
        default = False,
    )
    correlatedNoiseFile = lsst.pex.config.Field(
        doc = "File that contains list of k values and power spectrum",
        dtype = str,
        default = '',
    )
    maxArea = lsst.pex.config.Field(
        doc = "Maximum number of pixels in footprint to use",
        dtype = int,
        default = 10000,
    )
    minVariancePixels = lsst.pex.config.Field(
        doc = "Minimum number of pixels to use for the variance calculation",
        dtype = int,
        default = 20,
    )
    overWrite = lsst.pex.config.Field(
        doc = "Overwrite existing file?",
        dtype = bool,
        default = True,
    )

class MeasureCoaddTask(MeasureImageTask):
    """Specialization of MeasureImageTask for running on calexps
    """

    _DefaultName = "measureCoadd"
    ConfigClass = MeasureCoaddConfig

    def __init__(self, **kwargs):
        MeasureImageTask.__init__(self, **kwargs)
        self.dataPrefix = self.config.coaddName + "Coadd_"

        if self.config.calculateVariance:
            self.varKey = self.schema.addField("noise_variance", doc="variance of noise", type=float)
            self.meanKey = self.schema.addField("noise_mean", doc="mean of noise", type=float)

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
        
        if ref.getParent()==0 and ref.get('deblend_nchild')>0:
            source.set('bfd.flags.parent',True)
            source.set('bfd.flags',True)
            return False
        if ref.getFootprint().getArea() > self.config.maxArea:
            source.set('bfd.flags.too-big',True)
            source.set('bfd.flags',True)
            return False
        if ref.getFootprint().getArea() == 0:
            source.set('bfd.flags.footprint-empty',True)
            source.set('bfd.flags',True)
            return False
        if ref.get('flags_pixel_saturated_center'):
            source.set('bfd.flags.saturated.center',True)
            source.set('bfd.flags',True)
            return False
        if numpy.isnan(source.get('noise_variance')):
            return False
        return True

    def prepCatalog(self, inputs):
        """Prepare and return the output catalog
        """

        outCat =  lsst.afw.table.SourceCatalog(self.schema)
        metadata = dafBase.PropertyList()
        outCat.getTable().setMetadata(metadata)
        srcCat = inputs.sources

        exposurePsf = inputs.exposure.getPsf()
        exposureCalib = inputs.exposure.getCalib()

        # SchemaMapper will transfer ID, Coord, Footprint
        mapper = lsst.afw.table.SchemaMapper(srcCat.schema, self.schema)

        # Calculate the noise variance in the image
        sctrl = lsst.afw.math.StatisticsControl()
        sctrl.setNumIter(3)
        sctrl.setNumSigmaClip(5.0)
        mask = inputs.exposure.getMaskedImage().getMask()
        image =  inputs.exposure.getMaskedImage().getImage()
        #sctrl.setAndMask(map(lambda x, y: x | mask.getPlaneBitMask(y),
        #                       afwImage.Mask().getMaskPlaneDict().keys(), 0x0))
        sctrl.setNanSafe(True)
        stats = lsst.afw.math.makeStatistics(image,
                                             lsst.afw.math.VARIANCECLIP | lsst.afw.math.MEANCLIP)
        variance = stats.getValue(lsst.afw.math.VARIANCECLIP);
        mean = stats.getValue(lsst.afw.math.MEANCLIP); 
        outCat.getTable().getMetadata().set('noise_mean',mean)
        outCat.getTable().getMetadata().set('noise_variance',variance)

        if self.config.useCorrelatedNoise:
            data = numpy.genfromtxt(self.config.correlatedNoiseFile)
            outCat.getTable().getMetadata().set('kData',data[:,0])
            outCat.getTable().getMetadata().set('psData',data[:,1])

        if self.config.calculateVariance:
            x0 = inputs.exposure.getXY0().getX()
            y0 = inputs.exposure.getXY0().getY()
            self.xy0 = afwGeom.Extent2I(-x0,-y0)

        for srcRecord in srcCat:

            outRecord = outCat.addNew()
            #outRecord.set(self.psfSizeKey, psfSize)

            # Start by setting some miscellaneous calculated fields
            outRecord.assign(srcRecord, mapper)
            outRecord.setCoord(srcRecord.getCoord())
            # Increase footprint region?
            if self.config.nGrow > 0:
                outRecord.setFootprint(self.setupFitRegion(inputs.exposure, srcRecord, self.config.nGrow))

        return outCat

    def preMeasureImage(self, source, ref, exp, setFlags=True):

        if self.config.calculateVariance:
            # Calculate the noise around each object, maybe I should get this from the variance plane?
            image = exp.getMaskedImage().getImage()


            maskPlanes = ["BAD", "CR", "DETECTED", "DETECTED_NEGATIVE", "NO_DATA", "SAT",
                          "EDGE", "SUSPECT", "INTRP"]
            badBit = afwImage.Mask().getPlaneBitMask(maskPlanes)
            failed = True
            meanSrc = -1
            varSrc = -1
            # try multiples of the noiseGrow parameter until there are a sufficient number of pixels
            for scale in (1, 2, 4, 8):

                fp = afwDet.Footprint(source.getFootprint())
                fp.dilate(scale*self.config.noiseGrow, afwGeom.Stencil.MANHATTAN)

                if fp.getArea() == 0:
                    continue

                fp.clipTo(exp.getBBox(lsst.afw.image.PARENT))
                box = fp.getBBox()

                span = fp.spans.intersectNot(exp.getMaskedImage().getMask()[box], badBit)
                maskImage = afwImage.Mask(box, 0)
                span.setMask(maskImage, 0x1)
                mask = maskImage.array==0

                if span.getArea() <  self.config.minVariancePixels:
                    continue

                varSrc = numpy.nanvar(image[box].getArray()[mask])
                meanSrc = numpy.nanmean(image[box].getArray()[mask])

                failed = False
                break

            if setFlags and failed:
                source.set('bfd.flags.variance', True)
                source.set('bfd.flags', True)
                return False

            source.set('noise_variance', float(varSrc))
            source.set('noise_mean', float(meanSrc))
        return True


    def writeOutputs(self, dataRef, outCat):
        dataRef.put(outCat, self.dataPrefix+"moment", flags=afwTable.SOURCE_IO_NO_HEAVY_FOOTPRINTS)

    def getPreviousTaskClass(self):
        return MeasureCoaddTask

    def _getConfigName(self):
        """Return the name of the config dataset
        """
        return "%s_measureCoadd_config" % (self.config.coaddName,)

    def _getMetadataName(self):
        """Return the name of the metadata dataset
        """
        return "%s_measureCoadd_metadata" % (self.config.coaddName,)





