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
import lsst.daf.base as dafBase
import lsst.pex.config
import lsst.pipe.base
import lsst.afw.image
import lsst.afw.geom
import lsst.afw.table
from lsst.meas.base.noiseReplacer import NoiseReplacer, DummyNoiseReplacer

import numpy
from .measureCoadd import MeasureCoaddTask, MeasureCoaddConfig
import lsst.desc.bfd as bfd
import random
random.seed(11155)
from lsst.daf.persistence import Butler

import scipy.spatial
__all__ = ("MeasurePriorConfig", "MeasurePriorTask")


def matchCatalogs(ra_ref, dec_ref, ra_p, dec_p, dmin):
    ra_ref = ra_ref*numpy.pi/180
    dec_ref = dec_ref*numpy.pi/180
    ra_p = ra_p*numpy.pi/180
    dec_p = dec_p*numpy.pi/180
    posRef = numpy.dstack([numpy.sin(dec_ref)*numpy.cos(ra_ref), numpy.sin(dec_ref)*numpy.sin(ra_ref), numpy.sin(dec_ref)])[0]
    posP = numpy.dstack([numpy.sin(dec_p  )*numpy.cos(ra_p  ), numpy.sin(dec_p  )*numpy.sin(ra_p), numpy.sin(dec_p  )])[0]
    mytree = scipy.spatial.cKDTree(posRef)
    dist, index = mytree.query(posP)

    # convert to arsec
    dist*=3600.*(180/numpy.pi)
    close = dist < dmin
    closeIndices=index[close]
    return close,closeIndices


def matchTreeCatalogs(tree_ref, ra_p, dec_p, dmin):
    ra_p = ra_p*numpy.pi/180
    dec_p = dec_p*numpy.pi/180
    posP = numpy.dstack([numpy.sin(dec_p  )*numpy.cos(ra_p  ), numpy.sin(dec_p  )*numpy.sin(ra_p), numpy.sin(dec_p  )])[0]
    dist, index = tree_ref.query(posP)

    # convert to arsec
    dist*=3600.*(180/numpy.pi)

    close = dist < dmin
    closeIndices=index[close]
    return close,closeIndices

def matchXYTreeCatalogs(tree_ref, x, y, dmin):
    posP = numpy.dstack([x,y])[0]
    dist, index = tree_ref.query(posP)

    close = dist < dmin
    closeIndices=index[close]
    return close,closeIndices



class MeasurePriorConfig(MeasureCoaddConfig):
    snMin = lsst.pex.config.Field(
        dtype=float,
        default=5,
        optional=True,
        doc="Minimun flux S/N"
    )
    snMax = lsst.pex.config.Field(
        dtype=float,
        default=25,
        optional=True,
        doc="Maximum flux S/N"
    )
    fluxMin = lsst.pex.config.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Minimun flux"
    )
    fluxMax = lsst.pex.config.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Maximum flux"
    )
    magMin = lsst.pex.config.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Minimun mag"
    )
    magMax = lsst.pex.config.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Maximum mag"
    )
    noiseFactor = lsst.pex.config.Field(
        dtype=float,
        default=1,
        optional=True,
        doc="Noise boost factor for kernel smoothing"
    )
    priorSigmaCutoff = lsst.pex.config.Field(
        dtype=float,
        default=5.5,
        optional=True,
        doc="Maximum sigma range when sampling for prior"
    )
    priorSigmaStep = lsst.pex.config.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Step size when sampling for prior"
    )
    priorSigmaBuffer = lsst.pex.config.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Buffer width of KdTreePrior (in sigma)"
    )
    nSample = lsst.pex.config.Field(
        dtype=int,
        default=30000,
        optional=True,
        doc="Number of templates sampled per target"
    )
    maxXY = lsst.pex.config.Field(
        dtype=float,
        default=4.,
        optional=True,
        doc="Maximum translational displacement in sigma of the nominal covariance matrix"

    )
    sigma = lsst.pex.config.Field(
        dtype=float,
        default=0.8,
        optional=True,
        doc="Sigma used in k-space weight function"

    )
    wIndex = lsst.pex.config.Field(
        dtype=int,
        default=3,
        optional=True,
        doc="index used in k-space weight function"

    )
    centroidName = lsst.pex.config.Field(
        dtype=str,
        default='centroid.sdss',
        optional=True,
        doc="name of centroid to use from the catalog"

    )
    selectionOnly = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="only do selection in the prior"

    )
    covFile = lsst.pex.config.Field(
        dtype=str,
        default='./test.fits',
        optional=True,
        doc="file that contains the covariance matrices"

    )
    maxVar = lsst.pex.config.Field(
        dtype=float,
        default=0.15,
        optional=True,
        doc="Minimum Variance that will be considered"

    )
    minVar = lsst.pex.config.Field(
        dtype=float,
        default=1e-4,
        optional=True,
        doc="Minimum Variance that will be considered"

    )
    sample = lsst.pex.config.Field(
        dtype=float,
        default=0.2,
        optional=True,
        doc="Only use this fraction of the galaxies"

    )
    invariantCovariance = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"

    )
    reCentroid = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="recentroid galaxis"

    )
    reCentroidPsf = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="recentroid galaxis"
    )
    shallowRerun = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="rerun for the prior"
    )
    useBlends = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=True,
        doc="use blended galaxies in the prior"
    )
    label = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='',
        doc="additional label to add to data id"
    )
    matchFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file to match objects to"
    )
    zFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for redshift information"
    )
    readZFile = lsst.pex.config.Field(
        dtype=bool,
        optional=True,
        default=False,
        doc="read truth z file from butler"
    )
    zPdfFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for redshift pdf information"
    )
    zMaxCut = lsst.pex.config.Field(
        dtype=float,
        optional=True,
        default=None,
        doc=" maximum redshift cut"
    )
    zMinCut = lsst.pex.config.Field(
        dtype=float,
        optional=True,
        default=None,
        doc=" maximum redshift cut"
    )
    zQualityCut = lsst.pex.config.Field(
        dtype=float,
        optional=True,
        default=None,
        doc="redshift quality cut"
    )
    zField = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='frankenz_photoz_mode',
        doc="file for redshift information"
    )
    zRa = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='ira',
        doc="ra column from redshift file"
    )
    zDec = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='idec',
        doc="ra column from redshift file"
    )
    zId = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='object_id',
        doc="ra column from redshift file"
    )
    zQualityField = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='frankenz_photoz_zrisk_med',
        doc="field for quality cut"
    )
    zPdfField = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='pdf',
        doc="file for redshift information"
    )
    zPdfGridField = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default='zgrid',
        doc="file for redshift information"
    )
    zBin = lsst.pex.config.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting redshift bin"
    )
    noiseBin = lsst.pex.config.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting noise bin"
    )
    colorBin = lsst.pex.config.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting color bin"
    )
    useLabels = lsst.pex.config.ListField(
        dtype=str,
        default=[],
        optional=True,
        doc="List of labels for which to build the prior. "
    )
    ccFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for color-color cuts"
    )
    noClobber = lsst.pex.config.Field(
        dtype=bool,
        optional=True,
        default=False,
        doc="check if already exists"
    )
    extremeChi2 = lsst.pex.config.Field(
        dtype=float,
        optional=True,
        default=100000,
        doc="ignore objects with chi^2> than this"
    )
    verbose = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="print more info"

    )
    colorBinFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for binned color-color cuts"
    )
    useXY = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="us xy position instead of ra/dec"
    )
    zMatchDist = lsst.pex.config.Field(
        dtype=float,
        optional=True,
        default=0.5,
        doc="redshift quality cut"
    )
    mcSamples = lsst.pex.config.Field(
        dtype=int,
        default=500,
        optional=True,
        doc="number of Monte Carlo samples"
    )
    maxRatio = lsst.pex.config.Field(
        dtype=float,
        default=-1.,
        optional=True,
        doc="apply maximum ratio for selection"
    )

#    measureMomentCov = lsst.pex.config.ConfigField(
#        dtype=lsst.desd.bfd.bfd.measureMomentCov
#        doc="Procedure to measure the covariance"
#
class MeasurePriorTask(MeasureCoaddTask):
    ConfigClass = MeasurePriorConfig

    def __init__(self,  schema=None, **kwargs):
        """
        """
        MeasureCoaddTask.__init__(self, **kwargs)

        # make subtask to create covariance matrix
        # self.prior can create now set covariance later

    def getCovariances(self):
        cat = lsst.afw.table.SourceCatalog.readFits(self.config.covFile)

        covList = []

        for rec in cat:
            cov = rec.get('isoCov')
            label = rec.get('label')
            minVariance = rec.get('min')
            maxVariance = rec.get('max')
            cov = numpy.array(cov.reshape(6,6),dtype=numpy.float32)

            covList.append((cov,label,minVariance,maxVariance))


        return covList

    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
        of the pixel variance and covariance of the matrix as determined by the Psf
        """
        shallowSources = None
        shallowMoments = None
        shallowWcs = None
        if self.config.shallowRerun:
             shallowButler = Butler(self.config.shallowRerun)
             shallowSources = shallowButler.get(self.dataPrefix + "meas", dataId=dataRef.dataId, immediate=True)
             shallowCalexp = shallowButler.get(self.dataPrefix + "calexp", dataId=dataRef.dataId, immediate=True)
             shallowWcs = shallowCalexp.getWcs()
             #shallowMoments = shallowButler.get(self.dataPrefix + "moment", dataId=dataRef.dataId, immediate=True)

        exposure = dataRef.get(self.config.coaddName + "Coadd_calexp", immediate=True)

        return lsst.pipe.base.Struct(
            sources=dataRef.get(self.dataPrefix + "meas", immediate=True),
            moments=dataRef.get(self.dataPrefix + "moment", immediate=True),
            exposure=exposure,
            shallowSources=shallowSources,
            shallowWcs=shallowWcs
            #shallowMoments=shallowMoments
            )

    def run(self, dataRef):
        """Main driver
        """

        inputs = self.readInputs(dataRef)
        self.log.info("Preparing catalog")
        if self.config.maxObjects is not None:
            first = self.config.firstObject
            last = min(self.config.firstObject + self.config.maxObjects,len(inputs.sources))
            inputs.sources = inputs.sources[first:last]
        outCat = self.prepCatalog(inputs)

        self.runMeasure(inputs, outCat, dataRef)


        self.log.info("Writing outputs")
        self.writeOutputs(dataRef, outCat)
        return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

    def runMeasure(self, inputs, outCat, dataRef):
        # This should tell you how to select objects that apply and also
        # what the covariance matrix is
        covList = self.getCovariances()

        #centroidKey = inputs.sources.schema.find(self.config.centroidName).key
        noiseKey = outCat.schema.find('noise_variance').key

        if inputs.shallowSources:
            shallowCentroidKey = inputs.shallowSources.schema.find(self.config.centroidName).key
            shallow_ra = []
            shallow_dec = []
            for src in inputs.shallowSources:
                if not self.selection(None, src):
                    continue
                centroid = src.get(shallowCentroidKey)
                sky = inputs.shallowWcs.pixelToSky(centroid)
                shallow_ra.append(sky.toIcrs().getRa().asDegrees())
                shallow_dec.append(sky.toIcrs().getDec().asDegrees())

            shallow_ra = numpy.array(shallow_ra)*numpy.pi/180
            shallow_dec = numpy.array(shallow_dec)*numpy.pi/180
            posShallow = numpy.dstack([numpy.sin(shallow_dec)*numpy.cos(shallow_ra), numpy.sin(shallow_dec)*numpy.sin(shallow_ra), numpy.sin(shallow_dec)])[0]
            shallowTree = scipy.spatial.cKDTree(posShallow)

        if self.config.zFile or self.config.readZFile:
            if self.config.readZFile is False:
                import pyfits
                zFile = pyfits.open(self.config.zFile)[1].data
            else:
                zFile = dataRef.get('deepCoadd_truth', immediate=True)
            zRedshift = zFile[self.config.zField]
            zId = zFile[self.config.zId]
            if self.config.zQualityCut is not None:
                zQuality = zFile[self.config.zQualityField]

            if self.config.useXY is False:
                zRa = zFile[self.config.zRa]*numpy.pi/180
                zDec = zFile[self.config.zDec]*numpy.pi/180

                posRef = numpy.dstack([numpy.sin(zDec)*numpy.cos(zRa), numpy.sin(zDec)*numpy.sin(zRa), numpy.sin(zDec)])[0]
                zTree = scipy.spatial.cKDTree(posRef)
            else:
                zX = zFile[self.config.zRa]
                zY = zFile[self.config.zDec]

                posRef = numpy.dstack([zX, zY])[0]
                zTree = scipy.spatial.cKDTree(posRef)


            if self.config.zPdfFile:
                pdfId = pyfits.open(self.config.zPdfFile)[1].data['id']
                pdfData = pyfits.open(self.config.zPdfFile)[1].data[self.config.zPdfField]
                zPdfGrid = pyfits.open(self.config.zPdfFile)[2].data[self.config.zPdfGridField]

                zPdfIndex = (zPdfGrid < self.config.zMaxCut) & (zPdfGrid > self.config.zMinCut)
                print ('sum pdf index',numpy.sum(zPdfIndex))

            if self.config.colorBinFile:

                ccFile = pyfits.open(self.config.colorBinFile)[1].data
                colorBins = ccFile['color_bin']
                assert len(colorBins==len(zId))

        if self.config.ccFile:
            ccFile = numpy.load(self.config.ccFile)
            ccRa = ccFile['ra']*numpy.pi/180
            ccDec = ccFile['dec']*numpy.pi/180

            posRef = numpy.dstack([numpy.sin(ccDec)*numpy.cos(ccRa), numpy.sin(ccDec)*numpy.sin(ccRa), numpy.sin(ccDec)])[0]
            ccTree = scipy.spatial.cKDTree(posRef)

        for cov,label,varMin,varMax in covList:

            if (label not in self.config.useLabels) and len(self.config.useLabels) > 0:
                if self.config.verbose: self.log.info("Label %s not in %s"%(label,self.config.useLabels))
                continue

            # Build full label
            full_label = label + self.config.label
            if self.config.selectionOnly:
                full_label += '_selection'

            if self.config.noiseBin is not None:
                full_label += '_n%d'%self.config.noiseBin

            if self.config.zBin is not None:
                full_label += '_z%d'% (self.config.zBin)

            if self.config.colorBin is not None:
                full_label += '_c%d'%self.config.colorBin

            # if it already exists and we care don't run anything
            if self.config.noClobber:
                dataRef.dataId['label']=full_label
                if dataRef.datasetExists('deepCoadd_momentPrior'):
                    self.log.info('Prior already exists %s. skipping' % dataRef.dataId)
                    return

            self.log.info('Processing label %s' % label)
            sigmaFlux = numpy.sqrt(cov[0,0])

            minFlux = self.config.snMin*sigmaFlux
            if self.config.fluxMin is not None:
                minFlux = self.config.fluxMin
            elif self.config.magMin is not None:
                minFlux = 10**(-0.4*(self.config.magMin-27))

            maxFlux = self.config.snMax*sigmaFlux
            if self.config.fluxMax is not None:
                maxFlux = self.config.fluxMax
            elif self.config.magMax is not None:
                maxFlux = 10**(-0.4*(self.config.magMax-27))

            momentPrior = bfd.MomentPrior(minFlux, maxFlux,
                                          cov, self.config.invariantCovariance,#True,
                                          self.config.noiseFactor,
                                          self.config.priorSigmaCutoff,
                                          self.config.priorSigmaStep,
                                          self.config.nSample,
                                          self.config.priorSigmaBuffer,
                                          self.config.selectionOnly)
            if self.config.maxRatio > 0:
                momentPrior.setMaxRatio(self.config.maxRatio, self.config.mcSamples)

            momentPrior.setVarianceLimits(varMin, varMax)
            #noiseout = self.config.doReplaceWithNoise
            #if noiseout:
            #    self.replaceWithNoise.begin(inputs.exposure, inputs.sources)

            footprints = {measRecord.getId(): (measRecord.getParent(), measRecord.getFootprint())
                          for measRecord in inputs.sources}
            noiseImage=None
            exposureId = None
            beginOrder = None
            endOrder = None
            noiseReplacer = NoiseReplacer(self.config.noiseReplacer, inputs.exposure, footprints,
                                          noiseImage=noiseImage, log=self.log, exposureId=exposureId)

            ngood = 0
            iteration = 0
            for i, (source, ref, moment) in enumerate(zip(outCat, inputs.sources, inputs.moments)):

                if self.config.sample > 0:
                    if random.uniform(0,1) > self.config.sample:
                        continue

                noiseReplacer.insertSource(ref.getId())
                #with self.noiseContext(inputs, i, ref, noiseout) as inputs:
                if True:
                    if not self.selection(source, ref):
                        if self.config.verbose: self.log.info("Does not pass selection criteria")
                        continue

                    if moment.get('bfd.flags'):
                        continue
                    f = moment.get('bfd.moments')[0]
                    dflux = 0
                    dflux = max(dflux, minFlux - f)
                    dflux = max(dflux, f - maxFlux)
                    chisq = dflux*dflux / cov[0,0]

                    if chisq > self.config.extremeChi2:
                        if self.config.verbose: self.log.info("Chi2 too large")
                        continue

                    # This should happen in the prep catalog stake because we don't want to repeat
                    # I can fix this later...
                    if not self.preMeasure(source, ref, inputs.exposure):
                        continue

                    if (moment.get('bfd.momentsCov')[0] > self.config.maxVar  or
                        moment.get('bfd.momentsCov')[0] < self.config.minVar) :
                        if self.config.verbose: self.log.info("Does not pass variance cuts %f %0.2f/%0.2f"%(moment.get('bfd.momentsCov')[0],self.config.minVar,self.config.maxVar))
                        continue

                    sn = moment.get('bfd.moments')[0]/numpy.sqrt(moment.get('bfd.momentsCov')[0])
                    if self.config.verbose: self.log.info("Processing galaxy %d with S/N=%0.2f, centroid %0.2f,%0.2f"% (i,sn,moment.get('bfd.center_x'),moment.get('bfd.center_y')))

                    if inputs.shallowSources:
                        centroid = ref.getCentroid()#ref.get(centroidKey)
                        sky = inputs.exposure.getWcs().pixelToSky(centroid)
                        result = matchTreeCatalogs(shallowTree, numpy.array([sky.toIcrs().getRa().asDegrees()]),
                                                   numpy.array([sky.toIcrs().getDec().asDegrees()]), self.config.zMatchDist)
                        if numpy.sum(result[0]) < 1:
                            if self.config.verbose:
                                self.log.info("Failed Shallow match")
                                continue

                    # Weight is 1, unless using pdf
                    weight = 1./self.config.sample

                    if self.config.zFile or self.config.readZFile:
                        centroid = ref.get(centroidKey)
                        sky = inputs.exposure.getWcs().pixelToSky(centroid)
                        if self.config.useXY is False:
                            ra = numpy.array([sky.toIcrs().getRa().asDegrees()])
                            dec =  numpy.array([sky.toIcrs().getDec().asDegrees()])
                            result = matchTreeCatalogs(zTree, ra, dec, self.config.zMatchDist)
                        else:
                            x = numpy.array([centroid.getX()])
                            y = numpy.array([centroid.getY()])
                            result = matchXYTreeCatalogs(zTree, x, y, self.config.zMatchDist)

                        if numpy.sum(result[0]) < 1:
                            if self.config.verbose:
                                self.log.info("Failed to match redshift")
                            continue

                        if self.config.zQualityCut is not None:
                            if zQuality[result[1][0]] > self.config.zQualityCut:
                                if self.config.verbose:
                                    self.log.info("Failed redshift quality")
                                continue

                        if self.config.zMaxCut is not None and self.config.zMinCut is not None and self.config.zPdfFile is None:
                            redshift = zRedshift[result[1][0]]
                            if self.config.verbose: self.log.info('Redshift %f'%redshift)
                            if (redshift > self.config.zMaxCut) or (redshift < self.config.zMinCut):
                                if self.config.verbose:
                                    self.log.info("Redshift not range")
                                continue

                        if self.config.zPdfFile is not None:
                            idList = numpy.where(pdfId==zId[result[1][0]])[0]
                            if len(idList) == 0:
                                if self.config.verbose: self.log.info("Can't match to pdf")
                                continue

                            weight *= float(numpy.sum(pdfData[idList[0]][zPdfIndex]))
                            if self.config.verbose:
                                self.log.info("Using pdf weight %f"%weight)

                        if self.config.colorBinFile is not None:
                            colorBin = colorBins[result[1][0]]
                            if colorBin != self.config.colorBin:
                                continue
                            else:
                                if self.config.verbose:
                                    self.log.info("Selected color bin %d"%colorBin)

                    if self.config.ccFile:
                        centroid = ref.get(centroidKey)
                        sky = inputs.exposure.getWcs().pixelToSky(centroid)
                        result = matchTreeCatalogs(ccTree, numpy.array([sky.toIcrs().getRa().asDegrees()]),
                                                   numpy.array([sky.toIcrs().getDec().asDegrees()]), self.config.zMatchDist)

                        if numpy.sum(result[0]) < 1:
                                if self.config.verbose:
                                    self.log.info('Failed to match color-color')
                                    continue


                    if weight < 1e-6:
                        if self.config.verbose: self.log.info('Weight too low')
                        continue

                    bfd_control = bfd.BfdKMomentControl()
                    bfd_control.sigma = self.config.sigma
                    bfd_control.wIndex = self.config.wIndex
                    bfd_control.maxCentroid = self.config.maxXY
                    bfd_control.ignorePsf = False
                    bfd_control.shift = True
                    bfd_control.reCentroid = self.config.reCentroid
                    bfd_control.reCentroidPsf = self.config.reCentroidPsf

                    try:
                        if self.config.verbose:
                             self.log.info("Passed cuts")
                        pos = lsst.afw.geom.Point2D(moment.get('bfd.center_x'),moment.get('bfd.center_y'))
                        priorGalaxy = bfd.PriorGalaxy(bfd_control)
                        passed = priorGalaxy.addImage(ref, inputs.exposure, pos, source.get(noiseKey), True)
                        if passed is False:
                            if self.config.verbose: self.log.info('Failed add image')
                            continue

                        flip = True
                        momentPrior.addPriorGalaxy(priorGalaxy, self.config.maxXY, weight,
                                                   flip, source.getId())
                        ngood += 1
                    except MemoryError:
                        raise  # always let MemoryError propagate up, as continuing just causes
                               # more problems
                    except Exception as err:
                        self.log.warn("Error measuring source %s : %s"
                                      % (source.getId(), err))
                    iteration +=1
                noiseReplacer.removeSource(ref.getId())


            #self.replaceWithNoise.end(inputs.exposure, inputs.sources)
            noiseReplacer.end()


            selectionPqr = momentPrior.selectionProbability(cov)
            deselect = selectionPqr.copy()

            deselect[0] = 1 - selectionPqr[0]
            for i in range(1,6):
                deselect[i] *= -1.
            self.log.info('Used %d galaxies in prior'%ngood)
            catalog = momentPrior.getCatalog()
            metadata = dafBase.PropertyList()
            metadata.set('cov',numpy.array(cov.flatten(),dtype=float))
            metadata.set('selectionPqr',selectionPqr.astype(numpy.float))
            metadata.set('deselectPqr',deselect.astype(numpy.float))
            metadata.set('fluxMin',minFlux)
            metadata.set('fluxMax',maxFlux)
            metadata.set('varMin',varMin)
            metadata.set('varMax',varMax)
            if self.config.shallowRerun is not None:
                metadata.set('shallowRerun',self.config.shallowRerun)
            if self.config.zFile is not None:
                metadata.set('zFile',self.config.zFile)
                metadata.set('zField',self.config.zField)
            if self.config.zMaxCut is not None:
                metadata.set('zMaxCut',self.config.zMaxCut)
            if self.config.zMinCut is not None:
                metadata.set('zMinCut',self.config.zMinCut)

            if self.config.zPdfFile is not None:
                metadata.set('zPdfFile',self.config.zPdfFile)
            if self.config.colorBinFile is not None:
                metadata.set('colorBinFile',self.config.colorBinFile)

            if self.config.zBin is not None:
                metadata.set('zBin',self.config.zBin)
            if self.config.noiseBin is not None:
                metadata.set('noiseBin',self.config.noiseBin)
            if self.config.colorBin is not None:
                metadata.set('colorBin',self.config.colorBin)

            if self.config.maxRatio > 0:
                metadata.set('maxRatio',self.config.maxRatio)
                metadata.set('mcSamples',self.config.mcSamples)

            metadata.set('noiseFactor',self.config.noiseFactor)
            metadata.set('priorSigmaCutoff',self.config.priorSigmaCutoff)
            metadata.set('priorSigmaStep',self.config.priorSigmaStep)
            metadata.set('priorSigmaBuffer',self.config.priorSigmaBuffer)
            metadata.set('nsample',self.config.nSample)
            metadata.set('selectionOnly',self.config.selectionOnly)
            metadata.set('invariantCovariance',self.config.invariantCovariance)
            metadata.set('maxXY', self.config.maxXY)
            metadata.set('sigma', self.config.sigma)
            metadata.set('wIndex', self.config.wIndex)
            metadata.set('centroidName', self.config.centroidName)
            metadata.set('covFile', self.config.covFile)
            metadata.set('totalWeight', momentPrior.getTotalWeight())

            catalog.getTable().setMetadata(metadata)

            self.log.info('Created %d templates, total weight %f' % (len(catalog), momentPrior.getTotalWeight()))
            self.log.info('Selection %s' % selectionPqr)
            dataRef.dataId['label'] = full_label
#            dataRef.dataId['label'] = label + self.config.label
#
#            if self.config.selectionOnly:
#                dataRef.dataId['label']+='_selection'
#
#            if self.config.noiseBin is not None:
#                dataRef.dataId['label']+='_n%d'%self.config.noiseBin
#
#            if self.config.zPdfFile is not None:
#                dataRef.dataId['label']+='_z%d'% (self.config.zBin)
            dataRef.put(catalog, self.dataPrefix+"momentPrior")

    def prepCatalog(self, inputs):
        """Prepare the prior and return the output catalog
        """
        outCat =  lsst.afw.table.SourceCatalog(self.schema)
        srcCat = inputs.sources

        for srcRecord in srcCat:
            outRecord = outCat.addNew()
            #outRecord.setId(srcCat.get('id'))

        if self.config.calculateVariance:
            x0 = inputs.exposure.getXY0().getX()
            y0 = inputs.exposure.getXY0().getY()
            self.xy0 = lsst.afw.geom.Extent2I(-x0,-y0)

        return outCat

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        return

    def selection(self, source, ref):
        # Don't process blended parent objects
        if ref.getParent()==0 and ref.get('deblend_nChild')>0:
            return False
        # For now don't process objects in a blend
        if ref.getParent()!=0 and not self.config.useBlends:
            return False
        if ref.getFootprint().getArea() > self.config.maxArea:
            return False
        if ref.getFootprint().getArea() == 0:
            return False

        return True
