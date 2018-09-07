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
import lsst.meas.extensions.bfd
import numpy
from .measureCoadd import MeasureCoaddTask, MeasureCoaddConfig
from lsst.daf.persistence import Butler
import lsst.meas.extensions.bfd as bfd
import glob
import pyfits
import scipy.spatial


__all__ = ("MeasurePqrConfig", "MeasurePqrTask")

def matchTreeCatalogs(tree_ref, ra_p, dec_p, dmin):
    ra_p = ra_p*numpy.pi/180
    dec_p = dec_p*numpy.pi/180
    posP = numpy.dstack([numpy.sin(dec_p  )*numpy.cos(ra_p  ), numpy.sin(dec_p  )*numpy.sin(ra_p), numpy.sin(dec_p  )])[0]
    dist, index = tree_ref.query(posP)

    # convert to arsec
    dist*=3600.*(180/numpy.pi)

    close = dist < dmin
    closeIndices=index[close]
    return close,closeIndices,dist[close]
def matchXYTreeCatalogs(tree_ref, x, y, dmin):
    posP = numpy.dstack([x,y])[0]
    dist, index = tree_ref.query(posP)

    close = dist < dmin
    closeIndices=index[close]
    return close,closeIndices,dist[close]

class MeasurePqrZConfig(MeasureCoaddConfig):

    invariantCovariance = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"

    )
    numProc = lsst.pex.config.Field(
        dtype=int,
        default=1,
        optional=True,
        doc="number of processors to use"
    )
    chunk = lsst.pex.config.Field(
        dtype=int,
        default=100,
        optional=True,
        doc="number of chunks"
    )
    sampleSeed = lsst.pex.config.Field(
        dtype=int,
        default=0,
        optional=True,
        doc="number of processors to use"
    )
    sampleFraction = lsst.pex.config.Field(
        dtype=float,
        default=0.25,
        optional=True,
        doc="number of processors to use"
    )
    priorRerun = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        doc="rerun for the prior"
    )
    priorTracts = lsst.pex.config.ListField(
        dtype=int,
        default=[9813],
        optional=False,
        doc="tract for the prior"
    )
    priorLabel = lsst.pex.config.Field(
        dtype=str,
        default='b0',
        optional=False,
        doc="label for the prior"
    )
    priorFilter = lsst.pex.config.Field(
        dtype=str,
        default='HSC-I',
        optional=False,
        doc="filter for the prior"
    )
    maxPriorFiles = lsst.pex.config.Field(
        dtype=int,
        default=-1,
        optional=False,
        doc="filter for the prior"
    )
    priorPatches = lsst.pex.config.ListField(
        dtype=str,
        default=[],
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    checkExists = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    zFile = lsst.pex.config.Field(
        dtype=str,
        optional=True,
        default=None,
        doc="file for redshift information"
    )
    zField = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='mode',
        doc="file for redshift information"
    )
    noClobber = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="check if already exists"
    )
    ignoreZ = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="ignore redshift"
    )
    zBins = lsst.pex.config.ListField(
        dtype=float,
        default= [0, 0.55, 0.8, 1.1, 1.6, 3],
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    useAllZ = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=True,
        doc="attempt to use all available photo-z estimators to select appropriate prior"
    )
    zType = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='frankenz',
        doc="which photo-z estimator if not using all"
    )
    randomizeZ = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="Randomize redshift assignement"
    )
    readZFile = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="read truth z file from butler"
    )
    zField = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='frankenz_photoz_mode',
        doc="file for redshift information"
    )
    zRa = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='ira',
        doc="ra column from redshift file"
    )
    zDec = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='idec',
        doc="ra column from redshift file"
    )
    zId = lsst.pex.config.Field(
        dtype=str,
        optional=False,
        default='object_id',
        doc="ra column from redshift file"
    )
    useXY = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="us xy position instead of ra/dec"
    )
    zMatchDist = lsst.pex.config.Field(
        dtype=float,
        optional=False,
        default=0.5,
        doc="redshift quality cut"
    )
    writeEmpty = lsst.pex.config.Field(
        dtype=bool,
        default=True,
        optional=True,
        doc="write empty catalogs"
    )
    addTruth = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="add true shear values"
    )

    writeNonSelect = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="write out galaxies that were not selected by the flux cut with 1-P(s)"
    )



class MeasurePqrZTask(MeasureCoaddTask):
    ConfigClass = MeasurePqrZConfig

    def __init__(self,  schema=None, **kwargs):
        """
        """
        MeasureCoaddTask.__init__(self, **kwargs)

        self.schema =  lsst.afw.table.SourceTable.makeMinimalSchema()
        # Should move these into C++?
        self.noiseBinKey = self.schema.addField("bfd.noise.bin", doc="Noise bin", type=int)
        self.fluxBinKey = self.schema.addField("bfd.flux.bin", doc="Flux bin", type=int)
        self.zBinKey = self.schema.addField("bfd.redshift.bin", doc="Redshift", type=int)
        self.disZKey = self.schema.addField("bfd.redshift.distance", doc="distance to redshift", type=float)
        self.zKey = self.schema.addField("bfd.redshift", doc="redshift used in pqr", type=float)
        self.pqrKey = bfd.BfdPqrKey.addFields(self.schema, "bfd")
        self.flagKey = self.schema.find('bfd.flags').key
        if self.config.addTruth:
            self.g1Key = self.schema.addField("truth.g1", doc="true g1", type=float)
            self.g2Key = self.schema.addField("truth.g2", doc="true g2", type=float)
            self.xKey = self.schema.addField("x", doc="bfd x", type=float)
            self.yKey = self.schema.addField("y", doc="bfd y", type=float)
        if self.config.writeNonSelect:
            self.flagNSKey = self.schema.addField("bfd.flag.ns", doc="Failed flux selection", type="Flag")


    def readInputs(self, dataRef):
        """Return a lsst.pipe.base.Struct containing the Exposure to fit, catalog, measurement
        of the pixel variance and covariance of the matrix as determined by the Psf
        """

        return lsst.pipe.base.Struct(
            sources=dataRef.get(self.dataPrefix + "moment", immediate=True,
                                    flags=lsst.afw.table.SOURCE_IO_NO_FOOTPRINTS),
        )

    def prepCatalog(self, inputs):
        """Prepare the prior and return the output catalog
        """
        outCat =  lsst.afw.table.SourceCatalog(self.schema)
        srcCat = inputs.sources

        for srcRecord in srcCat:
            outRecord = outCat.addNew()
            #outRecord.setId(srcCat.get('id'))

        return outCat

    def prep(self):
        self.prior = bfd.MomentPrior()
        priorFiles=[]
        priorButler = Butler(self.config.priorRerun)
        prior_skyMap = priorButler.get('deepCoadd_skyMap')

        for tract in self.config.priorTracts:
            for patchInfo in prior_skyMap[tract]:
                patch = '%d,%d'%patchInfo.getIndex()

                if self.config.priorPatches:
                    if patch not in self.config.priorPatches:
                        continue

                if priorButler.datasetExists('deepCoadd_momentPrior', tract=tract, patch=patch,
                                             filter=self.config.priorFilter, label=self.config.priorLabel):
                    priorFiles.append(priorButler.get('deepCoadd_momentPrior_filename',
                                                       tract=tract, patch=patch,
                                                       filter=self.config.priorFilter,
                                                       label=self.config.priorLabel)[0])
        if len(priorFiles)==0:
            raise Exception('No Prior files found')
        max_file = len(priorFiles)
        if self.config.maxPriorFiles > 0:
            max_file = self.config.maxPriorFiles

        first = True
        self.zBin = None
        for file in priorFiles[:max_file]:
            if file.find('_parent') > 0:
                self.log.info("Skipping %s, from parent" % file)
                continue
            self.log.info("Adding prior %s" % file)
            try:
                cat = lsst.afw.table.BaseCatalog.readFits(file)

                self.prior.addCatalog(cat, self.config.invariantCovariance,
                                      self.config.sampleFraction, self.config.sampleSeed)
                # Should be same for all prior catalogs
                if first:
                    self.cov = numpy.array(cat.getTable().getMetadata().getArrayDouble('COV')).reshape(6,6)
                    self.zBin = cat.getTable().getMetadata().getInt('ZBIN')
                    self.fluxBin = cat.getTable().getMetadata().getInt('NOISEBIN')
                    first=False
                    self.log.info('Processing zbin: %d flux: %d'%(self.zBin,self.fluxBin))
            except  Exception as e:
                print 'Failed to read',e
                continue

        self.log.info("Building tree from prior")
        self.prior.prepare()
        self.fluxMin = self.prior.getFluxMin()
        self.fluxMax = self.prior.getFluxMax()
        self.varMin = self.prior.getVarMin()
        self.varMax = self.prior.getVarMax()

        if self.config.writeNonSelect:
            self.selectionPqr = self.prior.selectionProbability(self.cov.astype(numpy.float32))
            deselect = self.selectionPqr.copy()
            deselect[0] = 1 - self.selectionPqr[0]
            for i in range(1,6):
                deselect[i] *= -1.
            self.noSelectPqr = deselect
            self.log.info('Probability of not being selected: %s'%self.noSelectPqr)

        if self.config.readZFile is False:
            self.log.info("Reading redshift file")
            zFile = pyfits.open(self.config.zFile)[1].data

            self.zRedshift = numpy.zeros(len(zFile))
            self.zId = zFile[self.config.zId]
            if self.config.useXY is False:
                zRa = zFile[self.config.zRa]*numpy.pi/180
                zDec = zFile[self.config.zDec]*numpy.pi/180


                if self.config.useAllZ:
                    mask = zFile['frankenz_photoz_%s_isnull'%self.config.zField]==False
                    self.zRedshift[mask] = zFile['frankenz_photoz_%s'%self.config.zField][mask]
                    mask = (zFile['mizuki_photoz_%s_isnull'%self.config.zField]==False) & (self.zRedshift==0.0)
                    self.zRedshift[mask] = zFile['mizuki_photoz_%s'%self.config.zField][mask]
                    mask = (zFile['nnpz_photoz_%s_isnull'%self.config.zField]==False) & (self.zRedshift==0.0)
                    self.zRedshift[mask] = zFile['nnpz_photoz_%s'%self.config.zField][mask]
                    mask = (zFile['mlz_photoz_%s_isnull'%self.config.zField]==False) & (self.zRedshift==0.0)
                    self.zRedshift[mask] = zFile['mlz_photoz_%s'%self.config.zField][mask]
                else:
                    mask = zFile['%s_photoz_%s_isnull'%(self.config.zType,self.config.zField)]==False
                    self.zRedshift[mask] = zFile['%s_photoz_%s'%(self.config.zType,self.config.zField)][mask]



                    posRef = numpy.dstack([numpy.sin(zDec)*numpy.cos(zRa), numpy.sin(zDec)*numpy.sin(zRa), numpy.sin(zDec)])[0]
                self.log.info('Building redshift treee')
                self.zTree = scipy.spatial.cKDTree(posRef)
            else:
                zX = zFile[self.config.zRa]
                zY = zFile[self.config.zDec]

                posRef = numpy.dstack([zX, zY])[0]
                self.zTree = scipy.spatial.cKDTree(posRef)
                self.zRedshift = zFile[self.config.zField]

    def run(self, dataRef):
        """Main driver
        """
        self.log.info("Processing %s"% str(dataRef.dataId))

        if self.config.checkExists:
            dataRef.dataId['label'] = self.config.priorLabel
            if dataRef.datasetExists(self.dataPrefix+"pqr"):
                filename = dataRef.get(self.dataPrefix+"pqr_filename")[0]
                if filename.find('_parent') < 0 :
                    self.log.info("Skipping %s, file %s exists" % (str(dataRef.dataId), filename))
                    return

        if self.config.noClobber:
            dataRef.dataId['label'] = self.config.priorLabel
            if dataRef.datasetExists(self.dataPrefix+"pqr"):
                    self.log.info('Pqr already exists %s. skipping' % dataRef.dataId)
                    return

        if self.config.readZFile is True:
            self.log.info("Reading redshift file")
            zFile = dataRef.get('deepCoadd_truth', immediate=True)

            self.zRedshift = numpy.zeros(len(zFile))
            self.zId = zFile[self.config.zId]
            if self.config.useXY is False:
                zRa = zFile[self.config.zRa]*numpy.pi/180
                zDec = zFile[self.config.zDec]*numpy.pi/180


                if self.config.useAllZ:
                    mask = zFile['frankenz_photoz_%s_isnull'%self.config.zField]==False
                    self.zRedshift[mask] = zFile['frankenz_photoz_%s'%self.config.zField][mask]
                    mask = (zFile['mizuki_photoz_%s_isnull'%self.config.zField]==False) & (self.zRedshift==0.0)
                    self.zRedshift[mask] = zFile['mizuki_photoz_%s'%self.config.zField][mask]
                    mask = (zFile['nnpz_photoz_%s_isnull'%self.config.zField]==False) & (self.zRedshift==0.0)
                    self.zRedshift[mask] = zFile['nnpz_photoz_%s'%self.config.zField][mask]
                    mask = (zFile['mlz_photoz_%s_isnull'%self.config.zField]==False) & (self.zRedshift==0.0)
                    self.zRedshift[mask] = zFile['mlz_photoz_%s'%self.config.zField][mask]
                else:
                    mask = zFile['%s_photoz_%s_isnull'%(self.config.zType,self.config.zField)]==False
                    self.zRedshift[mask] = zFile['%s_photoz_%s'%(self.config.zType,self.config.zField)][mask]



                    posRef = numpy.dstack([numpy.sin(zDec)*numpy.cos(zRa), numpy.sin(zDec)*numpy.sin(zRa), numpy.sin(zDec)])[0]
                self.log.info('Building redshift treee')
                self.zTree = scipy.spatial.cKDTree(posRef)
            else:
                zX = zFile[self.config.zRa]
                zY = zFile[self.config.zDec]

                posRef = numpy.dstack([zX, zY])[0]
                self.zTree = scipy.spatial.cKDTree(posRef)
                self.zRedshift = zFile[self.config.zField]

            if self.config.addTruth:
                self.g1 = zFile['g1']
                self.g2 = zFile['g2']
                self.x = zFile['x']
                self.y = zFile['y']

        inputs = self.readInputs(dataRef)
        if self.config.maxObjects is not None:
            first = self.config.firstObject
            last = min(self.config.firstObject + self.config.maxObjects,len(inputs.sources))
            inputs.sources = inputs.sources[first:last]

        outCat = self.runMeasure(inputs, dataRef)
        if len(outCat) == 0:
            self.log.info("No objects processed")
            if self.config.writeEmpty is False:
                return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

        self.writeOutputs(dataRef, outCat)


        return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

    def runMeasureMulti(self, args):
        self.runMeasure(self, *args)

    def runMeasure(self, inputs, dataRef):

        sources = inputs.sources
        flags = sources.get('bfd.flags')
        flux = sources.get('bfd.moments')[:,0]
        noise = sources.get('bfd.momentsCov')[:,0]
        pqrKey = self.schema.find('bfd.pqr').key

        # Preslection cuts
        pre_sel = flags == False

        # redshift selection
        minimumDist = 0.5
        #zBins = [0, 0.55, 0.8, 1.1, 1.6, 4]

        #zBins = [0, 0.55, 0.8, 1.1, 1.6, 4]
        #zBins = [0,0.6,0.9, 1.25, 3]
        #zBins = [0,3]
        #self.zBin=1
        if self.config.useXY is False:
            ra = numpy.rad2deg(sources['coord.ra'])
            dec = numpy.rad2deg(sources['coord.dec'])
            result = matchTreeCatalogs(self.zTree, ra, dec, minimumDist)
        else:
            ra = numpy.array(sources['bfd.center.x'])
            dec = numpy.array(sources['bfd.center.y'])
            result = matchXYTreeCatalogs(self.zTree, ra, dec, self.config.zMatchDist)

        redshift = numpy.zeros(len(ra))
        distance = numpy.zeros(len(ra))

        redshift[result[0]==False] = -1.
        distance[result[0]==False] = -1.
        redshift[result[0]] = self.zRedshift[result[1]]
        distance[result[0]] = result[2]
        if self.config.addTruth:
            g1 = numpy.zeros(len(ra))
            g2 = numpy.zeros(len(ra))
            g1[result[0]==False] = -1.
            g1[result[0]] = self.g1[result[1]]
            g2[result[0]==False] = -1.
            g2[result[0]] = self.g2[result[1]]

            x = numpy.zeros(len(ra))
            y = numpy.zeros(len(ra))
            x[result[0]==False] = -1.
            x[result[0]] = self.x[result[1]]
            y[result[0]==False] = -1.
            y[result[0]] = self.y[result[1]]


        self.log.info("distance length %d,%d,%d"%(len(distance),len(redshift),len(sources)))

        self.log.info("Using redshift bins %s:"%self.config.zBins)
        redshiftBin = numpy.digitize(redshift, self.config.zBins)
        redshiftBin[redshiftBin==len(self.config.zBins)] = len(self.config.zBins)-1

        if self.config.ignoreZ:
            redshiftBin[:] = self.zBin

        if self.config.randomizeZ:
            numpy.random.shuffle(redshiftBin)

        noiseBins = numpy.arange(0.05,1.25,0.05)
        noiseBin = numpy.digitize(noise,noiseBins) - 1
        # Flux selection
        sel = numpy.logical_and.reduce((pre_sel,
                                        noise > self.varMin,
                                        noise < self.varMax,
                                        flux > self.fluxMin,
                                        flux < self.fluxMax,
                                        redshiftBin == self.zBin
        ))
        self.log.info("Surviving cuts:")
        self.log.info("   presel: %d" % numpy.sum(pre_sel))
        self.log.info("   noise: %d" % numpy.sum((noise > self.varMin)&(noise < self.varMax)))
        self.log.info("   flux: %d" % numpy.sum((flux > self.fluxMin)&(flux < self.fluxMax)))
        self.log.info("   redshift: %d" % numpy.sum(redshiftBin == self.zBin))
        self.log.info("  total:%d" % numpy.sum(sel))

        outCat =  lsst.afw.table.SourceCatalog(self.schema)

        for ii,(srcRecord,dis,zz,noi) in enumerate(zip(sources[sel],distance[sel],redshift[sel],noiseBin[sel])):
            print noi,self.varMin,self.varMax
            outRecord = outCat.addNew()
            outRecord.setId(srcRecord.getId())
            outRecord.setRa(srcRecord.getRa())
            outRecord.setDec(srcRecord.getDec())

            if self.config.addTruth:
                outRecord.set(self.g1Key, g1[sel][ii])
                outRecord.set(self.g2Key, g2[sel][ii])
                outRecord.set(self.xKey, x[sel][ii])
                outRecord.set(self.yKey, y[sel][ii])

        self.prior.getPqrCat(sources[sel], outCat, self.config.numProc, self.config.chunk)

        if self.config.writeNonSelect:

            not_sel = sel = (
                ((pre_sel)&(noise > self.varMin)&(noise < self.varMax)&(redshiftBin == self.zBin)) &
                ((flux < self.fluxMin)|(flux > self.fluxMax))
                )
            self.log.info("Adding non-selection values for %d objects" % numpy.sum(not_sel))
            for ii,(srcRecord,dis,zz,noi) in enumerate(zip(sources[not_sel],distance[not_sel],redshift[not_sel],noiseBin[not_sel])):
                outRecord = outCat.addNew()
                outRecord.setId(srcRecord.getId())
                outRecord.setRa(srcRecord.getRa())
                outRecord.setDec(srcRecord.getDec())
                outRecord.set(self.fluxBinKey, self.fluxBin)
                outRecord.set(self.zBinKey, self.zBin)
                outRecord.set(self.disZKey, dis)
                outRecord.set(self.zKey, zz)
                outRecord.set(self.noiseBinKey, noi)
                if self.config.addTruth:
                    outRecord.set(self.g1Key, g1[not_sel][ii])
                    outRecord.set(self.g2Key, g2[not_sel][ii])
                    outRecord.set(self.xKey, x[not_sel][ii])
                    outRecord.set(self.yKey, y[not_sel][ii])

                outRecord.set(pqrKey, numpy.array(self.noSelectPqr).astype(numpy.float32))
                outRecord.set(self.flagNSKey, True)



        return outCat

    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        dataRef.dataId['label'] = self.config.priorLabel
        dataRef.put(outCat, self.dataPrefix+"pqr")
        return

