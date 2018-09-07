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
import numpy
from .measureCoadd import MeasureCoaddTask, MeasureCoaddConfig
from lsst.daf.persistence import Butler
import lsst.desc.bfd as bfd
import glob

__all__ = ("MeasurePqrConfig", "MeasurePqrTask")

class MeasurePqrConfig(MeasureCoaddConfig):

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
        default=None,
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    checkExists = lsst.pex.config.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    noClobber = lsst.pex.config.Field(
        dtype=bool,
        optional=False,
        default=False,
        doc="check if already exists"
    )


class MeasurePqrTask(MeasureCoaddTask):
    ConfigClass = MeasurePqrConfig

    def __init__(self,  schema=None, **kwargs):
        """
        """
        MeasureCoaddTask.__init__(self, **kwargs)

        self.schema =  lsst.afw.table.SourceTable.makeMinimalSchema()
        # Should move these into C++?
        self.flagMomKey = self.schema.addField("bfd.flags.moment", doc="flag bad input on moments", type='Flag')
        self.notSelFluxKey = self.schema.addField("bfd.ns.flux", doc="not selected because of flux", type='Flag')
        self.notSelVarKey = self.schema.addField("bfd.ns.var", doc="not selected because of variance", type='Flag')
        self.pqrKey = bfd.BfdPqrKey.addFields(self.schema, "bfd")
        self.flagKey = self.schema.find('bfd.flags').key


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

        max_file = len(priorFiles)
        if self.config.maxPriorFiles > 0:
            max_file = self.config.maxPriorFiles

        first = True
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
                    first=False
            except  Exception as e:
                print ('Failed to read',e)
                continue

        self.prior.prepare()
        self.fluxMin = self.prior.getFluxMin()
        self.fluxMax = self.prior.getFluxMax()
        self.varMin = self.prior.getVarMin()
        self.varMax = self.prior.getVarMax()
        selectionPqr = self.prior.selectionProbability(self.cov.astype(numpy.float32))
        deselect = selectionPqr.copy()
        deselect[0] = 1 - selectionPqr[0]
        for i in range(1,6):
            deselect[i] *= -1.
        self.noSelectPqr = deselect

    def run(self, dataRef):
        """Main driver
        """
        self.log.info("Processing %s"% str(dataRef.dataId))

        if self.config.checkExists:
            dataRef.dataId['label'] = self.config.priorLabel
            if dataRef.datasetExists(self.dataPrefix+"pqr"):
                if self.config.noClobber:
                    self.log.info('Pqr already exists %s. skipping' % dataRef.dataId)
                    return
                filename = dataRef.get(self.dataPrefix+"pqr_filename")[0]
                if filename.find('_parent') < 0 :
                    self.log.info("Skipping %s, file %s exists" % (str(dataRef.dataId), filename))
                    return

        inputs = self.readInputs(dataRef)
        if self.config.maxObjects is not None:
            first = self.config.firstObject
            last = min(self.config.firstObject + self.config.maxObjects,len(inputs.sources))
            inputs.sources = inputs.sources[first:last]
        outCat = self.prepCatalog(inputs)
        print ('Number of sources',len(outCat))
        self.runMeasure(inputs.sources, outCat, dataRef)

        self.writeOutputs(dataRef, outCat)
        return lsst.pipe.base.Struct(outCat=outCat, inputs=inputs)

    def runMeasureMulti(self, args):
        self.runMeasure(self, *args)

    def runMeasure(self, sources, outCat, dataRef):

        flags = sources.get('bfd.flags')
        flux = sources.get('bfd.moments')[:,0]
        noise = sources.get('bfd.momentsCov')[:,0]
        pqrKey = self.schema.find('bfd.pqr').key

        # Preslection cuts
        fail_pre_sel = flags == True
        [rec.set(self.flagMomKey,True) for rec in outCat[fail_pre_sel]]
        [rec.set(self.flagKey,True) for rec in outCat[fail_pre_sel]]

        pre_sel = numpy.logical_not(fail_pre_sel)
        pre_frac =  1-1.*numpy.sum(fail_pre_sel)/len(sources)
        self.log.info('Fraction passing preselection %g (%d)' % (pre_frac, numpy.sum(numpy.logical_not(fail_pre_sel))))

        # Variance selection
        fail_var = numpy.logical_and.reduce((pre_sel,
                                             numpy.logical_or(noise < self.varMin, noise > self.varMax)
                                             ))
        if numpy.sum(pre_sel) >0:
            var_frac =  1-1.*numpy.sum(fail_var)/numpy.sum(pre_sel)
        else:
            var_frac = 0.
        self.log.info('Remaining fraction passing after variance cuts %g (%d)' % (var_frac, numpy.sum(numpy.logical_not(fail_var))))
        [rec.set(self.notSelVarKey, True) for rec in outCat[fail_var]]
        [rec.set(self.flagKey, True) for rec in outCat[fail_var]]

        # Flux selection
        flux_sel = numpy.logical_and.reduce((pre_sel,
                                             noise > self.varMin,
                                             noise < self.varMax,
                                             flux > self.fluxMin,
                                             flux < self.fluxMax
        ))
        print (self.varMin,self.varMax,self.fluxMin,self.fluxMax)
        not_flux_sel = numpy.logical_and.reduce((pre_sel,
                                                 noise > self.varMin,
                                                 noise < self.varMax,
                                                 numpy.logical_or(flux < self.fluxMin, flux > self.fluxMax)
        ))

        total = numpy.sum(flux_sel) + numpy.sum(not_flux_sel)

        if total == 0:
            sel_frac = 0
        else:
            sel_frac =  1.*numpy.sum(flux_sel)/(numpy.sum(flux_sel) + numpy.sum(not_flux_sel))

        [rec.set(pqrKey, numpy.array(self.noSelectPqr).astype(numpy.float32)) for rec in outCat[not_flux_sel]]

        # pqr will be set from no selection term
        [rec.set(self.notSelFluxKey, True) for rec in outCat[not_flux_sel]]
        [rec.set(self.flagKey, False) for rec in outCat[not_flux_sel]]
        self.log.info('Remaining fraction passing flux selection %g (%d) / %g (total)' % (sel_frac,numpy.sum(flux_sel),numpy.sum(flux_sel)/(1.*len(sources))))


        self.log.info("Surviving cuts:")
        self.log.info("   presel: %d" % numpy.sum(pre_sel))
        self.log.info("   noise: %d" % numpy.sum((noise > self.varMin)&(noise < self.varMax)))
        self.log.info("   flux: %d" % numpy.sum((flux > self.fluxMin)&(flux < self.fluxMax)))
        self.log.info("  total:%d" % numpy.sum(flux_sel))

        # Set flag key to false initially, a failure will be set in c++
        [rec.set(self.flagKey, False) for rec in outCat[flux_sel]]
        self.prior.getPqrCat(sources[flux_sel], outCat[flux_sel], self.config.numProc,  100)


            # Cut for failed moments or size too small
            #if (moment.FAILED or numpy.any(numpy.isnan(moment.moments)) or
            #    moment.moments[3]/moment.moments[0]>0.3):
            #    source.set(self.preFlagKey,True)
            #    continue
            #
            ## Cut on variance defined by nominal covariance matrix
            #if  ref.get('noise.variance')  < self.varMin or  ref.get('noise.variance')  > self.varMax:
            #    source.set(self.preFlagKey,True)
            #    continue
            #
            ## Cut out selection
            #if moment.moments[0] < self.fluxMin or moment.moments[0] > self.fluxMax:
            #    source.set(self.notSelKey,True)
            #    continue
            #print i
            #
            #try:
            #    pqr = self.prior.getPqr(moment)
            #except Exception as err:
            #    self.log.warn("Error measuring source %d" % i)
            #    source.set(self.flagKey,True)
            #    continue
            #print pqr.pqr
            #self.pqrKey.set(source,pqr)



    def writeOutputs(self, dataRef, outCat):
        """Write task outputs using the butler.
        """
        dataRef.dataId['label'] = self.config.priorLabel
        dataRef.put(outCat, self.dataPrefix+"pqr")
        return

    def selection(self, source, ref):
        # Don't process blended parent objects
        if ref.getParent()==0 and ref.get('deblend_nChild')>0:
            return False
        # For now don't process objects in a blend
        #if ref.getParent()!=0:
        #    return False
        if ref.getFootprint().getArea() > self.config.maxArea:
            return False
        if ref.getFootprint().getArea() == 0:
            return False
        #if ref.get('classification.extendedness') == 0:
        #    return False
        return True
