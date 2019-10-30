import time
import bfd
import numpy as np
from collections import defaultdict
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u
from astropy.table import Table
import lsst.daf.base as dafBase

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

from lsst.pipe.tasks.multiBand import getShortFilterName
from lsst.meas.base import NoiseReplacerConfig, NoiseReplacer
from lsst.afw.table import SourceCatalog, SchemaMapper
from lsst.pipe.base import Struct
from lsst.pex.config import Field, ListField, ConfigField, Config, ChoiceField

from .measureCoaddsTogether import ProcessCoaddsTogetherTask, ProcessCoaddsTogetherConfig
from .config import BFDConfig
from . import KSigmaWeightF, UniformDeviate
from .KColorGalaxy import KColorGalaxy


__all__ = ("MeasureCoaddsPriorConfig", "MeasureCoaddsPriorTask")


class MeasureCoaddsPriorConfig(ProcessCoaddsTogetherConfig):
    """
    basic config loads filters and misc stuff
    """
    filters = ListField(
        dtype=str,
        default=['r'],
        doc="List of expected bandpass filters."
    )
    start_index = Field(
        dtype=int,
        default=0,
        optional=True,
        doc='optional starting index for the processing',
    )
    num_to_process = Field(
        dtype=int,
        default=None,
        optional=True,
        doc='optional number to process',
    )
    use_mag = Field(dtype=bool, default=True, doc="include magnification")
    use_conc = Field(dtype=bool, default=True, doc="include concentration")
    weight_n = Field(dtype=int, default=4, doc="n for k-sigma weight function")
    weight_sigma = Field(dtype=float, default=0.5, doc="sigma for k-sigma weight function")
    max_shift = Field(dtype=float, default=4, doc="maximum centroid shift")
    coaddName = Field(dtype=str, default='deep', doc="name of coadd")
    grid_size = Field(dtype=int, default=48, doc="minimum size of pixel grid when computing moments")
    snMin = pexConfig.Field(
        dtype=float,
        default=5,
        optional=True,
        doc="Minimun flux S/N"
    )
    snMax = pexConfig.Field(
        dtype=float,
        default=25,
        optional=True,
        doc="Maximum flux S/N"
    )
    fluxMin = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Minimun flux"
    )
    fluxMax = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Maximum flux"
    )
    magMin = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Minimun mag"
    )
    magMax = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="Maximum mag"
    )
    noiseFactor = pexConfig.Field(
        dtype=float,
        default=1,
        optional=True,
        doc="Noise boost factor for kernel smoothing"
    )
    priorSigmaCutoff = pexConfig.Field(
        dtype=float,
        default=5.5,
        optional=True,
        doc="Maximum sigma range when sampling for prior"
    )
    priorSigmaStep = pexConfig.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Step size when sampling for prior"
    )
    priorSigmaBuffer = pexConfig.Field(
        dtype=float,
        default=1.,
        optional=True,
        doc="Buffer width of KdTreePrior (in sigma)"
    )
    nSample = pexConfig.Field(
        dtype=int,
        default=30000,
        optional=True,
        doc="Number of templates sampled per target"
    )
    maxXY = pexConfig.Field(
        dtype=float,
        default=4.,
        optional=True,
        doc="Maximum translational displacement in sigma of the nominal covariance matrix"
    )
    selectionOnly = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="only do selection in the prior"
    )
    covFile = pexConfig.Field(
        dtype=str,
        default='./run1.2_r_cov.fits',
        optional=True,
        doc="file that contains the covariance matrices"
    )
    maxVar = pexConfig.Field(
        dtype=float,
        default=30,
        optional=True,
        doc="Minimum Variance that will be considered"
    )
    minVar = pexConfig.Field(
        dtype=float,
        default=10,
        optional=True,
        doc="Minimum Variance that will be considered"
    )
    sample = pexConfig.Field(
        dtype=float,
        default=0.2,
        optional=True,
        doc="Only use this fraction of the galaxies"
    )
    invariantCovariance = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"
    )
    label = pexConfig.Field(
        dtype=str,
        optional=False,
        default='',
        doc="additional label to add to data id"
    )
    zFile = pexConfig.Field(
        dtype=str,
        optional=True,
        default='truth_cat.fits',
        doc="file for redshift information"
    )
    zDistCut = pexConfig.Field(
        dtype=float,
        optional=True,
        default=0.5,
        doc=" maximum distance for redshift match in arcsec"
    )
    zMaxCut = pexConfig.Field(
        dtype=float,
        optional=True,
        default=None,
        doc=" maximum redshift cut"
    )
    zMinCut = pexConfig.Field(
        dtype=float,
        optional=True,
        default=None,
        doc=" maximum redshift cut"
    )
    zField = pexConfig.Field(
        dtype=str,
        optional=True,
        default='redshift',
        doc="file for redshift information"
    )
    zRa = pexConfig.Field(
        dtype=str,
        optional=True,
        default='ra',
        doc="ra column from redshift file"
    )
    zDec = pexConfig.Field(
        dtype=str,
        optional=True,
        default='dec',
        doc="ra column from redshift file"
    )
    zId = pexConfig.Field(
        dtype=str,
        optional=True,
        default='id',
        doc="ra column from redshift file"
    )
    zBin = pexConfig.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting redshift bin"
    )
    noiseBin = pexConfig.Field(
        dtype=int,
        optional=True,
        default=None,
        doc="integer label represnting noise bin"
    )
    useLabels = pexConfig.ListField(
        dtype=str,
        default=[],
        optional=True,
        doc="List of labels for which to build the prior. "
    )
    zMatchDist = pexConfig.Field(
        dtype=float,
        optional=True,
        default=0.5,
        doc="redshift distance match"
    )
    mcSamples = pexConfig.Field(
        dtype=int,
        default=500,
        optional=True,
        doc="number of Monte Carlo samples"
    )
    maxRatio = pexConfig.Field(
        dtype=float,
        default=-1.,
        optional=True,
        doc="apply maximum ratio for selection"
    )
    randomSeed = pexConfig.Field(
        dtype=int,
        default=11155,
        optional=True,
        doc="apply maximum ratio for selection"
    )
    priorType = pexConfig.Field(
        doc="Name of prior type in policy",
        dtype=str,
        default='momentPrior',
    )

    def setDefaults(self):
        """
        prefix for the output file
        """
        self.output.name = "deepCoadd_prior"


class MeasureCoaddsPriorTask(ProcessCoaddsTogetherTask):
    """
    Base class for ngmix tasks
    """
    _DefaultName = "MeasureCoaddsPriorTask"
    ConfigClass = MeasureCoaddsPriorConfig

    def __init__(self, *, config=None, refSchema=None, butler=None, initInputs=None, **kwds):
        #import pdb;pdb.set_trace()
        ProcessCoaddsTogetherTask.__init__(self, config=config, refSchema=refSchema, butler=butler,
                                           initInputs=initInputs, **kwds)
        if refSchema is None:
            if butler is None:
                if initInputs is not None:
                    refSchema = initInputs.get("refSchema", None)
                if refSchema is None:
                    refSchema = SourceCatalog.Table.makeMinimalSchema()
            else:
                refSchema = butler.get(self.config.ref.name + "_schema").schema
        self.ncolors = len(self.config.filters) - 1
        self.bfd = BFDConfig(use_conc=self.config.use_conc, use_mag=self.config.use_mag,
                             ncolors=self.ncolors)
        self.n_even = self.bfd.BFDConfig.MSIZE
        self.n_odd = self.bfd.BFDConfig.XYSIZE
        self.weight = KSigmaWeightF(self.config.weight_sigma, self.config.weight_n)
        self.schema = self.defineSchema(refSchema)
        self.ud = UniformDeviate()

    def getCovariances(self):
        cat = afwTable.BaseCatalog.readFits(self.config.covFile)

        covList = []

        for rec in cat:
            cov_even = rec.get('isoCovEven')
            cov_odd = rec.get('isoCovOdd')
            label = rec.get('label')
            minVariance = rec.get('min')
            maxVariance = rec.get('max')
            full_cov_even = np.zeros((self.n_even, self.n_even), dtype=np.float32)
            full_cov_odd = np.zeros((self.n_odd, self.n_odd), dtype=np.float32)

            start = 0
            for i in range(self.n_even):
                full_cov_even[i][i:] = cov_even[start:start + self.n_even - i]
                start += self.n_even - i

            for i in range(self.n_even):
                for j in range(i):
                    full_cov_even[i, j] = full_cov_even[j, i]

            start = 0
            for i in range(self.n_odd):
                full_cov_odd[i][i:] = cov_odd[start:start + self.n_odd - i]
                start += self.n_odd - i

            for i in range(self.n_odd):
                for j in range(i):
                    full_cov_odd[i, j] = full_cov_odd[j, i]

            covList.append((full_cov_even, full_cov_odd, label, minVariance, maxVariance))

        return covList

    def defineSchema(self, refSchema):

        self.mapper = SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(SourceCatalog.Table.makeMinimalSchema(), True)
        schema = self.mapper.getOutputSchema()
        self.mKey = schema.addField("m", doc="template m", type="ArrayF",
                                    size=self.bfd.BFDConfig.MXYSIZE)
        self.dmKey = schema.addField("dm", doc="template m", type="ArrayF",
                                     size=self.bfd.BFDConfig.MSIZE*self.bfd.BFDConfig.DSIZE)
        self.dxyKey = schema.addField("dxy", doc="template m", type="ArrayF",
                                      size=self.bfd.BFDConfig.XYSIZE*self.bfd.BFDConfig.DSIZE)
        self.ndaKey = schema.addField("nda", doc="nda", type=np.float)
        self.idKey = schema.addField("bfd_id", doc="id", type=np.int64)
        if self.config.zFile:
            self.zKey = schema.addField("z", doc="redshift", type=np.float)
            # self.zIdKey = schema.addField("z_id", doc="redshift", type=np.int64)

        return schema

    def runDataRef(self, patchRefList):
        """Run this task via CmdLineTask and Gen2 Butler.
        Parameters
        ----------
        patchRefList : `list` of `lsst.daf.persistence.ButlerDataRef`
            A list of DataRefs for all filters in a single patch.
        """
        images = {}
        replacers = {}
        mergedDataId = {"tract": patchRefList[0].dataId["tract"],
                        "patch": patchRefList[0].dataId["patch"]}
        butler = patchRefList[0].butlerSubset.butler
        ref = butler.get("deepCoadd_ref", dataId=mergedDataId)
        imageId = butler.get("deepMergedCoaddId", dataId=mergedDataId)
        moments = butler.get('deepCoadd_moments', dataId=mergedDataId)

        for patchRef in patchRefList:
            filt = getShortFilterName(patchRef.dataId["filter"])
            images[filt] = patchRef.get(self.config.images.name)

            fpCat = patchRef.get(self.config.deblendCatalog.name)
            footprints = {rec.getId(): (rec.getParent(), rec.getFootprint()) for rec in fpCat}
            replacers[filt] = NoiseReplacer(self.config.deblendReplacer, exposure=images[filt],
                                            footprints=footprints, exposureId=imageId)
        results = self.run(butler, mergedDataId, images, ref, moments, imageId=imageId, replacers=replacers)

    def run(self, butler, dataId, images, ref, moments, imageId, replacers):
        """Process coadds from all bands for a single patch.
        This method should not add or modify self.
        So far all children are using this exact code so leaving
        it here for now. If we specialize a lot, might make a
        processor its own object
        Parameters
        ----------
        images : `dict` of `lsst.afw.image.ExposureF`
            Coadd images and associated metadata, keyed by filter name.
        ref : `lsst.afw.table.SourceCatalog`
            A catalog with one record for each object, containing "best"
            measurements across all bands.
        replacers : `dict` of `lsst.meas.base.NoiseReplacer`, optional
            A dictionary of `~lsst.meas.base.NoiseReplacer` objects that can
            be used to insert and remove deblended pixels for each object.
            When not `None`, all detected pixels in ``images`` will have
            *already* been replaced with noise, and this *must* be used
            to restore objects one at a time.
        imageId : `int`
            Unique ID for this unit of data.  Should be used (possibly
            indirectly) to seed random numbers.
        Returns
        -------
        results : `lsst.pipe.base.Struct`
            Struct with (at least) an `output` attribute that is a catalog
            to be written as ``self.config.output``.
        """

        if len(images) != len(self.config.filters) != len(moments):
            self.log.info('Number of filters does not match the list of images given.  Skipping')
            return None

        covList = self.getCovariances()

        if self.config.zFile:
            zFile = Table.read(self.config.zFile)

            z_redshift = zFile[self.config.zField]
            # z_id = zFile[self.config.zId]

            z_catalog = SkyCoord(ra=zFile[self.config.zRa]*u.deg,
                                 dec=zFile[self.config.zDec]*u.deg)
            self.z_list = {}
            # self.z_id_list = {}

            coord = SkyCoord(ra=ref['coord_ra']*u.rad, dec=ref['coord_dec']*u.rad)

            idx, d2d, d3d = match_coordinates_sky(coord, z_catalog)
            d2d = d2d.arcsec

        # we only need to compute the templates the first time through the loop
        first = True
        self.templates = []
        for cov_even, cov_odd, label, varMin, varMax in covList:
            tm0 = time.time()
            nproc = 0

            if (label not in self.config.useLabels) and len(self.config.useLabels) > 0:
                self.log.info("Label %s not in %s" % (label, self.config.useLabels))
                continue

            # Build full label
            full_label = label + self.config.label
            if self.config.selectionOnly:
                full_label += '_selection'

            if self.config.noiseBin is not None:
                full_label += '_n%d' % self.config.noiseBin

            if self.config.zBin is not None:
                full_label += '_z%d' % (self.config.zBin)

            self.log.info('Processing label %s' % label)
            sigmaFlux = np.sqrt(cov_even[0, 0])

            minFlux = self.config.snMin*sigmaFlux
            if self.config.fluxMin is not None:
                minFlux = self.config.fluxMin
            elif self.config.magMin is not None:

                # Current assumption is that coadd zeropoint is 27, true for HSC
                minFlux = 10**(-0.4*(self.config.magMin-27))

            maxFlux = self.config.snMax*sigmaFlux
            if self.config.fluxMax is not None:
                maxFlux = self.config.fluxMax
            elif self.config.magMax is not None:
                maxFlux = 10**(-0.4*(self.config.magMax-27))
            self.log.info('Min/Max %0.2f/%0.2f', minFlux, maxFlux)

            covMat = self.bfd.MomentCov(cov_even, cov_odd)
            momentPrior = self.bfd.KDTreePrior(minFlux, maxFlux, covMat, self.ud,
                                               self.config.nSample, self.config.selectionOnly,
                                               self.config.noiseFactor, self.config.priorSigmaStep,
                                               self.config.priorSigmaCutoff, self.config.priorSigmaBuffer,
                                               self.config.invariantCovariance)

            if first == True:
                last_index = self.config.num_to_process
                if last_index is None:
                    last_index = len(ref)

                for n, (refRecord, moment) in enumerate(zip(ref[:last_index], moments[:last_index])):

                    self.log.info(f'Processing {n}/{last_index}')
                    if refRecord.get('deblend_nChild') != 0:
                        continue

                    if moment.get('bfd_flag') is True:
                        continue

                    if self.config.zFile:
                        if d2d[n] < self.config.zDistCut:
                            match_z = z_redshift[idx[n]]
                            # id_z = z_id[idx]
                        else:
                            self.log.info('No matching redshift')
                            continue

                        if match_z < self.config.zMinCut or match_z > self.config.zMaxCut:
                            self.log.info(f'Does not match redshift range {match_z:0.2}')
                            continue
                        self.z_list[refRecord.getId()] = match_z

                    cov = moment.get(f'bfd_cov_even_{self.config.filters[0]}')
                    if cov[0] > self.config.maxVar and cov[0] < self.config.minVar:
                        continue

                    # Insert the deblended pixels for just this object into all images.
                    for r in replacers.values():
                        r.insertSource(refRecord.getId())

                    try:
                        kgals = self.buildKGalaxy(refRecord, images)
                        kc = KColorGalaxy(self.bfd, kgals, id=refRecord.getId())

                    except Exception as e:
                        kc = None
                        continue

                    dx, badcentering, msg = kc.recenter(self.config.weight_sigma)

                    if badcentering:
                        self.log.info('Bad centering %s', msg)
                        continue

                    template = kc.get_template(dx[0], dx[1])
                    self.templates.append(template)
                    nproc += 1

                    # Remove the deblended pixels for this object so we can process the next one.
                    for r in replacers.values():
                        r.removeSource(refRecord.getId())
                first = False

            for temp in self.templates:
                momentPrior.addTemplate(temp)

            tm = time.time()-tm0
            if nproc > 0:
                self.log.info(f'added {nproc} templates')
                self.log.info(f'time: {tm/60.0} min')
                self.log.info(f'time per: {tm/nproc} sec')

            outputCat = self.buildCatalog(momentPrior)

            dataId['label'] = full_label
            dataId['filter'] = self.config.filters[0]

            metadata = dafBase.PropertyList()
            metadata.set('cov_even', np.array(cov_even.flatten(), dtype=float))
            metadata.set('cov_odd', np.array(cov_odd.flatten(), dtype=float))
            metadata.set('fluxMin', minFlux)
            metadata.set('fluxMax', maxFlux)
            metadata.set('varMin', varMin)
            metadata.set('varMax', varMax)
            if self.config.zFile is not None:
                metadata.set('zFile', self.config.zFile)
                metadata.set('zField', self.config.zField)
            if self.config.zMaxCut is not None:
                metadata.set('zMaxCut', self.config.zMaxCut)
            if self.config.zMinCut is not None:
                metadata.set('zMinCut', self.config.zMinCut)
            if self.config.noiseBin is not None:
                metadata.set('noiseBin', self.config.noiseBin)
            metadata.set('noiseFactor', self.config.noiseFactor)
            metadata.set('priorSigmaCutoff', self.config.priorSigmaCutoff)
            metadata.set('priorSigmaStep', self.config.priorSigmaStep)
            metadata.set('priorSigmaBuffer', self.config.priorSigmaBuffer)
            metadata.set('nsample', self.config.nSample)
            metadata.set('selectionOnly', self.config.selectionOnly)
            metadata.set('invariantCovariance', self.config.invariantCovariance)
            metadata.set('maxXY', self.config.maxXY)
            metadata.set('sigma', self.config.weight_sigma)
            metadata.set('wIndex', self.config.weight_n)
            metadata.set('covFile', self.config.covFile)

            outputCat.getTable().setMetadata(metadata)

            butler.put(outputCat, self.config.output.name, dataId)

        # Restore all original pixels in the images.
        if replacers is not None:
            for r in replacers.values():
                r.end()

    def buildCatalog(self, momentPrior):

        outCat = afwTable.BaseCatalog(self.schema)

        for temp in momentPrior.templates:
            rec = outCat.addNew()
            rec.set(self.mKey, temp.m)
            rec.set(self.dmKey, temp.dm.flatten())
            rec.set(self.dxyKey, temp.dxy.flatten())
            rec.set(self.ndaKey, temp.nda)
            rec.set(self.idKey, temp.id)
            if self.config.zFile:
                rec.set(self.zKey, self.z_list[temp.id])
                # rec.set(self.zIdKey, self.z_id_list[temp.id])

        return outCat

    def buildKGalaxy(self, record, exposures):

        center = record.getCentroid()
        band = self.config.filters[0]
        local_lin_wcs = exposures[band].getWcs().linearizePixelToSky(center, afwGeom.arcseconds)

        jacobian = local_lin_wcs.getLinear().getMatrix()
        sky_pos = exposures[band].getWcs().pixelToSky(center)
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        box = record.getFootprint().getBBox()
        xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())

        bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)

        kgals = []
        for band in self.config.filters:
            exposure = exposures[band]
            factor = exposure.getMetadata().get('variance_scale')
            noise = np.sqrt(np.median(exposure.variance[box].array)/factor)
            image = exposure.image[box].array

            psf_image = exposure.getPsf().computeKernelImage(center).array

            kdata = bfd.generalImage(image, uvref, psf_image, wcs=bfd_wcs, pixel_noise=noise, size=self.config.grid_size)
            conjugate = set(np.where(kdata.conjugate.flatten()==False)[0])
            kgal = self.bfd.KGalaxy(self.weight, kdata.kval.flatten(), kdata.kx.flatten(), kdata.ky.flatten(),
                                    kdata.kvar.flatten(), kdata.d2k, conjugate)
            kgals.append(kgal)
        return kgals

    def selection(self, ref):
        childName = 'deblend_nChild'
        if ref.getParent() == 0 and ref.get(childName) > 0:
            return False
        return True