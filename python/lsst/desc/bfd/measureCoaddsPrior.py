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
import lsst.geom as geom

from lsst.pipe.tasks.multiBand import getShortFilterName
from lsst.meas.base import NoiseReplacerConfig, NoiseReplacer
from lsst.afw.table import SourceCatalog, SchemaMapper
from lsst.pipe.base import Struct
from lsst.pex.config import Field, ListField, ConfigField, Config, ChoiceField

from .measureCoaddsTogether import ProcessCoaddsTogetherTask, ProcessCoaddsTogetherConfig
from .config import BFDConfig
from . import KSigmaWeightF, UniformDeviate
from .KColorGalaxy import KColorGalaxy


__all__ = ("MeasureCoaddsPriorConfig", "MeasureCoaddsPriorTask", "make_templates")



def make_templates(rand,kc, sigma_xy, sigma_flux=1., sn_min=0., sigma_max=6.5, sigma_step=1., 
                   xy_max=2.,tid=0, weight_sigma=4, sample=1,
                   **kwargs):
    ''' Return a list of Template instances that move the object on a grid of
    coordinate origins that keep chisq contribution of flux and center below
    the allowed max.
    sigma_xy    Measurement error on target x & y moments (assumed equal, diagonal)
    sigma_flux  Measurement error on target flux moment
    sn_min      S/N for minimum flux cut applied to targets
    sigma_max   Maximum number of std deviations away from target that template will be used
    sigma_step  Max spacing between shifted templates, in units of measurement sigmas
    xy_max      Max allowed centroid shift, in sky units (prevents runaways)
    '''
#    xyshift, error, msg = self.recenter()
#    if error:
#       return None, "Center wandered too far from starting guess or failed to converge"

    dx, badcentering, msg = kc.recenter(6)
    #logging.debug(f'center {dx}')
    if badcentering:
        return []
    #logging.debug(f'moment: {kc.get_moment(dx[0], dx[1])}')
    # Determine derivatives of 1st moments on 2 principal axes,
    # and steps will be taken along these grid axes.

    #jacobian0 = self.xy_jacobian(np.zeros(2))
    jacobian0 = kc.xy_jacobian(dx)
    eval, evec = np.linalg.eigh(jacobian0)
    #logging.debug(f'eig {eval} jac0 {jacobian0}')
    if np.any(eval >= 0.):
        return []#None, "Template galaxy center is not at a flux maximum"

    detj0 = np.linalg.det(jacobian0)  # Determinant of Jacobian
    xy_step = np.abs(sigma_step * sigma_xy / eval)

    da = xy_step[0] * xy_step[1]
    #print('step',xy_step,da, detj0)
    # Offset the xy grid by random phase in the grid
    #xy_offset = np.random.random(2) - 0.5
    xy_offset = rand.uniform(0,1,size=2) - 0.5
    
    # Now explore contiguous region of xy grid that yields useful templates.
    result = []
    grid_try = set(((0, 0),))  # Set of all grid points remaining to try
    grid_done = set()           # Grid points already investigated

    flux_min = sn_min * sigma_flux
    #logging.debug(f'flux_min {flux_min}')
    while len(grid_try) > 0:
        # Try a new grid point
        mn = grid_try.pop()
        #logging.debug(f' try  {mn}')
        grid_done.add(mn)  # Mark it as already attempted
        # Offset and scale
        xy = np.dot(evec, xy_step*(np.array(mn) + xy_offset))
        # Ignore if we have wandered too far
        #logging.debug(f' xy, {xy}, {np.dot(xy, xy)}, {xy_max**2}')
        if np.dot(xy, xy) > xy_max*xy_max:
            continue

        #m = self.get_moment(xy[0], xy[1])
        mom, cov = kc.get_moment(dx[0]+xy[0], dx[1]+xy[1])
        #logging.debug(f' shift {dx[0]+xy[0]}, {dx[1]+xy[1]}, {mom.m}')
        even = mom.m
        odd = mom.xy

        #detj = 0.25 * ( e[m.MR]**2-e[m.M1]**2 - e[m.M2]**2)
        detj = 0.25 * (even[kc.bfd_config.MR]**2 -
                       even[kc.bfd_config.M1]**2 - even[kc.bfd_config.M2]**2)
        # Ignore if determinant of Jacobian has gone negative, meaning
        # we have crossed out of convex region for flux
        #logging.debug(f' detj {detj}')
        if detj <= 0.:
            #logging.debug(' det < 0')
            continue

        # Accumulate chisq that this template would have for a target
        # First: any target will have zero MX, MY

        #chisq = (m.odd[m.MX]**2 + m.odd[m.MY]**2) / sigma_xy**2
        chisq = (odd[kc.bfd_config.MX]**2 +
                 odd[kc.bfd_config.MY]**2) / sigma_xy**2
        #logging.debug(f' chi2 odd {chisq}')
        # Second: there is suppression by jacobian of determinant
        #logging.debug(f' jacobian supression {-2. * np.log(detj/detj0)}')
        chisq += -2. * np.log(detj/detj0)
        # Third: target flux will never be below flux_min

        # if (e[m.M0] < flux_min):
        #    chisq += ((flux_min -e[m.M0])/sigma_flux)**2

        if (even[kc.bfd_config.MF] < flux_min):
            #logging.debug(' below flux_min {even[kc.bfd_config.MF] }')
            chisq += ((flux_min - even[kc.bfd_config.MF])/sigma_flux)**2

        #logging.debug(f'chi2 , {chisq}, {sigma_max*sigma_max}')
        if chisq <= sigma_max*sigma_max:
            # This is a useful template!  Add it to output list

            #tmpl = self.get_template(xy[0],xy[1])
            tmpl = kc.get_template(dx[0]+xy[0], dx[1]+xy[1])
            tmpl.nda = tmpl.nda * da / sample
            tmpl.jSuppression = detj / detj0
            tmpl.id = tid
            result.append(tmpl)
            #logging.debug(f'deriv: {tmpl.realMDerivs()}')
            # Try all neighboring grid points not yet tried
            for mn_new in ((mn[0]+1, mn[1]),
                           (mn[0]-1, mn[1]),
                           (mn[0], mn[1]+1),
                           (mn[0], mn[1]-1)):
                if mn_new not in grid_done:
                    grid_try.add(mn_new)
    #logging.debug(grid_done)
    return result


class MeasureCoaddsPriorConfig(ProcessCoaddsTogetherConfig):
    """
    basic config loads filters and misc stuff
    """
    filters = ListField(
        dtype=str,
        default=['i'],
        doc="List of expected bandpass filters."
    )
    requirePeak = Field(
        dtype=bool,
        default=True,
        optional=True,
        doc='require peak detection in list of filters',
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
    use_mag = Field(
        dtype=bool,
        default=False,
        doc="include magnification"
    )
    use_conc = Field(
        dtype=bool,
        default=False,
        doc="include concentration"
    )
    weight_n = Field(
        dtype=int,
        default=4,
        doc="n for k-sigma weight function"
    )
    weight_sigma = Field(
        dtype=float,
        default=0.76,
        doc="sigma for k-sigma weight function"
    )
    max_shift = Field(
        dtype=float,
        default=4,
        doc="maximum centroid shift"
    )
    add_single_bands = Field(
        dtype=bool,
        default=False,
        doc="add single band measurements"
    )
    coaddName = Field(
        dtype=str,
        default='deep',
        doc="name of coadd"
    )
    grid_size = Field(
        dtype=int,
        default=-1,
        doc="minimum size of pixel grid when computing moments"
    )
    footprint_size = Field(
        dtype=int,
        default=64,
        doc="size of footprint, must be even"
    )
    snMin = pexConfig.Field(
        dtype=float,
        default=10,
        optional=True,
        doc="Minimun flux S/N"
    )
    snMax = pexConfig.Field(
        dtype=float,
        default=35,
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
        default=2.,
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
        default=10000,
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
        default=1.,
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
        default=None,
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
    fluxBin = pexConfig.Field(
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
        self.output = "deepCoadd_prior"


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
                refSchema = butler.get(self.config.ref_schema + "_schema").schema
        self.ncolors = len(self.config.filters) - 1
        self.bfd = BFDConfig(use_conc=self.config.use_conc, use_mag=self.config.use_mag,
                             ncolors=self.ncolors)
        self.n_even = self.bfd.BFDConfig.MSIZE
        self.n_odd = self.bfd.BFDConfig.XYSIZE
        self.weight = KSigmaWeightF(self.config.weight_sigma, self.config.weight_n)
        self.schema = self.defineSchema(refSchema)
        self.ud = UniformDeviate()

        if self.config.randomSeed:
            np.random.seed(self.config.randomSeed)
        self.rand = np.random.RandomState()
        self.templates = []

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
        self.mKey = schema.addField("m", doc="moments", type="ArrayF",
                                    size=self.bfd.BFDConfig.MXYSIZE)
        self.dmKey = schema.addField("dm", doc="moment derivatives", type="ArrayF",
                                     size=self.bfd.BFDConfig.MSIZE*self.bfd.BFDConfig.DSIZE)
        self.dxyKey = schema.addField("dxy", doc="xy moment derivaties", type="ArrayF",
                                      size=self.bfd.BFDConfig.XYSIZE*self.bfd.BFDConfig.DSIZE)
        self.ndaKey = schema.addField("nda", doc="nda", type=np.float)
        self.idKey = schema.addField("bfd_id", doc="id", type=np.int64)
        if self.config.zFile:
            self.zKey = schema.addField("z", doc="redshift", type=np.float)
            self.zIdKey = schema.addField("z_id", doc="redshift", type=np.int64)

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
            images[filt] = patchRef.get(self.config.images)

            fpCat = patchRef.get(self.config.deblendCatalog)
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
            self.z_id_list = {}

            coord = SkyCoord(ra=ref['coord_ra']*u.rad, dec=ref['coord_dec']*u.rad)

            idx, d2d, d3d = match_coordinates_sky(coord, z_catalog)
            d2d = d2d.arcsec


        last_index = self.config.num_to_process
        if last_index is None:
            last_index = len(ref)

        kgals = []
        for n, (refRecord, moment) in enumerate(zip(ref[:last_index], moments[:last_index])):

            if self.config.sample > 0:
                if np.random.rand() > self.config.sample:
                    continue

            if n%1000 == 0:
                self.log.info(f'Processing {n}/{last_index}')
            if refRecord.get('deblend_nChild') != 0:
                continue

            if refRecord.get('detect_isPrimary') == False:
                continue

            if self.config.requirePeak:
                hasPeak = False
                for filt in self.config.filters:
                    if refRecord.get(f'merge_peak_{filt}'):
                        hasPeak = True
                if hasPeak == False:
                    self.log.debug('No matching peak in filters')
                    continue

            if moment.get('bfd_flag') is True:
                continue

            if self.config.zFile:
                if d2d[n] < self.config.zDistCut:
                    match_z = z_redshift[idx[n]]
                    # id_z = z_id[idx]
                else:
                    self.log.debug('No matching redshift')
                    continue

                if match_z < self.config.zMinCut or match_z > self.config.zMaxCut:
                    self.log.debug(f'Does not match redshift range {match_z:0.2}')
                    continue
                self.z_list[refRecord.getId()] = match_z
                self.z_id_list[refRecord.getId()] = idx[n]

            cov = moment.get(f'bfd_cov_even')

            if cov[0] > self.config.maxVar and cov[0] < self.config.minVar:
                continue

            # Insert the deblended pixels for just this object into all images.
            for r in replacers.values():
                r.insertSource(refRecord.getId())

            try:
                kgal = self.buildKGalaxy(refRecord, images)
                kc = KColorGalaxy(self.bfd, kgal, id=refRecord.getId())

            except Exception:
                kc = None
                for r in replacers.values():
                    r.removeSource(refRecord.getId())
                continue

            dx, badcentering, msg = kc.recenter(self.config.weight_sigma)


            if badcentering:
                self.log.debug('Bad centering %s', msg)
                for r in replacers.values():
                    r.removeSource(refRecord.getId())
                continue

            kgals.append(kc)


            # Remove the deblended pixels for this object so we can process the next one.
            for r in replacers.values():
                r.removeSource(refRecord.getId())

        # Restore all original pixels in the images.
        if replacers is not None:
            for r in replacers.values():
                r.end()


        self.log.info('Measured %d moments from galaxies' % len(kgals))

        # we only need to compute the templates the first time through the loop
        for cov_even, cov_odd, label, varMin, varMax in covList:

            if (label not in self.config.useLabels) and len(self.config.useLabels) > 0:
                self.log.info("Label %s not in %s" % (label, self.config.useLabels))
                continue

            # Build full label
            full_label = label + self.config.label
            if self.config.selectionOnly:
                full_label += '_selection'

            if self.config.fluxBin is not None:
                full_label += '_f%d' % self.config.fluxBin

            if self.config.zBin is not None:
                full_label += '_z%d' % (self.config.zBin)

            self.log.info('Processing label %s' % label)
            sigmaFlux = np.sqrt(cov_even[0, 0])
            sigma_xy = np.sqrt(cov_odd[0, 0])

            minFlux = self.config.snMin*sigmaFlux
            if self.config.fluxMin is not None:
                minFlux = self.config.fluxMin
            elif self.config.magMin is not None:
                # Current assumption is that coadd zeropoint is 27, true for HSC/DC2
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

            atemplates = []
            for kc in kgals:
                ltemplates = make_templates(self.rand, kc, sigma_xy, sn_min=minFlux/sigmaFlux,
                                sigma_flux=sigmaFlux, sigma_step=self.config.priorSigmaStep,
                                sigma_max=self.config.priorSigmaCutoff, xy_max=self.config.maxXY, tid=kc.id,
                                weight_sigma=self.config.weight_sigma, sample=self.config.sample)

                atemplates.extend(ltemplates)

            for temp in atemplates:
                momentPrior.addTemplate(temp)

            self.log.info(f'added {len(atemplates)} templates')

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
            if self.config.zBin is not None:
                metadata.set('zBin', self.config.zBin)
            if self.config.fluxBin  is not None:
                metadata.set('fluxBin', self.config.fluxBin)
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

            butler.put(outputCat, self.config.output, dataId)



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
                rec.set(self.zIdKey, self.z_id_list[temp.id])

        return outCat

    def buildKGalaxy(self, record, exposures):

        center = record.getCentroid()

        # For wcs use the first band since a single tract should have the same wcs
        band = self.config.filters[0]
        local_lin_wcs = exposures[band].getWcs().linearizePixelToSky(center, afwGeom.arcseconds)

        jacobian = local_lin_wcs.getLinear().getMatrix()
        sky_pos = exposures[band].getWcs().pixelToSky(center)
        uvref = (sky_pos.getRa().asArcseconds(), sky_pos.getDec().asArcseconds())

        if self.config.footprint_size is None or self.config.footprint_size < 0:
            box = record.getFootprint().getBBox()
        else:
            box = geom.Box2I(geom.Point2I(int(center.getX()), int(center.getY())),
                             geom.Extent2I(1,1))
            box.grow(self.config.footprint_size//2)

        xy_pos = (center.getX() - box.getMinX(), center.getY() - box.getMinY())

        bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)

        kgals = []
        for band in self.config.filters:
            exposure = exposures[band]
            factor = exposure.getMetadata().get('variance_scale')
            exp_box = geom.Box2I(box)
            exp_box.clip(exposure.getBBox())
            noise = np.sqrt(np.median(exposure.variance[exp_box].array))
            image = exposure.image[exp_box].array

            psf_image = exposure.getPsf().computeKernelImage(center).array

            kdata = bfd.generalImage(image, uvref, psf_image, wcs=bfd_wcs, pixel_noise=noise, 
                                     size=self.config.grid_size)
            conjugate = set(np.where(kdata.conjugate.flatten()==False)[0])
            kgal = self.bfd.KGalaxy(self.weight, kdata.kval.flatten(), kdata.kx.flatten(),
                                    kdata.ky.flatten(), kdata.kvar.flatten(), kdata.d2k, conjugate)
            kgals.append(kgal)
        return kgals
