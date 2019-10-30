
import time
import numpy as np
import bfd
from collections import defaultdict
import lsst.afw.geom as afwGeom
from lsst.afw.table import SourceCatalog, SchemaMapper

from lsst.pipe.base import Struct
from lsst.meas.base import NoiseReplacerConfig, NoiseReplacer
from lsst.pex.config import Field, ListField, ConfigField, Config, ChoiceField
from lsst.pipe.tasks.multiBand import MergeSourcesRunner, getShortFilterName

from .measureCoaddsTogether import ProcessCoaddsTogetherTask, ProcessCoaddsTogetherConfig
from .config import BFDConfig
from .KColorGalaxy import KColorGalaxy
from . import KSigmaWeightF


__all__ = ("MeasureCoaddsBfdSingleConfig", "MeasureCoaddsBfdSingleTask")


class MeasureCoaddsBfdSingleConfig(ProcessCoaddsTogetherConfig):
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
    weight_sigma = Field(dtype=float, default=1.25, doc="sigma for k-sigma weight function")
    max_shift = Field(dtype=float, default=4, doc="maximum centroid shift")
    add_single_bands = Field(dtype=bool, default=True, doc="add single band measurements")
    coaddName = Field(dtype=str, default='deep', doc="name of coadd")
    grid_size = Field(dtype=int, default=48, doc="minimum size of pixel grid when computing moments")

    def setDefaults(self):
        """
        prefix for the output file
        """
        self.output.name = "deepCoadd_moments_single"


class MeasureCoaddsBfdSingleTask(ProcessCoaddsTogetherTask):
    """
    Base class for ngmix tasks
    """
    _DefaultName = "MeasureCoaddsBfdSingleTask"
    ConfigClass = MeasureCoaddsBfdSingleConfig

    def __init__(self, *, config=None, refSchema=None, butler=None, initInputs=None, **kwds):

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
        self.ncolors = 0
        self.bfd = BFDConfig(use_conc=self.config.use_conc, use_mag=self.config.use_mag,
                             ncolors=self.ncolors)
        self.n_even = self.bfd.BFDConfig.MSIZE
        self.n_odd = self.bfd.BFDConfig.XYSIZE
        self.weight = KSigmaWeightF(self.config.weight_sigma, self.config.weight_n)
        self.schema = self.defineSchema(refSchema)

    def defineSchema(self, refSchema):

        self.mapper = SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(SourceCatalog.Table.makeMinimalSchema(), True)
        schema = self.mapper.getOutputSchema()

        self.even = schema.addField('bfd_even', type="ArrayF",
                                    size=self.n_even, doc="Even Bfd moments")
        self.odd = schema.addField('bfd_odd', type="ArrayF",
                                   size=self.n_odd, doc="odd moments")
        self.shift = schema.addField('bfd_shift', type="ArrayF",
                                     size=2, doc="amount shifted to null moments")
        self.cov_even = schema.addField('bfd_cov_even', type="ArrayF",
                                        size=self.n_even*(self.n_even+1)//2,
                                        doc="even moment covariance matrix")
        self.cov_odd = schema.addField('bfd_cov_odd', type="ArrayF",
                                       size=self.n_odd*(self.n_odd+1)//2,
                                       doc="odd moment covariance matrix")
        self.flag = schema.addField('bfd_flag', type="Flag", doc="Set to 1 for any fatal failure")
        self.centroid_flag = schema.addField('bfd_flag_centroid', type="Flag",
                                             doc="Set to 1 for any fatal failure of centroid")
        self.parent_flag = schema.addField('bfd_flag_parent', type="Flag",
                                           doc="Set to 1 for parents")
        return schema

    def runDataRef(self, patchRefList):
        """Run this task via CmdLineTask and Gen2 Butler.
        Parameters
        ----------
        patchRefList : `list` of `lsst.daf.persistence.ButlerDataRef`
            A list of DataRefs for all filters in a single patch.
        """
        #import pdb;pdb.set_trace()
        images = {}
        replacers = {}
        mergedDataId = {"tract": patchRefList[0].dataId["tract"],
                        "patch": patchRefList[0].dataId["patch"]}
        butler = patchRefList[0].butlerSubset.butler
        ref = butler.get("deepCoadd_ref", dataId=mergedDataId)
        imageId = butler.get("deepMergedCoaddId", dataId=mergedDataId)
        for patchRef in patchRefList:
            filt = getShortFilterName(patchRef.dataId["filter"])
            images[filt] = patchRef.get(self.config.images.name)

            fpCat = patchRef.get(self.config.deblendCatalog.name)
            footprints = {rec.getId(): (rec.getParent(), rec.getFootprint()) for rec in fpCat}
            replacers[filt] = NoiseReplacer(self.config.deblendReplacer, exposure=images[filt],
                                            footprints=footprints, exposureId=imageId)
        results = self.run(images, ref, imageId=imageId, replacers=replacers)
        mergedDataId['filter'] = self.config.filters[0]
        butler.put(results.output, self.config.output.name, dataId=mergedDataId)

    def run(self, images, ref, replacers, imageId):
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

        if len(images) != len(self.config.filters):
            self.log.info('Number of filters does not match the list of images given.  Skipping')
            return None

        tm0 = time.time()
        nproc = 0
        #import pdb;pdb.set_trace()
        # Make an empty catalog
        output = SourceCatalog(self.schema)

        # Add mostly-empty rows to it, copying IDs from the ref catalog.
        output.extend(ref, mapper=self.mapper)

        min_index = self.config.start_index
        if self.config.num_to_process is None:
            max_index = len(ref)
        else:
            max_index = self.config.start_index + self.config.num_to_process

        #for n, (refRecord) in enumerate(zip(ref)):
        for n, (refRecord, outRecord) in enumerate(zip(ref, output)):
            if n < min_index or n >= max_index:
                continue

            if refRecord.get('deblend_nChild') != 0:
                outRecord.set(self.flag, 1)
                outRecord.set(self.parent_flag, 1)
                continue

            #outRecord = output.table.copyRecord(refRecord, self.mapper)
            #output._append(outRecord)

            self.log.info('index: %06d/%06d' % (n, max_index))
            nproc += 1

            outRecord.setFootprint(None)  # copied from ref; don't need to write these again

            # Insert the deblended pixels for just this object into all images.
            for r in replacers.values():
                r.insertSource(refRecord.getId())

            try:
                kgals = self.buildKGalaxy(refRecord, images)
                kc = KColorGalaxy(self.bfd, kgals)

            except Exception as e:
                kc = None

            if kc is None:
                outRecord.set(self.flag, 1)
                continue

            dx, badcentering, msg = kc.recenter(self.config.weight_sigma)

            if badcentering:
                self.log.info('Bad centering %s', msg)
                outRecord.set(self.flag, 1)
                outRecord.set(self.centroid_flag, 1)
                dx = [0, 0]

            mom, cov = kc.get_moment(dx[0], dx[1], True)
            mom_even = mom.m
            mom_odd = mom.xy
            cov_even = cov.m
            cov_odd = cov.xy

            cov_even_save = []
            cov_odd_save = []
            for ii in range(cov_even.shape[0]):
                cov_even_save.extend(cov_even[ii][ii:])
            for ii in range(cov_odd.shape[0]):
                cov_odd_save.extend(cov_odd[ii][ii:])

            outRecord.set(self.even, np.array(mom_even, dtype=np.float32))
            outRecord.set(self.odd, np.array(mom_odd, dtype=np.float32))
            outRecord.set(self.cov_even, np.array(cov_even_save, dtype=np.float32))
            outRecord.set(self.cov_odd, np.array(cov_odd_save, dtype=np.float32))
            outRecord.set(self.shift, np.array([dx[0], dx[1]], dtype=np.float32))

            # Remove the deblended pixels for this object so we can process the next one.
            for r in replacers.values():
                r.removeSource(refRecord.getId())
            del kgals
            del kc
            del mom
            del cov
        # Restore all original pixels in the images.
        if replacers is not None:
            for r in replacers.values():
                r.end()

        tm = time.time()-tm0
        self.log.info('time: %g min' % (tm/60.0))
        self.log.info('time per: %g sec' % (tm/nproc))

        return Struct(output=output[min_index:max_index])

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