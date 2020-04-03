
import time
import bfd

import numpy as np
import lsst.afw.geom as afwGeom
import lsst.geom as geom

from collections import defaultdict
from lsst.afw.table import SourceCatalog, SchemaMapper
from lsst.pipe.base import Struct
from lsst.pex.config import Field, ListField

from .measureCoaddsTogether import ProcessCoaddsTogetherTask, ProcessCoaddsTogetherConfig
from .config import BFDConfig
from . import KSigmaWeightF
from .KColorGalaxy import KColorGalaxy


__all__ = ("MeasureCoaddsBfdConfig", "MeasureCoaddsBfdTask")


class MeasureCoaddsBfdConfig(ProcessCoaddsTogetherConfig):
    """
    basic config loads filters and misc stuff
    """
    filters = ListField(
        dtype=str,
        #default=['u', 'g', 'r', 'i', 'z', 'y'],
        default=['i'],
        doc="List of expected bandpass filters."
    )
    weights = ListField(
        dtype=float,
        #default=[0, 0, 1, 1, 1, 0],
        default=[1],
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
    skip_parent = Field(
        dtype=bool,
        default=True,
        doc="skip parent objects"
    )

    def setDefaults(self):
        """
        prefix for the output file
        """
        self.output = "deepCoadd_moments"


class MeasureCoaddsBfdTask(ProcessCoaddsTogetherTask):
    """
    Base class for ngmix tasks
    """
    _DefaultName = "MeasureCoaddsBfdTask"
    ConfigClass = MeasureCoaddsBfdConfig

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
                refSchema = butler.get(self.config.ref_schema + "_schema").schema

        self.ncolors = len(self.config.filters) - 1
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
        self.no_peak_flag = schema.addField('bfd_no_peak', type="Flag",
                                            doc="Flag if object was not observed in given filters")
        if self.config.add_single_bands:
            self.filter_keys = defaultdict(dict)
            self.n_even_single = self.n_even - len(self.config.filters) + 1
            self.n_odd_single = self.n_odd
            for band in self.config.filters:
                self.filter_keys[band]['even'] = schema.addField(f'bfd_even_{band}', type="ArrayF",
                                                                 size=self.n_even_single,
                                                                 doc=f"Even Bfd moments for filter {band}")
                self.filter_keys[band]['odd'] = schema.addField(f'bfd_odd_{band}', type="ArrayF",
                                                                size=self.n_odd_single,
                                                                doc=f"Odd Bfd moments for filter {band}")
                self.filter_keys[band]['cov_even'] = schema.addField(f'bfd_cov_even_{band}', type="ArrayF",
                                                                     size=self.n_even_single*(self.n_even_single+1)//2,
                                                                     doc=f"even moment covariance matrix in filter {band}")
                self.filter_keys[band]['cov_odd'] = schema.addField(f'bfd_cov_odd_{band}', type="ArrayF",
                                                                    size=self.n_odd_single*(self.n_odd_single+1)//2,
                                                                    doc=f"odd moment covariance matrix in filter {band}")

        return schema

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

        # Make an empty catalog
        output = SourceCatalog(self.schema)

        # Add mostly-empty rows to it, copying IDs from the ref catalog.
        output.extend(ref, mapper=self.mapper)

        min_index = self.config.start_index
        if self.config.num_to_process is None:
            max_index = len(ref)
        else:
            max_index = self.config.start_index + self.config.num_to_process

        for n, (refRecord, outRecord) in enumerate(zip(ref, output)):
            if n < min_index or n >= max_index:
                continue


            if refRecord.get('deblend_nChild') != 0:
                outRecord.set(self.parent_flag, 1)
                if self.config.skip_parent:
                    outRecord.set(self.flag, 1)
                    continue

            if n % 1000 == 0:
                self.log.info('index: %06d/%06d' % (n, max_index))
            nproc += 1

            hasPeak = False
            for filt in self.config.filters:
                if refRecord.get(f'merge_peak_{filt}'):
                    hasPeak = True
            if hasPeak == False:
                outRecord.set(self.no_peak_flag, 1)

            outRecord.setFootprint(None)  # copied from ref; don't need to write these again

            # Insert the deblended pixels for just this object into all images.
            for r in replacers.values():
                r.insertSource(refRecord.getId())

            try:
                kgals = self.buildKGalaxy(refRecord, images)
                kc = KColorGalaxy(self.bfd, kgals, self.config.weights)

            except Exception:
                kc = None

            if kc is None:
                outRecord.set(self.flag, 1)
                for r in replacers.values():
                    r.removeSource(refRecord.getId())
                continue

            dx, badcentering, msg = kc.recenter(self.config.weight_sigma)

            # if centering failed, measure moments still without shifting
            if badcentering:
                self.log.debug('Bad centering %s', msg)
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

            if self.config.add_single_bands:
                for i, band in enumerate(self.config.filters):
                    seven, sodd, scov_even, scov_odd = kc.get_single_moment(i, dx[0], dx[1])
                    scov_even_save = []
                    scov_odd_save = []
                    for ii in range(scov_even.shape[0]):
                        scov_even_save.extend(scov_even[ii][ii:])
                    for ii in range(scov_odd.shape[0]):
                        scov_odd_save.extend(scov_odd[ii][ii:])

                    outRecord.set(self.filter_keys[band]['even'], np.array(seven, dtype=np.float32))
                    outRecord.set(self.filter_keys[band]['odd'], np.array(sodd, dtype=np.float32))
                    outRecord.set(self.filter_keys[band]['cov_even'], np.array(scov_even_save, dtype=np.float32))
                    outRecord.set(self.filter_keys[band]['cov_odd'], np.array(scov_odd_save, dtype=np.float32))

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


        kgals = []
        for band in self.config.filters:
            exposure = exposures[band]
            factor = exposure.getMetadata().get('variance_scale')
            exp_box = geom.Box2I(box)
            exp_box.clip(exposure.getBBox())

            xy_pos = (center.getX() - exp_box.getMinX(), center.getY() - exp_box.getMinY())
            bfd_wcs = bfd.WCS(jacobian, xyref=xy_pos, uvref=uvref)

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
