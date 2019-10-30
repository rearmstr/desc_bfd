import time
import bfd
import numpy as np
from collections import defaultdict
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
from astropy import units as u
from astropy.table import Table

import lsst.pex.config as pexConfig
import lsst.afw.geom as afwGeom
import lsst.afw.table as afwTable

from lsst.daf.persistence import Butler
from lsst.pipe.tasks.multiBand import getShortFilterName
from lsst.meas.base import NoiseReplacerConfig, NoiseReplacer
from lsst.afw.table import SourceCatalog, SchemaMapper
from lsst.pipe.base import Struct
from lsst.pex.config import Field, ListField, ConfigField, Config, ChoiceField

from .measureCoaddsTogether import ProcessCoaddsTogetherTask, ProcessCoaddsTogetherConfig
from .config import BFDConfig
from . import KSigmaWeightF, UniformDeviate
from .KColorGalaxy import KColorGalaxy
import glob
from lsst.pipe.base import TaskRunner
from lsst.pipe.base import (
    CmdLineTask, ArgumentParser,
    PipelineTask, InputDatasetField, OutputDatasetField, QuantumConfig,
)

__all__ = ("MeasureCoaddsPqrConfig", "MeasureCoaddsPqrTask")

class MyTaskRunner(TaskRunner):
    r"""A `TaskRunner` that combines
    """

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        files = glob.glob(parsedCmd.dir)

        return [(files, kwargs)]


class MeasureCoaddsPqrConfig(ProcessCoaddsTogetherConfig):
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
    zMatchDist = pexConfig.Field(
        dtype=float,
        optional=True,
        default=0.5,
        doc="redshift distance match"
    )
    invariantCovariance = pexConfig.Field(
        dtype=bool,
        default=False,
        optional=True,
        doc="will all the galaxies have the same covariance"

    )
    numProc = pexConfig.Field(
        dtype=int,
        default=1,
        optional=True,
        doc="number of processors to use"
    )
    chunk = pexConfig.Field(
        dtype=int,
        default=15,
        optional=True,
        doc="number of processors to use"
    )
    priorRerun = pexConfig.Field(
        dtype=str,
        optional=False,
        doc="rerun for the prior"
    )
    priorTracts = pexConfig.ListField(
        dtype=int,
        default=[5063, 4849, 4848],
        optional=False,
        doc="tract for the prior"
    )
    priorLabel = pexConfig.Field(
        dtype=str,
        default='b0',
        optional=False,
        doc="label for the prior"
    )
    priorFilter = pexConfig.Field(
        dtype=str,
        default='r',
        optional=False,
        doc="filter for the prior"
    )
    maxPriorFiles = pexConfig.Field(
        dtype=int,
        default=-1,
        optional=False,
        doc="filter for the prior"
    )
    priorPatches = pexConfig.ListField(
        dtype=str,
        default=None,
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    zBins = pexConfig.ListField(
        dtype=float,
        default=[0, 0.2, 0.4, 0.6, 0.8, 1.0],
        optional=True,
        doc="Dictionary for the fraction of the galaxies used for eac"
    )
    threads = pexConfig.Field(
        dtype=int,
        default=8,
        optional=False,
        doc="Number of threads"
    )
    chunk = pexConfig.Field(
        dtype=int,
        default=20,
        optional=False,
        doc="Number per chunk"
    )

    def setDefaults(self):
        """
        prefix for the output file
        """
        self.output.name = "deepCoadd_prior"

class MeasureCoaddsPqrTask(ProcessCoaddsTogetherTask):
    """
    Base class for ngmix tasks
    """
    _DefaultName = "MeasureCoaddsPqrTask"
    ConfigClass = MeasureCoaddsPqrConfig
    RunnerClass = MyTaskRunner

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
        print('Reinitialized')
        self.initialized = False
    @classmethod
    def _makeArgumentParser(cls):
        # Customize argument parsing for CmdLineTask.
        parser = ArgumentParser(name=cls._DefaultName)
        # This should be config.images.name, but there's no way to pass that
        # information in here in Gen2.
        parser.add_argument("--dir", dest="dir", default='./',
                            help="location of files")
        return parser

    def defineSchema(self, refSchema):

        self.mapper = SchemaMapper(refSchema)
        self.mapper.addMinimalSchema(SourceCatalog.Table.makeMinimalSchema(), True)
        schema = self.mapper.getOutputSchema()
        self.pqrKey = schema.addField("pqr", doc="pqr", type="ArrayF",
                                      size=self.bfd.BFDConfig.DSIZE)
        self.momKey = schema.addField("moment", doc="moment", type="ArrayF",
                                      size=self.n_even)
        self.momCovKey = schema.addField("moment_cov", doc="moment", type="ArrayF",
                                         size=self.n_even*(self.n_even+1)//2)
        self.numKey = schema.addField("n_templates", doc="number", type=np.int64)
        self.uniqKey = schema.addField("n_unique", doc="unique", type=np.int32)
        self.zKey = schema.addField("z", doc="redshift", type=np.float)
        self.g1Key = schema.addField("g1", doc="redshift", type=np.float)
        self.g2Key = schema.addField("g2", doc="redshift", type=np.float)
        self.kappaKey = schema.addField("kappa", doc="redshift", type=np.float)
        self.magKey = schema.addField("mag", doc="redshift", type=np.float)
        self.labelKey = schema.addField("label", doc="redshift", type=str, size=10)
            # self.zIdKey = schema.addField("z_id", doc="redshift", type=np.int64)

        return schema

    def runDataRef(self, files):
        """Run this task via CmdLineTask and Gen2 Butler.
        Parameters
        ----------
        patchRefList : `list` of `lsst.daf.persistence.ButlerDataRef`
            A list of DataRefs for all filters in a single patch.
        """
        self.prep()

        for file in files:
            cat = afwTable.BaseCatalog.readFits(file)

            mask = ((cat['moments_r_even'][:, 0] >= self.fluxMin) &
                    (cat['moments_r_even'][:, 0] < self.fluxMax) &
                    (cat['moments_r_cov_even'][:, 0] >= self.varMin) &
                    (cat['moments_r_cov_even'][:, 0] < self.varMax) &
                    (cat['z'] >= self.zMin) &
                    (cat['z'] < self.zMax))
            tgs = []
            for rec in cat[mask]:
                cov_even = rec.get('moments_r_cov_even')
                cov_odd = rec.get('moments_r_cov_odd')
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

                cov = self.bfd.MomentCov(full_cov_even, full_cov_odd)
                #mom_odd = np.array([0, 0]
                moment = self.bfd.Moment(rec.get('moments_r_even'), rec.get('moments_r_odd'))
                pos = np.array([rec.get('ra'), rec.get('dec')])
                tg = self.bfd.TargetGalaxy(moment, cov, pos, rec.get('id'))
                tgs.append(tg)

            results = self.prior.getPqrCatalog(tgs[:10], self.config.threads, self.config.chunk)

            outcat = afwTable.BaseCatalog(self.schema)
            raKey = self.schema['coord_ra'].asKey()
            decKey = self.schema['coord_ra'].asKey()
            idKey = self.schema['id'].asKey()

            for r, rec in zip(results, cat[mask][:10]):
                out = outcat.addNew()
                out.set(self.pqrKey, r[0]._pqr)
                out.set(self.numKey, r[1])
                out.set(self.uniqKey, r[2])
                out.set(self.momKey, rec.get('moments_r_even'))
                out.set(self.momCovKey, rec.get('moments_r_cov_even'))
                out.set(self.zKey, rec.get('z'))
                out.set(self.g1Key, rec.get('g1'))
                out.set(self.g1Key, rec.get('g2'))
                out.set(self.kappaKey, rec.get('kappa'))
                out.set(self.magKey, rec.get('mag'))
                out.set(self.labelKey, self.config.label)
                out.set(raKey, rec.get('ra')*afwGeom.degrees)
                out.set(decKey, rec.get('dec')*afwGeom.degrees)
                out.set(idKey, rec.get('id'))

            outfile = file.replace('moment', 'pqr')
            outcat.writeFits(outfile)

    def prep(self):

        if self.initialized:
            return

        self.prior = None
        priorFiles = []
        priorButler = Butler(self.config.priorRerun)
        prior_skyMap = priorButler.get('deepCoadd_skyMap')

        for tract in self.config.priorTracts:
            for patchInfo in prior_skyMap[tract]:
                patch = '%d,%d' % patchInfo.getIndex()

                if self.config.priorPatches:
                    if patch not in self.config.priorPatches:
                        continue

                if priorButler.datasetExists('deepCoadd_prior', tract=tract, patch=patch,
                                             filter=self.config.priorFilter, label=self.config.priorLabel):
                    priorFiles.append(priorButler.getUri('deepCoadd_prior',
                                                         tract=tract, patch=patch,
                                                         filter=self.config.priorFilter,
                                                         label=self.config.priorLabel))

        max_file = len(priorFiles)
        if self.config.maxPriorFiles > 0:
            max_file = self.config.maxPriorFiles

        self.zBin = None
        for file in priorFiles[:max_file]:
            if file.find('_parent') > 0:
                self.log.info("Skipping %s, from parent" % file)
                continue
            self.log.info("Adding prior %s" % file)
            try:
                cat = afwTable.BaseCatalog.readFits(file)
                md = cat.getTable().getMetadata().toDict()

                if self.prior is None:
                    self.fluxMin = md['FLUXMIN']
                    self.fluxMax = md['FLUXMAX']
                    self.varMin = md['VARMIN']
                    self.varMax = md['VARMAX']
                    cov_even = np.array(md['COV_EVEN'])
                    cov_odd = np.array(md['COV_ODD'])
                    self.zMax = md['ZMAXCUT']
                    self.zMin = md['ZMINCUT']
                    self.noiseFactor = md['noiseFactor']
                    self.priorSigmaCutoff = md['priorSigmaCutoff']
                    self.priorSigmaStep = md['priorSigmaStep']
                    self.priorSigmaBuffer = md['priorSigmaBuffer']
                    self.nSample = md['NSAMPLE']
                    self.selectionOnly = md['selectionOnly']
                    self.invariantCovariance = md['invariantCovariance']

                    covMat = self.bfd.MomentCov(cov_even.reshape(self.n_even, self.n_even),
                                                cov_odd.reshape(self.n_odd, self.n_odd))
                    self.prior = self.bfd.KDTreePrior(self.fluxMin, self.fluxMax, covMat, self.ud,
                                                      self.nSample, self.selectionOnly,
                                                      self.noiseFactor, self.priorSigmaStep,
                                                      self.priorSigmaCutoff, self.priorSigmaBuffer,
                                                      self.invariantCovariance)
                else:
                    fluxMin = md['FLUXMIN']
                    fluxMax = md['FLUXMAX']
                    varMin = md['VARMIN']
                    varMax = md['VARMAX']
                    cov_even = np.array(md['COV_EVEN'])
                    cov_odd = np.array(md['COV_ODD'])
                    zMax = md['ZMAXCUT']
                    zMin = md['ZMINCUT']
                    noiseFactor = md['noiseFactor']
                    priorSigmaCutoff = md['priorSigmaCutoff']
                    priorSigmaStep = md['priorSigmaStep']
                    priorSigmaBuffer = md['priorSigmaBuffer']
                    nSample = md['NSAMPLE']
                    selectionOnly = md['selectionOnly']
                    invariantCovariance = md['invariantCovariance']

                    mismatch = False
                    if fluxMin != self.fluxMin:
                        self.log.info('does not match fluxMin')
                        mismatch = True
                    if fluxMax != self.fluxMax:
                        self.log.info('does not match fluxMax')
                        mismatch = True
                    if varMin != self.varMin:
                        self.log.info('does not match varMin')
                        mismatch = True
                    if varMax != self.varMax:
                        self.log.info('does not match varMax')
                        mismatch = True
                    if zMin != self.zMin:
                        self.log.info('does not match zMin')
                        mismatch = True
                    if zMax != self.zMax:
                        self.log.info('does not match zMax')
                        mismatch = True
                    if noiseFactor != self.noiseFactor:
                        self.log.info('does not match fluxMin')
                        mismatch = True
                    if priorSigmaBuffer != self.priorSigmaBuffer:
                        self.log.info('does not match priorSigmaBuffer')
                        mismatch = True
                    if priorSigmaStep != self.priorSigmaStep:
                        self.log.info('does not match priorSigmaStep')
                        mismatch = True
                    if priorSigmaCutoff != self.priorSigmaCutoff:
                        self.log.info('does not match priorSigmaCutoff')
                        mismatch = True
                    if nSample != self.nSample:
                        self.log.info('does not match nSample')
                        mismatch = True
                    if selectionOnly != self.selectionOnly:
                        self.log.info('does not match selectionOnly')
                        mismatch = True
                    if invariantCovariance != self.invariantCovariance:
                        self.log.info('does not match invariantCovariance')
                        mismatch = True

                    if mismatch:
                        self.log.info('Skipping %s' % file)
                        continue

                for s in cat:
                    ti = self.bfd.TemplateInfo()
                    ti.m = s.get('m')
                    ti.dm = s.get('dm').reshape(self.bfd.BFDConfig.MSIZE, self.bfd.BFDConfig.DSIZE)
                    ti.dxy = s.get('dxy').reshape(self.bfd.BFDConfig.XYSIZE, self.bfd.BFDConfig.DSIZE)
                    ti.nda = s.get('nda')
                    ti.id = s.get('bfd_id')

                    self.prior.addTemplateInfo(ti)

            except Exception as e:
                print('Failed to read', e)
                continue
        self.prior.prepare()
        self.initialized = True

