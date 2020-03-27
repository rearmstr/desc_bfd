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

import numpy as np
import lsst.pex.config as pexConfig
import lsst.pipe.base as pipeBase
import lsst.afw.table as afwTable

from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, CoaddTaskRunner

import lsst.desc.bfd as dbfd


__all__ = ("BuildCovarianceConfig", "BuildCovarianceTask")


class PatchRunner(CoaddTaskRunner):
    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """Get bare butler into Task"""
        kwargs["butler"] = parsedCmd.butler
        return [(parsedCmd.id.refList, kwargs), ]


class BuildCovarianceConfig(pexConfig.Config):
    coaddName = pexConfig.Field(
        doc="coadd name: typically one of deep or goodSeeing",
        dtype=str,
        default="deep",
    )
    bins = pexConfig.Field(
        dtype=int,
        default=5,
        optional=True,
        doc="Split data set int this many bins and compute a covariance matrix for each one"
    )
    label = pexConfig.Field(
        dtype=str,
        default='nb',
        optional=True,
        doc="This will be the labe for each"
    )
    minPct = pexConfig.Field(
        dtype=float,
        default=10,
        optional=True,
        doc="The lower bound of noise variance to consider"
    )
    maxPct = pexConfig.Field(
        dtype=float,
        default=90,
        optional=True,
        doc="The upper bound of noise variance to consider"
    )
    minVal = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="The lower bound of noise variance to consider"
    )
    maxVal = pexConfig.Field(
        dtype=float,
        default=None,
        optional=True,
        doc="The upper bound of noise variance to consider"
    )
    outputFile = pexConfig.Field(
        dtype=str,
        default='./test.fits',
        optional=True,
        doc="Where to store the output"
    )
    momentFile = pexConfig.Field(
        dtype=str,
        default='',
        optional=False,
        doc="Write combined moment file?  If not empty it will be written to a file"
    )
    band = pexConfig.Field(
        dtype=str,
        default='i',
        doc="filter"
    )
    use_mag = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="include magnification"
    )
    use_conc = pexConfig.Field(
        dtype=bool,
        default=False,
        doc="include concentration"
    )


class BuildCovarianceTask(pipeBase.CmdLineTask):
    """
    """
    RunnerClass = PatchRunner
    ConfigClass = BuildCovarianceConfig
    _DefaultName = "buildCovariance"

    def __init__(self, *args, **kwargs):
        """
        """
        super(BuildCovarianceTask, self).__init__(*args, **kwargs)

        self.ncolors = 0
        self.bfd = dbfd.BFDConfig(use_conc=self.config.use_conc, use_mag=self.config.use_mag,
                                  ncolors=self.ncolors)
        self.n_even = self.bfd.BFDConfig.MSIZE
        self.n_odd = self.bfd.BFDConfig.XYSIZE
        self.size_even = self.n_even*(self.n_even+1)//2
        self.size_odd = self.n_odd*(self.n_odd+1)//2

        self.schema = afwTable.Schema()
        self.labelKey = self.schema.addField("label", type=str, doc="name of bin", size=10)
        self.minKey = self.schema.addField("min", type=float, doc="minimum value of the variance")
        self.maxKey = self.schema.addField("max", type=float, doc="maximum value of the variance")
        self.isoCovEvenKey = self.schema.addField("isoCovEven", doc="isotropized moment covariance matrix",
                                                  type="ArrayF", size=self.size_even)
        self.isoCovOddKey = self.schema.addField("isoCovOdd", doc="isotropized moment covariance matrix",
                                                 type="ArrayF", size=self.size_odd)
        self.covEvenKey = self.schema.addField("covEven", doc="moment covariance matrix", type="ArrayF",
                                               size=self.size_even)
        self.covOddKey = self.schema.addField("covOdd", doc="moment covariance matrix", type="ArrayF",
                                              size=self.size_odd)
        self.catalog = afwTable.BaseCatalog(self.schema)

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = pipeBase.ArgumentParser(name="buildCovariance", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def runDataRef(self, tractPatchRefList, butler):
        """
        """

        noiseVariance = []
        cov_even = []
        cov_odd = []
        for dataRef in tractPatchRefList:
            self.log.info('adding %s' % dataRef.dataId)
            try:
                cat = dataRef.get('deepCoadd_moments')
                maskGood = cat.get('bfd_flag') == False
                noiseVariance.extend(cat.get(f'bfd_cov_even_{self.config.band}')[:, 0][maskGood])
                cov_even.extend(cat.get(f'bfd_cov_even_{self.config.band}')[maskGood])
                cov_odd.extend(cat.get(f'bfd_cov_odd_{self.config.band}')[maskGood])
            except Exception as e:
                print('problem', e)
                continue

        cov_odd = np.array(cov_odd)
        cov_even = np.array(cov_even)
        noiseVariance = np.array(noiseVariance)
        if self.config.minVal is None:
            minVariance = np.percentile(noiseVariance, self.config.minPct)
        else:
            minVariance = self.config.minVal

        if self.config.maxVal is None:
            maxVariance = np.percentile(noiseVariance, self.config.maxPct)
        else:
            maxVariance = self.config.maxVal

        self.log.info(f'Using range: {minVariance},{maxVariance}')
        binWidth = (maxVariance - minVariance) / self.config.bins
        bins = np.arange(minVariance, maxVariance + binWidth, binWidth)

        for i in range(len(bins)-1):
            rec = self.catalog.addNew()
            mask = np.logical_and(noiseVariance >= bins[i], noiseVariance < bins[i+1])

            tot_even = np.zeros(self.size_even)
            for val in cov_even[mask]:
                tot_even += val
            tot_even /= np.sum(mask)
            rec.set('covEven', np.array(tot_even, dtype=np.float32))

            tot_odd = np.zeros(self.size_odd)
            for val in cov_odd[mask]:
                tot_odd += val
            tot_odd /= np.sum(mask)
            rec.set('covOdd', np.array(tot_odd, dtype=np.float32))
            #
            # isotropize matrix

            varE = 0.5*(tot_even[7] + tot_even[9])
            tot_even[2] = 0
            tot_even[3] = 0
            tot_even[5] = 0
            tot_even[6] = 0
            tot_even[8] = 0
            tot_even[7] = varE
            tot_even[9] = varE

            varX = 0.5*(tot_odd[0] + tot_odd[2])
            tot_odd[1] = 0
            tot_odd[0] = varX
            tot_odd[2] = varX

            rec.set('max', bins[i+1])
            rec.set('min', bins[i])
            rec.set('isoCovEven', np.array(tot_even, dtype=np.float32))
            rec.set('isoCovOdd', np.array(tot_odd, dtype=np.float32))

            rec.set('label', '%s%d' % (self.config.label, i))

        self.catalog.writeFits(self.config.outputFile)

    def write(self, cache, struct, focusMd=None):
        """Write the outputs.
        """

    def _getConfigName(self):
        return None

    def _getEupsVersionsName(self):
        return None

    def _getMetadataName(self):
        return None
