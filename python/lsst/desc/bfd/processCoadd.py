from __future__ import absolute_import, division, print_function
import os
from argparse import ArgumentError

from builtins import zip

from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.base import ArgumentParser, TaskRunner
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.ctrl.pool.pool import Pool, abortOnError, NODE
from lsst.meas.base.forcedPhotCoadd import ForcedPhotCoaddTask
from lsst.pipe.drivers.utils import getDataRef, TractDataIdContainer
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer
import lsst.afw.table as afwTable
import lsst.afw.geom as afwGeom
import argparse
from .measureCoadd import MeasureCoaddTask
import numpy as np
import scipy.spatial

from lsst.pipe.tasks.coaddBase import getSkyInfo

class AssociationDataIdContainer(CoaddDataIdContainer):

    def makeDataRefList(self, namespace):
        """Make self.refList from self.idList

        It's difficult to make a data reference that merely points to an entire
        tract: there is no data product solely at the tract level.  Instead, we
        generate a list of data references for patches within the tract.
        """
        datasetType = namespace.config.coaddName + "Coadd_calexp"
        validKeys = set(["tract", "filter", "patch"])

        def getPatchRefList(tract):
            return [namespace.butler.dataRef(datasetType=datasetType,
                                             tract=tract.getId(),
                                             filter=dataId["filter"],
                                             patch="%d,%d" % patch.getIndex())
                    for patch in tract]

        tractRefs = {}  # Data references for each tract
        for dataId in self.idList:
            for key in validKeys:
                if key in ("tract", "patch",):
                    # Will deal with these explicitly
                    continue
                if key not in dataId:
                    raise argparse.ArgumentError(
                        None, "--id must include " + key)

            skymap = self.getSkymap(namespace)

            if "tract" in dataId:
                tractId = dataId["tract"]
                if tractId not in tractRefs:
                    tractRefs[tractId] = []
                if "patch" in dataId:
                    tractRefs[tractId].append(namespace.butler.dataRef(datasetType=datasetType, tract=tractId,
                                                                       filter=dataId[
                                                                           'filter'],
                                                                       patch=dataId['patch']))
                else:
                    tractRefs[tractId] += getPatchRefList(skymap[tractId])
            else:
                tractRefs = dict((tract.getId(), tractRefs.get(tract.getId(), []) +
                                  getPatchRefList(tract)) for tract in skymap)

        self.refList = list(tractRefs.values())


class AssociationDriverConfig(Config):
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")
    measure = ConfigurableField(target=MeasureCoaddTask, doc="CCD processing task")


class AssociationDriverTaskRunner(TaskRunner):
    """TaskRunner for running MultiBandTask

    This is similar to the lsst.pipe.base.ButlerInitializedTaskRunner,
    except that we have a list of data references instead of a single
    data reference being passed to the Task.run, and we pass the results
    of the '--reuse-outputs-from' command option to the Task constructor.
    """

    def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
        TaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

    def makeTask(self, parsedCmd=None, args=None):
        """A variant of the base version that passes a butler argument to the task's constructor
        parsedCmd or args must be specified.
        """
        if parsedCmd is not None:
            butler = parsedCmd.butler
        elif args is not None:
            dataRefList, kwargs = args
            butler = dataRefList[0].butlerSubset.butler
        else:
            raise RuntimeError("parsedCmd or args must be specified")
        return self.TaskClass(config=self.config, log=self.log, butler=butler)


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)


class AssociationDriverTask(BatchPoolTask):
    """Multi-node driver for multiband processing"""
    ConfigClass = AssociationDriverConfig
    _DefaultName = "associationDriver"
    RunnerClass = AssociationDriverTaskRunner

    def __init__(self, butler=None, **kwargs):
        """!
        @param[in] butler: the butler can be used to retrieve schema or passed to the refObjLoader constructor
            in case it is needed.
        @param[in] schema: the schema of the source detection catalog used as input.
        @param[in] refObjLoader: an instance of LoadReferenceObjectsTasks that supplies an external reference
            catalog.  May be None if the butler argument is provided or all steps requiring a reference
            catalog are disabled.
        """
        BatchPoolTask.__init__(self, **kwargs)
        self.butler = butler
        self.makeSubtask("measure")


    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   butler=self.butler))

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2 filter=g^r^i",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numCpus):
        """!Return walltime request for batch job

        @param time: Requested time per iteration
        @param parsedCmd: Results of argument parsing
        @param numCores: Number of cores
        """
        return time

    @abortOnError
    def runDataRef(self, patchRefList):
        """!Run multiband processing on coadds

        Only the master node runs this method.

        No real MPI communication (scatter/gather) takes place: all I/O goes
        through the disk. We want the intermediate stages on disk, and the
        component Tasks are implemented around this, so we just follow suit.

        @param patchRefList:  Data references to run measurement
        """
        print (len(patchRefList))
        for patchRef in patchRefList:
            if patchRef:
                butler = patchRef.getButler()
                break
        else:
            raise RuntimeError("No valid patches")
        pool = Pool("all")
        pool.cacheClear()
        pool.storeSet(butler=butler)

        # Group by patch
        patches = {}
        tract = None
        for patchRef in patchRefList:
            dataId = patchRef.dataId
            if tract is None:
                tract = dataId["tract"]
            else:
                assert tract == dataId["tract"]

            patch = dataId["patch"]
            if patch not in patches:
                patches[patch] = []
            patches[patch].append(dataId)

        print (patches.values())
        dataRefList = [getDataRef(cache.butler, dataId, self.config.coaddName + "Coadd_calexp") for
                       dataId in patches.values()]
        pool.map(self.runAssociation, dataRefList)


    def runAssociation(self, cache, patchRef):
        """! Run detection on a patch"""


        with self.logOperation("processing %s " % (patchRef.dataId)):
            try:
                result = self.measure.run(patchRef)
            except Exception as e:
                self.log.warn("Failed to process %s: %s\n" % (patchRef.dataId, e))
                import traceback
                traceback.print_exc()
                return None


    def writeMetadata(self, dataRef):
        """We don't collect any metadata, so skip"""
        pass
