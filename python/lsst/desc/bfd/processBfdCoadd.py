
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
from lsst.pipe.base import Struct, ArgumentParser
from lsst.ctrl.pool.parallel import BatchPoolTask
from lsst.pipe.drivers.coaddDriver import CoaddDriverTaskRunner
from lsst.pipe.tasks.coaddBase import CoaddTaskRunner
from lsst.pipe.drivers.utils import getDataRef, TractDataIdContainer
from lsst.ctrl.pool.pool import Pool, abortOnError, NODE
from lsst.pex.config import Config, Field, ConfigurableField, ListField
from .measureCoadd import MeasureCoaddTask
from .processBfdPatch import PatchRunner



class ProcessBfdCoaddConfig(Config):
    measure = ConfigurableField(target=MeasureCoaddTask, doc="CCD processing task")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


class DriverTaskRunner(CoaddTaskRunner):

    def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
        CoaddTaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

    def makeTask(self, parsedCmd=None, args=None):
        return self.TaskClass(config=self.config, log=self.log)

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """!Get bare butler into Task
        @param parsedCmd results of parsing command input
        """
        kwargs["butler"] = parsedCmd.butler
        return [(parsedCmd.id.refList, kwargs), ]


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)

class ProcessBfdCoaddTask(BatchPoolTask):
    ConfigClass = ProcessBfdCoaddConfig
    _DefaultName = "processBfdCoadd"
    RunnerClass = DriverTaskRunner

    def __init__(self, **kwargs):
        BatchPoolTask.__init__(self, **kwargs)

        self.makeSubtask('measure')

    def __reduce__(self):
        """Pickler"""
        return unpickle, (self.__class__, [], dict(config=self.config, name=self._name,
                                                   parentTask=self._parentTask, log=self.log,
                                                   ))
    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name=cls._DefaultName, *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=TractDataIdContainer)
        return parser

    @classmethod
    def batchWallTime(cls, time, parsedCmd, numCores):
        return time
#    def run(self, tractPatchRefList, bmeautler, selectIdList=[]):
#        print "finished"
#    def runDataRef(self, tractPatchRefList, butler, selectIdList=[]):
#        """Determine which tracts are non-empty before processing"""
#        self.log.info('enter run')
#        pool = Pool("tracts")
#
#        return [self.runTract(patchRefList, butler) for patchRefList in tractPatchRefList]
#
    @abortOnError
    def runDataRef(self, tractPatchRefList, butler):
        """!Determine which tracts are non-empty before processing
        @param tractPatchRefList: List of tracts and patches to include in the coaddition
        @param butler: butler reference object
        @param selectIdList: List of data Ids (i.e. visit, ccd) to consider when making the coadd
        @return list of references to sel.runTract function evaluation for each tractPatchRefList member
        """
        pool = Pool("tracts")
        pool.storeSet(butler=butler, skymap=butler.get(self.config.coaddName + "Coadd_skyMap"))

        return [self.runTract(patchRefList, butler) for patchRefList in tractPatchRefList]

    def runTract(self, patchRefList, butler):

        pool = Pool("patches")
        pool.map(self.runBfd, patchRefList)

    def runBfd(self, cache, patchRef):

        with self.logOperation("processing %s " % (patchRef.dataId)):
            try:
                result = self.measure.run(patchRef)
            except Exception as e:
                self.log.warn("Failed to process %s: %s\n" % (patchRef.dataId, e))
                import traceback
                traceback.print_exc()
                return None

    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
