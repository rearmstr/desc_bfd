from lsst.pipe.base import ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.ctrl.pool.parallel import BatchParallelTask, BatchTaskRunner
from lsst.coadd.utils.coaddDataIdContainer import CoaddDataIdContainer

from .processBfdPatch import ProcessBfdPatchTask


class ProcessBfdCoaddConfig(Config):
    processPatch = ConfigurableField(target=ProcessBfdPatchTask, doc="Patch processing task")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


def unpickle(factory, args, kwargs):
    """Unpickle something by calling a factory"""
    return factory(*args, **kwargs)


class ProcessBfdCoaddTask(BatchParallelTask):
    """Process Patches in parallel
    """
    ConfigClass = ProcessBfdCoaddConfig
    _DefaultName = "processBfdCoadd"
    RunnerClass = BatchTaskRunner

    def __init__(self, butler=None, *args, **kwargs):
        """!
        Constructor
        """
        BatchParallelTask.__init__(self, *args, **kwargs)
        self.makeSubtask("processPatch")

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        kwargs.pop("doBatch", False)
        parser = ArgumentParser(name="processBfdCoadd", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def run(self, patchRef):
        """Process a patch, with scatter-gather-scatter using MPI.
        """
        with self.logOperation("processing %s" % (patchRef.dataId,)):
            self.processPatch.run(patchRef)

    def _getConfigName(self):
        return None

    def _getEupsVersionsName(self):
        return None

    def _getMetadataName(self):
        return None
