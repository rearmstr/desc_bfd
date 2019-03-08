import lsst.pipe.base as pipeBase
from lsst.pipe.base import ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer

from .measureMCShear import MeasureMCShearTask


class ProcessBfdPatchConfig(Config):
    measure = ConfigurableField(target=MeasureMCShearTask, doc="Patch Processing")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


class ProcessBfdPatchTask(pipeBase.CmdLineTask):
    """Process an tract at once.
    """

    ConfigClass = ProcessBfdPatchConfig
    _DefaultName = "processBfdPatch"

    def __init__(self, *args, **kwargs):
        """Constructor.

        All nodes execute this method.
        """
        super(ProcessBfdPatchTask, self).__init__(*args, **kwargs)
        self.makeSubtask('measure')

    @classmethod
    def _makeArgumentParser(cls, *args, **kwargs):
        parser = ArgumentParser(name="processBfdPatch", *args, **kwargs)
        parser.add_id_argument("--id", "deepCoadd", help="data ID, e.g. --id tract=12345 patch=1,2",
                               ContainerClass=CoaddDataIdContainer)
        return parser

    def run(self, patchRef):
        """
        """
        try:
            self.measure.prep()
        except Exception as e:
            self.log.warn("Failed to prepare %s: %s\n" % (patchRef.dataId, e))
            return

        try:
            self.measure.run(patchRef)
        except Exception as e:
            self.log.warn("Failed to process %s: %s\n" % (patchRef.dataId, e))
            import traceback
            traceback.print_exc()
            return

    def write(self, cache, struct, focusMd=None):
        """Write the outputs.
        """

    def _getConfigName(self):
        return None

    def _getEupsVersionsName(self):
        return None
        
    def _getMetadataName(self):
        return None
