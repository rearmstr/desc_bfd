import sys
import math
import collections

import lsst.daf.base as dafBase
import lsst.afw.table as afwTable
import lsst.afw.image as afwImage
import lsst.afw.cameraGeom as afwCg
from lsst.pipe.base import Struct, ArgumentParser
from lsst.pex.config import Config, Field, ConfigurableField, ListField
#from hsc.pipe.base.pool import abortOnError, NODE, Pool, Debugger
#from hsc.pipe.base.parallel import BatchPoolTask
from lsst.pipe.tasks.coaddBase import CoaddDataIdContainer, SelectDataIdContainer, CoaddTaskRunner
import lsst.pipe.base
from .measureCoadd import MeasureCoaddTask


from lsst.pipe.drivers.utils import TractDataIdContainer


#Debugger().enabled = True


class PatchRunner(CoaddTaskRunner):

    def __init__(self, TaskClass, parsedCmd, doReturnResults=False):
        CoaddTaskRunner.__init__(self, TaskClass, parsedCmd, doReturnResults)

    #def makeTask(self, parsedCmd=None, args=None):
    #    return self.TaskClass(schema=None, config=self.config, log=self.log)

    @staticmethod
    def getTargetList(parsedCmd, **kwargs):
        """!Get bare butler into Task
        @param parsedCmd results of parsing command input
        """
        kwargs["butler"] = parsedCmd.butler
        return [(parsedCmd.id.refList, kwargs), ]

# class PatchRunner(CoaddTaskRunner):
#     @staticmethod
#     def getTargetList(parsedCmd, **kwargs):
#         """Get bare butler into Task"""
#         kwargs["butler"] = parsedCmd.butler
#         return [(parsedCmd.id.refList, kwargs),]
# 

class ProcessBfdPatchConfig(Config):
    measure = ConfigurableField(target=MeasureCoaddTask, doc="Patch Processing")
    coaddName = Field(dtype=str, default="deep", doc="Name of coadd")


class ProcessBfdPatchTask(lsst.pipe.base.CmdLineTask):
    """Process an tract at once.
    """

    RunnerClass = PatchRunner
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
                               ContainerClass=TractDataIdContainer)
        return parser

    def run(self, tractPatchRefList, butler):
        pass
    #@abortOnError
    def runDataRef(self, tractPatchRefList, butler):
        """
        """
        print ('running prep')
        try:
            self.measure.prep()
        except Exception as e:
            print (e)
            return
        tractIdList = []
        for tract in tractPatchRefList:
            for patchRef in tract:
                result = self.process(patchRef)


    def process(self, dataRef):
        """Process a single Patch
        """
        try:
            result = self.measure.run(dataRef)
        except Exception as e:
            self.log.warn("Failed to process %s: %s\n" % (dataRef.dataId, e))
            import traceback
            traceback.print_exc()
            return None

        return result


    def write(self, cache, struct, focusMd=None):
        """Write the outputs.
        """

    def _getConfigName(self):
        return None
    def _getEupsVersionsName(self):
        return None
    def _getMetadataName(self):
        return None
