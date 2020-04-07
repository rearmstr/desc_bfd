#!/Users/armstrong46/opt/anaconda3/envs/mystack/bin/python  # noqa
from lsst.desc.bfd.processBfdCoadd import ProcessBfdCoaddTask
from lsst.desc.bfd.measurePrior import MeasurePriorTask

ProcessBfdCoaddTask.measure.retarget(MeasurePriorTask)
ProcessBfdCoaddTask.parseAndSubmit()

