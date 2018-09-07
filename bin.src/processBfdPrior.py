#!/usr/bin/env python
from lsst.desc.bfd.processBfdCoadd import ProcessBfdCoaddTask
from lsst.desc.bfd.measurePrior import MeasurePriorTask

ProcessBfdCoaddTask.measure.retarget(MeasurePriorTask)
ProcessBfdCoaddTask.parseAndSubmit()

