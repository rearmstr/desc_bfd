from lsst.desc.bfd.measurePrior import MeasurePriorTask
config.measure.retarget(MeasurePriorTask)

#config.measure.algorithms['bfd.kmoment'].sigma = 2
#config.measure.algorithms['bfd.kmoment'].shift = True
#config.measure.algorithms['bfd.kmoment'].wIndex = 4
#config.measure.algorithms['bfd.kmoment'].reCentroid = False
#config.measure.algorithms['bfd.kmoment'].reCentroidPsf = False
config.measure.sigma = 2
config.measure.wIndex = 4
config.measure.reCentroid = False
config.measure.reCentroidPsf = False
config.measure.priorSigmaStep = 2
config.measure.nSample = 50000
config.measure.sample = 0.2
config.measure.snMin = 5
config.measure.snMax = 40
