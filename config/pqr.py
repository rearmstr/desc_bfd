from lsst.desc.bfd.measurePqr import MeasurePqrTask
config.measure.retarget(MeasurePqrTask)

#root.measure.algorithms['bfd.kmoment'].sigma = 2
#root.measure.algorithms['bfd.kmoment'].shift = True
#root.measure.algorithms['bfd.kmoment'].wIndex = 4
#root.measure.algorithms['bfd.kmoment'].reCentroid = True
#root.measure.algorithms['bfd.kmoment'].reCentroidPsf = True

config.measure.priorTracts=[9813]#,8523]
#root.measure.patches=['4,4','5,5','7,4']

