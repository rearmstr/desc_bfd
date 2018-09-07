import lsst.desc.bfd
config.calibrate.measurement.plugins.names |= ["pythonBfdMoment", "bfdKMoment"]

config.calibrate.measurement.plugins['bfdKMoment'].sigma = 2
config.calibrate.measurement.plugins['bfdKMoment'].shift = True
config.calibrate.measurement.plugins['bfdKMoment'].wIndex = 4
config.calibrate.measurement.plugins['bfdKMoment'].reCentroid = True
config.calibrate.measurement.plugins['bfdKMoment'].reCentroidPsf = True
config.calibrate.measurement.plugins['bfdKMoment'].calculateVariance=True
config.calibrate.measurement.plugins['bfdKMoment'].useRecVariance=False
config.calibrate.measurement.plugins['bfdKMoment'].useTableVariance=False
#root.measure.algorithms['bfdKMoment'].useNoisePs = True



#config.calibrate.measurement.plugins.names |= ["base_Variance"]
#config.charImage.measurement.plugins.names |= ["pythonBfdMoment"]





