import lsst.desc.bfd
config.measure.plugins.names = ["base_PixelFlags",
                                "base_SdssCentroid",
                                "base_SdssShape",
                                "base_GaussianFlux",
                                "base_PsfFlux",
                                "base_CircularApertureFlux",
                                "base_SkyCoord",
                                "base_Variance",	
                                "bfdKMoment"]

config.measure.plugins['bfdKMoment'].sigma = 2
config.measure.plugins['bfdKMoment'].shift = True
config.measure.plugins['bfdKMoment'].wIndex = 4
config.measure.plugins['bfdKMoment'].reCentroid = True
config.measure.plugins['bfdKMoment'].reCentroidPsf = True
#config.measure.plugins['bfdKMoment'].calculateVariance=True
config.measure.plugins['bfdKMoment'].useRecVariance=True
config.measure.plugins['bfdKMoment'].useTableVariance=False
config.measure.algorithms['bfdKMoment'].useNoisePs = True

config.measure.oldHSC = True
config.measure.correlatedNoiseFile = "/tigress/rea3/bfd_update/desc_bfd/config/corr_ps.txt"
config.measure.useCorrelatedNoise = True
config.measure.minVariancePixels = 100




#config.calibrate.measurement.measure.plugins.names |= ["base_Variance"]
#config.charImage.measurement.measure.plugins.names |= ["pythonBfdMoment"]





