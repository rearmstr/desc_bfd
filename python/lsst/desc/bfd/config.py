import lsst.desc.bfd

class BFDConfig(object):
        """
        """
        def __init__(self, fix_center=False, use_conc=False, use_mag=False, ncolors=0, use_float=True):
                label = ""
                if fix_center:
                        label += "1"
                else:
                        label += "0"

                if use_conc:
                        label += "1"
                else:
                        label += "0"

                if use_mag:
                        label += "1"
                else:
                        label += "0"

                label += str(ncolors)
                if use_float:
                        label += "F"
                else:
                        label += "D"

                all_methods = dir(lsst.desc.bfd)
                for method in all_methods:
                        if method.endswith(label):
                                setattr(self, method.rsplit(label)[0], getattr(lsst.desc.bfd, method))