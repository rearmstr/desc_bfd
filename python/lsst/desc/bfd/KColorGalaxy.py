from .config import BFDConfig
from scipy.optimize import fsolve
import numpy as np


class KColorGalaxy:
    """Class to hold multi color information for several KGalaxies
    """

    def __init__(self, bfd_config, kgals, weights=None, nda=1., id=0):
        """
        bfd_config = BFD configuration
        id = galaxy id
        nda = sky density
        """
        self.bfd_config = bfd_config.BFDConfig
        self.bfd_build = bfd_config
        self.id = id
        self.nda = nda
        self.kgals = kgals
        if weights is None:
            self.weights = [1./len(kgals)]*len(kgals)
        else:
            self.weights = weights

    def get_moment(self, dx=0., dy=0., fill_cov=False):
        '''Return weighted moment vector with coordinate origin at dx, dy
         '''

        fluxes = np.zeros(len(self.kgals))
        single_cov = []
        for i, kgal in enumerate(self.kgals):
            shifted = kgal.getShifted(dx, dy)
            target = shifted.getTarget(True)
            fluxes[i] = target.mom.m[0]
            if fill_cov:
                single_cov.append(target.cov.m)

        even = np.zeros(self.bfd_config.MSIZE)
        odd = np.zeros(self.bfd_config.XYSIZE)
        even_cov = np.zeros((self.bfd_config.MSIZE, self.bfd_config.MSIZE))
        odd_cov = np.zeros((self.bfd_config.XYSIZE, self.bfd_config.XYSIZE))

        for kgal, weight in zip(self.kgals, self.weights):
            shifted = kgal.getShifted(dx, dy)
            target = shifted.getTarget(True)
            even += target.mom.m*weight
            odd += target.mom.xy*weight
            if fill_cov:
                even_cov += target.cov.m*weight**2
                odd_cov += target.cov.xy*weight**2

        # Replace fluxes with appropriate values from single band moments
        even[0] = fluxes[0]
        even[self.bfd_config.MC0:] = fluxes[1:]

        if fill_cov:
            # Flux covariances are same as the single band
            even_cov[0, 0] = single_cov[0][0, 0]
            for i in range(1, len(fluxes)):
                even_cov[self.bfd_config.MC0 + i - 1,
                         self.bfd_config.MC0 + i - 1] = single_cov[i][0, 0]

            # Set appropriate covariances between combined moments and fluxes
            for i in range(self.bfd_config.MR, self.bfd_config.MC0):
                even_cov[0, i] = self.weights[0]*single_cov[0][0, i]
                even_cov[0, i] = even_cov[i, 0]

            for j in range(1, len(fluxes)):
                for i in range(self.bfd_config.MR, self.bfd_config.MC0):
                    even_cov[self.bfd_config.MC0 + j - 1, i] = self.weights[j]*single_cov[j][0, i]
                    even_cov[i, self.bfd_config.MC0 + j - 1] = even_cov[self.bfd_config.MC0 + j - 1, i]

        moment = self.bfd_build.Moment(even, odd)
        moment_cov = self.bfd_build.MomentCov(even_cov, odd_cov)

        return moment, moment_cov

    def get_single_moment(self, index, dx=0., dy=0.):
        '''Return moments for single filter fiven by index
         '''
        msize = self.bfd_config.MSIZE - len(self.kgals) + 1
        xysize = self.bfd_config.XYSIZE
        shifted = self.kgals[index].getShifted(dx, dy)
        target = shifted.getTarget(True)

        even = target.mom.m
        odd = target.mom.xy
        even_cov = target.cov.m
        odd_cov = target.cov.xy

        return even[:msize], odd[:xysize], even_cov[:msize, :msize], odd_cov[:xysize, :xysize]

    def get_template(self, dx=0., dy=0.):

        t = self.bfd_build.TemplateGalaxy()
        even_out = np.zeros((t.MSIZE, t.DSIZE), dtype=np.complex64)
        odd_out = np.zeros((t.XYSIZE, t.DSIZE), dtype=np.complex64)

        single_templates = []
        for i, kgal in enumerate(self.kgals):
            shifted = kgal.getShifted(dx, dy)
            template = shifted.getTemplate()
            single_templates.append(template)

        for kgal, weight in zip(self.kgals, self.weights):
            shifted = kgal.getShifted(dx, dy)
            template = shifted.getTemplate()
            even_out += template.mDeriv
            odd_out += template.xyDeriv

        even_out[0, :] = single_templates[0].mDeriv[0, :]

        for j in range(1, t.MSIZE-t.MC0):
            even_out[t.MC0 + j, :] = single_templates[j].mDeriv[0, :]

        t = self.bfd_build.TemplateGalaxy(even_out, odd_out, self.nda, self.id)
        return t

    def xy_moment(self, dx):
        '''Interface to fsolve to return x,y moments given input origin shift
        '''
        mom, cov = self.get_moment(dx[0], dx[1])
        return mom.xy

    def xy_jacobian(self, dx):
        '''Function to return Jacobian of X & Y moments with respect to
        origin shift dx, for use in solver.  Use the fact that posn derivatives
        of the first moments are the second moments.
        '''
        mom, cov = self.get_moment(dx[0], dx[1])
        e = mom.m
        return -0.5 * np.array([[e[self.bfd_config.MR] + e[self.bfd_config.M1],
                                 e[self.bfd_config.M2]],
                                [e[self.bfd_config.M2],
                                 e[self.bfd_config.MR] - e[self.bfd_config.M1]]])

    def recenter(self, sigma):
        '''Find dx, dy that null the X and Y moments

        Returns:
        dx    Shift of origin applied
        error True if there was a failure to converge
        msg   message string on failure
        '''
        dx = np.zeros(2, dtype=float)
        dx, junk, ier, msg = fsolve(self.xy_moment, dx, fprime=self.xy_jacobian, full_output=True)

        threshold = np.sqrt(2.0) * np.sqrt(sigma**2)
        wandered_too_far = np.abs(dx) >= threshold
        badcentering = wandered_too_far[0] or wandered_too_far[1] or ier <= 0
        return dx, badcentering, msg
