from abc import ABC, abstractmethod
import numpy as np

class ParameterDomain(ABC):
    @abstractmethod
    def contains(self, parameters):
        """
        Given a Dict[str, float] of parameters, decide whether these parameters
        lie inside the parameter domain or not.
        """
        pass


class CubeDomain(ParameterDomain):

    def __init__(self, domain):
        """
        `domain` is a Dict[str, (float, float)].
        """
        self.domain = domain

    def contains(self, parameters):
        return all(self.domain[k][0] <= v <= self.domain[k][1] for k, v in parameters.items())


class EllipsoidDomain(ParameterDomain):

    # TODO don't hardcode
    # delta_chi2 of 5 for df = 9
    def __init__(self, fields, bestfit, covmat, delta_chi2=41.70122440330078):
        """
        `fields` is a list of strings giving the names (and their order) of the
        rows & columns of the `covmat` and of `bestfit`.
        """
        self.fields = fields
        self.bestfit = bestfit
        self.covmat = covmat
        self.inv_covmat = np.linalg.inv(covmat)
        self.delta_chi2 = delta_chi2

    def contains(self, parameters):
        pvec = np.array([parameters[f] for f in self.fields])
        d = pvec - self.bestfit
        return d.dot(self.inv_covmat.dot(d)) < self.delta_chi2



