import collections
from .utils import dotdict

class InputContainer:
    def __init__(self):
        self.tau_bg = None
        self.tau_th = None
        self.tau_source = None

        self.background = dotdict()
        self.thermo = dotdict()
        self.cosmos = dotdict()
        self.scalars = dotdict()

class TargetContainer:
    def __init__(self):
        self.k_source = None
        self.tau_source = None

        self.sources = collections.OrderedDict()
