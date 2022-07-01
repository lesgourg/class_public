from .T0_reco import Net_ST0_Reco
from .T0_reio import Net_ST0_Reio
from .T0_isw import Net_ST0_ISW
from .T1 import Net_ST1
from .T2_reco import Net_ST2_Reco
from .T2_reio import Net_ST2_Reio
from .phi_plus_psi import Net_phi_plus_psi

ALL_NETWORK_CLASSES = [
    Net_ST0_Reco,
    Net_ST0_Reio,
    Net_ST0_ISW,
    Net_ST1,
    Net_ST2_Reco,
    Net_ST2_Reio,
    Net_phi_plus_psi,
]

ALL_NETWORK_STRINGS = [
    'Net_ST0_Reco',
    'Net_ST0_Reio',
    'Net_ST0_ISW',
    'Net_ST1',
    'Net_ST2_Reco',
    'Net_ST2_Reio',
    'Net_phi_plus_psi',
]
