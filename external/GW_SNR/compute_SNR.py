from classy import Class
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
_c_ = 2.99792458e8


class Detector():
    detector_type = "CE_plus_ET"

    detector_num = 1
    detector_frequencies = []
    detector_psd = []
    detector_psd_omega = []


class Background():
    H0 = 0
    omega_freq_size = 1000
    frequencies = []
    Omega_cgwb = []


def calculate_background(pba: Background, omega_freq_size=1000, f_min=5.01, f_max=1000):
    M = Class()
    M.set({
        'output': 'OmGW, gwCl',
        'gwb_source_type':  'PBH_gwb',
        'f_pivot':          25,
        'f_min':            f_min,
        'f_max':            f_max,
        'A_star':           2e-5,
        'f_star':           100.,
        'f_NL':             1.,
        'tau_ini_gwb':      0.1,
    })

    M.compute()

    # OmGW = M.get_omega_gw()
    # pba.frequencies = OmGW['f [Hz]']
    # pba.Omega_cgwb = OmGW['Omega_GW(f)']

    pba.frequencies = np.linspace(f_min, f_max, omega_freq_size)
    pba.frequencies = np.geomspace(f_min, f_max, omega_freq_size)
    pba.Omega_cgwb = np.array([M.Omega_GW(f) for f in pba.frequencies])
    pba.omega_freq_size = len(pba.frequencies)

    # plt.figure()
    # plt.loglog(pba.frequencies, pba.Omega_cgwb)

    # pba.H0 = M.h() * 100
    back = M.get_background() 
    pba.H0 = back['H [1/Mpc]'][-1]

    return 0


def read_psds(pdet: Detector, det_filename, in_det):
    data = np.genfromtxt(det_filename).T
    pdet.detector_frequencies[in_det] = data[0]
    pdet.detector_psd[in_det] = data[1]
    return 0


def get_gw_detecotr_psd(pdet: Detector, pba: Background):
    # interpolation_type = 'linear'
    interpolation_type = 'quadratic'
    if (pdet.detector_type == "CE_plus_ET"):
        pdet.detector_num = 5

        pdet.detector_psd           = [[]]*pdet.detector_num
        pdet.detector_frequencies   = [[]]*pdet.detector_num

        read_psds(pdet,"./Detector_PSDs/ce1.txt",0)
        read_psds(pdet,"./Detector_PSDs/ce2.txt",1)
        read_psds(pdet,"./Detector_PSDs/ET.txt",2)
        read_psds(pdet,"./Detector_PSDs/ET.txt",3)
        read_psds(pdet,"./Detector_PSDs/ET.txt",4)

    pdet.detector_psd_omega = np.empty((pdet.detector_num, pba.omega_freq_size))

    for in_det in range(pdet.detector_num):
        f = interp1d(pdet.detector_frequencies[in_det], pdet.detector_psd[in_det], kind=interpolation_type)
        pdet.detector_psd_omega[in_det] = f(pba.frequencies)
        # plt.figure()
        # plt.loglog(pdet.detector_frequencies[in_det], pdet.detector_psd[in_det])
        # plt.loglog(pba.frequencies, pdet.detector_psd_omega[in_det])

    return 0


def compute_gw_SNR(pdet: Detector, pba: Background):
    f1 = 0
    f2 = 0
    dx = 0
    T_obs = 10
    snr = 0
    for in_det1 in range(pdet.detector_num):
        for in_det2 in range(pdet.detector_num):
            f1 = np.power(pba.Omega_cgwb[:-1], 2.) \
                    / np.power(pdet.detector_psd_omega[in_det1][:-1] * pdet.detector_psd_omega[in_det2][:-1], 2.) \
                    / np.power(pba.frequencies[:-1], 6.)
            f2 = np.power(pba.Omega_cgwb[1:], 2.) \
                    / np.power(pdet.detector_psd_omega[in_det1][1:] * pdet.detector_psd_omega[in_det2][1:], 2.) \
                    / np.power(pba.frequencies[1:], 6.)
            dx = (pba.frequencies[1:]-pba.frequencies[:-1])/2.
            
            x = np.power( np.power(_c_/1000*pba.H0*3.24*1e-20, 2.) / (4.*np.pi**2 * np.sqrt(4.*np.pi)), 2.) * (f1+f2) * dx
            snr += np.sum(x)

            # for in_freq in range(pba.omega_freq_size-1):
            #     f1 = np.power(pba.Omega_cgwb[in_freq], 2.) \
            #             / np.power(pdet.detector_psd_omega[in_det1][in_freq] * pdet.detector_psd_omega[in_det2][in_freq], 2.) \
            #             / np.power(pba.frequencies[in_freq], 6.)
            #     f2 = np.power(pba.Omega_cgwb[in_freq+1], 2.) \
            #             / np.power(pdet.detector_psd_omega[in_det1][in_freq+1] * pdet.detector_psd_omega[in_det2][in_freq+1], 2.) \
            #             / np.power(pba.frequencies[in_freq+1], 6.)
            #     dx = (pba.frequencies[in_freq+1]-pba.frequencies[in_freq])/2.
                
            #     snr += np.power( np.power(_c_/1000*pba.H0*3.24*1e-20, 2.) / (4.*np.pi**2 * np.sqrt(4.*np.pi)), 2.) * (f1+f2) * dx
    
    snr = np.sqrt(3.15e7*T_obs*snr)
    
    return snr


pba = Background()

calculate_background(pba)
# calculate_background(pba, 5)

pdet = Detector()

get_gw_detecotr_psd(pdet, pba)

snr = compute_gw_SNR(pdet, pba)
print('SNR = %g' % snr)

plt.show()