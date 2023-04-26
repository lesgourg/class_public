"""
Compute the SNR for the GWB energy density Omega_GW
author:: Florian Schulze <florian.tobias.schulze@rwth-aachen.de>

This script provides functions to calculate the SNR for the detection of a
CGWB energy density Omega_GW, given a detctor network.
For an example see in main().

You can import this script into your project using:
import os
import sys
sys.path.append(os.path.abspath('path/to/class/external/GW_SNR'))
import compute_GW_SNR

"""
import os
import numpy as np
from scipy.interpolate import interp1d

_folder_ = os.path.dirname(__file__)

def read_psds(det_filename, freqs):
    """Reads a PSD file and caluclates the PSD at the requested frequencies freqs

    Args:
        det_filename (str): detector PSD file
        freqs (array): requested frequencies for the PSD

    Returns:
        array: PSD of the dector at the frequencies freqs
    """

    data = np.genfromtxt(det_filename).T
    detector_frequencies = data[0]
    detector_psd = data[1]
    f = interp1d(detector_frequencies, detector_psd, kind='quadratic', fill_value='extrapolate')

    return f(freqs)


def get_gw_detector_psd(detector, freqs):
    """Calulates the PSD for a detector network

    Args:
        detector (str or list): Detector name or list with the PSD files of the detector (e.g. "CE+ET")
        freqs (array): requested frequencies for the PSD

    Returns:
        nd_array: detecor_psd
    """
    detector_num = 1
    psd_files = []

    if type(detector) == list:
        psd_files = detector
        detector_num = len(psd_files)

    elif detector == "LIGO":
        detector_num = 2
        psd_files = [_folder_+"/detector_PSD/aLIGO.txt", _folder_+"/detector_PSD/aLIGO.txt"]

    elif detector == "CE":
        detector_num = 2
        psd_files = [_folder_+"/detector_PSD/ce1.txt", _folder_+"/detector_PSD/ce2.txt"]

    elif detector == "ET":
        detector_num = 3
        psd_files = [_folder_+"/detector_PSD/ET.txt", _folder_+"/detector_PSD/ET.txt", _folder_+"/detector_PSD/ET.txt"]

    elif detector == "CE+ET":
        detector_num = 5
        psd_files = [_folder_+"/detector_PSD/ce1.txt", _folder_+"/detector_PSD/ce2.txt", _folder_+"/detector_PSD/ET.txt", _folder_+"/detector_PSD/ET.txt", _folder_+"/detector_PSD/ET.txt"]

    else:
        print("You entered an unknwon detector: %s." % detector)
        detector_num = 0

    detector_psd = np.empty((detector_num, len(freqs)))
    
    for in_det in range(detector_num):
        detector_psd[in_det] = read_psds(psd_files[in_det], freqs)
    
    return detector_psd


def compute_GW_SNR(freqs, Omega_GW, detector_psd, h=0.67, T_obs=10):
    """Calculates the SNR for the CGWB monopole Omega_GW given an detector detector_psd.

    Args:
        freqs (array): frequencies in Hz
        Omega_GW (array): CGWB energy density $\Omega_{GW}$
        detector_psd (nd_array): detector PSD
        h (int, optional): reduced Hubble rate h. Defaults to 0.67.
        T_obs (int, optional): observation time in years. Defaults to 10.

    Returns:
        float: Signal to Noise Ratio (SNR)
    """
    detector_num = len(detector_psd)
    snr = 0
    
    for in_det1 in range(detector_num):
        for in_det2 in range(detector_num):
            f1 = np.power(Omega_GW[:-1], 2.) \
                    / np.power(detector_psd[in_det1][:-1] * detector_psd[in_det2][:-1], 2.) \
                    / np.power(freqs[:-1], 6.)
            f2 = np.power(Omega_GW[1:], 2.) \
                    / np.power(detector_psd[in_det1][1:] * detector_psd[in_det2][1:], 2.) \
                    / np.power(freqs[1:], 6.)
            dx = (freqs[1:]-freqs[:-1])/2.
            
            x = np.power( np.power(h*3.24*1e-18, 2.) / (4.*np.pi**2 * np.sqrt(4.*np.pi)), 2.) * (f1+f2) * dx
            snr += np.sum(x)

    snr = np.sqrt(3.15e7*T_obs*snr)
    
    return snr


def main():
    from classy import Class

    M = Class()
    M.set({
        'output': 'OmGW, gwCl',
        'f_pivot': 25,
        'f_min':   5.01,
        'f_max':   1000,

        # 'gwb_source_type':  'analytic_gwb',
        # 'Omega_gwb':        3e-10,
        # 'n_gwb':            0.0,

        'gwb_source_type':  'PBH_gwb',
        'A_star':           2e-5,
        'f_star':           100.,
        'f_NL':             1.,
    })

    M.compute()

    OmGW = M.get_omega_gw()
    freqs = OmGW['f [Hz]']
    Omega_GW = OmGW['Omega_GW(f)']
    h = M.h()

    detecor_psd = get_gw_detector_psd("LIGO", freqs)
    #Alternanative:
    # detecor_psd = get_gw_detector_psd(["./detector_PSD/aLIGO.txt", "./detector_PSD/aLIGO.txt"], freqs)
    snr = compute_GW_SNR(freqs, Omega_GW, detecor_psd, h)
    print('SNR for LIGO:\t %g' % snr)

    detecor_psd = get_gw_detector_psd("CE", freqs)
    snr = compute_GW_SNR(freqs, Omega_GW, detecor_psd, h)
    print('SNR for CE:\t %g' % snr)

    detecor_psd = get_gw_detector_psd("ET", freqs)
    snr = compute_GW_SNR(freqs, Omega_GW, detecor_psd, h)
    print('SNR for ET:\t %g' % snr)

    detecor_psd = get_gw_detector_psd("CE+ET", freqs)
    snr = compute_GW_SNR(freqs, Omega_GW, detecor_psd, h)
    print('SNR for CE+ET:\t %g' % snr)

    return 0


if __name__ == '__main__':
    main()
