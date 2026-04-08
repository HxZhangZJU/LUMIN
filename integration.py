# !/usr/bin/env python3
# -*- coding=utf-8 -*-

import numpy as np
import pandas as pd


def integrate_peak_params(peak_key, peak_info):
    """
    Integrate one fitted peak analytically from fitted peak parameters.

    Parameters
    peak_key : str
        Peak identifier, e.g. 'peak0'.
    peak_info : dict
        Peak dictionary returned by fitting, with
        peak_info[peak_key]['peak_params'] =
        [offset, sigma, hwhm, amplitude, frac_gauss].

    Returns
    absolute : float
        Absolute integrated area of the mixed Gaussian/Lorentzian peak.
    """  
    # peak = peak_info[peak_key]
    # data, peak_params = peak
    peak_params = peak_info[peak_key]['peak_params']
    offset, sigma, hwhm, amplitude, frac_gauss = peak_params

    if frac_gauss > 1.0:
        frac_gauss = 1.0
    if frac_gauss < 0.0:
        frac_gauss = 0.0
        
    int_gauss = frac_gauss*(amplitude*np.sqrt(2.0*np.pi*sigma**2.0))
    int_lorentz = (1-frac_gauss)*amplitude*hwhm*np.pi
    absolute = int_gauss + int_lorentz
    return absolute


# def integrate_peak_params_(sigma, hwhm, amplitude, frac_gauss):
#     """

#     sigma: Gaussian sigma
#     hwhm: Lorentzian half-width-at-half-maximum
#     amplitude: height of peak
#     frac_gauss: fraction of peak to be Gaussian (Lorentzian fraction is 1-frac_gauss), 0-1
    
#     """
#     if frac_gauss > 1.0:
#         frac_gauss = 1.0
#     if frac_gauss < 0.0:
#         frac_gauss = 0.0
        
#     int_gauss = frac_gauss*(amplitude*np.sqrt(2.0*np.pi*sigma**2.0))
#     int_lorentz = (1-frac_gauss)*amplitude*hwhm*np.pi
#     absolute = int_gauss + int_lorentz
#     return absolute



###############################################################################

# peak_info = {'peak0':{'peak_y':[],'peak_params':[]},'peak1':{'peak_y':[],'peak_params':[]},.....,'peak n':{}}
# peak_key_ls = ['peak0', 'peak1', 'peak2']

def integrate_sum(peak_key, peak_info, x_ppm, r=None, algorithm='Trapezoidal'):
    """
    Numerically integrate one fitted peak over a local window.

    Parameters
    peak_key : str
        Peak identifier, e.g. 'peak0'.
    peak_info : dict
        Peak dictionary returned by fitting.
    x_ppm : np.ndarray
        Chemical-shift axis (ppm) used to estimate sampling interval.
    r : float, optional
        Half-width of integration window around peak center. If None,
        defaults to 3 * hwhm.
    algorithm : str, optional
        Integration method: 'Riemann' or 'Trapezoidal' (default).

    Returns
    integral_value : float
        Numerically integrated peak area (scaled by 10000).
    """
    # unpacking
    peak = peak_info[peak_key]
    data, peak_params = peak.values()
    offset, sigma, hwhm, amplitude, frac_gauss = peak_params

    # Define integration range
    if r is None:
        r = hwhm * 3
    integration_range = [offset - r, offset + r]

    dx = np.abs(x_ppm[1] - x_ppm[0])

    idx_start = max(int(integration_range[0]), 0)
    idx_end = min(int(integration_range[1]), len(data) - 1)

    if algorithm == None or algorithm =='Trapezoidal':
        integral_value = np.trapz(data[idx_start:idx_end+1], dx=dx)
    elif algorithm =='Riemann':
        integral_value = np.sum(data[idx_start:idx_end]) * dx
    else:
        import warnings
        warnings.warn("Unsupported algorithm '{}'. Using 'Trapezoidal' algorithm instead.".format(algorithm))
        integral_value = np.trapz(data[idx_start:idx_end+1], dx=dx)        

    integral_value = integral_value*10000

    return integral_value

