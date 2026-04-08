# !/usr/bin/env python3
# -*- coding=utf-8 -*-
import numpy as np
# {'peak0':{'peak_y':[],'peak_params':[]},'peak1':{'peak_y':[],'peak_params':[]},.....,'peak n': {},'residual':0.89,'residue_curve':[]}
# =>
#           [offset, sigma, hwhm, amplitude, frac_gauss]

def line_fitting(x_ppm, y_a, peaks, frac_gauss=None):
    """
    Fit NMR peaks with mixed Gaussian/Lorentzian components.

    Parameters
    x_ppm : np.ndarray
        Chemical-shift axis in ppm.
    y_a : np.ndarray
        Signal intensity values.
    peaks : list[int]
        Peak index list used as initial centers.
    frac_gauss : float, optional
        Initial Gaussian fraction in [0, 1].

    Returns
    fitting_result : dict or None
        Dictionary containing per-peak fitted parameters, reconstructed
        peak curves, total residual, and residue curve. Returns None if fitting fails.
    """
    import numpy as np
    import lmfit
    
    fitting_result = {}

    if frac_gauss is None:
        frac_gauss = 0.0

    p = _f_makep(y_a, peaks, frac_gauss)

    params = lmfit.Parameters()
    for i, peak_param in enumerate(p):
        for j, param in enumerate(peak_param):
            if j == 0:
                params.add(f'peak_{i}_{j}', value=param, vary=True, min=0.0, max=len(y_a) - 1)
            elif j == 3:
                params.add(f'peak_{i}_{j}', value=param, vary=True, min=0.0, max=np.max(y_a)*2)
            elif j == 4:  # frac_gauss
                # params.add(f'peak_{i}_{j}', value=param, vary=(frac_gauss is None), min=0, max=1)
                params.add(f'peak_{i}_{j}', value=param, vary=True, min=0.0, max=1)
            else:
                params.add(f'peak_{i}_{j}', value=param, vary=True, min=0.0)

    try:
        # result = lmfit.minimize(_f_res, params, args=(y_a,), calc_covar=False, method = 'least_squares',ftol=1.49012e-04, xtol=1.49012e-05)
        result = lmfit.minimize(_f_res, params, args=(y_a,),calc_covar=False, ftol=1.49012e-05, xtol=1.49012e-05)
        if result.success:
            fits = _parameters_to_list(result.params)

            for i, fit in enumerate(fits):
                peak_key = f"peak{i}"
                fitting_result[peak_key] = {
                    'peak_y': [],
                    'peak_params': fit
                }

            for peak_index, fit_info in fitting_result.items():
                peak_params = fit_info['peak_params']
                peak_y = _f_conv([peak_params], x_ppm)
                fit_info['peak_y'] = peak_y
            residual = np.sum((result.residual) ** 2)
            residue_curve = _f_conv(fits, x_ppm) - y_a
            
            fitting_result["residual"] = residual
            fitting_result["residue_curve"] = residue_curve
        else:
            fitting_result = None
    except Exception as e:
        print("line_fitting error:", e)
        fitting_result = None
        # return None, None

    return fitting_result

# Convert parameter list to summed fitted curve.
def _f_conv(parameterset_list, x_ppm):
    """
    Build summed fitted curve from multiple peak parameter sets.

    Parameters
    parameterset_list : list
        List of [offset, sigma, hwhm, amplitude, frac_gauss].
    x_ppm : np.ndarray
        Chemical-shift axis in ppm.

    Returns
    conv : np.ndarray
        Summed fitted curve.
    """
    conv = np.zeros_like(x_ppm)
    for params in parameterset_list:
        peak = _f_pk(offset=params[0], amplitude=params[3], gauss_sigma=params[1], 
                     lorentz_hwhm=params[2], x=np.arange(len(x_ppm)), frac_gauss=params[-1])
        # _f_pk(offset, amplitude, gauss_sigma, lorentz_hwhm, x, frac_gauss=0):
        conv += peak
    return conv


def _f_makep(data, peaks, frac_gauss=None):
    """
    Generate initial peak parameters for optimization.

    Parameters
    data : np.ndarray
        Input spectrum intensity.
    peaks : list[int]
        Initial peak indices.
    frac_gauss : float or None
        Initial Gaussian fraction.

    Returns
    p : np.ndarray
        Initial parameter array.
    """
    p = []
    for peak_index in peaks:
        pamp = 0.9 * abs(data[int(peak_index)])
        single_peak = [peak_index, 10, pamp/300, pamp, frac_gauss]
        p.append(single_peak)
    return np.array(p)

def _f_res(params, data):
    """
    Residual function used by lmfit minimization.

    Parameters
    params : lmfit.Parameters
        Optimization parameters.
    data : np.ndarray
        Observed spectrum intensity.

    Returns
    res : np.ndarray
        Residual vector (observed - fitted).
    """
    residual = np.zeros_like(data)
    for i_, param_name in enumerate(params.keys()):
        if i_ % 5 == 0:
            i = int(i_ / 5)
            offset = params[f"peak_{i}_0"].value
            sigma = params[f"peak_{i}_1"].value
            hwhm = params[f"peak_{i}_2"].value
            amplitude = params[f"peak_{i}_3"].value
            frac_gauss = params[f"peak_{i}_4"].value

            peak = _f_pk(offset, amplitude, sigma, hwhm, np.arange(len(data)), frac_gauss)
            residual += peak
    res = data - residual
    return res


def _f_gauss(offset, amplitude, gauss_sigma, x):
    return amplitude * np.exp(-((offset - x) ** 2.0) / (2.0 * gauss_sigma ** 2.0))

def _f_lorentz(offset, amplitude, lorentz_hwhm, x):
    return amplitude * lorentz_hwhm ** 2.0 / (lorentz_hwhm ** 2.0 + (x - offset) ** 2.0)

def _f_pk(offset, amplitude, gauss_sigma, lorentz_hwhm, x, frac_gauss=0):
    if frac_gauss > 1.0:
        frac_gauss = 1.0
    if frac_gauss < 0.0:
        frac_gauss = 0.0
    
    gauss_peak = _f_gauss(offset, amplitude, gauss_sigma, x)
    lorentz_peak = _f_lorentz(offset, amplitude, lorentz_hwhm, x)
    peak = frac_gauss * gauss_peak + (1 - frac_gauss) * lorentz_peak
    return peak


def _parameters_to_list(params):
    """
    Convert lmfit parameter container to list format.

    Parameters
    params : lmfit.Parameters
        Fitted parameter container.

    Returns
    p : list
        List of [offset, sigma, hwhm, amplitude, frac_gauss] per peak.
    """
    p = []
    peak_indices = set()
    for param_name in params.keys():
        parts = param_name.split('_')
        if len(parts) == 3:
            peak_index = int(parts[1])
            if peak_index not in peak_indices:
                peak_indices.add(peak_index)
                peak_params = [params[f'peak_{peak_index}_{j}'].value for j in range(5)]
                p.append(peak_params)
    return p

def get_peak_key(rs_fitting):
    """
    Extract ordered peak keys from fitting result.

    Parameters
    rs_fitting : dict
        Result dict returned by line_fitting().

    Returns
    peak_key_ls : list
        Peak keys, e.g. ['peak0', 'peak1', 'peak2'].
    """
    peak_key_ls = []
    for i in rs_fitting.keys():
        peak_key_ls.append(i)
        #     peak_key_ls.append(i)
        # else:
        #     pass
    del peak_key_ls[-2:]

    return peak_key_ls

