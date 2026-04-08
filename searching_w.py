
import pandas as pd
import numpy as np
from scipy import stats,interpolate,signal
from itertools import combinations
from line_fitting import line_fitting, get_peak_key
from integration import integrate_sum, integrate_peak_params

def data_preparing(clu_data, df_sample_, step_, bl_, rw_):
    '''
    Prepare cluster template and candidate sample window for matching.

    Parameters
    clu_data : pd.DataFrame
        Template rows for one cluster.
    df_sample_ : pd.DataFrame
        Sample spectrum dataframe with columns ['shift', 'intensity'].
    step_ : float
        Sampling step used for template interpolation.
    bl_ : float
        Baseline intensity scale.
    rw_ : float
        Search extension width around template window.

    Returns
    pw : int
        Number of detected peaks in sample window.
    w_info : list
        [clu_sy, px, window_aim, peak_indices], where px is expected peak count.
    '''
    # global df_sample,step,bl,rw    
 
    # step_ = step
    # df_sample_ = df_sample
    # bl_ = bl
    # rw_ = rw

    px = int(np.unique(clu_data['peaks'])[0])
    clu_start = clu_data.iloc[0,0]
    clu_end = clu_data.iloc[-1,0]
    clu_sx = np.arange(clu_start, clu_end, float(step_))
    clu_func = interpolate.UnivariateSpline(clu_data['shift'],clu_data['A'],s=0)
    clu_sy = clu_func(clu_sx)
    
    window_aim = df_sample_[(df_sample_['shift'] <= (clu_end + rw_)) & (df_sample_['shift'] >= (clu_start-rw_))]
    peak_indices,property = signal.find_peaks(window_aim['intensity'], height=bl_*3)

    if len(peak_indices) < px:
        for h_factor in [2.5, 2.0, 1.5, 1.2, 1.0, 0.8]:
            peak_indices_relaxed, _ = signal.find_peaks(window_aim['intensity'], height=bl_ * h_factor)
            if len(peak_indices_relaxed) >= px:
                peak_indices = peak_indices_relaxed
                break
    pw = len(peak_indices)

    w_info = [clu_sy, px, window_aim, peak_indices]

    return pw, w_info

 

def searching(rolling,clu_sy,thres):
    '''
    Compute Pearson similarity for each rolling window.

    Parameters
    rolling : pandas.core.window.Rolling
        Rolling windows over sample candidate region.
    clu_sy : np.ndarray
        Interpolated cluster template intensity.
    thres : float
        Similarity threshold.

    Returns
    out_r : list
        Matching windows as [pi, start_index, end_index].
    '''
    out_r = []

    for r_ in rolling:
        if len(r_) == len(clu_sy):
            pi = stats.pearsonr(clu_sy,r_.iloc[:,-1].values)[0] 

            if pi >= thres:  
                r_index = r_.index.to_list()
                out_r.append([pi,r_index[0],r_index[-1]]) 
            else:
                pass
        else:
            continue    
    return out_r      

def searching_re(rolling,clu_sy,thres,bl_):
    '''
    Secondary search with additional peak-count constraint.

    Parameters
    rolling : pandas.core.window.Rolling
        Rolling windows over candidate reconstructed curve.
    clu_sy : np.ndarray
        Interpolated cluster template intensity.
    thres : float
        Similarity threshold.
    bl_ : float
        Baseline intensity scale.

    Returns
    out_r : list
        Matching windows as [pi, start_index, end_index].
    '''
    out_r = []
    # out_debug_p = []
    # out_debug_s = []

    px_pick, property = signal.find_peaks(clu_sy,height=bl_*3)
    px = len(px_pick)
    for r_ in rolling:
        if len(r_) == len(clu_sy):
            peak_ind_rolling, property = signal.find_peaks(r_.iloc[:,-1].values,height=bl_*3)

            if len(peak_ind_rolling) >= (px-1):
                pi = stats.pearsonr(clu_sy,r_.iloc[:,-1].values)[0] 
                if pi > thres:  
                    r_index = r_.index.to_list()
                    out_r.append([pi,r_index[0],r_index[-1]]) 
                    # df = pd.DataFrame([clu_sy,r_.iloc[:,-1].values])
                    # df.to_excel('./val-c1-rolling.xlsx')
                    # out_debug_p.append([pi,r_.iloc[0,0],r_.iloc[-1,0]]) 
                else:
                    pass
            else:
                continue    
        else: 
            pass
    # df_debug_p = pd.DataFrame(out_debug_p,columns=['r','rolling-start','rolling-end'])
    # df_debug_s = pd.DataFrame(out_debug_s,columns=['spm','rolling-start','rolling-end'])
    # df_debug_p.to_excel('./0612-gaba-c1-pi.xlsx')
    # df_debug_s.to_excel('./0612-gaba-c1-spm.xlsx')
    return out_r 

def searching_continuing(cx,w_info,bl_,thres_,cpd_aim='7660-25-5',cpd_name=''):
    '''
    Continue searching for one cluster with fitting/recombination fallback.

    Parameters
    cx : str
        Cluster name.
    w_info : list
        Output from data_preparing().
    bl_ : float
        Baseline intensity scale.
    thres_ : float
        Similarity threshold.
    cpd_aim : str, optional
        Target CAS.
    cpd_name : str, optional
        Target compound name.

    Returns
    output_cx : list
        [fit_peaks, cx, Hs, abs_cx, pi, combination] or [] when not found.
    '''
    # global bl, threshold
    # bl_ = bl
    # thres_ = threshold

    window_aim = w_info[2]
    rolling = window_aim.rolling(len(w_info[0]), step = 3)

    out_r = searching(rolling=rolling,clu_sy=w_info[0],thres=thres_)
    
    if len(out_r) > 0: ###########################################################
        out_r_max = max(out_r, key=lambda x: x[0])

        fit_range = window_aim.loc[out_r_max[1]:out_r_max[-1]+1,:]

        # line fitting & integration
        y_a = fit_range['intensity'].values
        x_ppm = fit_range['shift'].values
        peak_indices,property = signal.find_peaks(y_a, height=bl_*10)

        fit_peaks = []
        for p_ind in peak_indices:
            peak_ppm = fit_range.iloc[p_ind,0]
            fit_peaks.append(peak_ppm)

        fitting_cx = line_fitting(x_ppm, y_a, peaks = peak_indices.tolist())
        peak_key_cx = get_peak_key(fitting_cx)

        abs_cx = 0
        for p_cx in peak_key_cx:
            area_cx = integrate_peak_params(peak_key=p_cx, peak_info=fitting_cx)
            abs_cx += area_cx

        output_cx = [fit_peaks,cx,2,abs_cx,out_r_max[0],'-']
    else: #####################################################
        fitting_w = line_fitting(x_ppm=window_aim['shift'].values, y_a=window_aim['intensity'].values, peaks=w_info[-1])
        peak_key_cx_re = get_peak_key(fitting_w)
        all_combinations = list(combinations(peak_key_cx_re, w_info[1]))

        for comb in all_combinations:
            comb_y = np.zeros_like(window_aim['shift'].values)
            for key in comb:
                comb_y += fitting_w[key]['peak_y']                    
            comb_df = window_aim.copy()
            comb_df['intensity'] = comb_y

            rolling_re = comb_df.rolling(len(w_info[0]),step=3)
            out_re = searching_re(rolling = rolling_re, thres=thres_*0.95, clu_sy=w_info[0],bl_=bl_)
            
            if len(out_re) > 0:
                out_re_max = max(out_re, key=lambda x: x[0])
                
                abs_cx_re = 0
                for p_cx_re in comb: 
                    area_cx_re = integrate_peak_params(peak_key=p_cx_re, peak_info=fitting_w)
                    abs_cx_re += area_cx_re                    
                ind0 = out_re_max[1]
                ind1 = out_re_max[2]
                out_range = comb_df.loc[ind0:ind1,:] 

                peak_ind_re, property = signal.find_peaks(out_range['intensity'].values, height=bl_*10)
                peak_ppm_re_ls = []
                for p_re_ind in peak_ind_re:
                    peak_ppm = out_range.iloc[p_re_ind,0]
                    peak_ppm_re_ls.append(peak_ppm)

                output_cx = [peak_ppm_re_ls,cx,2,abs_cx_re,out_re_max[0],comb]
            else:
                output_cx = []
    
    
    return output_cx            


