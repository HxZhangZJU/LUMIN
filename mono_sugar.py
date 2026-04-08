
import numpy as np
import pandas as pd
from scipy import stats,interpolate,signal
from itertools import combinations
from integration import integrate_peak_params, integrate_sum
from line_fitting import line_fitting, get_peak_key
from searching_w import data_preparing, searching, searching_re
import visualize



def sugar(cpd_aim,data_aim,cpd_name,cpd_m,thres_,rw_,df_sample_,step_,bl_):
    '''
    Quantify one sugar-like compound using staged cluster matching.

    Parameters
    cpd_aim : str
        Target CAS number.
    data_aim : pd.DataFrame
        Template rows of the target compound.
    cpd_name : str
        Compound display name.
    cpd_m : float
        Molecular weight.
    thres_, rw_, df_sample_, step_, bl_ :
        Search and fitting control parameters from upstream pipeline.

    Returns
    output_sugar : list or None
        [compound_name, cas, peak_positions, clusters, Hs, absolute_area,
         molecular_weight, pi_list, combination_info]

    '''   

    if cpd_aim == '7660-25-5':
        output_sugar=sugar_fru(data_aim,thres_,rw_,df_sample_,step_,bl_)
        return output_sugar
    else:
        pass
    

    w_dict = {}
    cx_config = {}
    configs = data_aim['config'].dropna().unique()
    for config_ in configs:
        data = data_aim[data_aim['config'] == config_]
        clusters = data['cluster'].dropna().unique()
        clusters.sort()

        count_cx = 0
        for cx in clusters:
            clu_data = data_aim[data_aim['cluster']==cx]
            cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 

            pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
            w_dict[cx] = w_info
            px = w_info[1]

            if pw < px-1:
                continue
            else:
                window_aim = w_info[2]
                rolling = window_aim.rolling(len(w_info[0]), step = 3)

                out_r = searching(rolling=rolling,clu_sy=w_info[0],thres=thres_)
                
                if len(out_r) > 0:
                    out_r_max = max(out_r, key=lambda x: x[0])

                    # line fitting & integration
                    fit_range = window_aim.loc[out_r_max[1]:out_r_max[-1]+1,:]
                    
                    y_a = fit_range['intensity'].values
                    x_ppm = fit_range['shift'].values
                    peak_indices,property = signal.find_peaks(y_a, height=bl_*10)

                    if len(peak_indices)==0:
                        print(f"{cpd_name}: no valid peaks found")
                        continue
                    else:
                        pass

                    fit_peaks = []
                    fitting_cx = line_fitting(x_ppm, y_a, peaks = peak_indices.tolist())
                    if fitting_cx is None:
                        for p_ind in peak_indices:
                            peak_ppm = fit_range.iloc[p_ind,0]
                            fit_peaks.append(peak_ppm)

                        print(f'{cpd_name} {cx} localized successfully, but fitting failed!')
                        print(f'Localization result: {cx} peak shifts={fit_peaks}, similarity={out_r_max[0]}')
                        continue
                    else:
                        info_name_ = f'./fitting_result_monosugar_{cpd_aim}_{cx}.hs'
                        visualize.save_fitting_results(info_name_, cx, window_aim, fit_range, fitting_cx)                                             
                    
                    peak_key_cx = get_peak_key(fitting_cx)
                    
                    if (px < 5)&(px >1):
                        if len(peak_indices) > px:
                            peak_heights = property['peak_heights'] # ndarray
                            peak_max = np.max(peak_heights)
                            i_ls = []
                            for i_ in range(0,len(peak_indices)):
                                if peak_heights[i_] <= peak_max/4:
                                    i_ls.append(i_)
                                else:
                                    pass
                                    
                            for i_pk in sorted(i_ls, reverse=True): 
                                del peak_key_cx[i_pk]  
                            peak_heights = np.delete(peak_heights, i_ls)
                            peak_indices = np.delete(peak_indices, i_ls)
                            
                            for p_ind in peak_indices:
                                peak_ppm = fit_range.iloc[p_ind,0]
                                fit_peaks.append(peak_ppm)      

                        else:
                            for p_ind in peak_indices:
                                peak_ppm = fit_range.iloc[p_ind,0]
                                fit_peaks.append(peak_ppm)                          
                    else:
                        for p_ind in peak_indices:
                            peak_ppm = fit_range.iloc[p_ind,0]
                            fit_peaks.append(peak_ppm)   
                    ###############################################

                    abs_cx = 0
                    # from integration import integrate_sum
                    for p_cx in peak_key_cx:
                        # area_cx = integrate_peak_params(peak_key=p_cx, peak_info=fitting_cx)
                        area_cx = integrate_sum(peak_key=p_cx, peak_info=fitting_cx, x_ppm=x_ppm, r=None, algorithm='Trapezoidal')
                        # area_cx = integrate_sum(peak_key=p_cx, peak_info=fitting_cx, x_ppm = x_ppm)
                        # abs_cx_ls_sum.append(area_cx)

                        abs_cx += area_cx

                    cx_config[config_] = [fit_peaks,cx,cx_h,abs_cx,out_r_max[0],'-']
                    count_cx += 1
                    break
                else:
                    continue

        if count_cx > 0:
            continue
        else:
            pass

        re_clu_data = data[data['peaks']<7]
        clusters_re = re_clu_data['cluster'].dropna().unique()
        clusters_re.sort()

        for cx_re in clusters_re:
            px_re = w_dict[cx_re][1]
            pw_re_ind = w_dict[cx_re][-1]

            if px_re > len(pw_re_ind):
                continue
            else:
                w_re = w_dict[cx_re][-2]
                clu_sy_re = w_dict[cx_re][0]

                fitting_re = line_fitting(w_re['shift'].values, w_re['intensity'].values, peaks =  w_dict[cx_re][-1].tolist())
                if fitting_re is None:
                    print(f'{cpd_name}: second-stage fitting failed for {cx_re} window W')
                    continue
                else:
                    pass                
                
                peak_key_cx_re = get_peak_key(fitting_re)
                all_combinations = list(combinations(peak_key_cx_re, w_dict[cx_re][1]))
                count=0

                for comb in all_combinations:
                    comb_y = np.zeros_like(w_re['shift'].values)
                    for key in comb:
                        comb_y += fitting_re[key]['peak_y']                    
                    comb_df = w_re.copy()
                    comb_df['intensity'] = comb_y

                    rolling_re = comb_df.rolling(len(clu_sy_re),step=3)
                    out_re = searching_re(rolling = rolling_re,thres=thres_*0.95,clu_sy=clu_sy_re,bl_=bl_)
                    
                    if len(out_re) > 0:
                        out_re_max = max(out_re, key=lambda x: x[0])
                        
                        abs_cx_re = 0
                        for p_cx_re in comb: 
                            # area_cx_re = integrate_peak_params(peak_key=p_cx_re, peak_info=fitting_re)
                            area_cx_re = integrate_sum(peak_key=p_cx_re, peak_info=fitting_re, x_ppm=x_ppm, r=None, algorithm='Trapezoidal')
                            abs_cx_re += area_cx_re                    
                        ind0 = out_re_max[1]
                        ind1 = out_re_max[2]
                        out_range = comb_df.loc[ind0:ind1,:] 

                        peak_ind_re, property = signal.find_peaks(out_range['intensity'].values, height=10*bl_)
                        peak_ppm_re_ls = []
                        for p_re_ind in peak_ind_re:
                            peak_ppm = out_range.iloc[p_re_ind,0]
                            peak_ppm_re_ls.append(peak_ppm)

                        cx_config[config_] = [peak_ppm_re_ls,cx_re,re_clu_data.loc[re_clu_data['cluster']==cx_re, 'Hs'].values[0],abs_cx_re,round(out_re_max[0],3),comb]    
                        count  += 1

                if count > 0:
                    break
                else:
                    continue 

    #############################################
    for c_ in configs:
        if c_ not in cx_config:
            print(f"configure {c_} not found!")
            break
        else:
            pass

    output_sugar = output_arrange(cx_config, cpd_name, cpd_aim, cpd_m)
    
    return output_sugar

def output_arrange(cx_config,cpd_name, cpd_aim, cpd_m):
    '''
    Aggregate per-cluster outputs into one result row.

    Parameters
    cx_config : dict
        Per-config matched output dictionary.
    cpd_name : str
        Compound name.
    cpd_aim : str
        CAS number.
    cpd_m : float
        Molecular weight.

    Returns
    output_arranged : list
        [compound_name, cas, peak_positions, clusters, Hs, absolute_area,
         molecular_weight, pi_list, combination_info]
    '''
    # global cpd_name, cpd_aim, cpd_m
    # cpd_name_=cpd_name
    # cpd_aim_ = cpd_aim
    # cpd_m_=cpd_m

    output_abs = 0
    output_peaks = []
    output_cx = []
    output_h = []
    output_pi = []
    output_comb = []    
    for c_key, c_output in cx_config.items():
        output_abs += c_output[3]
        output_peaks.append(c_output[0]) 
        output_cx.append(c_output[1])
        output_h.append(c_output[2])
        output_pi.append(c_output[-2])
        output_comb.append(c_output[-1])

    output_arranged = [cpd_name, cpd_aim, output_peaks, output_cx, list(set(output_h))[0],output_abs,cpd_m,output_pi,output_comb]

    return output_arranged


def sugar_fru(data_aim,thres_,rw_,df_sample_,step_,bl_,cpd_aim='7660-25-5',cpd_name=''):

    '''
    Specialized quantification path for fructose-like template behavior.

    Parameters
    data_aim : pd.DataFrame
        Template rows for fructose target.
    thres_, rw_, df_sample_, step_, bl_ :
        Search and fitting control parameters.
    cpd_aim : str, optional
        Target CAS number.
    cpd_name : str, optional
        Compound display name.

    Returns
    output_sugar : list or None
        Same output format as sugar().

    '''
    from searching_w import searching_continuing
    # global bl , threshold
    # bl_ = bl
    # thres_ = threshold 

    clusters = data_aim['cluster'].dropna().unique()
    clusters.sort() 
    # w_dict = {} 
    cx_config = {}

    


    
    sequence = {'c1':[31,47,64,79,96,128,140],'c2':[22,37],'c3':[19,49,61,97]}

    count_seq = 0
    for cx, peak_ind in sequence.items():
        clu_data = data_aim[data_aim['cluster']==cx]
        # cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 
        pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
        # w_dict[cx] = w_info
        px = w_info[1]
    
        window_aim = w_info[2]
        rolling = window_aim.rolling(len(w_info[0]), step = 3)
        out_r = searching(rolling=rolling,clu_sy=w_info[0],thres=thres_)


        if len(out_r) > 0: ###########################################################
            out_r_max = max(out_r, key=lambda x: x[0])

            fit_range = window_aim.loc[out_r_max[1]:out_r_max[-1]+1,:]

            # line fitting & integration
            y_a = fit_range['intensity'].values
            x_ppm = fit_range['shift'].values
            peak_ind_find,property = signal.find_peaks(y_a, height=bl_*5)
            if len(peak_ind_find) == len(peak_ind):
                peak_indices_cx = peak_ind_find
            else:
                peak_indices_cx = peak_ind

            fitting_cx = line_fitting(x_ppm, y_a, peaks = peak_indices_cx)
            if fitting_cx is None:
                peak_ppm_ls = []
                for p_ind in peak_ind_find:
                    peak_ppm = fit_range.iloc[p_ind,0]
                    peak_ppm_ls.append(peak_ppm)                 
                print(f'{cx} localized successfully, but fitting failed!')
                print(f'Localization result: {cx} peak shifts={peak_ppm_ls}, similarity={out_r_max[0]}')
                continue
            else:
                pass 
            peak_key_cx = get_peak_key(fitting_cx) 
            if cx == 'c1':
                abs_c11 = 0
                abs_c12 = 0
                peak_ppm_ls1 = []
                peak_ppm_ls2 = []

                for p_key in ['peak0','peak1','peak2','peak3','peak4']:
                    area_cx = integrate_peak_params(peak_key=p_key, peak_info=fitting_cx)
                    abs_c11 += area_cx

                    p_info = fitting_cx[p_key]
                    p_ind = p_info['peak_params'][0]
                    peak_ppm = fit_range.iloc[int(round(p_ind,0)),0]
                    peak_ppm_ls1.append(peak_ppm)
                cx_config['c1-1'] = [peak_ppm_ls1,cx,2,abs_c11,out_r_max[0],'-']

                for p_key in ['peak5','peak6']:
                    area_cx = integrate_peak_params(peak_key=p_key, peak_info=fitting_cx)
                    abs_c12 += area_cx

                    p_info = fitting_cx[p_key]
                    p_ind = p_info['peak_params'][0]
                    peak_ppm = fit_range.iloc[int(round(p_ind,0)),0]
                    peak_ppm_ls2.append(peak_ppm)
                cx_config['c1-2'] = [peak_ppm_ls2,cx,2,abs_c12,out_r_max[0],'-']
            else:
                abs_cx = 0
                peak_ppm_ls = []
                for p_key in peak_key_cx:
                    area_cx = integrate_peak_params(peak_key=p_key, peak_info=fitting_cx)
                    abs_cx += area_cx

                    p_info = fitting_cx[p_key]
                    p_ind = p_info['peak_params'][0]
                    peak_ppm = fit_range.iloc[int(round(p_ind,0)),0]
                    peak_ppm_ls.append(peak_ppm)
                cx_config[cx] = [peak_ppm_ls,cx,2,abs_cx,out_r_max[0],'-']
            count_seq += 1    
        else:
            if cx == 'c1':
                print('Failed at c1')
                return None # End ######################################## 
            else:
                pass

            if len(w_info[-1]) < px:
                print(f'Failed at {cx}: pw < px')
                continue
            else:
                pass

            fitting_w = line_fitting(x_ppm=window_aim['shift'].values, y_a=window_aim['intensity'].values, peaks=w_info[-1])
            peak_key_cx_re = get_peak_key(fitting_w)
            all_combinations = list(combinations(peak_key_cx_re, px)) 

            for comb in all_combinations:
                comb_y = np.zeros_like(window_aim['shift'].values)
                for key in comb:
                    comb_y += fitting_w[key]['peak_y']                    
                comb_df = window_aim.copy()
                comb_df['intensity'] = comb_y

                rolling_re = comb_df.rolling(len(w_info[0]),step=3)
                out_re = searching_re(rolling = rolling_re,thres=thres_*0.95,clu_sy=w_info[0],bl_=bl_)
                
                if len(out_re) > 0:
                    out_re_max = max(out_re, key=lambda x: x[0])
                    
                    abs_cx_re = 0
                    for p_cx_re in comb: 
                        area_cx_re = integrate_peak_params(peak_key=p_cx_re, peak_info=fitting_w)
                        abs_cx_re += area_cx_re                    
                    ind0 = out_re_max[1]
                    ind1 = out_re_max[2]
                    out_range = comb_df.loc[ind0:ind1,:] 

                    peak_ind_re, property = signal.find_peaks(out_range['intensity'].values, height=10*bl_)
                    peak_ppm_re_ls = []
                    for p_re_ind in peak_ind_re:
                        peak_ppm = out_range.iloc[p_re_ind,0]
                        peak_ppm_re_ls.append(peak_ppm)

                    cx_config[cx] = [peak_ppm_re_ls,cx,2,abs_cx_re,out_re_max[0],comb]
                    count_seq += 1
                else:
                    pass

    ########################################################################
    if count_seq == 3:
        output_sugar = output_arrange_fru(cx_config,cpd_name, cpd_aim, cpd_m=180.1559)
        output_sugar[4] = 2
        return  output_sugar
    else:
        print('Not all required clusters were resolved')
        
        if 'c1-1' not in cx_config:
            print('Missing c1-1 after processing')
            return None
        else:
            if 'c2' in cx_config:
                c_ = 'c6'
                clu_data = data_aim[data_aim['cluster']==c_]
                # cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 
                pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
                # w_dict[c_] = w_info
                px = w_info[1]

                if pw >= px:
                    c_output = searching_continuing(cx=c_,w_info=w_info,bl_=bl_,thres_=thres_)
                    if c_output != []:
                        output_sugar = output_arrange_fru(cx_config,cpd_name, cpd_aim, cpd_m=180.1559)
                        output_sugar[4] = 2
                        return  output_sugar # End #################################################################
                    else:
                        pass
                else:
                    pass

                for c_ in ['c4','c5','c7']:
                    clu_data = data_aim[data_aim['cluster']==c_]
                    # cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 
                    pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
                    # w_dict[c_] = w_info
                    px = w_info[1]

                    if pw < px:#######################################
                        continue
                    else:
                        pass
                    
                    c_output = searching_continuing(cx=c_,w_info=w_info,bl_=bl_,thres_=thres_)

                    if c_output == []:
                        print(f'c3 branch failed at {c_}')
                        break
                    else:
                        cx_config[cx] = c_output
                 
                    for c in ['c2','c4','c5','c7']:
                        if c not in cx_config:
                            print(f'c3 branch missing {c}')
                            return None # End ####################################
                        else:
                            pass                         
                    output_sugar = output_arrange_fru(cx_config,cpd_name, cpd_aim, cpd_m=180.1559)
                    output_sugar[4] = 3
                    return  output_sugar # End 
            else:
                c_ = 'c4'
                clu_data = data_aim[data_aim['cluster']==c_]
                # cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 
                pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
                # w_dict[c_] = w_info
                px = w_info[1]

                if pw < px:
                    print('c2 & c4 branch failed: pw < px')
                    return None # End ####################################
                else:
                    c_output = searching_continuing(cx=c_,w_info=w_info,bl_=bl_,thres_=thres_)
                    if c_output == []:
                        print('c2 & c4 branch failed during continuing search')
                        return None # End
                    else:
                        cx_config[cx] = c_output


                c_ = 'c3'
                clu_data = data_aim[data_aim['cluster']==c_]
                # cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 
                pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
                # w_dict[c_] = w_info
                px = w_info[1]
                if pw >= px:
                    c_output = searching_continuing(cx=c_,w_info=w_info)
                    if c_output != []:
                        output_sugar = output_arrange_fru(cx_config,cpd_name, cpd_aim, cpd_m=180.1559)
                        output_sugar[4] = 2
                        return  output_sugar # End
                    else:
                        pass
                else:
                    pass

                c_ = 'c6'
                clu_data = data_aim[data_aim['cluster']==c_]
                # cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 
                pw, w_info = data_preparing(clu_data,df_sample_,step_,bl_,rw_)
                # w_dict[c_] = w_info
                px = w_info[1]

                if pw >= px:
                    c_output = searching_continuing(cx=c_,w_info=w_info,bl_=bl_,thres_=thres_)
                    if c_output != []:
                        output_sugar = output_arrange_fru(cx_config,cpd_name, cpd_aim, cpd_m=180.1559)
                        output_sugar[4] = 2
                        return  output_sugar # End
                    else:
                        pass
                else:
                    print('c2 & c3 & c6 branch failed')
                    return None # End 

                    

def output_arrange_fru(cx_config,cpd_name, cpd_aim, cpd_m):
    '''
    Aggregate fructose-path cluster outputs into one result row.

    Parameters
    cx_config : dict
        Cluster outputs from sugar_fru().
    cpd_name : str
        Compound name.
    cpd_aim : str
        CAS number.
    cpd_m : float
        Molecular weight.

    Returns
    output_arrange_frud : list
        [compound_name, cas, peak_positions, clusters, Hs, absolute_area,
         molecular_weight, pi_list, combination_info]
    '''
    # global cpd_name, cpd_aim, cpd_m
    # cpd_name_=cpd_name
    # cpd_aim_ = cpd_aim
    # cpd_m_=cpd_m

    output_abs = 0
    output_peaks = []
    output_cx = []
    output_h = []
    output_pi = []
    output_comb = []
    # cx_config['c1-1'] = [peak_ppm_ls1,cx,2,abs_c11,out_r_max[0],'-']
    if 'c2' not in cx_config.keys():
        # output_abs += cx_config['c1-1'][3]
        # output_peaks.append(cx_config['c1-1'][0])
        # output_abs += cx_config['c1-2'][3]
        # output_peaks.append(cx_config['c1-2'][0])
        # output_cx.append('c1')
        # output_h.append(2)
        # output_pi.append(cx_config['c1-1'][-2])
        # output_comb.append('-')
                           
        for c_key, c_output in cx_config.items():
            if c_key == 'c1-2':
                pass
            else:
                output_abs += c_output[3]
                output_peaks.append(c_output[0]) 
                output_cx.append(c_output[1])
                output_h.append(c_output[2])
                output_pi.append(round(c_output[-2],3))
                output_comb.append(c_output[-1])

            del output_cx[0]
            output_cx.insert(0,'c1')            
                 
    else:
        for c_key, c_output in cx_config.items():
            output_abs += c_output[3]
            output_peaks.append(c_output[0]) 
            output_cx.append(c_output[1])
            output_h.append(c_output[2])
            output_pi.append(round(c_output[-2],3))
            output_comb.append(c_output[-1])
        del output_cx[0:2]
        del output_pi[0]
        output_cx.insert(0,'c1')

        # output_abs += cx_config['c1-1'][3]
        # output_peaks.append(cx_config['c1-1'][0])
        # output_cx.append('c1')
        # output_h.append(2)
        # output_pi.append(cx_config['c1-1'][-2])
        # output_comb.append('-')
        # for c_key, c_output in cx_config.items():
        #     if c_key[0:2] == 'c1':
        #         output_abs += c_output[3]
        #         output_peaks.append(c_output[0]) 
        #         output_cx.append(c_output[1])
        #         output_h.append(c_output[2])
        #         output_pi.append(round(c_output[-2],3))
        #         output_comb.append(c_output[-1])
        #     else:
        #         output_abs += c_output[3]
        #         output_peaks.append(c_output[0]) 
        #         output_cx.append(c_output[1])
        #         output_h.append(c_output[2])
        #         output_pi.append(round(c_output[-2],3))
        #         output_comb.append(c_output[-1])

                
        
    output_arrange_frud = [cpd_name, cpd_aim, output_peaks, output_cx, list(set(output_h))[0],
                           output_abs,cpd_m,output_pi,output_comb]

    return output_arrange_frud

