import numpy as np
import pandas as pd
from scipy import stats,signal,spatial
import re
import os
from ast import literal_eval
from glob import glob
from itertools import combinations
from datetime import datetime
import plotly.graph_objects as go
import plotly.express as pex
from integration import integrate_sum, integrate_peak_params
from line_fitting import line_fitting, get_peak_key
from searching_w import data_preparing, searching,searching_re
from mono_sugar import sugar
import visualize


rw = 0.02
threshold = 0.80
export_html = False
html_output_dir = './fitting_html'
export_lumin_style = True

if export_html:
    os.makedirs(html_output_dir, exist_ok=True)


def export_lumin_style_plot(df_spectrum, df_quant, out_dir, name_sample_):
    """Export Lumin-style overview plot: full spectrum background + compound peak markers."""
    fig = go.Figure()
    cas_to_eng = dict(zip(cpd_info['cas'].astype(str), cpd_info['ENG'].astype(str)))

    def normalize_eng_name(name):
        return str(name)

    fig.add_trace(go.Scatter(
        x=df_spectrum['shift'],
        y=df_spectrum['intensity'],
        mode='lines',
        name='Spectrum',
        line=dict(color='#696969', width=1.5),
        showlegend=False
    ))

    if df_quant is not None and (not df_quant.empty):
        df_results = df_quant[df_quant['Compound Name'] != 'TSP'].copy()

        def parse_peaks_robust(peak_data):
            if isinstance(peak_data, list):
                return peak_data
            if isinstance(peak_data, str):
                try:
                    return literal_eval(peak_data)
                except (ValueError, SyntaxError):
                    try:
                        return [float(peak_data)]
                    except ValueError:
                        return []
            if isinstance(peak_data, (int, float, np.floating, np.integer)):
                return [float(peak_data)]
            return []

        df_results['peak_list'] = df_results['Peak Positions'].apply(parse_peaks_robust)

        df_results['compound_eng'] = df_results.apply(
            lambda r: normalize_eng_name(cas_to_eng.get(str(r['cas']), str(r['Compound Name']))),
            axis=1
        )

        compounds = df_results['compound_eng'].astype(str).unique()
        colormap = pex.colors.qualitative.Plotly + pex.colors.qualitative.Alphabet
        color_mapping = {cpd: colormap[i % len(colormap)] for i, cpd in enumerate(compounds)}

        for cpd_name, df_cpd in df_results.groupby('compound_eng', sort=False):
            peak_values = []
            for peaks_raw in df_cpd['peak_list']:
                if peaks_raw and isinstance(peaks_raw[0], list):
                    for sub in peaks_raw:
                        peak_values.extend(sub)
                else:
                    peak_values.extend(peaks_raw if isinstance(peaks_raw, list) else [])

            if len(peak_values) == 0:
                continue

            peak_idx = pd.Index(df_spectrum['shift']).get_indexer(peak_values, method='nearest')
            valid_idx = [i for i in peak_idx if i != -1]
            if len(valid_idx) == 0:
                continue

            peak_shifts = df_spectrum.iloc[valid_idx]['shift']
            peak_intensities = df_spectrum.iloc[valid_idx]['intensity']
            rel_vals = df_cpd['Relative'].dropna().values if 'Relative' in df_cpd.columns else []
            rel_text = f"{float(rel_vals[0]):.4f}" if len(rel_vals) > 0 else 'N/A'

            hover_texts = [
                f"<b>{cpd_name}</b><br>Peak: {s:.4f} ppm<br>Relative: {rel_text}"
                for s in peak_shifts
            ]

            fig.add_trace(go.Scatter(
                x=peak_shifts,
                y=peak_intensities,
                mode='markers',
                name=str(cpd_name),
                marker=dict(
                    color=color_mapping.get(str(cpd_name), '#000000'),
                    size=10,
                    line=dict(width=2, color='#c80000')
                ),
                legendgroup=str(cpd_name),
                showlegend=True,
                hoverinfo='text',
                text=hover_texts
            ))

    fig.update_layout(
        title_text=f"Lumin Quantitative Visualization - {name_sample_}",
        xaxis_title='Chemical Shift (ppm)',
        yaxis_title='Intensity',
        xaxis=dict(range=[10.0, -0.5], zeroline=False, ticks='inside', linewidth=1.5, tickwidth=1.5),
        yaxis=dict(zeroline=False, ticks='', showticklabels=False, linewidth=1.5),
        template='simple_white',
        height=900,
        width=1500,
        showlegend=True
    )

    os.makedirs(out_dir, exist_ok=True)
    html_path = os.path.join(out_dir, f'visualization-{name_sample_}.html')
    json_path = os.path.join(out_dir, f'visualization-{name_sample_}.json')

    html_content = fig.to_html(full_html=True, include_plotlyjs='cdn')
    with open(html_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

    json_content = fig.to_json()
    with open(json_path, 'w', encoding='utf-8') as f:
        f.write(json_content)

    print(f'Lumin-style visualization saved: {html_path}')
    print(f'Lumin-style JSON saved: {json_path}')

p_cas = re.compile(r"./templates\\(\S*).csv")
p_sample = re.compile(r"./sample\\(.*).csv")
files_template = glob('./templates/*.csv')
files_sample = glob('./sample/*.csv')
cpd_info = pd.read_excel('./sum_cpdinfo.xlsx')

# G10-01 Input the list of CAS numbers for compounds to analyze. The program will automatically check if the corresponding template files exist in the ./templates folder, and skip the compound if the template is not found.
cpd_aim_ls = ['7004-03-7', '56-41-7', '98-79-3', '70-26-8', '50-21-5', '56-12-2', 'anomeric-7660-25-5', '62-49-7', '58-61-7']

file_s = files_sample[0]

name_sample = p_sample.match(file_s).group(1)
output_ls = []

df_sample = pd.read_csv(file_s, names=['shift','intensity']) * 1
df_sample = df_sample.dropna()
df_sample = df_sample.reset_index(drop=True)

mask = (0.2 <= df_sample.iloc[:,0]) & (df_sample.iloc[:,0] <= 0.6)
data_bl = df_sample.loc[mask]
max_ = np.max(data_bl.iloc[:, -1])
min_ = np.min(data_bl.iloc[:, -1])
bl = max_-min_
step = (data_bl['shift'].max() - data_bl['shift'].min())/(data_bl.shape[0]-1)

######################
data_tsp = df_sample[(df_sample['shift'] < 0.05) & (df_sample['shift'] > -0.05)]
y_a = data_tsp['intensity'].values
x_ppm = data_tsp['shift'].values
peak_indices,property = signal.find_peaks(y_a, height=bl*5)
fitting_tsp = line_fitting(x_ppm, y_a, peaks = peak_indices.tolist())
peak_key_tsp = get_peak_key(fitting_tsp)

abs_tsp = 0
for p_ in peak_key_tsp:
    # area_ = integrate_peak_params(peak_key=p_, peak_info=fitting_tsp)
    area_ = integrate_sum(peak_key=p_, peak_info=fitting_tsp, x_ppm=x_ppm, r=None, algorithm='Trapezoidal')
    abs_tsp = abs_tsp + area_

output_is = ['TSP','24493-21-8','0.00','-',9, abs_tsp, 172.263,'-','-']
output_ls.append(output_is) 
######################

for cpd_aim in cpd_aim_ls:
    filepath_aim = './templates/' + cpd_aim + '.csv'
    data_aim = pd.read_csv(filepath_aim, header=0)

    if data_aim.shape[1] == 6:
        cpd_aim_ = cpd_aim[9:]
        cpd_name = cpd_info.loc[cpd_info['cas']==cpd_aim_,'ENG'].values[0]
        cpd_m = cpd_info.loc[cpd_info['cas']==cpd_aim_,'M'].values[0]

        output_sugar = sugar(cpd_aim=cpd_aim_, cpd_name=cpd_name,cpd_m=cpd_m,data_aim=data_aim, thres_= threshold, rw_=rw,df_sample_ =df_sample,step_=step,bl_=bl)
        if output_sugar == None:
            print(f'{cpd_name}: automatic quantification failed')
            pass
        else:
            output_ls.append(output_sugar)

        continue
    else:
        pass

    cpd_name = cpd_info.loc[cpd_info['cas']==cpd_aim,'ENG'].values[0]
    cpd_m = cpd_info.loc[cpd_info['cas']==cpd_aim,'M'].values[0]

    w_dict = {}
    clusters = data_aim['cluster'].dropna().unique()
    clusters.sort()        
    count_cx = 0
    for cx in clusters:
        clu_data = data_aim[data_aim['cluster']==cx]
        cx_h = data_aim.loc[data_aim['cluster']==cx, 'Hs'].values[0] 

        pw, w_info = data_preparing(clu_data,df_sample_ =df_sample,step_=step,bl_=bl,rw_=rw)
        w_dict[cx] = w_info
        px = w_info[1]

        if pw < px-1:
            continue
        else:
            window_aim = w_info[2]

            if px > 1:
                rolling = window_aim.rolling(len(w_info[0]), step = 3)
                out_r = searching(rolling=rolling,clu_sy=w_info[0],thres=threshold)
            else:
                if pw == 1:
                    amplitude_pw = window_aim.iloc[w_info[-1][0],-1]
                    if  amplitude_pw >= bl*10:
                        peak_indices = w_info[-1]
                        rolling = window_aim.rolling(len(w_info[0]), step = 3)
                        out_r = searching(rolling=rolling,clu_sy=w_info[0],thres=threshold*0.95)
                    else:
                        continue
                else:
                    continue

            if len(out_r) > 0:
                out_r_max = max(out_r, key=lambda x: x[0])
                pi_out = out_r_max[0]
                
                fit_range = window_aim.loc[out_r_max[1]:out_r_max[-1]+1,:]

                # line fitting & integration
                y_a = fit_range['intensity'].values
                x_ppm = fit_range['shift'].values
                peak_indices,property = signal.find_peaks(y_a, height=bl*5)

                if len(peak_indices)==0:
                    for h_factor in [4, 3, 2, 1.5, 1.2, 1.0, 0.8]:
                        peak_indices_relaxed, property_relaxed = signal.find_peaks(y_a, height=bl*h_factor)
                        if len(peak_indices_relaxed) > 0:
                            peak_indices = peak_indices_relaxed
                            property = property_relaxed
                            break
                
                if len(peak_indices)==0:
                    print(f"{cpd_name} {cx} localized successfully ({[np.min(x_ppm), np.max(x_ppm)]}, similarity={out_r_max[0]}), but no valid peak was found in fit range.")
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
                    # output_cpd = [cpd_name ,cpd_aim, fit_peaks,cx,cx_h,0,cpd_m,pi_out,'-']
                    continue
                else:                   
                    ###################################
                    pass
                    # info_name_ = f'./fitting_result_{cpd_aim}_{name_sample}_frac00.hs'
                    # visualize.save_fitting_results(info_name_, cx, window_aim, fit_range, fitting_cx)

                    # if export_html:
                    #     fig_fit = visualize.plot_fitted_curve(fit_range, fitting_cx, w_range=window_aim)
                    #     html_name = f'{html_output_dir}/fitting_result_{cpd_aim}_{name_sample}_{cx}.html'
                    #     fig_fit.write_html(html_name)

                    ###################################             

                peak_key_cx = get_peak_key(fitting_cx)
                # for p_ind in peak_indices:
                #     peak_ppm = fit_range.iloc[p_ind,0]
                #     fit_peaks.append(peak_ppm)   

                #########################     
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
                             
                        if len(peak_key_cx)==px:
                            for p_ind in peak_indices:
                                peak_ppm = fit_range.iloc[p_ind,0]
                                fit_peaks.append(peak_ppm)

                            pass
                        else:
                            peak_ind_t,property_cx = signal.find_peaks(clu_data.iloc[:,1],height=3)
                            peak_heights_t = property_cx['peak_heights']
                            i_comb = []
                            for i_c in range(0,len(peak_key_cx)):
                                i_comb.append(i_c)

                            comb = list(combinations(i_comb, px))
                            peak_max_index = np.argmax(peak_heights)

                            rs = []
                            peak_heights_norm = peak_heights / np.linalg.norm(peak_heights)
                            peak_heights_t_norm = peak_heights_t / np.linalg.norm(peak_heights_t)
                            
                            for comb_ in comb:                                
                                if peak_max_index in comb_:
                                    # ratio_peaks = stats.pearsonr(peak_heights[list(comb_)], peak_heights_t)[0]
                                    cos_peaks = 1 - spatial.distance.cosine(peak_heights_norm[list(comb_)], peak_heights_t_norm)
                                    rs.append([comb_, cos_peaks])
                                else:
                                    continue
                            
                            rs_df = pd.DataFrame(rs,columns=['peak combination','cos'])
                            comb_final = rs_df.iloc[rs_df['cos'].idxmax(),0]

                            peak_key = []
                            for i_ in comb_final:
                                peak_key.append(peak_key_cx[i_])
                                peak_ppm = fit_range.iloc[peak_indices[i_],0]
                                fit_peaks.append(peak_ppm) 
                            peak_key_cx = peak_key                                                   
                    else:
                        for p_ind in peak_indices:
                            peak_ppm = fit_range.iloc[p_ind,0]
                            fit_peaks.append(peak_ppm)                          
                else:
                    for p_ind in peak_indices:
                        peak_ppm = fit_range.iloc[p_ind,0]
                        fit_peaks.append(peak_ppm)         


                ###############################
                abs_cx = 0
                for p_cx in peak_key_cx:
                    # area_cx = integrate_peak_params(peak_key=p_cx, peak_info=fitting_cx) 

                    area_cx = integrate_sum(peak_key=p_cx, peak_info=fitting_cx, x_ppm=x_ppm, r=None, algorithm='Trapezoidal')

                    abs_cx += area_cx

                output_cpd = [cpd_name ,cpd_aim, fit_peaks,cx,cx_h,abs_cx,cpd_m,pi_out,'-']
                output_ls.append(output_cpd)
                count_cx += 1
                break
            else:
                continue

    if count_cx > 0:
        continue
    else:
        pass  

    # ##############################################################        
    re_clu_data = data_aim[(data_aim['peaks']<7)&(data_aim['peaks']>1)]
    clusters_re = re_clu_data['cluster'].dropna().unique()
    clusters_re.sort()
    
    count=0
    for cx_re in clusters_re:
        px_re = w_dict[cx_re][1]
        pw_re_ind = w_dict[cx_re][-1]
        w_re = w_dict[cx_re][-2]

        if px_re > len(pw_re_ind):
            pw_re_ind_relaxed, _ = signal.find_peaks(w_re['intensity'].values, height=bl*0.8)
            if len(pw_re_ind_relaxed) >= px_re:
                pw_re_ind = pw_re_ind_relaxed
            else:
                print(f'{cpd_name}: second-stage {cx_re} has pw < px')
                continue

        clu_sy_re = w_dict[cx_re][0]

        fitting_re = line_fitting(x_ppm=w_re['shift'].values, y_a = w_re['intensity'].values, peaks = pw_re_ind.tolist())
        if fitting_re is None:
            print(f'{cpd_name}: second-stage fitting failed for {cx_re} window W')
            continue
        else:
            pass
                    
        peak_key_cx_re = get_peak_key(fitting_re)
        all_combinations = list(combinations(peak_key_cx_re, w_dict[cx_re][1]))
        
        for comb in all_combinations:
            comb_y = np.zeros_like(w_re['shift'].values)
            for key in comb:
                comb_y += fitting_re[key]['peak_y']                    
            comb_df = w_re.copy()
            comb_df['intensity'] = comb_y

            rolling_re = comb_df.rolling(len(clu_sy_re),step=3)
            out_re = searching_re(rolling = rolling_re,thres=threshold*0.95,clu_sy=clu_sy_re, bl_=bl)
            
            if len(out_re) > 0:
                out_re_max = max(out_re, key=lambda x: x[0])
                
                abs_cx_re = 0
                for p_cx_re in comb: 
                    area_cx_re = integrate_peak_params(peak_key=p_cx_re, peak_info=fitting_re)
                    abs_cx_re += area_cx_re                    
                ind0 = out_re_max[1]
                ind1 = out_re_max[2]
                out_range = comb_df.loc[ind0:ind1,:] 

                peak_ind_re, property = signal.find_peaks(out_range['intensity'].values, height=5*bl)
                peak_ppm_re_ls = []
                for p_re_ind in peak_ind_re:
                    peak_ppm = out_range.iloc[p_re_ind,0]
                    peak_ppm_re_ls.append(peak_ppm)                         

                # cpd_normalize_re = round(cpd_normalize_re,2)
                # output_re = [cpd_name, cpd_aim, peak_ppm_re_ls, cx_re, data_aim.loc[data_aim['cluster']==cx_re, 'Hs'].values[0],
                #               abs_cx_re, cpd_normalize_re, cpd_m, out_re_max[0], comb] 
                output_re = [cpd_name, cpd_aim, peak_ppm_re_ls, cx_re, data_aim.loc[data_aim['cluster']==cx_re, 'Hs'].values[0],
                                abs_cx_re, cpd_m, out_re_max[0], comb]                         
                output_ls.append(output_re)
                count  += 1
            else:
                pass

        if count > 0:
            break
        else:
            continue
    
    if count == 0:
        print(f'{cpd_name}: failed to localize after two-stage processing; automatic quantification failed')
    else:
        pass        

col_names = ['Compound Name','cas','Peak Positions','Cluster Name','Hs','Absolute','M','pi','Comb-re']
output_df = pd.DataFrame(output_ls,columns=col_names)

normalized = output_df.iloc[:,5] / abs_tsp
output_nor = round(normalized,2)
output_df.insert(9,column='Normalized',value=output_nor)
# Cx = Ax x Ns x Mx x Cs /(As x Nx x Ms)
cpd_relative = (output_df.iloc[:,5] * 9 * output_df.iloc[:, 6]) / (abs_tsp * output_df.iloc[:,4] * 172.263)
output_df.insert(10,column='Relative',value=cpd_relative)

# output_df.to_excel(f'./quantified-thre.xlsx') 
output_path = f'./0624quantified-{name_sample}.xlsx'
try:
    output_df.to_excel(output_path)
except PermissionError:
    ts = datetime.now().strftime('%Y%m%d-%H%M%S')
    output_path = f'./0624quantified-{name_sample}-{ts}.xlsx'
    output_df.to_excel(output_path)
    print(f'Output file is locked; saved as: {output_path}')
else:
    print(f'Result saved: {output_path}')

if export_lumin_style:
    export_lumin_style_plot(df_sample, output_df, out_dir=html_output_dir, name_sample_=name_sample)
                  


