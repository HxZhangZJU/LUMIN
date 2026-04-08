import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import pandas as pd
import numpy as np
import os
import plotly.graph_objs as go
import plotly.express as px
from searching_w import data_preparing, searching

cpd_info = pd.read_excel('sum_cpdinfo.xlsx')
template_ids = {
    os.path.splitext(name)[0]
    for name in os.listdir('./templates')
    if name.lower().endswith('.csv')
}

def resolve_template_name(cas_value):
    cas_str = str(cas_value)
    if cas_str in template_ids:
        return cas_str
    anomeric_name = f'anomeric-{cas_str}'
    if anomeric_name in template_ids:
        return anomeric_name
    return None

cpd_info['template_name'] = cpd_info['cas'].apply(resolve_template_name)
cpd_info = cpd_info[cpd_info['template_name'].notna()].copy()
compound_list = cpd_info['ENG'].tolist()
clusters_selected = []
file_sample = None

window = tk.Tk()
window.title("Similarity Threshold Evaluation Tool")
window.geometry("550x500")

window.grid_rowconfigure(0, weight=1)
window.grid_rowconfigure(1, weight=1)
window.grid_rowconfigure(2, weight=1)
window.grid_rowconfigure(3, weight=1)
window.grid_rowconfigure(4, weight=1)
window.grid_columnconfigure(0, weight=1)

title_label = ttk.Label(window, text="Similarity Threshold Evaluation Tool", font=("Arial", 14, "bold"))
title_label.grid(row=0, column=0, columnspan=3, pady=(5, 5), sticky='n')

compound_frame = ttk.LabelFrame(window,text="Compound for evaluation", padding=10)
compound_frame.grid(row=1, column=0, padx=15, pady=(5, 10), sticky='ew')

ttk.Label(compound_frame, text="Search Compound:").grid(row=0, column=0, sticky='w')

compound_search = ttk.Entry(compound_frame)
compound_search.grid(row=0, column=1, padx=5)
compound_dropdown = ttk.Combobox(compound_frame, values=compound_list)
compound_dropdown.grid(row=0, column=2, padx=5)

def update_compound_list(event=None):
    search_term = compound_search.get()
    matched_compounds = [cpd for cpd in compound_list if search_term.lower() in cpd.lower()]
    compound_dropdown['values'] = matched_compounds

compound_search.bind("<KeyRelease>", update_compound_list)


cluster_frame = ttk.LabelFrame(window,text="Select Cluster", padding=10)
cluster_frame.grid(row=2, column=0, padx=15, pady=(10, 10), sticky='ew')

ttk.Label(cluster_frame, text="Clusters:").grid(row=0, column=0, padx=5, sticky='w')
load_button = ttk.Button(cluster_frame, text="Load Clusters", command=lambda: load_clusters())
load_button.grid(row=0, column=1, padx=5, sticky='e')

checkbox_frame = ttk.Frame(cluster_frame, padding=5)
checkbox_frame.grid(row=1, column=0, columnspan=3, padx=10, pady=5, sticky='ew')

def load_clusters():
    compound_name = compound_dropdown.get()
    if not compound_name:
        messagebox.showwarning("Warning", "Please select a compound.")
        return

    template_name = cpd_info[cpd_info['ENG'] == compound_name].iloc[0]['template_name']
    filepath_aim = f'./templates/{template_name}.csv'
    
    if not os.path.exists(filepath_aim):
        messagebox.showerror("Error", f"Template file not found for {compound_name}.")
        return

    global cpd_aim
    cpd_aim = pd.read_csv(filepath_aim)
    clusters = cpd_aim['cluster'].dropna().unique()
    clusters.sort()

    for widget in checkbox_frame.winfo_children():
        widget.destroy()

    global clusters_selected
    clusters_selected = []

    columns_per_row = 3  
    row, col = 0, 0

    for cx in clusters:
        var = tk.BooleanVar(value=(cx in ['c1', 'c2']))
        chk = ttk.Checkbutton(checkbox_frame, text=cx, variable=var)
        chk.grid(row=row, column=col, padx=5, pady=5, sticky='w')
        clusters_selected.append((cx, var))
        
        col += 1
        if col >= columns_per_row:
            col = 0
            row += 1


file_frame = ttk.LabelFrame(window,text="Select Sample File", padding=10)
file_frame.grid(row=3, column=0, padx=15, pady=(10, 10), sticky='ew')
ttk.Label(file_frame, text="Sample File:").grid(row=0, column=0, sticky='w')

file_label = ttk.Label(file_frame, text="")
file_label.grid(row=0, column=1, padx=5)  # Show selected file name

def select_file():
    global file_sample
    file_sample = filedialog.askopenfilename(initialdir='./', filetypes=[("CSV files", "*.csv")])
    if file_sample:
        file_label.config(text=os.path.basename(file_sample))

file_button = ttk.Button(file_frame, text="Browse", command=select_file)
file_button.grid(row=0, column=2, padx=5)

def start_analysis():
    compound_name = compound_dropdown.get()
    if not compound_name:
        messagebox.showwarning("Warning", "Please select a compound.")
        return

    if not file_sample:
        messagebox.showwarning("Warning", "Please select a sample file.")
        return

    selected_clusters = [cx for cx, var in clusters_selected if var.get()]
    if not selected_clusters:
        messagebox.showwarning("Warning", "Please select at least one cluster.")
        return

    df_sample = pd.read_csv(file_sample, names=['shift', 'intensity']).dropna().reset_index(drop=True)
    run_analysis(compound_name, df_sample, selected_clusters)

button_frame = ttk.Frame(window,padding=10)
button_frame.grid(row=4, column=0, padx=15, pady=(10, 15), sticky='ew')
button_frame.columnconfigure(0, weight=1)
button_frame.rowconfigure(0, weight=1)

start_button = ttk.Button(button_frame, text="Start Analysis", command=start_analysis)
start_button.grid(row=0, column=0, padx=5, pady=5, sticky='nsew')

def run_analysis(compound_name, df_sample, clusters):

    mask = (0.2 <= df_sample.iloc[:, 0]) & (df_sample.iloc[:, 0] <= 0.6)
    data_bl = df_sample.loc[mask]
    max_ = np.max(data_bl.iloc[:, -1])
    min_ = np.min(data_bl.iloc[:, -1])
    bl = max_ - min_

    rw = 0.02  # Adjustable as needed
    threshold = -1

    marker_records = []

    for cx in clusters:
        clu_data = cpd_aim[cpd_aim['cluster'] == cx]
        step = (df_sample['shift'].max() - df_sample['shift'].min()) / (df_sample.shape[0] - 1)
        pw, w_info = data_preparing(clu_data, df_sample_=df_sample, step_=step, bl_=bl, rw_=rw)
        window_aim = w_info[2]

        rolling = window_aim.rolling(len(w_info[0]), step=3)
        out_r = searching(rolling=rolling, clu_sy=w_info[0], thres=threshold)
        if len(out_r) == 0:
            continue

        out_r_max = max(out_r, key=lambda x: x[0])
        fit_range = window_aim.loc[out_r_max[1]:out_r_max[-1]+1, :]

        peak_idx = fit_range['intensity'].idxmax()
        marker_records.append({
            'cluster': cx,
            'shift': float(fit_range.loc[peak_idx, 'shift']),
            'intensity': float(fit_range.loc[peak_idx, 'intensity']),
            'pi': float(out_r_max[0])
        })

    plot_lumin_style_figure(compound_name, df_sample, marker_records)


def plot_lumin_style_figure(compound_name, df_sample, marker_records):
    """Lumin-style visualization: full-spectrum background + cluster markers."""
    fig = go.Figure()

    fig.add_trace(go.Scatter(
        x=df_sample['shift'],
        y=df_sample['intensity'],
        mode='lines',
        name='Spectrum',
        line=dict(color='#696969', width=1.5),
        showlegend=False
    ))

    if marker_records:
        colors = px.colors.qualitative.Plotly + px.colors.qualitative.Alphabet
        color_map = {}

        for i, item in enumerate(marker_records):
            cx = item['cluster']
            if cx not in color_map:
                color_map[cx] = colors[len(color_map) % len(colors)]

            fig.add_trace(go.Scatter(
                x=[item['shift']],
                y=[item['intensity']],
                mode='markers',
                name=f"{compound_name} - {cx}",
                marker=dict(
                    color=color_map[cx],
                    size=10,
                    line=dict(width=2, color='#c80000')
                ),
                hoverinfo='text',
                text=[f"<b>{compound_name}</b><br>Cluster: {cx}<br>Peak: {item['shift']:.4f} ppm<br>pi: {item['pi']:.4f}"]
            ))

    fig.update_layout(
        title_text=f'Lumin-style Visualization - {compound_name}',
        xaxis_title='Chemical Shift (ppm)',
        yaxis_title='Intensity',
        xaxis=dict(range=[10.0, -0.5], zeroline=False, ticks='inside', linewidth=1.5, tickwidth=1.5),
        yaxis=dict(zeroline=False, ticks='', showticklabels=False, linewidth=1.5),
        template='simple_white',
        height=700,
        width=900,
        showlegend=True
    )

    fig.show(renderer='browser')

# def plot_combined_figure(compound_name, df_out_r, df_sample, window_aim, w_info, out_r_max, cx):
#     fig = make_subplots(rows=2, cols=1, subplot_titles=[f'Rolling Window Similarity Curve of {cx}', 
#                                                         f'Sample Spectrum with Template Overlay of {cx}'])

#     x_values = df_sample['shift'].iloc[df_out_r['rolling index [-1]']]
#     pi_values = df_out_r['pi']

#     fig.add_trace(
#         go.Scatter(x=x_values, y=pi_values, mode='lines+markers', name='Similarity Index'),
#         row=1, col=1
#     )
#     fig.update_yaxes(title_text='Similarity Index', row=1, col=1)
#     fig.update_xaxes(title_text='Chemical Shift of W (ppm)', row=1, col=1)

#     start_idx_r = int(out_r_max[1])
#     end_idx_r = int(out_r_max[2])
#     window_rolling = window_aim.loc[start_idx_r:end_idx_r]

#     template_curve = w_info[0]
#     template_curve_scaled = template_curve * (window_rolling['intensity'].max() / template_curve.max())

#     fig.add_trace(
#         go.Scatter(x=window_aim['shift'], y=window_aim['intensity'], mode='lines', name='Sample Spectrum'),
#         row=2, col=1
#     )
#     fig.add_trace(
#         go.Scatter(x=window_rolling['shift'], y=template_curve_scaled, mode='lines', 
#                    name='Template Curve (Scaled)', line=dict(dash='dash')),
#         row=2, col=1
#     )
#     fig.update_yaxes(showticklabels=False, title_text='Intensity', row=2, col=1)
#     fig.update_xaxes(title_text='Chemical Shift of W (ppm)', row=2, col=1)

#     fig.update_layout(
#         height=800, width=1000, 
#         title_text=f'{compound_name} - {cx} Analysis',
#         template='simple_white',
#         showlegend=True,
#         legend=dict(
#             orientation='h',
#             y=-0.2,
#             x=0.5,
#             xanchor='center',
#             yanchor='top'
#         )
#     )

#     fig.show(renderer='browser')

###############################################
# def plot_combined_figure(compound_name, df_out_r, df_sample, window_aim, w_info, out_r_max, cx):
    # fig = make_subplots(
    #     rows=2,
    #     cols=1,
    #     subplot_titles=[f'Rolling Window Similarity Curve of {cx}',
    #                     f'Sample Spectrum with Template Overlay of {cx}']
    # )

    # x_values = df_sample['shift'].iloc[df_out_r['rolling index [-1]']]
    # pi_values = df_out_r['pi']

    # fig.add_trace(
    #     go.Scatter(x=x_values, y=pi_values, mode='lines+markers', name='Similarity Index'),
    #     row=1, col=1
    # )
    # fig.update_yaxes(
    #     title_text='Similarity Index',
    #     row=1,
    #     col=1,
    #     title_font=dict(size=12),
    #     tickfont=dict(size=12)
    # )
    # fig.update_xaxes(
    #     title_text='Chemical Shift of W (ppm)',
    #     row=1,
    #     col=1,
    #     title_font=dict(size=12),
    #     tickfont=dict(size=12)
    # )

    # start_idx_r = int(out_r_max[1])
    # end_idx_r = int(out_r_max[2])
    # window_rolling = window_aim.loc[start_idx_r:end_idx_r]
    # template_curve_scaled = w_info[0] * (window_rolling['intensity'].max() / w_info[0].max())

    # fig.add_trace(
    #     go.Scatter(x=window_aim['shift'], y=window_aim['intensity'], mode='lines', name='Sample Spectrum'),
    #     row=2, col=1
    # )
    # fig.add_trace(
    #     go.Scatter(x=window_rolling['shift'], y=template_curve_scaled, mode='lines',
    #              name='Template Curve (Scaled)', line=dict(dash='dash')),
    #     row=2, col=1
    # )
    # fig.update_yaxes(
    #     showticklabels=False,
    #     title_text='Intensity',
    #     row=2,
    #     col=1,
    #     title_font=dict(size=12)
    # )
    # fig.update_xaxes(
    #     title_text='Chemical Shift of W (ppm)',
    #     row=2,
    #     col=1,
    #     title_font=dict(size=12),
    #     tickfont=dict(size=12)
    # )

    # fig.update_layout(
    #     height=700,
    #     width=550,
    #     title_text=f'{compound_name} - {cx} Analysis',
    #     title_font=dict(size=14),
    #     font=dict(size=12),
    #     template='simple_white',
    #     showlegend=True,
    #     legend=dict(
    #         font=dict(size=12),
    #         orientation='h',
    #         y=-0.2,
    #         x=0.5,
    #         xanchor='center',
    #         yanchor='top'
    #     ),
    #      margin=dict(l=50, r=50, t=80, b=50)
    # )

    # fig.update_annotations(font_size=14)

    # fig.show(renderer='browser')

##############################################

window.mainloop()

