# This script visualizes fitting results.
import plotly.graph_objects as go
import numpy as np
import pandas as pd

import visualize

# cpd_aim_ls = ['121521-90-2','76822-21-4','1135-24-6','139-85-5']  # Newly added in GXNI
cpdname = '7004-03-7'
samplename = 'G10-01-2S_frac00'

filename = f'./fitting_result_{cpdname}_{samplename}.hs'

#############

cx, window_aim, fit_range, fitting_cx = visualize.load_fitting_results(filename)

fig = visualize.plot_fitted_curve(fit_range,fitting_cx,w_range=window_aim)
fig.show()

fig.write_html(f'./fitting_result_{cpdname}_{samplename}.html')

for i in range(len(fitting_cx)-2):
    peak_key = f'peak{i}'
    peak_fit = fitting_cx[peak_key]['peak_params']
    print(f'Peak {i} Parameters: {peak_fit}')

