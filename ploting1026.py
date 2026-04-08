import plotly.graph_objects as go
import numpy as np
import pandas as pd

import visualize

cpdname = '50-21-5'
samplename = 'CHO-bl_frac00'

filename = f'./fitting_result_{cpdname}_{samplename}.hs'


#############

cx, window_aim, fit_range, fitting_cx = visualize.load_fitting_results(filename)

fig = visualize.plot_fitted_curve(fit_range,fitting_cx,w_range=window_aim)
fig.show()

# fig.write_html(f'./visual-GXNI/frac005-GXNI/fitting_result_{cpdname}_GXN-KB-bl.html')
fig.write_html(f'./fitting_result_{cpdname}_{samplename}.html')
###############
for i in range(len(fitting_cx)-2):
    peak_key = f'peak{i}'
    peak_fit = fitting_cx[peak_key]['peak_params']
    print(f'Peak {i} Parameters: {peak_fit}')

