import plotly.graph_objects as go
import numpy as np
import pandas as pd
import h5py


def save_fitting_results(filename, cx, window_aim, fit_range, fitting_cx):
    """
    Save one fitting result package to an HDF5 file.

    Parameters
    filename : str
        Output .hs/.h5 file path.
    cx : str
        Cluster name.
    window_aim : pd.DataFrame
        Candidate sample window used in searching.
    fit_range : pd.DataFrame
        Final fitting range.
    fitting_cx : dict
        Fitting result dict from line_fitting().

    Stored structure
    {'cx': cx, 'Window aim': window_aim, 'Fit range': fit_range, 'Fit result': fitting_cx}
    """
    with h5py.File(filename, 'w') as f:
        f.attrs['cx'] = cx
        # Save DataFrames
        window_grp = f.create_group('Window aim')
        for col in window_aim.columns:
            window_grp.create_dataset(col, data=window_aim[col].values)
        
        range_grp = f.create_group('Fit range')
        for col in fit_range.columns:
            range_grp.create_dataset(col, data=fit_range[col].values)

        # Save fitting_cx dictionary
        fit_grp = f.create_group('Fit result')
        for key, data in fitting_cx.items():
            if key.startswith('peak'):
                peak_grp = fit_grp.create_group(key)
                peak_grp.create_dataset('peak_y', data=data['peak_y'])
                peak_grp.create_dataset('peak_params', data=data['peak_params'])
            elif key == 'residual':
                fit_grp.attrs['residual'] = data
            elif key == 'residue_curve':
                fit_grp.create_dataset('residue_curve', data=data)



def load_fitting_results(filename):
    """
    Load one fitting result package from an HDF5 file.
    """    
    with h5py.File(filename, 'r') as f:
        cx = f.attrs['cx']
        
        window_aim = pd.DataFrame({col: f['Window aim'][col][:] for col in f['Window aim'].keys()})
        fit_range = pd.DataFrame({col: f['Fit range'][col][:] for col in f['Fit range'].keys()})

        fitting_cx = {}
        fit_grp = f['Fit result']
        for key in fit_grp.keys():
            if key.startswith('peak'):
                fitting_cx[key] = {
                    'peak_y': fit_grp[key]['peak_y'][:],
                    'peak_params': fit_grp[key]['peak_params'][:]
                }
        fitting_cx['residual'] = fit_grp.attrs['residual']
        fitting_cx['residue_curve'] = fit_grp['residue_curve'][:]

    return cx, window_aim, fit_range, fitting_cx



def plot_fitted_curve(fit_range,rs_fitting,w_range=None):
    """
    Plot fitted peaks, residue curve, and summed fit.

    Parameters
    fit_range : pd.DataFrame
        DataFrame containing fitting-range spectrum.
    rs_fitting : dict
        Result dictionary from line_fitting().
    w_range : pd.DataFrame, optional
        Original window spectrum for background display.
    """

    fig = go.Figure()

    if w_range is None:
        fig.add_trace(go.Scatter(x=fit_range['shift'].values, y=fit_range['intensity'].values,
                                 mode='lines', name='Original spectrum of fitting range'))
    else:
        fig.add_trace(go.Scatter(x=w_range['shift'].values, y=w_range['intensity'].values, 
                                 mode='lines', name='Original spectrum of Window aim'))

    x_ppm = fit_range['shift'].values
    residue_curve = rs_fitting['residue_curve']
    sum_fitted = np.zeros_like(fit_range['intensity'].values)
    for i in range(len(rs_fitting)-2):
        peak_key = f'peak{i}'
        peak_fit = rs_fitting[peak_key]['peak_y']
        sum_fitted += peak_fit
        fig.add_trace(go.Scatter(x=x_ppm, y=peak_fit, mode='lines', name=f'Peak {i} Fit'))

    fig.add_trace(go.Scatter(x=x_ppm, y=residue_curve, mode='lines', name='Residue Curve'))
    fig.add_trace(go.Scatter(x=x_ppm, y=sum_fitted, mode='lines', name = 'Sum of fitted'))

    fig.update_layout(title='Fitting Results',
                    template='simple_white',
                    xaxis_title='Chemical Shift (ppm)',
                    yaxis_title='Intensity',
                        legend=dict(
                        orientation='h',
                        y=-0.2,
                        x=0.5,
                        xanchor='center',
                        yanchor='top' )
                    )

    return fig




