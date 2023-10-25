""" Fitting the Peroxisome data from Sick Kids."""

import pandas as pd
import os
import matplotlib.cm as cmx
import pylab as pl
from scipy.integrate import quad
from scipy.optimize import curve_fit
from mpl_toolkits import mplot3d
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def double_exponential(time, rate_1, rate_2, a1):
    """ Double exponential in time, takes time, then the rate parameters (1 and
    2) then the amplitude parameters(1 and 2)"""
    #(rate_1, rate_2, a1, a2) = params
    return a1 / rate_1 * pl.exp(- time / rate_1) +  ( 
        (1-a1) / rate_2 * pl.exp(- time / rate_2))

def fitting_function(time, r1, r2, a1):
    """Calculates the function for the real time sequence and integrates the
    last bin to infinity"""
    fitted_dist = pl.empty(len(time))
    #print(type(params), 'fitted function')

    for i, t in enumerate(time):
        a = t
        if i < len(time) - 1:
            b = time[i+1]
        else:
            b = pl.inf
        integrated_bin = quad(func = double_exponential, a = a, b = b,
            args = (r1, r2, a1))[0]
        fitted_dist[i] = integrated_bin
    return fitted_dist

def integrated_exponential(time, rate, amplitude):
    """integrate up the contributions from either exponential"""
    fitted_dist = pl.empty(len(time))
    #print(type(params), 'fitted function')

    for i, t in enumerate(time):
        a = t
        if i < len(time) - 1:
            b = time[i+1]
        else:
            b = pl.inf
        integrated_bin = quad(func = my_exp, a = a, b = b,
            args = (rate, amplitude))[0]
        fitted_dist[i] = integrated_bin
    return fitted_dist

def my_exp(x, r, a):
    """exponential with magnitude"""
    return a / r * pl.exp(- x / r)

def format_sigfigs(value, error):
    """ return strings of the value and error with correct sig figs"""
    
    error = '{:.1g}'.format(error)
    n_sig = len(error)
    if pl.log10(value) > 1:
        n_sig += int(pl.log10(value))
    value = str(value)[:n_sig]

    return '{0} $\pm$ {1}'.format(value, error)

def import_data(file_name):
    """import the data from an excel file with filename 'file_name'"""
    if '.xlsx' not in file_name:
        file_name += '.xlsx'

    data_df = pd.read_excel(file_name, header=None)
    data = data_df.values
    return data

def plot_distributions(distributions, time):
    """ Plot the distributions for a set of trials """
    
    N = len(distributions)
    dt = time[1]/(N+1)
    
    pl.figure(figsize = (8,6))
    for i, d in enumerate(distributions):
        pl.bar(time + (i + 1/2)*dt, d, width = dt, color = 'C0')    


def plot_3d(data_dicts, sampling_types, pair = False):
    """ Function to plot the parameter spaces of various fitting thingies.
    Pair boolean is to put the controls together and the AT trials together"""
    def onpick(event):
        # on the pick event, find the orig line corresponding to the
        # legend proxy line, and toggle the visibility
        legline = event.artist
        origline = lined[legline]
        vis = not origline.get_visible()
        origline.set_visible(vis)
        # Change the alpha on the line in the legend so we can see what lines
        # have been toggled
        if vis:
            legline.set_alpha(1.0)
        else:
            legline.set_alpha(0.2)
        fig.canvas.draw()


    """
    values = pl.linspace(0, 1, len(data_dicts[0].keys()))
    cm = mpl.cm.viridis(values)
    cm = mpl.colors.ListedColormap(cm)
    c_norm = pl.Normalize(vmin=0, vmax=1)
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap = cm)
    colors = [scalar_map.to_rgba(v) for v in values]
    """
    colors = ['C{}'.format(i%10) for i in range(len(data_dicts[0].keys()))]
    lines = []
    mean_lines = []
    i = 0
    fig = pl.figure(figsize = (8,6))
    fig.canvas.mpl_connect('pick_event', onpick)
    ax = pl.axes(projection = '3d')
    ax.set_xlabel('Slow Timescale')
    ax.set_ylabel('Slow Proportion')
    ax.set_zlabel('Fast Timescale')

    for k, v in data_dicts[0].items():
        ax.set_title(k)
        c_index = i
        tao_slow = v[0]
        p_slow = v[1]
        tao_fast = v[2]
        """
        mean_line = ax.plot3D(pl.average(tao_slow), pl.average(p_slow),
            pl.average(tao_fast), marker = '*', c = colours[c_index],
            markersize = 100, lw = 0)
        """
        line, = ax.plot3D(tao_slow, p_slow, tao_fast, label = k,
            c = colors[c_index], lw = 0, marker = 'o')
        
        lines.append(line)
        #lines.append(mean_line)
        i += 1

    leg = ax.legend(loc='upper left', fancybox=False, shadow=False)
    leg.get_frame().set_alpha(0.4)

    # we will set up a dict mapping legend line to orig line, and enable
    # picking on the legend line
    lined = dict()
    for legline, origline in zip(leg.get_lines(), lines):
        legline.set_picker(5)  # 5 pts tolerance
        lined[legline] = origline

def print_ranked_params(all_params, order, trial_name, output_dir):
    """ just spew parameter values for cells ranked by their tao slow,
    largest to smallest"""

    alphabet = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    ranked_params = pl.empty([5, len(order)], dtype = '<U5')
    columns = ['Cell Number', 'Excel Column', 'Slow Timescale',
        'Slow Proportion', 'Fast Timescale']
    for i in range(len(order)):
        
        letter_index = int(order[i]%26)
        letter = alphabet[letter_index]
        
        if order[i] >= 26:
            letter = 'A' + letter

        print("Cell {0} (Excel {4}):\tSlow Timescale {1:.1f}, Slow Proportion\
 {2:.1f}%, Fast Timescale {3:.1f}".format(order[i],
                all_params[0][i],
                all_params[1][i]*100,
                all_params[2][i],
                letter
                ))

        ranked_params[:,i] = [order[i],
                letter,
                all_params[0][i],
                all_params[1][i]*100,
                all_params[2][i]
                ]
    
    data_frame = pd.DataFrame(ranked_params.T, columns = columns)
    data_frame.to_csv('{0}/{1}.csv'.format(output_dir, trial_name))



if __name__ == '__main__':
    pl.close('all')
    seed = 2
    pl.seed(seed)
    time_scale = 5.203
    time = pl.arange(10) * time_scale
    bounds = ([0.001, 0.001, 0.00], [10000, 100, 1])
    initial = [20, 5, 0.5]
 
    # Parameter holding stuff
    params_dicts = []
    individual_fit_dict = {}
    sampling_types = []
    sampling_types.append('Individual Trials')


    plot_all_3d = True
    plot_fits = True
    plot_trial_distributions = True

    # Data is pulled from the Data folder
    data_dir = 'Data/'
    files = sorted(os.listdir(data_dir))

    # Output directory for rank ordered parameters
    output_dir = 'Parameters/'

    # could also just pull specific ones by name from that folder (uncomment)
    #files = ['1_siCNTRL.xslx']

    # Add directory for loading.
    files = [data_dir + f for f in files]
    names = [f.strip('.xlsx').strip(data_dir) for f in files]

    for i, f in enumerate(files):
        print("loading: ", f)
        data = import_data(f)#pl.loadtxt(f)
        name = names[i]
        n_samp = data.shape[1]
        params =  [] #pl.empty([n_samp], dtype = pl.ndarray)
        fast_params = [] #pl.empty([n_samp], dtype = pl.ndarray)
        fast_counts = pl.empty([n_samp], dtype = pl.ndarray)
        slow_params = [] #pl.empty([n_samp], dtype = pl.ndarray)
        slow_counts = pl.empty([n_samp], dtype = pl.ndarray)
        fit_counts = pl.empty([n_samp], dtype = pl.ndarray)

        distributions = pl.empty([n_samp], dtype = pl.ndarray)

        for n in range(n_samp):
            """
            dwell_times = data[n,:]
            dwell_times = dwell_times[~pl.isnan(dwell_times)] * time_scale
            counts, bins = pl.histogram(dwell_times)
            dist = counts/sum(counts)
            """
            dist = data[:,n]
            dist = dist/pl.sum(dist)
            distributions[n] = dist
            res = curve_fit(fitting_function, time, dist,
                p0 = initial, bounds = bounds)
            r1, r2, a1 = res[0]
            a2 = 1 - a1
            fit_counts[n] = fitting_function(time, r1, r2, a1)
            if r1 > r2:
                slow_params.append([r1, a1])
                fast_params.append([r2, a2])
            else:
                fast_params.append([r1, a1])
                slow_params.append([r2, a2])
            fast_counts[n] = integrated_exponential(time, fast_params[n][0],
                fast_params[n][1])
            slow_counts[n] = integrated_exponential(time, slow_params[n][0],
                slow_params[n][1])
            params.append(res[0])
        slow_param_means = pl.average(slow_params, axis = 0)
        slow_param_medians = pl.median(slow_params, axis = 0)
        slow_param_errors = pl.std(slow_params, axis = 0)/pl.sqrt(n_samp - 1)
        
        print(n_samp)
        print(slow_param_means, ' +- ', slow_param_errors)
        fast_param_means = pl.average(fast_params, axis = 0)
        fast_param_medians = pl.median(fast_params, axis = 0)
        fast_param_errors = pl.std(fast_params, axis = 0)/pl.sqrt(n_samp - 1)
        print(fast_param_means, ' +- ', fast_param_errors)
        
        slow_means = pl.average(slow_counts, axis = 0)
        slow_errors = pl.std(slow_counts, axis = 0)/pl.sqrt(n_samp-1)
        
        fast_means = pl.average(fast_counts, axis = 0)
        fast_errors = pl.std(fast_counts, axis = 0)/pl.sqrt(n_samp-1)
        
        fit_means = pl.average(fit_counts, axis = 0)
        fit_errors = pl.std(fit_counts, axis = 0)/pl.sqrt(n_samp-1)
 
        #print(res)
        mean_dist = pl.average(distributions, axis = 0)
        dist_errors = pl.std(distributions, axis = 0)/(n_samp-1)
        if plot_fits:
            pl.figure(figsize = (8,6))
            pl.bar(time + time_scale/2, mean_dist, width = time_scale,
                label = 'Binned Data', yerr = dist_errors, capsize = 5,
                ecolor = 'C0')
            pl.errorbar(time + time_scale/2, fit_means, fmt = 'ko',
                label = 'Fitted Curve', yerr = fit_errors, capsize = 3.5)
            pl.bar(time + time_scale/4, fast_means, color = 'C2',
                width = time_scale/2, yerr = fast_errors,
                label = 'Fast Rate Contribution')
            pl.bar(time + 3*time_scale/4, slow_means, color = 'C1',
                width = time_scale/2, yerr = slow_errors,
                label = 'Slow Rate Contribution')

            pl.annotate(('Slow Timescale: ' + 
                format_sigfigs(slow_param_means[0], slow_param_errors[0]) + 
                's') + ', Median: {:.1f}'.format(slow_param_medians[0]),
                xy = (20, 0.2))
            pl.annotate(('Slow Proportion: ' + 
                format_sigfigs(slow_param_means[1]*100,
                slow_param_errors[1]*100) + '%') + ', Median: {:.1f}'.format(
                slow_param_medians[1]*100),
                xy = (20, 0.175))
            pl.annotate(('Fast Timescale: ' + 
                format_sigfigs(fast_param_means[0], fast_param_errors[0]) + 
                's') + ', Median: {:.1f}'.format(fast_param_medians[0]),
                xy = (20, 0.15))
            pl.annotate(('Fast Proportion: ' + 
                format_sigfigs(fast_param_means[1]*100,
                fast_param_errors[1]*100) + '%') + ', Median: {:.1f}'.format(
                fast_param_medians[1]*100),
                xy = (20, 0.125))


            pl.legend()
            pl.xlabel('time (s)')
            pl.ylabel('Average % Contact Events')
            pl.title('{0} Individual Trial Fits'.format(names[i]))
            p_dir = 'WeirdPlots/'
            pl.savefig(p_dir+'{0}ByTrial.pdf'.format(names[i]))

        tao_slow = pl.asarray([l[0] for l in slow_params])
        p_slow = pl.asarray([l[1] for l in slow_params])
        tao_fast = pl.asarray([l[0] for l in fast_params])

        order = (-tao_slow).argsort()

        all_params = [tao_slow, p_slow, tao_fast]
        all_params = [a[order] for a in all_params]
        individual_fit_dict[names[i]]  = all_params
        print_ranked_params(all_params, order, name, output_dir)

        if plot_trial_distributions:
            plot_distributions(distributions[order], time)
            pl.title(names[i])
            pl.xlabel('Time (s)')
            pl.ylabel('Occupation')


        #params_array = 

    params_dicts.append(individual_fit_dict)

    if plot_all_3d:
        plot_3d(params_dicts, sampling_types)

    

    pl.show()
