# Copyright (C) 2026 Yangyao Chen (yangyaochen.astro@foxmail.com) - All Rights 
# Reserved
# 
# You may use, distribute and modify this code under the MIT license. We kindly
# request you to give credit to the original author(s) of this code, and cite 
# the following paper(s) if you use this code in your research: 
# - Chen, Y. and Wang K. 2023, ASCL:2301.030, NASA/ADS entry: 
#   https://ui.adsabs.harvard.edu/abs/2023ascl.soft01030C/abstract 
#   (the pyhipp package).
# - Zhang Z. et al. 2025. Nature 642, 47-52 (for the computation of relative 
#   bias).
# - Wang H. et al. 2016. ApJ, 831, 164 (for the reconstructed density field 
#   and halo catalogs).
#
# Usage:
# python bias_curve.py
#
# This takes input files associated with the this script, computes the bias 
# curve, saves the result to a JSON file, and makes a plot comparing with 
# observational data from Zhang et al. 2025.

import numpy as np
from pyhipp.io import h5, json
from pyhipp import plot
from pyhipp.astro.stats import ccf
from pathlib import Path

### Options ###
l_box = 500.0                                           # [cMpc/h]
data_dir = Path('./data/model/')                        # where input files are stored
path_ref = data_dir/'z0_ref.hdf5'                       # reference sample
path_dst = data_dir/'z0_massive_dwarfs.hdf5'            # target sample
path_Sigma = data_dir/'z0_subhalos_fountain_sigma.csv'  # Sigma file (containing stellar surface density)
skiprows = 1                                            # number of rows to skip in the Sigma file
usecol = 23                                             # index of the column to use in the Sigma file

ccf_kwargs = {
    'n_threads': 4,                                     # number of threads for parallel computation
    'n_bootstrap': 25,                                  # number of bootstrap resamplings for error estimation
    'pi_max': 40.,                                      # maximum line-of-sight separation for projected correlation function [in cMpc/h]
}
bias_kwargs = {
    'bin_edges': [0., 7., 15., 25., 10000],             # bin edges for Sigma [in the same units as the input Sigma file]
    'ref_bin': -1,                                      # index of the bin of Sigma to normalize the bias curve (-1 means the last bin)
    'r_min': 2., 'r_max': 10.                           # radial range to define the bias [in cMpc/h]
}

path_obs_data = data_dir/'rel_bias_Zhang25_Fig1b.json'  # observation, for plot
out_dir = Path("./output/")
path_out_bias_curve = out_dir/'bias_curve.json'
path_out_plot = out_dir/'bias_curve.pdf'


### Execute the analysis ###

# Load datasets
ref_sample = ccf.SimSample(l_box, h5.File.load_from(path_ref, key='objs'))
dst_sample = ccf.SimSample(l_box, h5.File.load_from(path_dst, key='objs'))
dst_sample.data['Sigma'] = np.loadtxt(
    path_Sigma, skiprows=skiprows, delimiter=',', usecols=usecol)
print(f'Data loaded.')

# Compute and save the bias curve.
executor = ccf.SimCCFProjected(ref_sample, **ccf_kwargs)
print(f'CCF executor initialized: {executor}')

bias = executor.relative_bias_curve(
    dst_sample, bin_by_key='Sigma', **bias_kwargs)
print(f'Bias curve computed: {bias}')
json.File.dump_file(bias, path_out_bias_curve, flag='w', indent=2)


# Make a plot showing the bias curve and compare with Zhang+25.
fig, ax = plot.subplots(
    1, figsize=(5.5, 5.25), margin=[0.02, 0.02, 0.11, 0.135], layout='none')

x, (x_lo, x_hi) = bias['x/median'], bias['x']['sigma_1'].T
err_x_lo, err_x_hi = x - x_lo, x_hi - x

y, (y_lo, y_hi) = bias['y/median'], bias['y']['sigma_1'].T
err_y_lo, err_y_hi = y - y_lo, y_hi - y

ax.errorbar(x, y, yerr=[err_y_lo, err_y_hi], xerr=[err_x_lo, err_x_hi],
            label=r'$\rm Model\,(This\ work)$')

obs_data = json.File.load_file(path_obs_data)
x, y, el, eh = np.array(obs_data['Fig1b']['main']).T
ax.c('orange').fmt_marker('o',).errorbar(x, y, yerr=[el, eh],
    label=r'$\rm SDSS\ (Zhang+25)$', capsize=0, lw=2)

ax.lim([-2., 52.], [0.7, 2.75])\
    .label(r'\Sigma_*\,[{\rm M}_\odot\,{\rm pc}^{-2}]', r'\text{Relative bias}')\
    .leg(loc='ur', numpoints=1, handlelength=1.2)
plot.savefig(path_out_plot)
print('Plot saved.')