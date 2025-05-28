# Copyright (C) 2025 Yangyao Chen (yangyaochen.astro@foxmail.com) - All Rights 
# Reserved
# 
# This file is part of the Dwarf Assembly Bias project. You may use, distribute 
# and modify this code under the MIT license. We kindly request you to give 
# credit to the original author(s) of this code, and cite the following paper(s) 
# if you use this code in your research: 
# - Zhang Z. et al. 2025. Nature ???, ??? (for the galaxy-cosmic web cross 
#   correlation as an environmental indicator).
# - Wang H. et al. 2016. ApJ, 831, 164 (for the reconstructed field).
# - Chen, et al. 2020. ApJ, 899, 81 (for the numerical method of field 
#   analysis).
#
# This file is similar to "galaxy_web_cross.py" but takes a observation sample
# and correct its redshift-space distortion (RSD) by matching it with the 
# reconstruction sample of galaxies.
#
# Usage:
# Suppose the data layout is as follows:
# ./
#    |- data/
#       |- gal_web_cross/
#          |- diffuse.txt                  - Galaxy sample for which the galaxy-cosmic-web 2PCCFs are computed
#          |- galaxy_recon.hdf5            - Reconstruction sample of galaxies (used to correct RSD)
#          |- domain.tidal.s99.sm1.hdf5    - Tidal field
#    |- output/                            - Directory for the output 
#
# Run this script by:
# $ python galaxy_web_cross_from_obs_sample.py
#
# The 2PCCFs and a plot will be saved to the output directory upon successful
# execution.

from __future__ import annotations
import numpy as np
from pathlib import Path
from scipy.spatial import KDTree
from pyhipp import plot
from pyhipp.io import h5
from pyhipp.core import DataDict
from dwarf_assembly_bias import samples, statistics

# Parameters
file_target_sample = './data/gal_web_cross/diffuse.txt'
file_recon_sample = './data/gal_web_cross/galaxy_recon.hdf5'
file_tidal_field = './data/gal_web_cross/domain.tidal.s99.sm1.hdf5'
dir_output = './output'

# Maximal distance tolerance in the matching of target galaxies with the 
# reconstruction sample (to correct RSD).
ra_match_max = .1 / 3600.0       # RA in degree
dec_match_max = .1 / 3600.0      # Dec in degree
z_match_max = 1.0e-4             # redshift

# The following parameters are used to control the computation of 2PCCFs
lr_range, n_bins =[0.1, 1.2], 8 # log10(r) [Mpc/h], range of radial bins, and total number of bins
rng = 10086                     # Random seed            
n_bootstrap = 100               # Number of bootstrapping
mass_weighted = True            # Whether to weight each field point by its density


def correct_rsd(target_sample: DataDict[str, np.ndarray]):
    recon_sample = h5.File.load_from(file_recon_sample)
    mask = recon_sample['is_recon_safe']
    recon_sample = DataDict({k: v[mask] for k, v in recon_sample.items()})
    
    x_recon = np.column_stack(recon_sample['ra', 'dec', 'z_obs'])
    x_target = np.column_stack(target_sample['ra', 'dec', 'z_obs'])

    scale = (ra_match_max, dec_match_max, z_match_max)
    x_recon /= scale
    x_target /= scale

    _, ids = KDTree(x_recon).query(x_target)
    ds = np.abs(x_recon[ids] - x_target)
    x_sim_cor = recon_sample['x_sim_cor'][ids]
    matched = (ds < 1.).all(1).nonzero()[0]

    print(f'Matched {len(matched)} out of {len(x_target)} galaxies.')
    print(f'Maximal distance {ds[matched].max(0)}')

    matched_sample = DataDict({k: v[matched] for k, v in target_sample.items()}) | {
        'row_id': matched,
        'x_sim_cor': x_sim_cor[matched],
    }

    return matched_sample


def find_2pccf(target_sample: samples.galaxy.GalaxySample, 
             tidal_field: samples.field.TidalField):
    gwc = statistics.cosmic_web.GalaxyWebCross(
        target_sample, tidal_field, lr_range=lr_range, n_bins=n_bins, rng=rng)
    
    lrs = gwc.lrs
    web_types = 'void', 'sheet', 'filament', 'knot'
    ccfs = gwc.corrs(mass_weighted=mass_weighted)['web_typed']
    out = {}
    for ccf, web_type in zip(ccfs, web_types):
        _ccf = gwc.bootstrap(ccf, n_bootstrap)
        med, (lo, hi) = _ccf.median, _ccf.sigma_1
        out[web_type] = np.column_stack([lrs, med, lo, hi])
    return out

# Load tidal field and target sample of galaxies; RSD is corrected by matching
# the target sample with the reconstruction sample.
tidal_field = samples.field.TidalField.from_file(file_tidal_field)
keys = 'ra', 'dec', 'z_obs', 'weight'
vals = np.loadtxt(file_target_sample, dtype=float).T
target_sample = DataDict(dict(zip(keys, vals, strict=True)))

target_sample = correct_rsd(target_sample)
target_sample = samples.galaxy.GalaxySample(target_sample)

# Find 2PCCF and save to file
header = '''Columns:
1: log(r): [h^{-1} Mpc] center of radial bin
2: xi(r): galaxy-cosmic-web cross-correlation function (median among bootstrap results)
3. 16th percentile of xi(r)
4. 84th percentile of xi(r)'''
out_2pccf = find_2pccf(target_sample, tidal_field)
for key, val in out_2pccf.items():
    file_out = f'{dir_output}/galaxy_{key}.txt'
    np.savetxt(file_out, val, header=header)    

# Make a plot and save it
fig, axs = plot.subplots((2,2), share=False, space=.3, subsize=5.5, layout='none')
axs_f = axs.flat

web_types = 'void', 'sheet', 'filament', 'knot'
for i_t, web_type in enumerate(web_types):
    lrs, y, yl, yh  = out_2pccf[web_type].T
    e_lo, e_hi = y-yl, yh-y
        
    ax = axs_f[i_t]
    ax.errorbar(lrs, y, yerr=[e_lo, e_hi], label=r'Galaxy\text{-}'+web_type)

axs.leg(loc='ur').label(r'\log\,r\,(h^{-1}{\rm Mpc})', r'\xi(r)')
plot.savefig(f'{dir_output}/galaxy_web_cross.png')