# Copyright (C) 2025 Yangyao Chen (yangyaochen.astro@foxmail.com) - All Rights 
# Reserved
# 
# You may use, distribute and modify this code under the MIT license. We kindly
# request you to give credit to the original author(s) of this code, and cite 
# the following paper(s) if you use this code in your research: 
# - Zhang Z. et al. 2025. Nature ???, ??? (for the galaxy-cosmic web cross 
#   correlation as an environmental indicator).
# - Wang H. et al. 2016. ApJ, 831, 164 (for the reconstructed field).
# - Chen, et al. 2020. ApJ, 899, 81 (for the numerical method of field 
#   analysis).

from __future__ import annotations
from pyhipp import plot
import numpy as np
from pathlib import Path
from dwarf_assembly_bias import samples, statistics

# Load data
#
# The galaxy sample here consists of compact dwarf galaxies, stored as a HDF5 
# file. Alternatively, you may pass a dict to construct `GalaxySample`directly.
#
# The tidal field is obtained from the reconstruction pipeline of ELUCID.
# This is too large to be hosted, and you may contact the author for to get it.
data_dir =  Path('data/gal_web_cross')
g_samp = samples.galaxy.GalaxySample.from_file(data_dir / 'compact.hdf5')              # Galaxy sample
tf = samples.field.TidalField.from_file(data_dir / f'domain.tidal.s99.sm1.hdf5')      # Tidal field

# Compute cross correlation
lr_range, n_bins =[0.1, 1.2], 8
rng = 10086
n_bootstrap = 100
mass_weighted = True

gwc = statistics.cosmic_web.GalaxyWebCross(
    g_samp, tf, lr_range=lr_range, n_bins=n_bins, rng=rng)
ccf = gwc.corrs(mass_weighted=mass_weighted)
ccf = [gwc.bootstrap(_ccf, n_bootstrap).as_dict() for _ccf in ccf['web_typed']]
output = {
    'lrs': gwc.lrs,
    'ccf': ccf,
}

# Make a plot. This is shown in Fig. 2 (compact dwarfs) of Zhang Z. et al. 2025,
# with some potential difference due to the platform and the random number 
# implementation.
fig, axs = plot.subplots((2,2), share=False, space=.3, subsize=5.5, layout='none')
axs_f = axs.flat

web_types = ['void', 'sheet', 'filament', 'knot']
lrs = output['lrs']
for i_t, web_type in enumerate(web_types):
    _ccf = output['ccf'][i_t] 
    y, (yl, yh) = np.array(_ccf['median']), np.array(_ccf['sigma_1'])
    e_lo, e_hi = y-yl, yh-y
        
    ax = axs_f[i_t]
    ax.errorbar(lrs, y, yerr=[e_lo, e_hi])
    
axs.label(r'\log\,r\,(h^{-1}{\rm Mpc})', r'\xi(r)')
plot.savefig('output/galaxy_web_cross.png')