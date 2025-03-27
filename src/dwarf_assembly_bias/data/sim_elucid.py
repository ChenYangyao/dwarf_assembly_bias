from __future__ import annotations
import typing
from typing import Self
from pyhipp.core import DataTable
from pyhipp_sims import sims
from pyhipp.io import h5
from elucid.geometry.frame import SimFrameSdssL500
from elucid.geometry.mask import ReconstAreaMaskSdssL500
from .. import sample
from ..sample_dumper import SampleTngDark, ExtraPropTngDark
from ..config import ProjPaths
import gc
import numpy as np

def create_for_cell(info: sims.SimInfo, cell_id, m_lb_in_sol = 0.):
    info_s = info.get_sub(cell_id)
    tr_ld = sims.TreeLoaderElucidExtV2(info_s)

    u_m = 1.0e10 / info_s.cosmology.hubble
    m_lb = m_lb_in_sol / u_m

    spin_dset = SampleTngDark(tr_ld, only_c=False, m_sub_lb=m_lb)
    spin_dset.run()
    path = ProjPaths.sim_dir_of(info) / f'cells/z0_cell{cell_id}.hdf5'
    with h5.File(path, 'w') as f:
        spin_dset.dump_to(f)
            
    with h5.File(path, 'a') as f:
        spin_extr = ExtraPropTngDark(tr_ld, f, f['objs'])
        spin_extr.load_tree_props(['m_crit200', 'm_mean200',
                        'spin', 'm_sub', 'm_tophat', 'x', 'v', 'subhalo_id'])
        spin_extr.find_ftimes()
        spin_extr.find_peak()
        spin_extr.find_infall()
        tr_ld.cache.clear(); gc.collect()
        
# Combined Datasets

def take_cell(sim_info: sims.SimInfo, cell_id, m_lb):
    
    s_all = sample.SimSample.load(sim_info, f'cells/z0_cell{cell_id}.hdf5')\
        .requires_spin().requires_assembly()
    sel_c = (s_all['m_mean200'] >= m_lb) & s_all['is_c']
    sel_s = (s_all['m_peak'] >= m_lb)
    sel = sel_c | sel_s

    x_sim, v_sim = s_all['x'][sel], s_all['v'][sel]
    l_box = sim_info.full_box_size
    x_sim[x_sim < 0.] += l_box
    x_sim[x_sim >= l_box] -= l_box
    assert np.all((x_sim >= 0.) & (x_sim < l_box))

    mask = ReconstAreaMaskSdssL500()
    is_recon = mask.is_in_reconst_area(x_sim)

    df = SimFrameSdssL500()
    x_j2k, v_j2k = df.pos_sim_to_j2k(x_sim), df.vel_sim_to_j2k(v_sim)
    x_norm = np.linalg.norm(x_j2k, axis=1, keepdims=True)
    x_unit = x_j2k / x_norm
    v_los = (v_j2k * x_unit).sum(axis=1, keepdims=True)

    H0 = 100.
    dx_los = v_los / H0                         # Mpc/h
    x_norm_rsd = x_norm + dx_los
    x_j2k_rsd = x_norm_rsd * x_unit
    ra, dec, z = df.pos_j2k_to_ra_dec_z(x_j2k_rsd)

    keys = (
        'm_peak', 'm_mean200', 'last_sat_z', 'is_c',
        'spin', 'v_max', 'x', 'v', 'z_half')
    d = DataTable({
        k: s_all[k][sel] for k in keys
    })
    d |= {'ra': ra, 'dec': dec, 'z': z, 'd_c': x_norm_rsd[:, 0],
          'is_recon': is_recon
          }

    return d


def take_cells(sim_info, n_cell, m_lb):
    ds = [take_cell(sim_info, c, m_lb)
          for c in range(n_cell)]
    d = ds[0].concat(*ds[1:])
    return d