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

from pyhipp.core import DataDict
from pyhipp.stats.random import Rng
from pyhipp.stats.summary import Summary
from pyhipp.stats.binning import _BinnedData
import numpy as np
from ..samples.field import TidalField
from ..samples.galaxy import GalaxySample

class GalaxyWebCross:
    '''
    Compute the cross-correlation function of the galaxies with the cosmic web
    (i.e. field points classified by the tidal field as given types).
    Two fields: 'x_sim_cor' (coordinates in cMpc/h in the simulation box) and 
    'weight' (weights of the galaxies) are used. If 'weight' is not provided, 
    the galaxies are equally weighted. 
    
    @lr_range: range of log10(r) [Mpc/h]. 
    @n_bins: number of radial bins.
    @rng: random number generator seed
    '''

    def __init__(self, g_samp: GalaxySample,
                 tf: TidalField, *,
                 rng: Rng, lr_range=[-0.1, 1.2], n_bins=8,
                 ):
        if 'weight' in g_samp:
            ws = g_samp['weight']
            ws = ws / ws.sum()
        else:
            n_xs = g_samp.n_objs
            ws = np.ones(n_xs) / n_xs

        self.xs = g_samp['x_sim_cor']
        self.ws = ws
        self.tf = tf
        self.rng = Rng(rng)
        self.set_bin(lr_range, n_bins)

    def profiles(self, n_bootstrap=100, lam_off=0.):
        tf = self.tf
        lams = tf.lams[tf.recon_mask]
        n_lam = (lams >= lam_off).sum(1)
        f_webs = []
        for i_web in range(4):
            vals = (n_lam == i_web).astype(float)
            p = self._v_profiles(vals)
            f_web = self._bootstrap(p, n_bootstrap)
            f_webs.append(f_web)

        rho = tf.delta[tf.recon_mask] + 1.
        p = self._v_profiles(rho)
        den_prof = self._bootstrap(p, n_bootstrap)

        return DataDict({
            'f_webs': f_webs, 'den_prof': den_prof, 'lr': self.lrs,
        })

    def _v_profiles(self, vals: np.ndarray):
        assert len(vals) == len(self.tf.xs_recon)
        cnts, nn_ids, bin_ids = self.cnts, self.nn_ids, self.bin_ids
        b_data = self.b_data
        tot_wgts = np.zeros_like(cnts)
        for i, (nn_id, bin_id) in enumerate(zip(nn_ids, bin_ids)):
            b_data.reset()
            b_data.add_n_chked(bin_id, vals[nn_id])
            tot_wgts[i] = b_data.data
        avg_wgts = tot_wgts / cnts.clip(1.0e-2)
        return avg_wgts

    def corrs(self, lam_off=0., mass_weighted=True):
        '''
        Return the correlation functions for every points with four types of
        cosmic-web points.
        
        @mass_weighted: whether to weight the field point by mass.
        '''
        tf = self.tf
        n_lam = (tf.lams[tf.recon_mask] >= lam_off).sum(1)
        rho = tf.delta[tf.recon_mask] + 1.
        corrs = []
        for i_web in range(4):
            wgts = (n_lam == i_web).astype(float)
            if mass_weighted:
                wgts *= rho
            corr = self._w_corr(wgts)
            corrs.append(corr)
        return DataDict({
            'web_typed': corrs
        })

    def _w_corr(self, wgts: np.ndarray):
        '''
        @wgts: e.g. mass or density.
        Return the correlation functions for every points.
        '''
        assert len(wgts) == len(self.tf.xs_recon)
        nn_ids, bin_ids = self.nn_ids, self.bin_ids
        b_data = self.b_data
        tot_wgts = np.zeros_like(self.cnts)
        for i, (nn_id, bin_id) in enumerate(zip(nn_ids, bin_ids)):
            b_data.reset()
            b_data.add_n_chked(bin_id, wgts[nn_id])
            tot_wgts[i] = b_data.data

        mesh = self.tf.mesh
        V_cell = mesh.l_grid**3
        rho_mean = wgts.sum() / (len(wgts) * V_cell)
        wgt_exp = rho_mean * self.dVs
        corr = tot_wgts / wgt_exp.clip(1.0e-6) - 1.
        return corr

    def bootstrap(self, outs: np.ndarray, n: int):
        ws = self.ws
        means = np.zeros((n, outs.shape[1]))
        for i in range(n):
            ids = self.rng.choice(len(outs), len(outs), p=ws)
            means[i] = outs[ids].mean(0)
        return Summary.on(means)

    def _bootstrap(self, outs: np.ndarray, n: int):
        ws = self.ws
        if n == 0:
            mean = (outs * ws[:, None]).sum(0)
            sd = np.zeros_like(mean)
            return mean, sd

        means = np.zeros((n, outs.shape[1]))
        for i in range(n):
            ids = self.rng.choice(len(outs), len(outs), p=ws)
            means[i] = outs[ids].mean(0)
        mean, sd = means.mean(0), means.std(0)
        return mean, sd

    def set_bin(self, lr_range, n_bins):
        lr_min, lr_max = lr_range
        dlr = (lr_max - lr_min) / n_bins
        b_data = _BinnedData(n_bins + 1)

        kdt = self.tf.kdt_recon
        xs = self.xs
        nn_ids, nn_ds = kdt.query_radius(
            xs, r=10**lr_max, return_distance=True)

        bin_ids = []
        for nn_d in nn_ds:
            bin_id: np.ndarray = np.floor(
                (np.log10(nn_d + 1.0e-10) - lr_min) / dlr).astype(int)
            bin_id += 1                  # fill lower-OOB to 0-th bin
            bin_id.clip(0, out=bin_id)
            bin_ids.append(bin_id)

        cnts = np.zeros((len(xs), n_bins+1), dtype=float)
        for i, bin_id in enumerate(bin_ids):
            b_data.reset()
            b_data.cnt_n_chked(bin_id)
            cnts[i] = b_data.data

        lr_es = np.linspace(lr_min, lr_max, n_bins + 1)
        lr_es = np.concatenate([[-10.], lr_es])
        r_es = 10**lr_es
        V_es = 4.0 / 3.0 * np.pi * r_es**3
        dVs = np.diff(V_es)
        lrs = 0.5 * (lr_es[:-1] + lr_es[1:])
        lrs[0] = lr_min - dlr

        self.lr_min, self.lr_max = lr_min, lr_max
        self.dlr = dlr
        self.lrs = lrs
        self.lr_es = lr_es
        self.dVs = dVs
        self.n_bins = n_bins
        self.b_data = b_data

        self.nn_ids = nn_ids
        self.nn_ds = nn_ds
        self.bin_ids = bin_ids
        self.cnts = cnts