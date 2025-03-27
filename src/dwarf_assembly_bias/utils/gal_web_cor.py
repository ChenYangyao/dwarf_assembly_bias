from __future__ import annotations
import typing
from typing import Self
from pyhipp.field.cubic_box import TidalClassifier
from pyhipp.core import Num, DataDict, abc
from pyhipp.stats import Rng
from elucid.geometry.mask import ReconstAreaMaskSdssL500
import numpy as np
from functools import cached_property
from scipy.spatial import KDTree
from .weights import ZWeighting
from ..obs_catalogs import ReconSample, TidalField
from ..sample import SimSample
from ..utils.stats import _BinnedData, Summary
from dataclasses import dataclass

def _DD_nonperiodic(x1, x2, rs, w1, w2, n_threads):
    
    import Corrfunc
    
    '''
    Return total pair counts and weights.
    '''
    X1, Y1, Z1 = x1.T
    X2, Y2, Z2 = x2.T
    kw = {
        'autocorr': False, 'nthreads': n_threads, 'periodic': False,
        'binfile': rs, 'weight_type': 'pair_product',
        'X1': X1, 'Y1': Y1, 'Z1': Z1,
        'X2': X2, 'Y2': Y2, 'Z2': Z2,
        'weights1': w1, 'weights2': w2,

    }
    out = Corrfunc.theory.DD(**kw)
    n, w_avg = out['npairs'].astype(float), out['weightavg'].astype(float)
    w = n * w_avg
    return n, w

class GalWebCor:
    def __init__(self,
                 samp: ReconSample,
                 tidal_clsf: TidalClassifier,
                 mask: ReconstAreaMaskSdssL500 = None,
                 zweighting: ZWeighting = None,
                 rng=0,
                 n_threads=1):

        xs = samp.data['x_sim_cor']
        lams = tidal_clsf.lam_at(xs)

        if zweighting is not None:
            ws = zweighting.weights_of(samp.data['z_cor'])
        else:
            n_xs = len(xs)
            ws = np.ones(n_xs) / n_xs
        if mask is not None:
            assert mask.mask_reconst_area.shape == tidal_clsf.lam.shape[:3]

        self.samp = samp
        self.tidal_clsf = tidal_clsf
        self.mask = mask
        self.zweighting = zweighting
        self.rng = Rng(rng)
        self.n_threads = n_threads

        self.xs = xs
        self.lams = lams
        self.ws = ws

    def d_to_webt(self, n_lam=2, lam_off=0.):
        x1, w1 = self.xs, self.ws
        x2 = self._xs_webt(n_lam, lam_off)
        d, _ = KDTree(x2).query(x1)
        return DataDict({
            'd': d, 'w': w1,
        })

    def corr_webt(self, rs=np.logspace(-.5, 1.2, 12),
                  n_lam=2, lam_off=0.,
                  n_bootstrap=0):
        rs_c = 0.5 * (rs[:-1] + rs[1:])
        lg_rs_c = np.log10(rs_c)
        out = DataDict({'lg_rs_c': lg_rs_c, 'rs': rs})
        x1, w1 = self.xs, self.ws
        x2 = self._xs_webt(n_lam, lam_off)
        x2_rand = self._xs_recon_grids
        if n_bootstrap == 0:
            xi = self._ref_cnt(rs, x1, w1, x2, x2_rand)
            return out | {'xi': xi}
        xis = []
        for _ in range(n_bootstrap):
            _x1, _w1 = self._rand(x1, w1)
            _x2 = self._rand(x2)
            _x2_rand = self._rand(x2_rand)
            xi = self._ref_cnt(rs, _x1, _w1, _x2, _x2_rand)
            xis.append(xi)
        xis = np.array(xis)
        xi, xi_med, xi_sd = xis.mean(
            axis=0), np.median(
            xis, axis=0), np.std(
            xis, axis=0)
        return out | {'xi': xi, 'xi_med': xi_med, 'xi_sd': xi_sd, 'xis': xis}

    def _ref_cnt(self, rs, x1, w1, x2, x2_rand):
        n2, n2_rand = len(x2), len(x2_rand)
        _, sum_w = _DD_nonperiodic(x1, x2, rs, w1, None, self.n_threads)
        _, sum_w_rand = _DD_nonperiodic(x1, x2_rand, rs, w1, None, self.n_threads)
        xi = (sum_w / sum_w_rand.clip(1.0e-10)) * (n2_rand / n2) - 1.
        return xi

    def _rand(self, x, *xs):
        n = len(x)
        ids = self.rng.choice(n, n)
        if len(xs) == 0:
            return x[ids]
        return tuple(_x[ids] for _x in (x,) + xs)

    @cached_property
    def _xs_all_grids(self):
        mesh = self.tidal_clsf.mesh
        n, h = mesh.n_grids, mesh.l_grid

        idx_1d = np.arange(n)
        idx_3d = np.meshgrid(idx_1d, idx_1d, idx_1d, indexing='ij')
        x = np.stack([idx * h for idx in idx_3d], axis=-1)
        return x

    @cached_property
    def _xs_recon_grids(self):
        xs = self._xs_all_grids[self.mask.mask_reconst_area]
        return xs.reshape(-1, 3)

    def _xs_webt(self, n_lam, lam_off):
        sel = (self.tidal_clsf.lam >= lam_off).sum(-1) == n_lam
        sel &= self.mask.mask_reconst_area
        xs = self._xs_all_grids[sel]
        return xs.reshape(-1, 3)

    def frac_webt(self, lam_off: float | np.ndarray):
        if not np.isscalar(lam_off):
            return np.array([self.frac_webt(o) for o in lam_off]).T
        ws = self.ws
        n_lam = (self.lams >= lam_off).sum(1)
        fracs = []
        for n in (0, 1, 2, 3):
            sel = n_lam == n
            frac = ws[sel].sum() / ws.sum()
            fracs.append(frac)
        return fracs

    def G_f(self, lam_off=0., norm=True):
        signs = np.array([[-1., 1., 1.]])
        return self._G(lam_off, signs, norm)

    def G_s(self, lam_off=0., norm=True):
        signs = np.array([[-1., -1., 1.]])
        return self._G(lam_off, signs, norm)

    def _G(self, lam_off, signs, norm):
        lams = self.lams - lam_off
        G = (lams * signs).sum(1)
        if norm:
            G_max = np.abs(lams).sum(1)
            G /= G_max
        return G

class GalEnv:

    def __init__(self, samp: ReconSample,
                 tf: TidalField, *,
                 rng: Rng,
                 zwgt: ZWeighting = None,
                 lr_range=[-0.1, 1.2], n_bins=8,
                 ):
        
        zs = samp.data['z_cor']
        if zwgt is not None:
            ws = zwgt.weights_of(zs)
        else:
            n_xs = len(zs)
            ws = np.ones(n_xs) / n_xs
        
        self.samp = samp
        self.tf = tf
        self.zwgt = zwgt
        self.ws = ws
        self.rng = Rng(rng)
        self.set_bin(lr_range, n_bins)

    def profiles(self, n_bootstrap=100, lam_off = 0.):
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

    def corrs(self, lam_off = 0., mass_weighted=True):
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
        xs = self.samp.data['x_sim_cor']
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

class GalEnvSimed:
    def __init__(self, samp: SimSample,
                 tf: TidalField, *,
                 rng: Rng,
                 lr_range=[-0.1, 1.2], n_bins=8,
                 ):
        
        self.samp = samp
        self.tf = tf
        self.rng = Rng(rng)
        self.set_bin(lr_range, n_bins)

    def profiles(self, n_bootstrap=100, lam_off = 0.):
        tf = self.tf
        lams = tf.lams.reshape(-1, 3)
        n_lam = (lams >= lam_off).sum(1)
        f_webs = []
        for i_web in range(4):
            vals = (n_lam == i_web).astype(float)
            p = self._v_profiles(vals)
            f_web = self._bootstrap(p, n_bootstrap)
            f_webs.append(f_web)
        
        rho = tf.delta.reshape(-1) + 1.
        p = self._v_profiles(rho)
        den_prof = self._bootstrap(p, n_bootstrap) 
        
        return DataDict({
            'f_webs': f_webs, 'den_prof': den_prof, 'lr': self.lrs,
        })

    def _v_profiles(self, vals: np.ndarray):
        assert len(vals) == self.n_grids
        cnts, nn_ids, bin_ids = self.cnts, self.nn_ids, self.bin_ids
        b_data = self.b_data
        tot_wgts = np.zeros_like(cnts)
        for i, (nn_id, bin_id) in enumerate(zip(nn_ids, bin_ids)):
            b_data.reset()
            b_data.add_n_chked(bin_id, vals[nn_id])
            tot_wgts[i] = b_data.data
        avg_wgts = tot_wgts / cnts.clip(1.0e-2)
        return avg_wgts

    def _bootstrap(self, outs: np.ndarray, n: int):
        if n == 0:
            mean = outs.mean(0)
            sd = np.zeros_like(mean)
            return mean, sd
        
        means = np.zeros((n, outs.shape[1]))
        for i in range(n):
            ids = self.rng.choice(len(outs), len(outs))
            means[i] = outs[ids].mean(0)
        mean, sd = means.mean(0), means.std(0)
        return mean, sd
        
    def set_bin(self, lr_range, n_bins):
        lr_min, lr_max = lr_range
        dlr = (lr_max - lr_min) / n_bins
        b_data = _BinnedData(n_bins)
        
        kdt = self.tf.kdt_grids
        xs = self.samp['x']
        nn_ids, nn_ds = kdt.query_radius(
            xs, r=10**lr_max, return_distance=True)

        bin_ids = []
        for nn_d in nn_ds:
            bin_id = np.floor(
                (np.log10(nn_d + 1.0e-10) - lr_min) / dlr).astype(int)
            bin_ids.append(bin_id)

        cnts = np.zeros((len(xs), n_bins), dtype=float)
        for i, bin_id in enumerate(bin_ids):
            b_data.reset()
            b_data.cnt_n_chked(bin_id)
            cnts[i] = b_data.data

        self.lr_min, self.lr_max = lr_min, lr_max
        self.dlr = dlr
        self.lrs = (np.arange(n_bins) + 0.5) * dlr + lr_min
        self.n_bins = n_bins
        self.b_data = b_data

        self.n_grids = self.tf.mesh.total_n_grids
        self.nn_ids = nn_ids
        self.nn_ds = nn_ds
        self.bin_ids = bin_ids
        self.cnts = cnts