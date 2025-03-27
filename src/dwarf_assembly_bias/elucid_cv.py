from __future__ import annotations
import typing
from pyhipp.core import DataDict
from pyhipp_sims import sims
from pyhipp.stats.random import Rng
from typing import Self
from scipy.interpolate import interp1d
from functools import cached_property
from .sample import SubhaloSet, SimSample
from .clustering import SubSample
import numpy as np

class ReconSample:
    
    def __init__(self, sim_info: sims.SimInfo, data_file, 
        weight_curve_file,
        n_interp = 32,
        m_h_range_in_sol = (10**10.5, 10**11.),
        z_lc_lb = 15., 
        rng = 10086,
        rc_mask='is_recon') -> None:
        
        u_m = 1.0e10 / sim_info.cosmology.hubble
        m_lb, m_ub = m_h_range_in_sol
        m_lb, m_ub = m_lb / u_m, m_ub / u_m
        
        self.sim_info = sim_info
        self.data_file = data_file
        self.weight_curve_file = weight_curve_file
        self.n_interp = n_interp
        self.m_lb = m_lb
        self.m_ub = m_ub
        self.z_lc_lb = z_lc_lb
        self.rng = Rng(rng)
        self.rc_mask = rc_mask
        
        self.__load_weight_curve()
        self.__load_data()
        
    def find_weight(self, z: np.ndarray):
        z_range = 0., max(z.max(), 0.12) 
        hs, es = np.histogram(z, range=z_range, bins=self.n_interp, density=True)
        zs = 0.5 * (es[:-1] + es[1:])
        f_at_z = interp1d(zs, hs, kind='slinear', 
            bounds_error=False, fill_value=(0., 0.))
        hs_pred = f_at_z(z)
        hs_dst = self.f_at_z(z)
        w = hs_dst / hs_pred.clip(1.0e-10)
        w /= w.sum()
        return w
        
    @cached_property
    def s_c(self):
        return self.s_all.subset_by_value('is_c', eq=True)
    
    @cached_property
    def s_dwarf(self):
        s_dwarf = self.s_c\
            .subset_by_value('m_mean200', lo=self.m_lb, hi=self.m_ub)\
            .subset_by_value('last_sat_z', lo=15.)\
            .subset_by_value(self.rc_mask, eq=True)
        return s_dwarf
    
    def ramdom_subs_dwarf(self):
        s = self.s_dwarf
        n = s.n_objs
        ids = self.rng.choice(n, n)
        s_rand = s.subset(ids)
        return s_rand
    
    @cached_property
    def s_halo(self):
        s_halo = self.s_c\
            .subset_by_value(self.rc_mask, eq=True)
        return s_halo
    
    @cached_property
    def s_ref(self):
        return self.s_all.subset_by_value('m_peak', lo=self.m_lb)
    
    @cached_property
    def subs_ref(self):
        return SubSample(self.s_ref)
        
    def __load_weight_curve(self):
        z, frac = np.loadtxt(self.weight_curve_file, dtype=float).T
        f_at_z = interp1d(z, frac, kind='slinear', 
                            bounds_error=False, fill_value=(0.,0.))
        
        self.f_at_z = f_at_z
        
    def __load_data(self):
        s_all = SimSample.load(self.sim_info, self.data_file)
        self.s_all = s_all
        
coarse_rs = np.concatenate([
    np.logspace(-2., 0., 5),
    np.logspace(0.1, 1.2, 12)
])
default_ps = [0.00, 0.02, 0.33, 0.67, 0.98, 1.00]

class ObsProjCross:
    def __init__(self, rc_samp: ReconSample) -> None:
        self.rc_samp = rc_samp

    def find_wp(self, subs: SubSample, pimax=40., rs=coarse_rs,
                n_threads=4, **corr_kw):
        out = self.rc_samp.subs_ref.obs_proj_cross_with(
            subs, pimax=pimax, rs=rs, n_threads=n_threads, n_repeat=1,
            **corr_kw)
        return out

    def find_wp_in_simframe(self, subs: SubSample, pimax=10., n_threads=4,
            n_repeat=20, **corr_kw):
        out = self.rc_samp.subs_ref.proj_cross_with(
            subs, pimax=pimax, n_threads=n_threads, n_repeat=n_repeat,
            **corr_kw)
        return out

    def wp_of(self, samp_name='s_dwarf',
              sub_key=None, ps=None, qs=None,
              n_z_match=5000, n_repeat=5, in_simframe=False,
              **corr_kw):
        
        s = self.sub_sample_of(samp_name=samp_name, 
            sub_key=sub_key, ps=ps, qs=qs)
        
        outs = {}
        wp_samples = []
        for i in range(n_repeat):
            print(f'resample {i}', end=' ', flush=True)
            _s = self.resample_match_z(s, n_z_match=n_z_match)
            subs = SubSample(_s)
            if in_simframe:
                out = self.find_wp_in_simframe(subs, **corr_kw)
            else:
                out = self.find_wp(subs, **corr_kw)
            wp_samples.append(out['wp'])
            if i == 0:
                outs['lg_rs_c'] = out['lg_rs_c']
        wp_samples = np.array(wp_samples)
        wp, wp_sd = wp_samples.mean(axis=0), np.std(wp_samples, axis=0)
        outs |= {'wp': wp, 'wp_sd': wp_sd, 'wp_samples': wp_samples}
        print('Done')
        return DataDict(outs)

    def wps_of(self, samp_name='s_dwarf',
            sub_key = 'z_half', ps=default_ps, qs=None,
            n_z_match=5000, n_repeat=5, in_simframe=False,
            **corr_kw):
        if ps is not None:
            assert qs is None
            ps = list(zip(ps[:-1], ps[1:]))
            qs = (None,)*len(ps)
        else:
            assert qs is not None
            qs = list(zip(qs[:-1], qs[1:]))
            ps = (None,)*len(qs)
        outs = DataDict()
        for i in range(len(ps)):
            outs[str(i)] = self.wp_of(samp_name=samp_name,
                sub_key=sub_key, ps=ps[i], qs=qs[i],
                n_z_match=n_z_match, n_repeat=n_repeat, 
                in_simframe=in_simframe, **corr_kw)
        return outs
    
    
    def sub_sample_of(self, samp_name='s_dwarf',
        sub_key=None, ps=None, qs=None):
        s = getattr(self.rc_samp, samp_name)
        s = self.__sub_sample(s, sub_key=sub_key, ps=ps, qs=qs)
        return s
    
    def resample_match_z(self, s: SubhaloSet, n_z_match=5000):
        n = s.n_objs
        ids = self.rc_samp.rng.choice(n, n)
        _s = s.subset(ids)
        _s = self.__re_sample(_s, n_z_match=n_z_match) 
        return _s

    def __sub_sample(self, s: SubhaloSet,
                     sub_key=None, ps=None, qs=None):
        if sub_key is not None:
            if ps is not None:
                assert qs is None
                qs = np.quantile(s[sub_key], ps)
            q_lo, q_hi = qs
            s = s.subset_by_value(sub_key, lo=q_lo, hi=q_hi)
        return s
    
    def __re_sample(self, s: SubhaloSet, n_z_match=5000):
        rc_samp = self.rc_samp
        wgt = rc_samp.find_weight(s['z'])
        if n_z_match is not None:
            ids = rc_samp.rng.choice(len(wgt), n_z_match, p=wgt)
            s = s.subset(ids)
            wgt = np.ones(n_z_match)
        s.objs['wgt'] = wgt
        return s