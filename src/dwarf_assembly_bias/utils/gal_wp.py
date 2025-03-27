from __future__ import annotations
import typing
from pyhipp.core import DataDict
from typing import Self
from ..clustering import SubSampleManager, SubSample, SimSample, Rng, pair_count_proj_periodic, n2wp_periodic
import numpy as np


class WpVsProp:
    def __init__(
            self, sman: SubSampleManager, samp_name='s_dwarf', samp_kw=None,
            prop_key='z_half', ps=np.linspace(0., 1., 201),
            rs_range=(2., 10.),
            pimax=40., n_threads=1, n_reps=5, rng=10086):

        sim_info = sman.sim_info
        wp_kw = dict(
            l_box=sim_info.full_box_size,
            rs=np.linspace(*rs_range, 2),
            pimax=pimax,
            n_threads=n_threads,)

        if samp_kw is None:
            s1_parent = getattr(sman, samp_name)
        else:
            s1_parent = getattr(sman, samp_name)(**samp_kw)
        s1_parent: SimSample
        s2 = SubSample(sman.s_ref, add_rsd=sman.add_rsd)

        self.sman = sman
        self.sim_info = sim_info
        self.s1_parent = s1_parent
        self.s2 = s2
        self.prop_key = prop_key
        self.ps = np.asarray(ps)
        self.wp_kw = wp_kw
        self.n_reps = n_reps
        self.rng = Rng(rng)

    def run(self):
        pcs = [self._get_pc_i(i) for i in range(len(self.ps)-1)]
        self.pcs = pcs

    def _get_pc_i(self, i):
        p_lo, p_hi = self.ps[i], self.ps[i+1]
        s1 = self.s1_parent.subset_by_percentile(self.prop_key, p_lo, p_hi)
        s1 = SubSample(s1, add_rsd=self.sman.add_rsd)
        x1, x2 = s1.x_boxed, self.s2.x_boxed
        n1 = len(x1)
        pcs = []
        for _ in range(self.n_reps):
            _x1 = x1[self.rng.choice(n1, n1)]
            pc = pair_count_proj_periodic(_x1, x2, **self.wp_kw)
            pcs.append(pc)
        return pcs

    def summary(self, inv_cum):
        n_reps, pcs = self.n_reps, self.pcs
        wps_all = []
        for i_rep in range(n_reps):
            ds = [d[i_rep] for d in pcs]
            wps = self._summary_1(ds, inv_cum)
            wps_all.append(wps)
        wps_all = np.array(wps_all)

        ps = self.ps
        vals = self.s1_parent[self.prop_key]
        if inv_cum:
            cum_ps = 1.0 - ps[:-1][::-1]
            qs = -np.quantile(-vals, cum_ps)
        else:
            cum_ps = ps[1:]
            qs = np.quantile(vals, cum_ps)

        out = DataDict({
            'ps': ps,
            'cum_ps': cum_ps,
            'qs': qs,
            'wps': wps_all,
            'n_objs': self.s1_parent.n_objs,
        })
        self.wps = out

    @staticmethod
    def _summary_1(ds, inv_cum):
        n_pairs, n1 = 0., 0.
        wps = []
        if inv_cum:
            ds = ds[::-1]
        for d in ds:
            n_pairs = n_pairs + np.array(d['n_pairs'])
            n1 = n1 + d['n1']
            n2 = d['n2']
            rs = np.array(d['rs'])
            l_box = d['l_box']
            wp = n2wp_periodic(n_pairs=n_pairs, n1=n1, n2=n2,
                               rs=rs, l_box=l_box)['wp']
            wps.append(wp)
        wps = np.concatenate(wps)
        return wps
