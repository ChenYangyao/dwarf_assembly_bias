from __future__ import annotations
import typing
from typing import Self
import numpy as np
from scipy.interpolate import interp1d
from pyhipp.stats.binning import BiSearchBins

class _HistInterp:
    def __init__(self, xs, range, bins):
        ys, es = np.histogram(xs, range=range, bins=bins)
        ys = np.array(ys, dtype=np.float64)
        xs = 0.5 * (es[:-1] + es[1:])
        n_at_x = interp1d(
            xs, ys, kind='slinear', bounds_error=False, fill_value=(0., 0.))
        
        self._range = range
        self._bins = bins
        self._es = es
        self._ys = ys
        self._xs = xs
        self._n_at_x = n_at_x
        self._bin = BiSearchBins(es)
        
    def __call__(self, xs, by='step'):
        assert by in ('interp', 'step')
        if by == 'interp':
            return self._n_at_x(xs)
        ys_out = np.zeros_like(xs)
        ids = self._bin.locate(xs)
        sel = (ids >= 0) & (ids < self._bins)
        ys_out[sel] = self._ys[ids[sel]]
        return ys_out

class ZWeighting:
    def __init__(self, zs_dst: np.ndarray, range=(0., 0.12), bins=16):
        self.zs_dst = zs_dst
        self.range = range
        self.bins = bins
        self.f_n_at_z = _HistInterp(zs_dst, range, bins)
        
    def weights_of(self, zs):
        ns = _HistInterp(zs, self.range, self.bins)(zs).clip(0.1)
        ns_dst = self.f_n_at_z(zs)
        w = ns_dst / ns
        w /= w.sum()
        return w
