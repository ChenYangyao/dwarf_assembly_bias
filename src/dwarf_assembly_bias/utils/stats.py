from __future__ import annotations
import typing
from typing import Self
from pyhipp.stats.sampling import RandomNoise, Rng
import numba
from numba.experimental import jitclass
import numpy as np
from scipy.interpolate import interp1d
from dataclasses import dataclass
from pyhipp.core import DataDict
from pyhipp.stats.prob_dist import ProbTransToNorm

@jitclass
class _BinnedData:

    data: numba.float64[:]

    def __init__(self, n_bins: int) -> None:

        self.data = np.zeros(n_bins, dtype=np.float64)

    def reset(self) -> None:
        self.data[:] = 0.0

    def add_n_chked(self, inds, vals):
        data = self.data
        n = len(data)
        for ind, val in zip(inds, vals):
            if 0 <= ind < n:
                data[ind] += val

    def cnt_n_chked(self, inds):
        data = self.data
        n = len(data)
        for ind in inds:
            if 0 <= ind < n:
                data[ind] += 1


def select_val(val: np.ndarray, lo=None, hi=None, p_lo=None, p_hi=None):
    sel = np.ones(len(val), dtype=bool)
    if lo is not None:
        sel &= val >= lo
    if hi is not None:
        sel &= val < hi
    if p_lo is not None:
        sel &= val >= np.quantile(val, p_lo)
    if p_hi is not None:
        sel &= val < np.quantile(val, p_hi)
    return sel


class RelativeBias:
    def __init__(self, wps: DataDict, wp_ref: int | DataDict | None):
        n_wps = len(wps)
        if wp_ref is None:
            wp_ref = n_wps - 1
        if isinstance(wp_ref, int):
            wp_ref = wps[str(wp_ref)]

        data = []
        for i in range(n_wps):
            data.append(self.__find_one(wps[str(i)], wp_ref))
        x, xel, xeh, y, yel, yeh = np.array(data).T
        self.data = DataDict({
            'x': x, 'xe': (xel, xeh),
            'y': y, 'ye': (yel, yeh),
        })

    @staticmethod
    def __find_one(d_in: DataDict, d_rel: DataDict):
        wp_rel = d_in['wp_samples'] / d_rel['wp_samples']
        assert wp_rel.shape[1] == 1
        wp_rel = wp_rel[:, 0]
        y, yl, yh = np.median(wp_rel), np.quantile(
            wp_rel, 0.16), np.quantile(wp_rel, 0.84)
        yel, yeh = y - yl, yh - y
        x, (xl, xh) = d_in['sub_sample']['median', '1sigma']
        xel, xeh = x - xl, xh - x
        return x, xel, xeh, y, yel, yeh


class Summary:

    @dataclass
    class FullResult:
        mean: np.ndarray
        sd: np.ndarray
        median: np.ndarray
        sigma_1: np.ndarray
        sigma_2: np.ndarray
        sigma_3: np.ndarray

    @staticmethod
    def on(vals: np.ndarray):
        vals = np.asarray(vals)
        n_vals = len(vals)
        assert n_vals > 0
        if n_vals == 1:
            mean = vals.mean(axis=0)
            sd = np.zeros_like(mean)
            median = mean.copy()
            sigma_1 = np.array([mean, mean])
            sigma_2 = np.array([mean, mean])
            sigma_3 = np.array([mean, mean])
        else:
            mean = np.mean(vals, axis=0)
            sd = np.std(vals, axis=0)
            median = np.median(vals, axis=0)
            sigma_1 = np.quantile(vals, [0.16, 0.84], axis=0)
            sigma_2 = np.quantile(vals, [0.025, 0.975], axis=0)
            sigma_3 = np.quantile(vals, [0.005, 0.995], axis=0)
        return Summary.FullResult(mean, sd, median, sigma_1, sigma_2, sigma_3)


def abundance_match(x_src, x_dst, noise: RandomNoise = None):
    i_src, i_dst = np.argsort(x_src), np.argsort(x_dst)
    x_src, x_dst = x_src[i_src], x_dst[i_dst]
    n_src, n_dst = len(x_src), len(x_dst)
    p_src = np.arange(n_src) / (n_src-1.0)
    p_dst = np.arange(n_dst) / (n_dst-1.0)

    pred = interp1d(p_dst, x_dst)(p_src)
    out = np.zeros_like(pred)
    out[i_src] = pred

    if noise is not None:
        out = noise.add_to(out)

    return out

def abudance_match_with_corr_coef(xs_src, xs_dst, rho=1., rng: Rng = 0):
    '''
    @rho: can be negativa.
    '''
    rng = Rng(rng)
    p_src, p_dst = ProbTransToNorm(xs_src), ProbTransToNorm(xs_dst)

    ys = p_src.forw(xs_src)
    eps = rng.normal(size=len(ys))
    ys = rho * ys + np.sqrt(1.0 - rho**2) * eps
    xs_pred = p_dst.back(ys)

    return xs_pred
