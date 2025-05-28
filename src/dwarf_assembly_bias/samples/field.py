# Copyright (C) 2025 Yangyao Chen (yangyaochen.astro@foxmail.com) - All Rights 
# Reserved
# 
# You may use, distribute and modify this code under the MIT license. We kindly
# request you to give credit to the original author(s) of this code, and cite 
# the following paper(s) if you use this code in your research: 
# - Zhang Z. et al. 2025. Nature ???, ???.
# - Wang H. et al. 2016. ApJ, 831, 164 (for the reconstructed field).
# - Chen, et al. 2020. ApJ, 899, 81 (for the numerical method of field 
#   analysis).

from __future__ import annotations
from pyhipp.io import h5
from pyhipp.core import abc
from pyhipp.field.cubic_box import Mesh
import numpy as np
from pathlib import Path
from sklearn.neighbors import KDTree
from functools import cached_property

class TidalField(abc.HasLog):
    def __init__(self, lams: np.ndarray, delta: np.ndarray,
                 mesh: Mesh, lam_off: float = 0.,
                 recon_mask: np.ndarray = None, 
                 verbose = True):

        super().__init__(verbose=verbose)

        self.log(f'TidalField: threshold={lam_off}')
        if recon_mask is None:
            self.log('All field volume are used.')
            recon_mask = np.ones_like(delta, dtype=bool)

        self.lams = lams
        self.delta = delta
        self.mesh = mesh
        self.lam_off = lam_off
        self.recon_mask = recon_mask
        self._xs_web_t_cache = {}
        self.n_lams = (lams >= lam_off).sum(-1)

    @classmethod
    def from_file(cls, path: Path | str, recon_only=True, **init_kw):
        with h5.File(path) as f:
            lams, delta, n_grids, l_box = f.datasets[
                'lam', 'delta_sm_x', 'n_grids', 'l_box']
            if recon_only:
                recon_mask = f.datasets['reconstruction_mask']
            else:
                recon_mask = None
            mesh = Mesh.new(n_grids, l_box)
            lams = np.array(lams, dtype=np.float32)
        return cls(lams, delta, mesh, recon_mask=recon_mask, **init_kw)

    @cached_property
    def xs_grids(self):
        mesh = self.mesh
        n, h = mesh.n_grids, mesh.l_grid

        idx_1d = np.arange(n)
        x_1d = idx_1d * h
        x, y, z = np.meshgrid(x_1d, x_1d, x_1d, indexing='ij')
        X = np.stack([x, y, z], axis=-1)

        return X

    @cached_property
    def xs_recon(self):
        xs = self.xs_grids[self.recon_mask]
        return xs.reshape(-1, 3)

    @cached_property
    def kdt_recon(self):
        return KDTree(self.xs_recon)

    @cached_property
    def kdt_grids(self):
        xs = self.xs_grids.reshape(-1, 3)
        return KDTree(xs)

    def xs_web_t(self, n_lam):
        cache = self._xs_web_t_cache
        if n_lam in cache:
            return cache[n_lam]
        mask = self._mask_web_t(n_lam)
        xs = self.xs_grids[mask].reshape(-1, 3)
        cache[n_lam] = xs
        return xs

    def wet_t_at(self, xs: np.ndarray):
        xs = np.asarray(xs)
        n_xs = len(xs)
        assert xs.shape == (n_xs, 3)
        xis = self.mesh._impl.x_to_xi_nd(xs)
        xi0, xi1, xi2 = xis.T
        return self.n_lams[xi0, xi1, xi2]

    def _mask_web_t(self, n_lam):
        mask = (self.lams >= self.lam_off).sum(-1) == n_lam
        mask &= self.recon_mask
        return mask