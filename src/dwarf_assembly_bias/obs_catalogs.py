from __future__ import annotations
import typing
from typing import Self
from pyhipp.io import h5
from pyhipp.field.cubic_box import Mesh
from pyhipp.core import abc, DataTable, DataDict
from pyhipp.astro.cosmology.model import LambdaCDM
from functools import cached_property
from elucid.geometry.mask import ReconstAreaMaskSdssL500, reconst_area_masks
import numpy as np
from pyhipp_sims import sims
from sklearn.neighbors import KDTree
from .config import ProjPaths
from .sample import SubhaloSet
from .utils.stats import select_val

default_cosm = sims.predefined['elucid'].cosmology


class Dwarfs(abc.HasDictRepr):

    data_dir = ProjPaths.obs_dir
    dtype = [
        ('ra', 'f8'), ('dec', 'f8'), ('z', 'f8'),
        ('id0', 'i8'), ('log_m_s', 'f8'), ('id1', 'i8'),
    ]

    def __init__(self, file_name: str, cosm: LambdaCDM = default_cosm) -> None:
        self.file_name = file_name
        self.cosm = cosm

    @cached_property
    def data(self) -> DataTable:
        path, dtype = self.data_dir / self.file_name, self.dtype
        d_in = np.loadtxt(path, dtype=dtype)
        d_out = DataTable({name: d_in[name] for name, _ in dtype})

        cosm = self.cosm
        z0, z1 = max(d_out['z'].min()-0.01, 0.0), d_out['z'].max()+0.01

        zs = np.linspace(z0, z1, 256)
        d_cs = cosm.distances.comoving_at(zs)
        d_cs = np.interp(d_out['z'], zs, d_cs)
        d_out['d_c'] = d_cs

        return d_out


class TotalDwarfs(abc.HasDictRepr):

    data_dir = ProjPaths.obs_dir
    dtype = [
        ('z', 'f8'), ('lm_s', 'f8'), ('r_50', 'f8'),
        ('r_50_in_arcsec', 'f8'), ('Mag_r', 'f8'), ('color_gr', 'f8'),
        ('n_sersic', 'f8'), ('Sigma_star', 'f8'),
    ]

    def __init__(self, file_name, cosm: LambdaCDM = default_cosm) -> None:
        self.file_name = file_name
        self.cosm = cosm

    @cached_property
    def data(self) -> DataTable:
        path, dtype = self.data_dir / self.file_name, self.dtype
        d_in = np.loadtxt(path, dtype=dtype)
        d_out = DataTable({name: d_in[name] for name, _ in dtype})

        return d_out

    @cached_property
    def s_massive(self) -> DataTable:
        data = self.data
        z, lm, color, n = data['z', 'lm_s', 'color_gr', 'n_sersic']
        sel = (z < 0.04) & (lm >= 8.5) & (color <= 0.6) & (n <= 1.6)
        return data.subset(sel)
    
    @cached_property
    def s_massive_z003(self) -> DataTable:
        data = self.data
        z, lm, color, n = data['z', 'lm_s', 'color_gr', 'n_sersic']
        sel = (z < 0.03) & (lm >= 8.5) & (color <= 0.6) & (n <= 1.6)
        return data.subset(sel)


diffuse_dwarfs = Dwarfs('samp_0')
compact_dwarfs = Dwarfs('samp_3')
tot_dwarfs = TotalDwarfs('total_dwarfs')


class ReconSample(abc.HasDictRepr):

    data_dir = ProjPaths.obs_dir / 'recon_samples'

    repr_attr_keys = ('cat_name', 'n_objs', 'keys')

    def __init__(
            self, cat_name: str, data: DataDict, sim_info: sims.SimInfo) -> None:

        keys = tuple(data.keys())
        n_objs = len(data[keys[0]])

        self.cat_name = cat_name
        self.data = data
        self.sim_info = sim_info
        self.n_objs = n_objs
        self.keys = keys

    @classmethod
    def new_from_file(cls, cat_name: str, sim_info=sims.predefined['elucid']):
        data = h5.File.load_from(cls.data_dir / (cat_name+'.hdf5'))
        return cls(cat_name, data, sim_info)

    def subset_by(self, sel):
        data = DataDict({k: v[sel] for k, v in self.data.items()})
        return type(self)(self.cat_name, data, self.sim_info)

    def subset_by_val(self, key, lo=None, hi=None, p_lo=None, p_hi=None):
        val = self.data[key]
        sel = select_val(val, lo, hi, p_lo, p_hi)
        return self.subset_by(sel)


class _ReconSamples:
    @cached_property
    def diffuse_dwarfs(self):
        return ReconSample.new_from_file('diffuse_dwarfs')

    @cached_property
    def diffuse_dwarfs_2(self):
        return ReconSample.new_from_file('diffuse_dwarfs.2')

    @cached_property
    def compact_dwarfs(self):
        return ReconSample.new_from_file('compact_dwarfs')

    @cached_property
    def compact_dwarfs_2(self):
        return ReconSample.new_from_file('compact_dwarfs.2')

    @cached_property
    def groups(self):
        return ReconSample.new_from_file('groups')


recons_samples = _ReconSamples()


class TidalField:
    def __init__(self, lams: np.ndarray, delta: np.ndarray,
                 mesh: Mesh, lam_off: float,
                 recon_mask: np.ndarray):
        self.lams = lams
        self.delta = delta
        self.mesh = mesh
        self.lam_off = lam_off
        self.recon_mask = recon_mask
        self._xs_web_t_cache = {}

        self.n_lams = (lams >= lam_off).sum(-1)

    @classmethod
    def from_file(cls, fname: str, lam_off: float, recon_mask: np.ndarray):
        with h5.File(fname) as f:
            lams, delta, n_grids, l_box = f.datasets['lam',
                                                     'delta_sm_x', 'n_grids', 'l_box']
            mesh = Mesh.new(n_grids, l_box)
            lams = np.array(lams, dtype=np.float32)
        return cls(lams, delta, mesh, lam_off, recon_mask)

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


class _TidalFields:
    def __init__(self, sim_info=sims.predefined['elucid'],
                 mask_name='large',
                 lam_off: float = 0.) -> None:

        self.sim_info = sim_info
        self.field_dir = ProjPaths.sim_dir_of(sim_info) / 'fields'
        self.mask_name = mask_name
        self.lam_off = lam_off

    @cached_property
    def recon_mask_n500(self):
        return getattr(reconst_area_masks, self.mask_name).mask_reconst_area

    @cached_property
    def recon_mask_n512(self):
        return getattr(
            reconst_area_masks, self.mask_name + '_n512').mask_reconst_area

    def new_by_file(self, fname: str, n_grids):
        return TidalField.from_file(
            fname,
            self.lam_off,
            getattr(self, f'recon_mask_n{n_grids}'))

    @cached_property
    def dom_z0_sm1(self):
        return self.new_by_file(
            self.field_dir / 'domain.tidal.s99.sm1.hdf5', 500)

    @cached_property
    def dom_z0_sm2(self):
        return self.new_by_file(
            self.field_dir / 'domain.tidal.s99.sm2.hdf5', 500)

    @cached_property
    def dom_z0_sm4(self):
        return self.new_by_file(
            self.field_dir / 'domain.tidal.s99.sm4.hdf5', 500)

    @cached_property
    def resim_z0_sm1(self):
        return self.new_by_file(self.field_dir / 'tidal.s99.sm1.hdf5', 512)

    @cached_property
    def resim_z0_sm2(self):
        return self.new_by_file(self.field_dir / 'tidal.s99.sm2.hdf5', 512)


tidal_fields = _TidalFields()
