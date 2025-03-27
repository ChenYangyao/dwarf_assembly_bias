from __future__ import annotations
import typing
from typing import Self
import numba
import gc
import numpy as np, pandas as pd
from pyhipp.core import abc, DataDict
from pyhipp.io import h5
from pyhipp_sims import sims
from functools import cached_property
from pathlib import Path


@numba.njit
def _find_one_ftime(snaps, masses, is_cs, m_dst) -> None:
    n_snaps = len(snaps)
    i_form = 0
    for i in range(n_snaps):
        mass, is_c = masses[i], is_cs[i]
        if not is_c:
            continue
        if mass < m_dst:
            i_form = i
            break
    i_afterform = 0
    for i in range(i_form):
        if is_cs[i]:
            i_afterform = i
    return i_form, snaps[i_form], snaps[i_afterform]


@numba.njit
def find_ftimes(root_offs, leaf_offs, snaps, masses, is_cs, f_m=0.5):
    offs_form = np.empty_like(root_offs)
    snaps_form = np.empty_like(root_offs)
    snaps_afterform = np.empty_like(root_offs)

    n_roots = len(root_offs)
    for i in range(n_roots):
        b, e = root_offs[i], leaf_offs[i]
        m_root = masses[b]
        m_dst = f_m * m_root
        i_form, snap_form, snap_afterform = _find_one_ftime(
            snaps[b:e], masses[b:e], is_cs[b:e], m_dst)
        offs_form[i] = b + i_form
        snaps_form[i] = snap_form
        snaps_afterform[i] = snap_afterform
    return offs_form, snaps_form, snaps_afterform


@numba.njit
def _find_one_peak(snaps, masses, v_maxs, is_cs):
    mp, mpv, mps = masses[0], v_maxs[0], snaps[0]
    for s, m, v, is_c in zip(snaps, masses, v_maxs, is_cs):
        if not is_c:
            continue
        if m > mp:
            mp, mpv, mps = m, v, s
    vp, vpm, vps = v_maxs[0], masses[0], snaps[0]
    for s, m, v, is_c in zip(snaps, masses, v_maxs, is_cs):
        if not is_c:
            continue
        if v > vp:
            vp, vpm, vps = v, m, s
    return mp, mpv, mps, vp, vpm, vps


@numba.njit
def find_peak(root_offs, leaf_offs, snaps, masses, v_maxs, is_cs):
    n_roots = len(root_offs)
    m_peak = np.zeros(n_roots, dtype=np.float32)
    m_peak_v = np.zeros(n_roots, dtype=np.float32)
    m_peak_snap = np.zeros(n_roots, dtype=np.int32)
    v_peak = np.zeros(n_roots, dtype=np.float32)
    v_peak_m = np.zeros(n_roots, dtype=np.float32)
    v_peak_snap = np.zeros(n_roots, dtype=np.int32)
    for i in range(n_roots):
        b, e = root_offs[i], leaf_offs[i]
        mp, mpv, mps, vp, vpm, vps = _find_one_peak(snaps[b:e], masses[b:e],
                                                    v_maxs[b:e], is_cs[b:e])
        m_peak[i] = mp
        m_peak_v[i] = mpv
        m_peak_snap[i] = mps
        v_peak[i] = vp
        v_peak_m[i] = vpm
        v_peak_snap[i] = vps
    return m_peak, m_peak_v, m_peak_snap, v_peak, v_peak_m, v_peak_snap


@staticmethod
@numba.njit
def find_last_sat(root_offs, leaf_offs, is_cs, snaps):
    last_sat_snaps = np.empty_like(root_offs)
    for i, (b, e) in enumerate(zip(root_offs, leaf_offs)):
        out = -1
        for j in range(b, e):
            snap, is_s = snaps[j], not is_cs[j]
            if is_s:
                out = snap
                break
        last_sat_snaps[i] = out
    return last_sat_snaps


def dump_sample_to(path: Path, objs: dict, sim_info: sims.SimInfo,
                   snap_dst: int):
    h5.File.dump_to(path, {
        'objs': objs,
        'header': {
            'm_sub_lb': 0.,
            'z_dst': sim_info.redshifts[snap_dst],
            'snap_dst': snap_dst,
        },
        'ctx': {
            'sim_name': sim_info.name
        }
    }, f_flag='w')


class SampleTngDark(abc.HasLog):
    def __init__(self, tr_ld: sims.TreeLoaderTngDark, z_dst: float = 0.0,
                 only_c=True, m_sub_lb=5.0e-1, verbose=True) -> None:
        super().__init__(verbose=verbose)

        sim_info = tr_ld.sim_info
        snap_dst = sim_info.z_to_snap(z_dst)

        self.tr_ld = tr_ld
        self.z_dst = z_dst
        self.snap_dst = snap_dst
        self.sim_info = sim_info
        self.only_c = only_c
        self.m_sub_lb = m_sub_lb

    def run(self):
        self.objs = self.__select_branches()

    def dump_to(self, g: h5.Group):
        g.dump({
            'header': {
                'z_dst': self.z_dst, 'snap_dst': self.snap_dst,
                'm_sub_lb': self.m_sub_lb,
            },
            'objs': self.objs,
            'ctx': {
                'sim_name': self.sim_info.name,
            }
        })

    def __select_branches(self) -> None:
        tr_ld, snap_dst, m_sub_lb = self.tr_ld, self.snap_dst, self.m_sub_lb
        snaps = tr_ld['snap']
        hids, leafs, cids, m_sub = tr_ld['subhalo_id',
                                         'main_leaf', 'f_in_grp', 'm_sub']
        is_cs = cids == hids
        sel = (snaps == snap_dst) & (m_sub >= m_sub_lb)
        if self.only_c:
            sel &= is_cs
        root_offs = sel.nonzero()[0]
        leaf_offs = root_offs + (leafs[root_offs] - hids[root_offs]) + 1
        n_objs = len(root_offs)
        self.log(f'Found {n_objs} galaxies')

        return DataDict({
            'is_c': is_cs[root_offs],
            'root_offs': root_offs,
            'leaf_offs': leaf_offs,
        })


class ExtraPropTngDark(abc.HasLog):
    def __init__(self, tr_ld: sims.TreeLoaderTngDark,
                 in_group: h5.Group,
                 out_group: h5.Group,
                 verbose=True, **kw) -> None:
        super().__init__(verbose=verbose, **kw)

        self.tr_ld = tr_ld
        self.in_group = in_group
        self.out_group = out_group

    @cached_property
    def sample_info(self) -> DataDict:
        dsets = self.in_group['objs'].datasets
        keys = 'root_offs', 'leaf_offs'
        return dsets.load(keys=keys)

    def dump(self, d: DataDict, flag='ac'):
        self.out_group.dump(d, flag=flag)

    def load_tree_props(self, keys: tuple[str] | None = None):
        tr_ld = self.tr_ld
        if keys is None:
            keys = ('m_crit200', 'm_mean200', 'r_v_max',
                    'spin', 'm_sub', 'm_tophat', 'r_tophat', 'x')
            if tr_ld.sim_info.is_baryon:
                keys += ('r_half_mass_star', 'sfr', 'm_star')
        root_offs = self.sample_info['root_offs']
        for k in keys:
            self.log(f'Loading {k}...')
            self.dump({k: tr_ld[k][root_offs]})
            tr_ld.cache.pop(k)
            gc.collect()

    def find_ftimes(self):
        root_offs, leaf_offs = self.sample_info['root_offs', 'leaf_offs']
        snaps, m_hs, cids, hids = self.tr_ld[
            'snap', 'm_crit200', 'f_in_grp', 'subhalo_id']
        is_cs = cids == hids
        # ftime only valid for central
        off_form, snap_form, snap_afterform = find_ftimes(
            root_offs, leaf_offs, snaps, m_hs, is_cs)
        self.dump({
            'snap_form': snap_form,
            'snap_afterform': snap_afterform,
        })

        keys = 'm_crit200', 'm_mean200', 'spin'
        for k in keys:
            k_out = f'{k}_form'
            self.log(f'Loading {k_out}...')
            self.dump({k_out: self.tr_ld[k][off_form]})

    def find_peak(self):
        root_offs, leaf_offs = self.sample_info['root_offs', 'leaf_offs']
        snaps, m_hs, cids, hids, v_maxs = self.tr_ld[
            'snap', 'm_crit200', 'f_in_grp', 'subhalo_id', 'v_max']
        is_cs = cids == hids
        # peak property only valid for central
        m_peak, m_peak_v, m_peak_snap, v_peak, v_peak_m, v_peak_snap = \
            find_peak(root_offs, leaf_offs, snaps, m_hs, v_maxs, is_cs)
        self.dump({
            'm_peak': m_peak,
            'm_peak_v': m_peak_v,
            'm_peak_snap': m_peak_snap,
            'v_peak': v_peak,
            'v_peak_m': v_peak_m,
            'v_peak_snap': v_peak_snap,
            'v_max': v_maxs[root_offs],
        })

    def find_infall(self) -> Self:
        root_offs, leaf_offs = self.sample_info['root_offs', 'leaf_offs']
        snaps, cids, hids = self.tr_ld['snap', 'f_in_grp', 'subhalo_id']
        is_cs = cids == hids
        last_sat_snap = find_last_sat(root_offs, leaf_offs, is_cs, snaps)
        self.dump({
            'last_sat_snap': last_sat_snap,
        })


class LGalaxiesCat:

    default_model_keys = [
        'StellarBulgeMass',     # 10^10 M_sun/h
        'StellarDiskMass',
        'StellarDiskRadius',    # ckpc/h
        'Type',                 # 0: central, 1: satellite, 2 orphan
        'Central_M_Crit200',
        'M_Crit200',
        'StellarDiskSpin',      # (kpc/h)(km/s)
        'StarFormationRate',    # Msun/yr
        'InfallSnap',
    ]

    def __init__(self, dst_path: Path, mod_path: Path):

        dst_hs = h5.File.load_from(dst_path, 'objs')
        dst_subf_id = dst_hs['subfind_id']
        with h5.File(mod_path) as f:
            tps, subf_ids = f['Galaxy'].datasets[
                'Type', 'SubhaloIndex_TNG-Dark']
            dst_offs = pd.Index(dst_subf_id).get_indexer(subf_ids)
            mod_sel = (tps < 2) & (subf_ids >= 0) & (dst_offs >= 0)

            mod_offs = mod_sel.nonzero()[0]
            dst_offs = dst_offs[mod_sel]

        self.mod_offs = mod_offs
        self.dst_offs = dst_offs
        self.mod_path = mod_path
        self.dst_path = dst_path
        self.n_dsts = len(dst_subf_id)

        print(f'Found {len(self.dst_offs)} / {self.n_dsts}')

    def dump(self, suffix='.LGalaxies-H15', mod_keys=default_model_keys):
        with h5.File(self.mod_path) as f_mod, \
                h5.File(str(self.dst_path) + suffix, 'x') as f_dst:
            g_dst = f_dst.create_group('objs')
            for mod_key in mod_keys:
                mod_v = f_mod['Galaxy'].datasets[mod_key]
                out_shape = (self.n_dsts,) + mod_v.shape[1:]
                out_dt = mod_v.dtype
                out_v = np.full(out_shape, -1, dtype=out_dt)
                out_v[self.dst_offs] = mod_v[self.mod_offs]
                g_dst.datasets.create(mod_key, out_v)
