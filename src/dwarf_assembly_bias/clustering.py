from __future__ import annotations
import typing
from typing import Tuple
from pyhipp.stats import Rng, reduction, Bootstrap
from pyhipp.core import Num, DataDict, abc
from .sample import SubhaloSet, SimSample
from dataclasses import dataclass
from pyhipp_sims import sims
import numpy as np
from sklearn.neighbors import KDTree
from scipy import spatial
from functools import cached_property


def r2rc(rs: np.ndarray):
    lg_rs = Num.safe_lg(rs)
    lg_rs_c = .5 * (lg_rs[1:] + lg_rs[:-1])
    return DataDict({
        'rs': rs.copy(), 'lg_rs': lg_rs, 'lg_rs_c': lg_rs_c
    })


def n2wp_periodic(n_pairs: np.ndarray,
                  n1: int, n2: int, rs: np.ndarray, l_box: float):
    vol = np.pi * rs**2 * 2.0
    dvol = np.diff(vol)
    exp_n = dvol * n2 / l_box**3
    xi = n_pairs / n1 / exp_n[:, None] - 1.
    wp = np.sum(xi * 2., axis=1)
    return DataDict({
        'exp_n': exp_n, 'dvol': dvol, 'xi': xi, 'wp': wp
    }) | r2rc(rs)


def pair_count_proj_periodic(
        x1: np.ndarray, x2: np.ndarray | None, l_box: float, rs: np.ndarray,
        n_threads=1, pimax=10.0):

    import Corrfunc
    
    kw = {'nthreads': n_threads, 'pimax': pimax, 'binfile': rs,
          'periodic': True, 'boxsize': l_box}
    kw['X1'], kw['Y1'], kw['Z1'] = x1.T
    if x2 is not None:
        kw['autocorr'] = False
        kw['X2'], kw['Y2'], kw['Z2'] = x2.T
    else:
        kw['autocorr'] = True

    n = Corrfunc.theory.DDrppi(**kw)['npairs']
    n = n.reshape(len(rs)-1, -1).astype(np.float64)
    n1 = len(x1)
    n2 = len(x2) if x2 is not None else n1

    return DataDict({
        'n_pairs': n,
        'n1': n1, 'n2': n2, 'rs': rs, 'pimax': pimax, 'l_box': l_box
    })
    
def wp_with_ratio(out: DataDict, out_ref: DataDict):
    ys, lx = out['wp_samples', 'lg_rs_c']
    ys_ref, lx_ref = out_ref['wp_samples', 'lg_rs_c']
    assert (np.abs(lx-lx_ref) < 1.0e-5).all()
    assert len(ys) == len(ys_ref)
    
    x = 10.**lx
    y_mean = np.mean(ys, axis=0)
    y_median = np.median(ys, axis=0)
    y_sd = np.std(ys, axis=0)
    y_lo, y_hi = np.percentile(ys, [16., 84.], axis=0)
    y_ebar = y_median, (y_median - y_lo, y_hi - y_median)
    y = DataDict({
        'mean': y_mean, 'sd': y_sd, 'median': y_median, '1sigma': (y_lo, y_hi),
        'ebar': y_ebar
    })
    
    yrs = Num.safe_div(ys, ys_ref)
    yr_mean = np.mean(yrs, axis=0)
    yr_median = np.median(yrs, axis=0)
    yr_sd = np.std(yrs, axis=0)
    yr_lo, yr_hi = np.percentile(yrs, [16., 84.], axis=0)
    yr_ebar = yr_median, (yr_median - yr_lo, yr_hi - yr_median)
    yr = DataDict({'mean': yr_mean, 'sd': yr_sd, 'median': yr_median,
                   '1sigma': (yr_lo, yr_hi), 'ebar': yr_ebar})
    
    return DataDict({'x': x, 'y': y, 
                     'y_ratio': yr})

class _Clustering:

    Result = DataDict[str, np.ndarray]

    def __init__(self, sim_info: sims.SimInfo) -> None:
        self.sim_info = sim_info
        self.l_box = sim_info.full_box_size

    def cross_in_rad(self, x1, x2, r_min: float, r_max: float, n_threads=1):
        import Corrfunc
        
        l_box = self.l_box
        n = Corrfunc.theory.DD(
            False, n_threads, [r_min, r_max],
            x1[:, 0],
            x1[:, 1],
            x1[:, 2],
            periodic=True, boxsize=l_box, X2=x2[:, 0],
            Y2=x2[:, 1],
            Z2=x2[:, 2])['npairs'][0]

        n1, n2 = len(x1), len(x2)
        exp_n = 4./3. * np.pi * (r_max**3 - r_min**3) * n2 / l_box**3
        xi = float(n) / n1 / exp_n - 1.

        return self.Result({
            'r_min': r_min, 'r_max': r_max,
            'n1': n1, 'n2': n2, 'n': n, 'xi': xi
        })

    def cross(self, x1, x2, rs: np.ndarray, n_threads=1):
        
        import Corrfunc
        
        l_box = self.l_box
        rs = np.array(rs)
        n = Corrfunc.theory.DD(
            False, n_threads, rs,
            x1[:, 0],
            x1[:, 1],
            x1[:, 2],
            periodic=True, boxsize=l_box, X2=x2[:, 0],
            Y2=x2[:, 1],
            Z2=x2[:, 2])['npairs'].astype(float)

        n1, n2 = len(x1), len(x2)
        vol = 4./3. * np.pi * rs**3
        dvol = np.diff(vol)
        exp_n = dvol * n2 / l_box**3
        xi = n / n1 / exp_n - 1.
        lg_rs, lg_rs_c = self.__cvt_r(rs)
        return self.Result({
            'rs': rs, 'lg_rs': lg_rs, 'lg_rs_c': lg_rs_c, 'xi': xi
        })

    def proj_cross(self, x1, x2, rs: np.ndarray, pimax=10., n_threads=1):
        
        import Corrfunc
        
        l_box = self.l_box
        rs = np.array(rs)
        n = Corrfunc.theory.DDrppi(
            False, n_threads, pimax=pimax, binfile=rs, X1=x1[:, 0],
            Y1=x1[:, 1],
            Z1=x1[:, 2],
            periodic=True, boxsize=l_box, X2=x2[:, 0],
            Y2=x2[:, 1],
            Z2=x2[:, 2])['npairs']
        n = n.reshape(len(rs)-1, -1).astype(np.float64)
        n1, n2 = len(x1), len(x2)
        vol = np.pi * rs**2 * 2.0
        dvol = np.diff(vol)
        exp_n = dvol * n2 / l_box**3
        xi = n / n1 / exp_n[:, None] - 1.
        wp = np.sum(xi * 2., axis=1)
        lg_rs, lg_rs_c = self.__cvt_r(rs)
        return self.Result({
            'rs': rs, 'lg_rs': lg_rs, 'lg_rs_c': lg_rs_c, 'pimax': pimax,
            'xi': xi, 'wp': wp
        })

    def min_dist_2d(self, x1, x2, dz_max: float, root=True):
        l_box = self.l_box
        id = KDTree(x2[:, [2]]).query_radius(
            x1[:, [2]], dz_max, return_distance=False)
        x1_2d, x2_2d = x1[:, [0, 1]], x2[:, [0, 1]]
        dist = []
        for _x1_2d, _id in zip(x1_2d, id):
            d = l_box if _id.size == 0 else np.linalg.norm(
                _x1_2d - x2_2d[_id], axis=1).min()
            dist.append(d)
        if not root:
            return np.array(dist)
        x2_shift = x2.copy()
        x2_shift[:, 2] += l_box
        d_shift = self.min_dist_2d(x1, x2_shift, dz_max, root=False)
        dist = np.where(dist < d_shift, dist, d_shift)

        x2_shift[:, 2] -= 2*l_box
        d_shift = self.min_dist_2d(x1, x2_shift, dz_max, root=False)
        dist = np.where(dist < d_shift, dist, d_shift)

        return dist

    def iter_periodic_3d(self, x):
        '''
        Do not modify returned value.
        '''
        l_box = self.l_box
        id_list = 0, -1, 1
        y = x.copy()
        for off0 in id_list:
            y[:, 0] = x[:, 0] + off0 * l_box
            for off1 in id_list:
                y[:, 1] = x[:, 1] + off1 * l_box
                for off2 in id_list:
                    y[:, 2] = x[:, 2] + off2 * l_box
                    yield y

    def min_dist_3d(self, x1, x2, root=True):
        l_box = self.l_box

        if not root:
            dist, _ = KDTree(x2).query(x1, return_distance=True)
            return dist[:, 0]

        dist = np.full(len(x1), l_box)
        for di0 in (-1, 0, 1):
            y0 = x2[:, 0] + di0 * l_box
            for di1 in (-1, 0, 1):
                y1 = x2[:, 1] + di1 * l_box
                for di2 in (-1, 0, 1):
                    y2 = x2[:, 2] + di2 * l_box
                    y = np.column_stack((y0, y1, y2))
                    dist_shift = self.min_dist_3d(x1, y, root=False)
                    dist = np.where(dist < dist_shift, dist, dist_shift)
        return dist

    def min_dist_3d_by_kdt(self, kdt: KDTree, x2):
        out = np.full(len(x2), self.l_box)
        for y2 in self.iter_periodic_3d(x2):
            d, _ = kdt.query(y2, return_distance=True)
            d = d[:, 0]
            sel = d < out
            out[sel] = d[sel]
        return out

    def __cvt_r(self, rs: np.ndarray):
        lg_rs = Num.safe_lg(rs)
        lg_rs_c = .5 * (lg_rs[1:] + lg_rs[:-1])
        return lg_rs, lg_rs_c


class Clustering:
    '''
    For periodic box.
    '''

    def __init__(self, hset_ref: SubhaloSet) -> None:

        self.impl = _Clustering(hset_ref.sim_info)
        self.x_ref = hset_ref['x']
        self.rng = hset_ref.ctx.rng

    def proj_cross(self, hset: SubhaloSet,
                   rs=np.logspace(-2., 1.2, 21),
                   pimax=10., n_threads=1, n_repeat=50,
                   n_max_rand=None, bootstrap_kw={}):

        x, x_ref = hset['x'], self.x_ref
        return self.proj_cross_of(x, x_ref, rs, pimax, n_threads,
                                  n_repeat, n_max_rand, bootstrap_kw)

    def proj_cross_of(self, x: np.ndarray, x_ref: np.ndarray,
                      rs=np.logspace(-2., 1.2, 21),
                      pimax=10., n_threads=1, n_repeat=50, n_max_rand=None,
                      bootstrap_kw={}):
        '''
        @x: gets bootstrapped.
        '''
        def fn(_x: np.ndarray):
            return self.impl.proj_cross(_x, x_ref, rs, pimax, n_threads)
        if n_repeat <= 1:
            return fn(x)

        if n_max_rand is not None:
            n_max_rand = [n_max_rand]

        return Bootstrap.resampled_call(
            fn, dsets_in=[{'_x': x}],
            keys_out=['wp'], n_resample=n_repeat, rng=self.rng,
            dsets_max_sizes=n_max_rand,
            **bootstrap_kw
        )

    def cross(self, hset: SubhaloSet,
              rs: np.ndarray = np.logspace(-2., 1.2, 21),
              n_threads=1, n_repeat=50,
              n_max_rand=None):
        x, x_ref = hset['x'], self.x_ref

        def fn(_x: np.ndarray):
            return self.impl.cross(_x, x_ref, rs, n_threads=n_threads)
        if n_repeat <= 1:
            return fn(x)

        if n_max_rand is not None:
            n_max_rand = [n_max_rand]

        return Bootstrap.resampled_call(
            fn, dsets_in=[{'_x': x}],
            keys_out=['xi'],
            n_resample=n_repeat, rng=self.rng, dsets_max_sizes=n_max_rand)

    @cached_property
    def kdt(self):
        return KDTree(self.x_ref, leaf_size=16)

    def counts_in_r(self, x: np.ndarray, r: float, period=True):
        assert x.ndim == 2
        kdt = self.kdt
        counts = kdt.query_radius(x, r=r, count_only=True)

        if not period:
            return counts

        l_box = self.impl.l_box
        for d0 in -1, 0, 1:
            x0 = x[:, 0] + d0 * l_box
            sel = (x0 >= -r) & (x0 <= l_box + r)
            ids0 = sel.nonzero()[0]
            if ids0.size == 0:
                continue
            for d1 in -1, 0, 1:
                x1 = x[ids0, 1] + d1 * l_box
                sel = (x1 >= -r) & (x1 <= l_box + r)
                ids1 = ids0[sel]
                if ids1.size == 0:
                    continue
                for d2 in -1, 0, 1:
                    if d0 == d1 == d2 == 0:
                        continue
                    x2 = x[ids1, 2] + d2 * l_box
                    sel = (x2 >= -r) & (x2 <= l_box + r)
                    ids2 = ids1[sel]
                    if ids2.size == 0:
                        continue
                    x_new = x[ids2] + np.array([d0, d1, d2]) * l_box
                    counts_new = kdt.query_radius(x_new, r=r, count_only=True)
                    counts[ids2] += counts_new
        return counts

    def min_dist_3d(self, hset: SubhaloSet):
        return self.impl.min_dist_3d_by_kdt(self.kdt, hset['x'])


class ObsClustering:

    Result = DataDict[str, np.ndarray]

    def __init__(self, hset_ref: SubhaloSet) -> None:
        ra, dec, d_c = hset_ref['ra', 'dec', 'd_c']

        self.x_ref = np.column_stack([ra, dec, d_c]).astype(np.float32)
        self.l_box = hset_ref.sim_info.full_box_size
        self.rng = hset_ref.ctx.rng

    def proj_cross(self, hset: SubhaloSet,
                   rs: np.ndarray = np.logspace(-2., 1.2, 21),
                   pimax=10., n_threads=1, n_repeat=50,
                   bootstrap_kw={}):
        '''
        @hset: 
            get bootstraped. 
            ra, dec, d_c, wgt will be used.
        '''
        x = np.column_stack(hset['ra', 'dec', 'd_c', 'wgt']).astype(np.float32)
        x_ref = self.x_ref

        def fn(_x: np.ndarray):
            return self.__proj_cross_impl(_x, x_ref, rs, pimax,
                                          n_threads, self.l_box)
        if n_repeat <= 1:
            return fn(x)

        return Bootstrap.resampled_call(
            fn, dsets_in=[{'_x': x}],
            keys_out=['wp'],
            n_resample=n_repeat, rng=self.rng,
            **bootstrap_kw)

    def __proj_cross_impl(self, x1, x2, rs, pimax, n_threads, l_box):
        
        import Corrfunc
        
        kw = {
            'autocorr': False,
            'cosmology': 2,
            'nthreads': n_threads,
            'pimax': pimax,
            'binfile': np.asarray(rs),
            'RA1': x1[:, 0],
            'DEC1': x1[:, 1],
            'CZ1': x1[:, 2],
            'weights1': x1[:, 3],
            'RA2': x2[:, 0],
            'DEC2': x2[:, 1],
            'CZ2': x2[:, 2],
            'is_comoving_dist': True,
            'weight_type': 'pair_product',
        }
        out = Corrfunc.mocks.DDrppi_mocks(**kw)
        n = out['weightavg'] * out['npairs']
        n = n.reshape(len(rs)-1, -1)
        n1, n2 = kw['weights1'].sum(), len(kw['RA2'])

        vol = np.pi * rs**2 * 2.0
        dvol = np.diff(vol)
        exp_n = dvol * n2 / l_box**3
        xi = n / n1 / exp_n[:, None] - 1.
        wp = np.sum(xi * 2., axis=1)
        lg_rs, lg_rs_c = self.__cvt_r(rs)

        return self.Result({
            'rs': rs, 'lg_rs': lg_rs, 'lg_rs_c': lg_rs_c, 'pimax': pimax,
            'xi': xi, 'wp': wp
        })

    def __cvt_r(self, rs: np.ndarray):
        lg_rs = Num.safe_lg(rs)
        lg_rs_c = .5 * (lg_rs[1:] + lg_rs[:-1])
        return lg_rs, lg_rs_c


class ClusteringOfSubhaloSets:
    def __init__(self, samp: SubhaloSet,
                 ref_samp: SubhaloSet, r_range=[2., 10.], rng=0) -> None:

        x, x_ref = samp['x'], ref_samp['x']
        cls = _Clustering(samp.sim_info)
        xi_all_ref = cls.cross_in_rad(x, x_ref, *r_range).xi

        self.samp = samp
        self.x, self.x_ref = x, x_ref
        self.r_range = r_range
        self.rng = Rng(rng)
        self.cls = cls
        self.xi_all_ref = xi_all_ref
        self.red = reduction.Errorbar('median+1sigma')

    def find_bias_with_binnings(self, key='S', es=[0., 1.], n_workers=8,
                                n_resample=50):
        n = len(es) - 1

        def f(lo, hi):
            return self.find_bias_with_binning(
                key, lo, hi, n_resample=n_resample)
        outs = [
            f(es[i], es[i+1]) for i in range(n)]
        return np.array(outs)

    def find_bias_with_binning(self, key='S', lo=0., hi=1., n_resample=50):
        val, x = self.samp[key], self.x
        n = len(val)
        id_list = self.rng.choice(n, (n_resample, n))
        bs = []
        for id in id_list:
            val_re = val[id]
            sel = (val_re >= lo) & (val_re < hi)
            id_sel = id[sel]
            x_sel = x[id_sel]
            bs.append(self.find_bias(x_sel))
        b_out = self.red(bs)
        sel = (val >= lo) & (val < hi)
        val_out = np.median(val[sel])
        return np.concatenate([[val_out], b_out])

    def find_bias(self, x, **kw):
        return self.cls.cross_in_rad(
            x, self.x_ref, *self.r_range, **kw).xi / self.xi_all_ref


class SubSample(abc.HasDictRepr):
    '''
    Wrap a subsample for clustering analysis.
    '''

    repr_attr_keys = ('samp', 'l_box', 'add_rsd', 'header')

    def __init__(self, samp_parent: SubhaloSet, key: str = None,
                 lo: float = None, hi: float = None,
                 add_rsd=False) -> None:

        if lo is None and hi is None:
            samp = samp_parent
        else:
            samp = samp_parent.subset_by_value(key, lo=lo, hi=hi)

        f_n = samp.n_objs / samp_parent.n_objs

        self.samp = samp
        self.l_box = samp.sim_info.full_box_size
        self.add_rsd = add_rsd

        self.header = DataDict({
            'f_n': f_n, 'n_objs': samp.n_objs
        })

    @classmethod
    def from_percent(cls, samp_parent: SubhaloSet, key: str = None,
                     p_lo=None, p_hi=None, **kw):

        q_lo = None if p_lo is None else \
            np.quantile(samp_parent[key], p_lo)
        q_hi = None if p_hi is None else \
            np.quantile(samp_parent[key], p_hi)
        return cls(samp_parent, key=key, lo=q_lo, hi=q_hi, **kw)

    @classmethod
    def list_from_percents(cls, samp_parent: SubhaloSet,
                           key: str, ps: np.ndarray, **kw):
        p_los, p_his = ps[:-1], ps[1:]
        return [cls.from_percent(samp_parent, key=key, p_lo=p_lo, p_hi=p_hi, **kw)
                for p_lo, p_hi in zip(p_los, p_his)]

    @classmethod
    def list_from_bins(cls, samp_parent: SubhaloSet,
                       key: str, edges: np.ndarray, **kw):
        los, his = edges[:-1], edges[1:]
        return [cls(samp_parent, key=key, lo=lo, hi=hi, **kw)
                for lo, hi in zip(los, his)]

    @cached_property
    def cls_ref(self):
        return Clustering(self.samp)

    @cached_property
    def obs_cls_ref(self):
        return ObsClustering(self.samp)

    @cached_property
    def x_boxed(self):
        l_box = self.l_box
        x = self.x_pos
        x[x < 0.] += l_box
        x[x >= l_box] -= l_box
        assert np.all(0. <= x) and np.all(x < l_box)
        return x

    @cached_property
    def x_pos(self):
        return self.x_added_rsd if self.add_rsd else self.samp['x']

    @cached_property
    def x_added_rsd(self):
        samp = self.samp

        H0 = 100.0                   # (km/s)/(Mpc/h)
        dim = 2
        dx = samp['v'][:, dim] / H0  # Mpc/h

        x = samp['x'].copy()
        x[:, dim] += dx

        return x

    def iter_periodic(self, x):
        l_box = self.l_box
        id_list = 0, -1, 1
        y = x.copy()
        for off0 in id_list:
            y[:, 0] = x[:, 0] + off0 * l_box
            for off1 in id_list:
                y[:, 1] = x[:, 1] + off1 * l_box
                for off2 in id_list:
                    y[:, 2] = x[:, 2] + off2 * l_box
                    yield y

    def x_added_periodic(self, dx_max: np.ndarray):
        l_box = self.l_box
        dx_max = np.array(dx_max)
        x_lb, x_ub = 0.0 - dx_max, l_box + dx_max

        xs, ids = [], []
        for x in self.iter_periodic(self.x_boxed):
            sel = ((x >= x_lb) & (x < x_ub)).all(1)
            id = sel.nonzero()[0]
            ids.append(id)
            xs.append(x[id])
        xs = np.concatenate(xs, axis=0)
        ids = np.concatenate(ids)
        return ids, xs

    def obs_proj_cross_with(self, subs: SubSample, pimax=10.,
                            n_repeat=5, n_threads=4,
                            **corr_kw):
        '''
        @self: reference sample
        @subs: get bootstraped
        '''
        return self.obs_cls_ref.proj_cross(
            subs.samp, pimax=pimax, n_threads=n_threads,
            n_repeat=n_repeat, **corr_kw) | self.header

    def proj_cross_with(self, subs: SubSample, pimax=10.,
                        n_repeat=5, n_threads=4, **corr_kw):
        '''
        @self: reference sample
        @subs: get bootstraped
        '''
        x, x_ref = subs.x_boxed, self.x_boxed
        return self.cls_ref.proj_cross_of(
            x, x_ref, pimax=pimax, n_threads=n_threads,
            n_repeat=n_repeat, **corr_kw) | self.header

    def min_dist_3d_with(self, subs: SubSample):
        '''
        self = reference sample
        '''
        return self.cls_ref.min_dist_3d(subs.samp)

    def min_dist_2d_with(self, subs: SubSample,
                         max_dx_proj=15.,
                         max_dv_h=3.0):
        ref_samp = self.samp
        m_h = ref_samp['m_mean200']
        ht = ref_samp.sim_info.cosmology.halo_theory
        u_v = ref_samp.sim_info.cosmology.unit_system.u_v_to_kmps
        vir = ht.vir_props_mean(m_h, z=0., f=200.)
        v_h = vir.v * u_v   # [physical km/s]
        r_h = vir.r         # [comoving Mpc/h]

        H0 = 100.0   # km/s/(Mpc/h)

        max_dx_los = v_h * max_dv_h / H0
        max_dx = np.array([max_dx_proj, max_dx_proj, max_dx_los.max()])
        print(max_dx)

        id_ref, x_ref = self.x_added_periodic(max_dx)
        max_dx_los_ref = max_dx_los[id_ref]
        r_h_ref = r_h[id_ref]
        x_ref_scaled = x_ref / max_dx
        kdt = KDTree(x_ref_scaled)

        x = subs.x_boxed
        x_scaled = x / max_dx

        nn_ids_list = kdt.query_radius(x_scaled, r=np.sqrt(3))
        print('Found nn list')

        min_ds = np.zeros(len(x), dtype=float) + 1.0e6
        r_hs = min_ds.copy()
        for i, nn_ids in enumerate(nn_ids_list):
            if nn_ids.size == 0:
                continue
            _x, _x_ref, _max_dx_los_ref = x[i], x_ref[nn_ids], max_dx_los_ref[nn_ids]
            _r_h_ref = r_h_ref[nn_ids]

            sel = np.abs((_x[2] - _x_ref)[:, 2]) <= _max_dx_los_ref
            sel = sel.nonzero()[0]
            if sel.size == 0:
                continue
            _x_sel, _r_h_sel = _x_ref[sel], _r_h_ref[sel]
            _ds = np.linalg.norm(_x[:2] - _x_sel[:, :2], axis=1)
            i_found = (_ds / _r_h_sel).argmin()
            min_ds[i] = _ds[i_found]
            r_hs[i] = _r_h_sel[i_found]

        return min_ds, r_hs


class SubSampleManager:
    def __init__(self, sim_info: sims.SimInfo, data_file,
                 m_s_range_in_sol=(10**8.0, 10**9.),
                 m_h_range_in_sol=(10**10.5, 10**11.0),
                 z_lc_lb=15.,
                 ssfr_lb_in_sol=10**-11.0,
                 ssfr_ub_in_sol=1.0e10,
                 baryon_excl_ej=False,
                 m_h_ref_lb = None,
                 m_s_ref_lb = None,
                 extra_mask_dwarfs = {},
                 rng=10086, add_rsd=False, dmo=False) -> None:

        h = sim_info.cosmology.hubble
        u_m = 1.0e10 / h
        u_ssfr = h / 1.0e9
        m_s_lb, m_s_ub = m_s_range_in_sol
        m_s_lb, m_s_ub = m_s_lb / u_m, m_s_ub / u_m
        m_h_lb, m_h_ub = m_h_range_in_sol
        m_h_lb, m_h_ub = m_h_lb / u_m, m_h_ub / u_m
        ssfr_lb = ssfr_lb_in_sol / u_ssfr
        ssfr_ub = ssfr_ub_in_sol / u_ssfr
        
        m_h_ref_lb = m_h_ref_lb / u_m if m_h_ref_lb is not None else m_h_lb
        m_s_ref_lb = m_s_ref_lb / u_m if m_s_ref_lb is not None else m_s_lb

        self.sim_info = sim_info
        self.data_file = data_file
        self.m_s_lb = m_s_lb
        self.m_s_ub = m_s_ub
        self.m_h_lb = m_h_lb
        self.m_h_ub = m_h_ub
        self.z_lc_lb = z_lc_lb
        self.ssfr_lb = ssfr_lb
        self.ssfr_ub = ssfr_ub
        self.bar_excl_ej = baryon_excl_ej
        self.m_h_ref_lb = m_h_ref_lb
        self.m_s_ref_lb = m_s_ref_lb
        self.extra_mask_dwarfs: dict = extra_mask_dwarfs
        self.rng = Rng(rng)
        self.add_rsd = add_rsd
        self.dmo = dmo

    @cached_property
    def s_all(self):
        return SimSample.load(self.sim_info, self.data_file)
    
    @cached_property
    def s_c(self):
        return self.s_all.subset_by_value('is_c', eq=True)

    @cached_property
    def s_dwarf(self):
        s_c = self.s_c
        if self.dmo:
            s = s_c\
                .subset_by_value('m_mean200', lo=self.m_h_lb, hi=self.m_h_ub)\
                .subset_by_value('last_sat_z', lo=self.z_lc_lb)
        else:
            s = s_c\
                .subset_by_value('m_star', lo=self.m_s_lb, hi=self.m_s_ub)\
                .subset_by_value('ssfr', lo=self.ssfr_lb, hi=self.ssfr_ub)
            if self.bar_excl_ej:
                s = s.subset_by_value('last_sat_z', lo=self.z_lc_lb)
        for k, v in self.extra_mask_dwarfs.items():
            s = s.subset_by_value(k, **v)
        return s

    def s_diffuse_dwarf(self, S_lb=None, S_ub=None, S_key='Sigma_star'):
        assert not self.dmo
        s = self.s_dwarf.subset_by_value(S_key, lo=S_lb, hi=S_ub)
        return s

    @cached_property
    def s_halo(self):
        return self.s_c

    @cached_property
    def s_ref(self):
        s_all = self.s_all
        if self.dmo:
            s = s_all.subset_by_value('m_peak', lo=self.m_h_ref_lb)
        else:
            s = s_all.subset_by_value('m_star', lo=self.m_s_ref_lb)
        return s

    @cached_property
    def subs_ref(self):
        return SubSample(self.s_ref, add_rsd=self.add_rsd)

    


default_ps = [0.00, 0.02, 0.33, 0.67, 0.98, 1.00]


class SimCross:

    def __init__(self, sman: SubSampleManager) -> None:

        self.sman = sman

    def find_wp(self, subs: SubSample, pimax=10., n_threads=4,
                n_repeat=20, **corr_kw):
        out = self.sman.subs_ref.proj_cross_with(
            subs, pimax=pimax, n_threads=n_threads, n_repeat=n_repeat,
            **corr_kw)
        return out

    def wp_of(self, samp_name='s_dwarf', sub_key=None, ps=None, qs=None,
              subsamp_sel = [], **corr_kw):
        sman = self.sman

        s: SimSample = getattr(sman, samp_name)
        for sel in subsamp_sel:
            s = s.subset_by_value(**sel)
        s = self.__sub_sample(s, sub_key=sub_key, ps=ps, qs=qs)
        s = SubSample(s, add_rsd=sman.add_rsd)
        n_objs = s.samp.n_objs
        print(f'Start {samp_name=}, {sub_key=}, {ps=}, {qs=}, {n_objs=}')
        outs = self.find_wp(s, **corr_kw)
        print('Done')

        if sub_key is not None:
            val = s.samp[sub_key]
            outs['sub_sample'] = {
                'key': sub_key,
                'min': val.min(), 'max': val.max(),
                'median': np.median(val),
                'mean': val.mean(),
                'stddev': val.std(),
                '1sigma': np.percentile(val, [16., 84.]),
                'n_subhalos': len(val),
            }

        return outs

    def wps_of(self, samp_name='s_dwarf', sub_key='z_half',
               ps=default_ps, qs=None, **corr_kw):
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
            outs[str(i)] = self.wp_of(
                samp_name=samp_name,
                sub_key=sub_key, ps=ps[i], qs=qs[i],
                **corr_kw)
        return outs

    def __sub_sample(self, s: SubhaloSet,
                     sub_key=None, ps=None, qs=None):
        if sub_key is not None:
            if ps is not None:
                assert qs is None
                qs = np.quantile(s[sub_key], ps)
            q_lo, q_hi = qs
            s = s.subset_by_value(sub_key, lo=q_lo, hi=q_hi)
        return s


'''
Deprecated - old implementations.
'''


def get_cls_binned(s_parent: SubhaloSet,
                   cls_ref: Clustering,
                   key='S',
                   ranges=[(0., 1.0e6),
                           (0., 7.),
                           (7., 15.),
                           (15., 25.),
                           (25., 1.0e6)],
                   n_rep=5, n_threads=4):
    outs = {}
    for i, (vlo, vhi) in enumerate(ranges):
        s = s_parent.subset_by_value(key, lo=vlo, hi=vhi)
        out = cls_ref.proj_cross(s, n_repeat=n_rep, n_threads=n_threads) | {
            'n_objs': s.n_objs, 'f_n': s.n_objs / s_parent.n_objs
        }
        outs[str(i)] = out
    return outs


def get_cls_by_quantiles(s_parent: SubhaloSet,
                         cls_ref: Clustering,
                         key='S',
                         q_ranges=[
                             (0.0, 1.0),
                             (0.0, 0.02),
                             (0.02, 0.33),
                             (0.33, 0.67),
                             (0.67, 0.98),
                             (0.98, 1.00),
                         ],
                         n_rep=5, n_threads=4):
    val = s_parent[key]
    ranges = [
        np.quantile(val, q) for q in q_ranges
    ]
    return get_cls_binned(s_parent, cls_ref, key=key,
                          ranges=ranges, n_rep=n_rep, n_threads=n_threads)


'''
@dataclass
class _CountRResult:
    rs: np.ndarray
    counts: np.ndarray
    expected_counts: np.ndarray

    def xis_betweens(self, inds=[0, 1]):
        cnts, exp_cnts = self.counts, self.expected_counts
        i0, i1 = inds
        return (cnts[:, i1] - cnts[:, i0]) / (exp_cnts[i1] - exp_cnts[i0]) - 1.


class Clustering(abc.HasDictRepr):

    repr_attr_keys = ('samp', 'l_box')

    def __init__(self, samp: SubhaloSet) -> None:
        super().__init__()

        sim_info = samp.sim_info
        l_box = sim_info.box_size

        self.samp = samp
        self.sim_info = sim_info
        self.l_box = l_box

    def acf(self, rs=np.linspace(.3, 1.2, 10),
            n_resample: int | None = 100, impl_kw={}):
        x, l_box, rng = self.samp['x'], self.l_box, self.samp.rng
        rs = np.array(rs)

        def f(x):
            return self._acf(x, l_box=l_box, rs=rs, impl_kw=impl_kw)
        if n_resample is None:
            out = f(x) | {
                'xi_sd': 0.,
            }
        else:
            out = Bootstrap.resampled_call(
                stats_fn=f, dsets_in=[{'x': x},], keys_out=['xi', 'lg_xi'],
                n_resample=n_resample,
                rng=rng
            )
        rs_c = 0.5*(rs[1:] + rs[:-1])
        lg_rs_c = Num.safe_lg(rs_c)
        out |= {
            'rs': rs, 'rs_c': rs_c, 'lg_rs_c': lg_rs_c,
        }
        return out

    @staticmethod
    def _acf(
            x: np.ndarray, l_box: float, rs=np.linspace(0.3, 1.2, 10),
            impl_kw={}):
        impl_kw = {'nthreads': 1, **impl_kw}
        cf = Corrfunc.theory.xi(l_box, binfile=rs,
                                X=x[:, 0], Y=x[:, 1], Z=x[:, 2], **impl_kw)
        xi = cf['xi']
        lg_xi = Num.safe_lg(xi)
        return DataDict({
            'xi': xi, 'lg_xi': lg_xi
        })

    def count_r(self, rs: np.ndarray):
        samp, l_box = self.samp, self.l_box
        rs = np.array(rs)
        r_max = rs.max()
        assert ((rs >= 0.) & (rs < l_box)).all()

        x = samp['x']
        x_p = self._x_periodic(x, l_box, r_max)
        kdt = KDTree(x_p, leaf_size=32)
        cnts = []
        for r in rs:
            cnt = kdt.query_radius(x, r, count_only=True).astype(float)
            cnts.append(cnt)
        cnts = np.column_stack(cnts)

        vols = 4./3. * np.pi * rs**3
        vol_box = l_box**3
        exp_cnts = len(x) / vol_box * vols

        return _CountRResult(rs, cnts, exp_cnts)

    def count_r_by_ref(self, rs: np.ndarray, ref: Clustering):
        assert self.sim_info is ref.sim_info

        samp, l_box = self.samp, self.l_box
        rs = np.array(rs)
        r_max = rs.max()
        assert ((rs >= 0.) & (rs < l_box)).all()

        x, x_ref = samp['x'], ref.samp['x']
        n_refs = len(x_ref)
        kdt = KDTree(self._x_periodic(x_ref, l_box, r_max), leaf_size=32)
        cnts = []
        for r in rs:
            cnt = kdt.query_radius(x, r, count_only=True).astype(float)
            cnts.append(cnt)
        cnts = np.column_stack(cnts)

        vols = 4./3. * np.pi * rs**3
        vol_box = l_box**3
        exp_cnts = n_refs / vol_box * vols

        return _CountRResult(rs, cnts, exp_cnts)

    @staticmethod
    def _x_periodic(x: np.ndarray, l_box: float, r: float):
        x_lb, x_ub = -r, l_box + r
        x_out = []
        for di0 in range(-1, 2):
            x0 = x[:, 0] + l_box * di0
            sel0 = (x0 >= x_lb) & (x0 < x_ub)
            for di1 in range(-1, 2):
                x1 = x[:, 1] + l_box * di1
                sel1 = (x1 >= x_lb) & (x1 < x_ub)
                for di2 in range(-1, 2):
                    x2 = x[:, 2] + l_box * di2
                    sel2 = (x2 >= x_lb) & (x2 < x_ub)
                    sel = sel0 & sel1 & sel2
                    x_new = np.column_stack((x0, x1, x2))[sel]
                    x_out.append(x_new)
        return np.concatenate(x_out, axis=0)


class Bias(abc.HasDictRepr):

    repr_attr_keys = ('xi', 'xi_sd', 'b', 'b_sd')

    def __init__(self, samp: SubhaloSet,
                 ref_bias: Bias = None,
                 r_range: tuple[float, float] = [2., 10.],
                 n_resample: int = 100) -> None:
        super().__init__()

        acf = Clustering(samp).acf(r_range, n_resample=n_resample)
        xi, xi_sd = acf['xi'][0], acf['xi_sd'][0]

        if ref_bias is not None:
            b, b_sd = xi / ref_bias.xi, xi_sd / ref_bias.xi
        else:
            b, b_sd = None, None

        self.xi = xi
        self.xi_sd = xi_sd
        self.b = b
        self.b_sd = b_sd
'''
