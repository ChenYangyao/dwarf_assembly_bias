from __future__ import annotations
import typing
from typing import Tuple, Any, Self
if typing.TYPE_CHECKING:
    from .clustering import Clustering
from .config import ProjPaths
from pyhipp.core import abc, DataDict, Num, DataTable
from pyhipp.io import h5
from pyhipp.stats import Rng, Bootstrap
from pyhipp.field.cubic_box import TidalClassifier
import numpy as np
from pyhipp_sims import sims
from scipy.spatial import KDTree
from dataclasses import dataclass


class SubhaloSet(sims.abc.SimObjectCatalog):
    pass

    @classmethod
    def from_data(cls, objs: DataTable, header: DataDict,
                  sim_info: sims.SimInfo, **kw) -> Self:
        ctx = sims.abc.SimContext(sim_info)
        return cls(objs=objs, header=header, ctx=ctx, **kw)

    @classmethod
    def from_cosmic_web(cls, clsf: TidalClassifier,
                        sim_info: sims.SimInfo, web_type='filament'):
        x = clsf.grid_points_of_web_type(web_type)
        return cls.from_data(DataTable({
            'x': x}), DataDict(), sim_info)


class SimSample(SubhaloSet):

    @classmethod
    def load(cls, sim_info: sims.SimInfo, in_path='z0.hdf5', **init_kw) -> Self:
        in_path = ProjPaths.sim_dir_of(sim_info) / in_path
        out = cls.from_h5_file(in_path, **init_kw)
        assert out.sim_info is sim_info

        subhs = out.objs
        if 'x' in subhs:
            x = subhs['x']
            l_box = sim_info.full_box_size
            x[x < 0.] += l_box
            x[x >= l_box] -= l_box
            assert np.all((x >= 0.) & (x < l_box))

        out.is_root_dataset = True
        out.in_path = in_path

        return out

    @property
    def z_dst(self) -> float:
        return self.header['z_dst']

    def requires_cosmic_type(self, clsf: TidalClassifier) -> Self:
        assert np.abs(clsf.mesh.l_box - self.sim_info.full_box_size) < 1.0e-6
        x = self['x']
        cosmic_type = clsf.web_type_at(x)

        for k, v in clsf.web_types.items():
            self.objs[f'is_{k}'] = cosmic_type == v
        self.objs |= {'cosmic_type': cosmic_type}

        return self

    def requires_dist_to(self, stored_key: str,
                         x_dst: np.ndarray,
                         k_th=1) -> Self:
        x = self['x']
        l_box = self.sim_info.full_box_size
        d, _ = KDTree(x_dst, leafsize=16, boxsize=l_box).query(
            x, k=k_th, workers=8)
        if k_th > 1:
            assert d.shape == (len(x), k_th)
            d = d[:, k_th-1]
        self.objs[stored_key] = d
        return self

    def requires_dist_to_halos(
            self, m_h_lbs: float = 10 ** np.array([2., 2.5, 3.]),
            keys: str = ['h12', 'h125', 'h13']):
        assert (np.diff(m_h_lbs) >= 0.).all()
        x = self['x']
        samp = self.subset_by_value('is_c', eq=True)
        out = DataDict()
        for m_h_lb, key in zip(m_h_lbs, keys):
            samp = samp.subset_by_value('m_h', lo=m_h_lb)
            x_h = samp['x']
            d, _ = KDTree(x_h).query(x)
            out[f'd_{key}'] = d

        self.objs |= out

        return self

    def requires_spin(self) -> Self:
        subhs, z_dst = self.objs, self.z_dst
        ht = self.sim_info.cosmology.halo_theory
        a = 1. / (1. + z_dst)
        suf_map = {
            'crit200': (ht.rho_vir_crit, 'h'),
            'mean200': (ht.rho_vir_mean, 'mean200')
        }

        for k_in, (fn_rho, k_out) in suf_map.items():
            m_h = subhs[f'm_{k_in}']
            spin, spin_form = subhs['spin', 'spin_form']
            m_h = m_h.clip(1.0e-3)                      # 1.0e10 Msun/h
            rho_h = fn_rho(z=z_dst)
            r_h = ht.r_vir(m_h, rho_h)
            r_h_phy = r_h*a
            v_h = ht.v_vir(m_h, r_h_phy, to_kmps=True)  # physical km/s
            j_h = np.sqrt(2.) * v_h * r_h_phy
            spin = Num.norm(spin, 1) / j_h
            spin_form = Num.norm(spin_form, 1) / j_h

            subhs |= {
                f'spin_{k_out}': spin,
                f'spin_form_{k_out}': spin_form,
                f'v_{k_out}': v_h,
                f'm_{k_out}': m_h,
                f'r_{k_out}': r_h}
        return self

    def require_c(self):
        ctx, subhs = self.ctx, self.objs
        pf = ctx.sim_info.cosmology.halo_theory.nfw_profile
        for suf in 'h', 'mean200':
            v_max, v_h = subhs['v_max', f'v_{suf}']
            v_max2h = Num.safe_div(v_max, v_h)
            v_max2h = v_max2h.clip(1.01)
            c = pf.v_max2vir_to_c(v_max2h)
            subhs |= {f'c_{suf}': c}
        return self

    def requires_disk(self, sham_key='v_peak_m', spin_key='spin') -> Self:
        
        from dynamic_hotness.observations import Behroozi2019MsToMh
        
        subhs, z_dst = self.objs, self.z_dst
        b19 = Behroozi2019MsToMh()

        m_h = subhs[sham_key].clip(1.0e-3)
        m_s = b19.find_m_s(m_h, z=z_dst)

        spin, m_h, v_h = subhs[spin_key, 'm_h', 'v_h']
        m_h = m_h.clip(1.0e-3)
        spin = spin.clip(1.0e-3)
        hubble = self.sim_info.cosmology.hubble
        f_d = Num.safe_div(m_s, m_h)
        Dz = 2.0**0.5
        S = 207. * hubble * (f_d / 0.05) * (spin / 0.05)**(-2.0) * (v_h / 200.) * Dz

        r_v_max, r_h = subhs['r_v_max', 'r_h']
        c = 2.162582 * r_h / r_v_max
        p1 = -0.06 + 2.71 * f_d + 0.0047 / spin
        p2 = 1. - 3. * f_d + 5.2 * f_d**2
        p3 = 1. - 0.019 * c + 0.00025 * c**2 + 0.52 / c
        fR = (spin / 0.1)**p1 * p2 * p3
        fR = 1. / fR
        S_fR = S * fR**2

        subhs |= {
            'm_s_pred': m_s,
            'S': S, 'S_fR': S_fR,
        }
        return self

    def requires_assembly(self) -> Self:
        subhs, sim_info, rng = self.objs, self.sim_info, self.ctx.rng
        is_b = sim_info.is_baryon

        s_f, s_af = subhs['snap_form', 'snap_afterform']
        z_f, z_af = sim_info.redshifts[s_f], sim_info.redshifts[s_af]
        z_half = rng.uniform(low=z_af, high=z_f, size=len(z_f))
        x_t_half = np.log10(1.+z_half)
        subhs |= {
            'z_half': z_half, 'x_t_half': x_t_half,
        }

        s_lt = subhs['last_sat_snap'].clip(0)
        rv = rng.uniform(0., 0.001, size=len(s_lt))
        lt_z = self.sim_info.redshifts[s_lt] + rv
        subhs['last_sat_z'] = lt_z

        v_max, v_peak, v_h = subhs['v_max', 'v_peak', 'v_h']
        rv = self.ctx.rng.uniform(0., 0.01, size=len(v_max))
        v_max2h = Num.safe_div(v_max, v_h)
        v_peak2max = Num.safe_div(v_peak, v_max) + rv
        subhs |= {
            'v_max2h': v_max2h, 'v_peak2max': v_peak2max,
        }

        if is_b:
            m_s, r_s, m_h = subhs['m_star', 'r_half_mass_star', 'm_h']
            a = 1. / (1. + self.z_dst)
            h = self.sim_info.cosmology.hubble

            _r = (r_s * a).clip(1.0e-5) * 1.0e6 / h       # pc
            _m = m_s * 1.0e10 / h                               # Msun
            Sigma_s = (0.5 * _m) / (np.pi * _r**2)
            subhs['Sigma_star'] = Sigma_s

            m_s2h = Num.safe_div(m_s, m_h)
            subhs['m_s2h'] = m_s2h

            ssfr = Num.safe_div(subhs['sfr'], m_s)
            subhs |= {'ssfr': ssfr}

        return self

    def require_LGalaxies_H15(self) -> Self:
        mod_path = str(self.in_path) + '.LGalaxies-H15'
        mod = h5.File.load_from(mod_path, 'objs')
        subhs = self.objs
        assert mod['StellarBulgeMass'].size == subhs['is_c'].size
        m_s_b, m_s_d, r_d, sfr = mod[
            'StellarBulgeMass', 'StellarDiskMass', 'StellarDiskRadius',
            'StarFormationRate']
        m_s = m_s_b + m_s_d
        invalid = m_s < 1.0e-10
        m_s.clip(1.0e-10, out=m_s)
        r_d.clip(1.0e-10, out=r_d)
        sfr.clip(1.0e-10, out=sfr)
        
        r_scale = r_d / 3.
        r_eff =  r_scale * 1.678
        A_d = np.pi * r_eff**2

        h = self.sim_info.cosmology.hubble
        u_m = 1.0e10 / h
        u_l = 1.0e3 / h
        u_S = u_m / u_l**2

        S = (0.5 * m_s) / A_d * u_S       # Msun / pc^2
        f_d = m_s_d / m_s
        ssfr = sfr / m_s * .1             # h/Gyr

        f_d[invalid] = -1.
        S[invalid] = -1.
        m_s[invalid] = -1.
        ssfr[invalid] = -1.

        subhs |= {
            'f_disk': f_d,
            'Sigma_star': S,
            'm_star': m_s,
            'ssfr': ssfr,
            'is_valid': ~invalid,
        }
        return self

    def require_rsd(self, at_axis=2) -> Self:
        self['v']
        # TBD


class ElucidLGalaxiesGuo13Sample(SubhaloSet):
    @classmethod
    def load(cls,
             sim_info=sims.predefined['elucid'],
             fin='z0_m1e8_recon.hdf5',
             **init_kw) -> Self:

        fin = ProjPaths.sim_dir_of(sim_info) / 'LGalaxiesGuo13' / fin
        out = cls.from_h5_file(fin, **init_kw)
        assert out.sim_info is sim_info

        subhs = out.objs
        x = subhs['x']
        l_box = sim_info.full_box_size
        x[x < 0.] += l_box
        x[x >= l_box] -= l_box

        return out


class LGalaxiesGuo13Sample(SubhaloSet, abc.HasName):

    @classmethod
    def load(cls,
             sim_info=sims.predefined['millennium_2_scaled_wmap7'],
             in_path='MS2_z0_all.hdf5', **init_kw) -> Self:

        in_dir = ProjPaths.models_dir / 'lgalaxies_guo13' / sim_info.name
        in_file = in_dir / in_path
        out = cls.from_h5_file(in_file, **init_kw)
        assert out.sim_info is sim_info

        subhs = out.objs
        x = subhs['x']
        l_box = sim_info.full_box_size
        x[x < 0.] += l_box
        x[x >= l_box] -= l_box

        return out


class LGalaxiesGuo13SampleRaw(LGalaxiesGuo13Sample):

    @classmethod
    def load(cls, sim_info=sims.predefined['millennium_2_scaled_wmap7'],
             in_path='MS2_z0_all.hdf5', **init_kw) -> Self:
        out = super().load(sim_info=sim_info, in_path=in_path, **init_kw)

        subhs = out.objs

        m_s, m_s_b, r_d = subhs['stellarMass', 'bulgeMass', 'stellarDiskRadius']
        m_s_d = m_s - m_s_b
        m_h = subhs['mvir']

        h = out.sim_info.cosmology.hubble
        u_m = 1.0e10 / h
        u_l = 1.0e6 / h
        u_S = u_m / u_l**2

        A_d = np.pi * r_d**2
        S = (0.5 * m_s) / A_d * u_S
        S_d = (0.5 * m_s_d) / A_d * u_S
        f_d = m_s_d / m_s

        sfr = subhs['sfr']
        ssfr = sfr / m_s * .1      # h/Gyr

        subhs |= {
            'm_s': m_s, 'm_s_b': m_s_b, 'm_s_d': m_s_d, 'f_d': f_d,
            'S': S, 'S_d': S_d, 'm_h': m_h, 'ssfr': ssfr,
        }

        for k in 'stellarMass', 'bulgeMass', 'stellarDiskRadius', 'mvir':
            del subhs[k]

        return out


class LGalaxiesGuo13SampleUdg(LGalaxiesGuo13Sample):

    @classmethod
    def load(cls, sim_info=sims.predefined['millennium_2_scaled_wmap7'],
             in_path='z0_UDG_sample.hdf5', **init_kw) -> Self:
        out = super().load(sim_info=sim_info, in_path=in_path, **init_kw)

        subhs = out.objs
        sfr, m_s = subhs['sfr', 'stellarmass']
        m_s = 10.0**(m_s - 10.0)   # 10^10 Msun/h
        ssfr = sfr / m_s * .1      # h/Gyr

        subhs |= {
            'm_s': m_s, 'ssfr': ssfr,
        }

        for k in ['stellarmass', ]:
            del subhs[k]

        return out

    def find_cen_flag(self, ref_samp: LGalaxiesGuo13Sample):
        x, is_c = ref_samp['x', 'is_c']
        self_x = self['x']

        d, id = KDTree(x).query(self_x)
        is_c = is_c[id]

        self.objs |= {'is_c': is_c, 'd_to_ref': d}
