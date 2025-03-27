from __future__ import annotations
import typing
from typing import Self
from ..sample import SimSample
from . import profiles, models
from pyhipp.core import DataDict
import numpy as np
import numba
from astropy import units as U
from pyhipp.core import abc


class SIDMDriverJiang22(abc.HasDictRepr):
    repr_attr_keys = ('lg_M_b', 'r0', 't_age', 'lg_M_h', 'c',
                      'sigma_per_mX', 'Delta_h', 'redshift',
                      'profile', 'fR')

    def __init__(
        self,
        lg_M_b,                 # [Msun]
        r0,                     # [kpc]
        lg_M_h,                 # [Msun]
        t_age,                  # [Gyr]
        c,                      # concentration
        sigma_per_mX=0.5,       # [cm^2/g]
        Delta_h=200.0,          # halo overdensity in rho_crit
        redshift=0.,
        N_nodes=256,
        output_prof=True,
        enable_ada=True,
    ):
        self.lg_M_b = lg_M_b
        self.r0 = r0
        self.t_age = t_age
        self.lg_M_h = lg_M_h
        self.c = c

        self.sigma_per_mX = sigma_per_mX
        self.Delta_h = Delta_h
        self.redshift = redshift
        self.N_nodes = N_nodes
        self.output_prof = output_prof
        self.enable_ada = enable_ada

        self.__run()

    def __run(self):
        from .jiang22 import profiles as pr, galhalo as gh

        halo = pr.NFW(10.0**self.lg_M_h, c=self.c,
                      Delta=self.Delta_h, z=self.redshift)
        gal = pr.Hernquist(10.0**self.lg_M_b, a=self.r0)

        R_h = halo.rh
        V_h = halo.Vcirc(R_h)
        r_full = self.__rs(1.0e-3, 10.0 * R_h)
        if self.enable_ada:
            halo_contracted = gh.contra(r_full, halo, gal)[0]
        else:
            halo_contracted = halo

        r1 = pr.r1(halo_contracted, sigmamx=self.sigma_per_mX, tage=self.t_age)
        self.r1 = r1
        if not self.output_prof:
            return

        rho_h0, sigma_0, rho, V_c, r = pr.stitchSIDMcore(
            r1, halo_contracted, gal, N=self.N_nodes)

        r_out = self.__rs(1.01*r1, 1.0 * R_h)
        V_c_out = halo_contracted.Vcirc(r_out)
        rho_out = halo_contracted.rho(r_out)

        r = np.concatenate([r, r_out])
        V_c_dm = np.concatenate([V_c, V_c_out])
        rho_dm = np.concatenate([rho, rho_out])
        V_c_b = gal.Vcirc(r)
        rho_b = gal.rho(r)
        V_c = np.sqrt(V_c_dm**2 + V_c_b**2)
        rho = rho_dm + rho_b

        self.halo = halo
        self.halo_contracted = halo_contracted
        self.gal = gal
        self.R_h = R_h
        self.V_h = V_h
        self.profile = DataDict({
            'r': r, 'V_c': V_c, 'V_c_b': V_c_b, 'V_c_dm': V_c_dm,
            'rho': rho, 'rho_b': rho_b, 'rho_dm': rho_dm,
            'rho_h0': rho_h0, 'sigma_0': sigma_0,
        })
        self.J_integ = self.__find_fR()

    def __find_fR(self):
        out = DataDict()

        rs, V_cs = self.profile['r', 'V_c']
        us = rs / self.r0
        vs = V_cs / self.V_h

        du = us[1:] - us[:-1]
        us = 0.5 * (us[1:] + us[:-1])
        vs = 0.5 * (vs[1:] + vs[:-1])

        integ = us**2 / (1. + us)**3 * vs * du
        fR = 1.0 / integ.sum()
        out['H'] = {
            'fR': fR, 'integ': integ,
        }

        integ = np.exp(-us) * us**2 * vs * du * 0.5
        fR = 1.0 / integ.sum()
        out['exp'] = {
            'fR': fR, 'integ': integ,
        }
        return out

    def __rs(self, r_lb, r_ub, N=None):
        if N is None:
            N = self.N_nodes
        lr_lb, lr_ub = np.log10(r_lb), np.log10(r_ub)
        return np.logspace(lr_lb, lr_ub, N)


class SIDMDriver:
    def __init__(self, samp: SimSample,
                 M_s_in_sol=10**8.8,
                 S_key='Sigma_star_AM_lm10p5-11_rho0.85',
                 R_eff_key=None,
                 ada_model='raw',
                 sidm_with_b=True,
                 sigma_in_cm2perg=0.5,
                 keep_halo=False,
                 ):

        sim_info = samp.sim_info
        us = sim_info.cosmology.unit_system
        u_m = us.u_m_to_sol
        u_l = us.u_l_to_pc
        u_sig = us.u_l**2 / us.u_m
        sigma = (sigma_in_cm2perg * U.cm**2 / U.g).to(u_sig).value

        if S_key is not None:
            assert R_eff_key is None
            S = samp[S_key]                                     # Msun / pc^2
            R_eff = np.sqrt(M_s_in_sol / S / (2.0 * np.pi))     # pc
            R_eff = R_eff / u_l
        else:
            R_eff = samp[R_eff_key]
            R_eff_in_pc = R_eff * u_l
            S = M_s_in_sol / (2.0 * np.pi * R_eff_in_pc**2)     # Msun / pc^2
        R_scl = R_eff / 1.678
        M_s = M_s_in_sol / u_m

        M_h = samp['m_mean200']
        R_h = samp['r_mean200']
        c_h = samp['c_mean200']
        z_f = samp['z_half']
        t_age = sim_info.z_to_lookback_time(z_f)

        self.data = DataDict({
            'S': S, 'R_eff': R_eff,
            'M_s': M_s, 'R_scl': R_scl,
            'M_h': M_h, 'R_h': R_h, 'c_h': c_h, 't_age': t_age,
        })
        self.m_ada = profiles.AdiabaticModel(ada_model,
                                             M_soft=0.01*M_s)
        self.sidm_with_b = sidm_with_b
        self.sigma = sigma
        self.u_l = u_l
        self.u_m = u_m
        self.u_sig = u_sig
        self.keep_halo = keep_halo

    def require_r1(self):
        r1s = []
        for h in self.__iter_halos():
            dm_st = self.__get_one(*h)
            r1s.append(10.0**dm_st.lr1)
        self.data |= {'r1': np.array(r1s)}

    def require_rho0(self):
        outs, sidm_halos = [], []
        for h in self.__iter_halos():
            dm_st = self.__get_one(*h)
            r1 = 10.0 ** dm_st.lr1
            rho0 = dm_st.iso_core._rho_0
            rho1 = dm_st.rho_r1

            pf = dm_st.profile._impl
            lrs, rhos = pf._lr_cents, pf._rho_cents
            rhalf = self.__find_r(lrs, rhos, rho0 * 0.5)
            rh4 = self.__find_r(lrs, rhos, rho0 * 0.25)

            outs.append((rhalf, rh4, r1, rho0, rho1))
            if self.keep_halo:
                sidm_halos.append(dm_st)

        rhalf, rh4, r1, rho0, rho1 = np.array(outs).T
        self.data |= {'rhalf': rhalf, 'rh4': rh4, 'r1': r1,
                      'rho0': rho0, 'rho1': rho1, 'sidm_halos': sidm_halos}

    def __iter_halos(self):
        M_s, R_scl, M_h, R_h, c_h, t_age = self.data[
            'M_s', 'R_scl', 'M_h', 'R_h', 'c_h', 't_age']
        for i in range(len(M_h)):
            yield M_s, R_scl[i], M_h[i], R_h[i], c_h[i], t_age[i]

    def __get_one(self, M_s, R_scl, M_h, R_h, c_h, t_age):
        b = profiles.ExpSph(M_s, R_scl, f_in=1.0e-3)
        dm = profiles.NFW(M_h - M_s, R_h, c_h)
        dm_ada = self.m_ada(dm, b)
        if not self.sidm_with_b:
            b = profiles.ExpSph(M_s*1.0e-3, R_scl)
        dm_st = models.StitchedHalo(dm_ada, b, t_age, self.sigma)
        return dm_st

    def __find_r(self, lrs, rhos, rho_dst):
        args = (rhos >= rho_dst).nonzero()[0]
        if args.size == 0:
            r_dst = 0.
        else:
            arg = args[-1]
            r_dst = 10.0**(lrs[arg])
        return r_dst
