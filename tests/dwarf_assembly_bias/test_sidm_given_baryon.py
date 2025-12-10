from pyhipp import plot
from dwarf_assembly_bias.sidm import profiles, models
from astropy import units as U
from functools import cached_property
import numpy as np

cosm = profiles.default_cosm                # Planck 15
ctx = profiles.Context.from_cosm(cosm)


class Halo:
    def __init__(self,
                 cdm: profiles.Halo,
                 age: float,
                 bary: profiles.Baryon):

        self.cdm = cdm
        self.age = age
        self.bary = bary

        self.r1: float = None
        self.rho0: float = None
        self.rho1: float = None
        self.sidm: profiles.SphSymmHalo = None


class SIDMDriverFixedBaryonProfile:
    def __init__(self,
                 ada_model='raw',
                 ada_M_soft=1.0e-4,
                 sigma_in_cm2perg=0.5,
                 ):

        us = cosm.unit_system
        u_sig = us.u_l**2 / us.u_m
        sigma = (sigma_in_cm2perg * U.cm**2 / U.g).to(u_sig).value

        self.u_sig = u_sig
        self.m_ada = profiles.AdiabaticModel(ada_model, M_soft=ada_M_soft)
        self.sigma = sigma

    def on_halo(self, halo: Halo):
        '''
        @halo [in|out].
        '''
        dm = halo.cdm
        b = halo.bary

        dm_ada = self.m_ada(dm, b)
        dm_st = models.StitchedHalo(dm_ada, b, halo.age, self.sigma)

        halo.r1 = 10.0**dm_st.lr1
        halo.rho0 = dm_st.iso_core._rho_0
        halo.rho1 = dm_st.rho_r1
        halo.sidm = profiles.SphSymmHalo(dm.R_v, dm_st.profile)


def run_sidm(
    m_h_in_sol,
    c,
    age_in_yr,
    z,
    r_nodes_in_pc,
    dm_bary_bins_in_sol,
    r_outs_in_pc,
    sigma_in_cm2perg,
    halo_def=('mean', {'f': 200.0}),
):
    us = cosm.unit_system

    u_l = us.u_l_to_pc
    u_m = us.u_m_to_sol
    u_t = us.u_t_to_yr
    u_rho = u_m / u_l**3
    u_v2kmps = us.u_v_to_kmps

    lr_nodes = np.log10(r_nodes_in_pc / u_l)
    dm_baryon_bins = dm_bary_bins_in_sol / u_m
    bary = profiles.SphSymmBaryon.from_bins(lr_nodes, dm_baryon_bins)
    m_bary = bary.M

    m_h = m_h_in_sol / u_m
    # find halo radius
    def_name, def_kw = halo_def
    ht = cosm.halo_theory
    if def_name == 'mean':
        r_h = ht.vir_props_mean(m_h, **def_kw, z=z).r_phy
    elif def_name == 'crit':
        r_h = ht.vir_props_crit(m_h, **def_kw, z=z).r_phy
    else:
        raise ValueError(f'Unknown halo_def: {halo_def}')
    age = age_in_yr / u_t
    cdm = profiles.NFW(m_h - m_bary, r_h, c=c)

    halo = Halo(cdm, age, bary)
    sidm_driver = SIDMDriverFixedBaryonProfile(
        sigma_in_cm2perg=sigma_in_cm2perg)
    sidm_driver.on_halo(halo)

    sidm = halo.sidm
    lr_outs = np.log10(r_outs_in_pc / u_l)
    Vcs_sidm = np.sqrt(sidm.rotation.Vc_sqs(lr_outs)) * u_v2kmps
    Vcs_cdm = np.sqrt(cdm.rotation.Vc_sqs(lr_outs)) * u_v2kmps
    Vcs_bary = np.sqrt(bary.rotation.Vc_sqs(lr_outs)) * u_v2kmps

    rhos_sidm = sidm.profile.rhos(lr_outs) * u_rho
    rhos_cdm = cdm.profile.rhos(lr_outs) * u_rho
    rhos_bary = bary.profile.rhos(lr_outs) * u_rho

    r1_in_pc = halo.r1 * u_l
    rho0_in_solperpc3 = halo.rho0 * u_rho
    rho1_in_solperpc3 = halo.rho1 * u_rho

    return {
        'r_outs_in_pc': r_outs_in_pc,
        'r1_in_pc': r1_in_pc,                         # one-scattering radius
        'rho0_in_solperpc3': rho0_in_solperpc3,       # central density
        'rho1_in_solperpc3': rho1_in_solperpc3,       # density at r1
        'Vcs_in_kmps': {
            'sidm': Vcs_sidm,
            'cdm': Vcs_cdm,
            'bary': Vcs_bary,
        },
        'rho_in_solperpc3': {
            'sidm': rhos_sidm,
            'cdm': rhos_cdm,
            'bary': rhos_bary,
        },
    }


def test_run():

    m_h_in_sol = 1.0e12                      # halo mass (DM + baryon) [Msun]
    c = 10.0                                 # halo concentration
    age_in_yr = 5.0e9                        # halo age [yr]
    z = 0.                                   # redshift of the halo
    sigma_in_cm2perg = 0.5                   # SIDM cross section [cm^2/g]

    # the input profile of baryon, specified by the edges (i.e. nodes) [pc]
    # of continous bins and the masses [Msun] in the bins.
    # For non-spherical baryon distribution, one should give the masses in
    # spherical shells.
    # Here I take an exponential spheroid as an example.
    r_nodes_in_pc = np.array([
        1.000e+00, 1.258e+00, 1.584e+00, 1.995e+00,
        2.511e+00, 3.162e+00, 3.981e+00, 5.011e+00,
        6.309e+00, 7.943e+00, 1.000e+01, 1.258e+01,
        1.584e+01, 1.995e+01, 2.511e+01, 3.162e+01,
        3.981e+01, 5.011e+01, 6.309e+01, 7.943e+01,
        1.000e+02, 1.258e+02, 1.584e+02, 1.995e+02,
        2.511e+02, 3.162e+02, 3.981e+02, 5.011e+02,
        6.309e+02, 7.943e+02, 1.000e+03, 1.258e+03,
        1.584e+03, 1.995e+03, 2.511e+03, 3.162e+03,
        3.981e+03, 5.011e+03, 6.309e+03, 7.943e+03,
        1.000e+04, 1.258e+04, 1.584e+04, 1.995e+04,
        2.511e+04, 3.162e+04, 3.981e+04, 5.011e+04,
        6.309e+04, 7.943e+04, 1.000e+05])
    dm_bary_bins_in_sol = np.array([
        8.322e+00, 1.681e+01, 3.042e+01, 4.750e+01,
        7.677e+01, 1.201e+02, 1.914e+02, 3.034e+02,
        4.768e+02, 7.659e+02, 1.197e+03, 1.919e+03,
        3.024e+03, 4.782e+03, 7.633e+03, 1.191e+04,
        1.921e+04, 3.009e+04, 4.784e+04, 7.586e+04,
        1.188e+05, 1.908e+05, 2.977e+05, 4.751e+05,
        7.474e+05, 1.174e+06, 1.867e+06, 2.896e+06,
        4.617e+06, 7.166e+06, 1.120e+07, 1.749e+07,
        2.675e+07, 4.183e+07, 6.308e+07, 9.601e+07,
        1.431e+08, 2.093e+08, 3.045e+08, 4.243e+08,
        5.842e+08, 7.636e+08, 9.553e+08, 1.126e+09,
        1.227e+09, 1.215e+09, 1.072e+09, 8.097e+08,
        5.085e+08, 2.565e+08])

    # at which radii to output the density and Vc profiles [pc]
    r_outs_in_pc = np.logspace(0., 6., 61)

    # run SIDM model
    out = run_sidm(m_h_in_sol, c, age_in_yr, z,
                r_nodes_in_pc, dm_bary_bins_in_sol,
                r_outs_in_pc, sigma_in_cm2perg)


    # Make a plot
    fig, axs = plot.subplots((1,2), share=(True, False), 
                            space=0.3, 
                            subsize=(5.0, 5.0), 
                            margin=[0.1, 0.1, 0.15, 0.1], layout='none')
    axs_f = axs.flat

    ax = axs_f[0]

    rs = out['r_outs_in_pc']
    Vcs = out['Vcs_in_kmps']
    rhos = out['rho_in_solperpc3']

    ax.plot(rs, Vcs['sidm'], c='purple', label=r'$\rm SIDM$')
    ax.plot(rs, Vcs['cdm'], c='purple', lw=1, ls='--', label=r'$\rm CDM$')
    ax.plot(rs, Vcs['bary'], c='green', label=r'$\rm Baryon$')

    ax.scale('log', 'log').lim([1.0, 1.0e6], [1.0e-2, 1.0e3])\
        .label(r'r\,[{\rm pc}]', r'V_c\,[{\rm km\,s^{-1}}]')\
        .leg()

    ax = axs_f[1]

    ax.plot(rs, rhos['sidm'], c='purple')
    ax.plot(rs, rhos['cdm'], c='purple', lw=1, ls='--')
    ax.plot(rs, rhos['bary'], c='green')

    rho0, rho1 = out['rho0_in_solperpc3'], out['rho1_in_solperpc3']
    r1 = out['r1_in_pc']
    ax.plot([1.0e-8, 1.0e6], [rho0, rho0], c='black', lw=1, ls='-', label=r'$\rho_0$')
    ax.plot([1.0e-8, 1.0e6], [rho1, rho1], c='red', lw=1, ls='-', label=r'$\rho_1$')
    ax.plot([r1, r1], [1.0e-8, 1.0e2], c='blue', lw=1, ls='-', label=r'$r_1$')

    ax.scale('log', 'log').lim(y=[1.0e-8, 1.0e2])\
        .label(r'r\,[{\rm pc}]', r'\rho\,[{\rm M_\odot\,pc^{-3}}]')\
        .leg(loc='ur')

    plot.savefig('out.pdf')
