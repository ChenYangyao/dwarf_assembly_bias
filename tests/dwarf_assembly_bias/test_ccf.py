import pytest
from pyhipp.stats import Rng
from pyhipp.core import DataDict, DataTable
from pyhipp.astro.stats import ccf
import numpy as np

rng = Rng(10086)
l_boxs = 500.


@pytest.fixture
def samp1():
    n = 1000
    x = rng.uniform(0, l_boxs, size=(n, 3))
    return ccf.SimSample(l_boxs, DataTable({'x': x}))

@pytest.fixture
def samp2():
    n = 1000
    x = rng.uniform(0, l_boxs, size=(n, 3))
    q = rng.uniform(0, 1, size=n)
    return ccf.SimSample(l_boxs, DataTable({'x': x, 'q': q}))


@pytest.fixture
def samp_ref():
    n = 2000
    x = rng.uniform(0, l_boxs, size=(n, 3))
    return ccf.SimSample(l_boxs, DataTable({'x': x}))


@pytest.fixture
def sim_ccf_projected(samp_ref):
    return ccf.SimCCFProjected(samp_ref, rng=rng, n_threads=1,
                               pi_max=20., n_bootstrap=10)


def test_wp(sim_ccf_projected: ccf.SimCCFProjected, samp1: ccf.SimSample):
    return sim_ccf_projected.wp(samp1, rs=np.array([1., 2., 5., 15.]))


def test_relative_bias_curve(sim_ccf_projected: ccf.SimCCFProjected,
                             samp2: ccf.SimSample):

    return sim_ccf_projected.relative_bias_curve(
        samp2, bin_by_key='q', bin_edges=[-0.1, 0.25, 0.5, 0.75, 1.1],
        ref_bin=-1, r_min=2., r_max=10.)
