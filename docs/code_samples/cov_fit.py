import matplotlib.pyplot as plt
from scipy.optimize import minimize
import numpy as np
import emcee
import sys

#---------------------------------------

def log_likelihood(b, mean_corr, inv_cov):
    residual = mean_corr - b
    return -0.5 * residual.T @ inv_cov @ residual

def log_prior(b):
    if 0.0 < b < 10.0:
        return 0.0  # 平坦先验
    return -np.inf

def log_posterior(b, mean_corr, inv_cov):
    lp = log_prior(b)
    if not np.isfinite(lp):
        return -np.inf
    return lp + log_likelihood(b, mean_corr, inv_cov)

def cal_cov(boots):
    #6, 100

    num0    =len(boots)
    num1    =len(boots[0])

    cov     =np.zeros((num0, num0))

    for i in range(num0):
        for j in range(num0):
            m0  =np.mean(boots[i])
            m1  =np.mean(boots[j])

            tem =0
            for k in range(num1):
                tem +=(boots[i][k] - m0) * (boots[j][k] - m1)

            cov[i][j]   =tem / (num1 - 1)

    return cov

def fit(mean_two_point_corr, two_point_corr, tlab):
    #two_point_corr, row:bootstrap; colume:rp

    cov_matrix = np.cov(two_point_corr.T)
    inv_cov_matrix = np.linalg.inv(cov_matrix)

    """
    tcov        =cal_cov(two_point_corr.T)
    #inv_tcov    =cal_inv(tcov)

    for i in range(len(tcov)):
        for j in range(len(tcov[i])):
            print(round(tcov[i][j]  - cov_matrix[i][j], 3))
    sys.exit(0)
    """

    ndim        =1  # 参数维度
    nwalkers    =600  # MCMC 行走者数量
    nsteps      =1500  # MCMC 采样步数

    initial_b   =1.0  # 初始拟合值猜测
    initial_positions = np.random.normal(initial_b, 0.1, (nwalkers, ndim))  # 初始分布

    sampler = emcee.EnsembleSampler(
        nwalkers,
        ndim,
        log_posterior,
        args=(mean_two_point_corr, inv_cov_matrix),
    )

    sampler.run_mcmc(initial_positions, nsteps, progress=False)

    burnin      =1400
    b_samples   =sampler.chain[:,burnin:,:].reshape((-1,ndim))

    #flat_samples = sampler.get_chain(discard=100, thin=15, flat=True)
    #b_samples   = flat_samples[:, 0]

    import seaborn as sns

    plt.figure(figsize=(6, 5))
    sns.heatmap(cov_matrix, annot=True, fmt='.2f', cmap='coolwarm', cbar=True, square=True)
    plt.title("Covariance Matrix, "+tlab)

    save= '/home/zwzhang/data/Dats/box/fig_cov_'+tlab+'.png'
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)

    return b_samples
