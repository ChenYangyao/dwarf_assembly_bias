import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)

from multiprocessing import Pool
import numpy as np
import myfuncs
import cov_fit
import random
import emcee
import math
import sys
import os

#-----------------------------------------------

def lnprior(b):
    if 0.0001 < b < 100:
        return 0.0

    return -np.inf

def lnlike(b, ys,er):
    ys  =np.array(ys)
    er  =np.array(er)

    bs  =np.zeros(len(ys))+b
 
    invers= 1.0/er/er
    diff  = -0.5*((ys - bs)**2*invers - np.log(invers))

    return diff.sum()   

def lnprob(b, ys,er):
    lp = lnprior(b)
    if not np.isfinite(lp):
        return -np.inf

    return lp+lnlike(b, ys,er)

def read_wrp(rou, name):
    with open(r''+rou+name+'.wrp', 'r') as data:
        fp=data.readlines()

    xs,ys,er,lis    =[],[],[],[]
    for row in fp:
        i=row.split()
        i=list(map(lambda j:float(j), i))

        if i[0] < 1.1:
            continue

        xs.append(i[0])
        ys.append(i[1]*2)
        er.append(i[2]*2 + 0.000001)

        tem =list(np.array(i[3:])*2)
        tem1=random.sample(tem, 100)

        #print(len(tem), len(i[3:]), len(tem1))

        lis.append(np.array(tem1))

    return np.array(xs),np.array(ys),np.array(er),lis

def cal_ratio0(rou1,name1, rou2,name2):
    xs,ys1,er1,lis1 =read_wrp(rou1, name1)
    xs,ys2,er2,lis2 =read_wrp(rou2, name2)

    if Mode == 0:
        ys          =ys1/ys2

    else:
        ys          =[]

    er          =[]
    boots       =[]

    for i in range(len(lis1)):
        tem     =lis1[i]/lis2[i]
        boots.append(tem)

        if Mode == 0:
            er.append(np.std(tem))

        elif Mode == 1:
            ys.append(np.mean(tem))
            er.append([np.percentile(tem, 50) - np.percentile(tem, 16), np.percentile(tem, 84) - np.percentile(tem, 50)])

    myfuncs.output3(xs,ys,er, rou1, name1+name2+'.2pccf_ratio')

    return np.array(xs),np.array(ys),np.array(er),np.array(ys1),np.array(er1),np.array(ys2),np.array(er2), np.array(boots)

def cal_ratio1(rou1,name1, rou2,name2):
    xs,ys1,er1  =myfuncs.read3(rou1, name1)
    xs,ys2,er2  =myfuncs.read3(rou2, name2)

    ys          =np.array(ys1)/np.array(ys2)
    er          =((np.array(er1)/np.array(ys1))**2 + (np.array(er2)/np.array(ys2))**2)**0.5

    return np.array(xs),ys,er,ys1,er1,ys2,er2,[]

def main():
    xs,ys,er,ys1,er1,ys2,er2,boots  =cal_ratio0(rou1,name1, rou2,name2)

    ids     =np.where((2 < xs) & (xs < 10))[0]
    nxs     =xs[ids]
    nys     =ys[ids]
    ner     =er[ids]

    if Mode == 0:
        pars            =np.array([1.0])
        ndim,nwalkers   =1,600
        pos             =[pars+1e-4*np.random.randn(ndim) for i in range(nwalkers)]

        with Pool(5) as pool:
            sampler = emcee.EnsembleSampler(nwalkers,ndim,lnprob,args=(nys,ner), pool=pool)
            sampler.run_mcmc(pos,1500, progress=False)

        burnin  =1400
        samples =sampler.chain[:,burnin:,:].reshape((-1,ndim))  

    elif Mode == 1:
        nboots  =boots[ids]

        nboots1 =[]
        for i in range(len(nboots[0])):
            tem =[]
            for j in range(len(nboots)):
                tem.append(nboots[j][i])

            nboots1.append(np.array(tem))

        samples =cov_fit.fit(np.array(nys), np.array(nboots1), name1+name2)

    #-----------------------------------------------

    dat0    =open(rou1+'Samples_'+name1+name2, 'w+')
    for i in samples:
        print(*i, file=dat0)
    dat0.close()

    perc    =[50, 16, 84]
    tb      =[]
    for i in perc:
        tb.append(np.percentile(samples, i))

    dat     =open(rou1+'fit_b_'+name1+name2, 'w+')
    print(tb[0], tb[0]-tb[1], tb[2]-tb[0], file=dat)
    dat.close()

    #print(name1,name2, tb)

    #-----------------------------------------------

    if Draw == 1:
        fig,axes=plt.subplots(2, 1, sharex=True,gridspec_kw = {'hspace':0,  'height_ratios': [2, 1]})
        fig.set_size_inches(10, 10)

        ax1,ax2 =axes[0],axes[1]

        ax1.errorbar(xs,ys1, yerr=er1, fmt='o',ms=20, color='blue',ecolor='blue',elinewidth=7,capsize=20, zorder=1, linestyle='-',linewidth=4, label=name1)
        ax1.errorbar(xs,ys2, yerr=er2, fmt='o',ms=20, color='red', ecolor='red', elinewidth=7,capsize=20, zorder=2, linestyle='-',linewidth=4, label=name2)

        if Mode == 1:
            ter1,ter2   =[],[]
            for n in range(len(er)):
                ter1.append(er[n][0])
                ter2.append(er[n][1])
        else:
            ter1,ter2   =er,er

        ax2.errorbar(xs,ys, yerr=(ter1,ter2), fmt='o',ms=20, color='grey', ecolor='grey', elinewidth=7,capsize=20, zorder=2, linestyle='-',linewidth=4, label=name1+'/'+name2)

        ax2.axhline(y=tb[0],zorder=100,linewidth=3,color='grey',linestyle='--', label='fit relative bias, $'+str(round(tb[0], 2))+'^{+'+str(round(tb[2]-tb[0], 2))+'}_{-'+str(round(tb[0]-tb[1], 2))+'}$')
        ax2.fill_between(nxs, [tb[1]]*len(nxs), [tb[2]]*len(nxs), zorder=100, color='grey', alpha=0.4)

        ax1.legend(fontsize=21,loc='upper right', frameon=False)
        ax2.legend(fontsize=18,loc='upper left', frameon=False)
        ax2.axhline(y=1,linewidth=2,color='black')

        ax1.set_xlim(1, 20)
        ax2.set_xlim(1, 20)

        ax1.set_xscale("log")
        ax2.set_xscale("log")
        ax1.set_yscale("log")

        ax1.set_ylabel(r'$w_{\rm p}(r_{\rm p})[h^{-1}\rm Mpc]$',fontsize=30,labelpad=2)
        ax2.set_ylabel(r'$\rm Ratio$',fontsize=30,labelpad=1.5)
        ax2.set_xlabel(r'$r_{\rm p} [h^{-1}\rm Mpc]$', fontsize=30)

        save= '/home/zwzhang/data/Dats/box/fig_2pccf_'+name1+name2+'_'+str(Mode)+'.png'
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

#------------------------------------------

Mode    =1
Draw    =1
rou1    =sys.argv[1]
name1   =sys.argv[2]
rou2    =sys.argv[3]
name2   =sys.argv[4]

main()
