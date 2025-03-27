import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import matplotlib.font_manager as fm

from pdf2image import convert_from_path
from PIL import Image
import pandas as pd
import numpy as np
import warnings
import time
import os

#--------------------------------

f_size  =16 #figure size
t_size  =35 #tick size
l_size  =40 #label size
le_size =40 #legend size
s_size  =60 #scatter size
e_size  =30 #errorbar size
c_size  =20 #capsize
l_pad   =20 #label pad
l_wid   =5  #linewidth
el_wid  =7  #errorbar linewidth
b_wid   =4  #the width of the frame

#--------------------------------

def fig1():
    fig,axes= plt.subplots(1, 3)
    fig.set_size_inches(f_size*3+5, f_size)

    alabs   =['a', 'b', 'c']
    for i in range(len(alabs)):
        ax  =axes[i]
        ax.text(0.07, 0.93, alabs[i], transform=ax.transAxes,fontsize=60, ha='center', va='center', color='black', weight='bold')

    names   =['samp_0', 'samp_3', 'halo_1', 'halo_2', 'halo_3']
    labs    =['Diffuse dwarfs', 'Compact dwarfs', r'Groups, 11.0$\leq$log $M_{\rm h}$ (M$_{\odot}$)<12.0', r'Groups, 12.0$\leq$log $M_{\rm h}$ (M$_{\odot}$)<13.0', r'Groups, 13.0$\leq$log $M_{\rm h}$ (M$_{\odot}$)']
    colos   =['blue', 'red', 'dimgray', 'darkred', 'turquoise']
    lines   =['-', '-', '--', '--', '--']

    for i in range(len(names)):
        if i==0:
            tz  =-100
        else:
            tz  =-i

        xs,ys,er1,er2   =funcs().read_data('Fig1.xlsx', 'a_'+names[i], ['x', 'y', 'y_er_low', 'y_er_up'])

        if 'samp' in names[i]:
            axes[0].errorbar(xs,ys,yerr=(er1,er2), zorder=tz,fmt='o',ms=e_size,color=colos[i],ecolor=colos[i], elinewidth=el_wid, linestyle=lines[i], linewidth=l_wid, capsize=c_size, label=labs[i])

        elif 'halo' in names[i]:
            axes[0].errorbar(xs,ys,yerr=(er1,er2), zorder=tz,fmt='D',ms=e_size,mew=3, mfc='none', mec=colos[i],color=colos[i],ecolor=colos[i], elinewidth=el_wid, linestyle=lines[i], linewidth=l_wid, capsize=c_size, label=labs[i])

    axes[0].fill_betweenx(np.linspace(3, 100, 20), [2]*20, [10]*20, color='grey', alpha=0.2, zorder=-1000)
    axes[0].text(1.21, 105, "The fitting range for\n       relative bias", fontsize=35, color="black")
    style().sub_draw_2pccf(axes[0], 0)

    #--------------------------------

    colos           =['blue', 'green', 'orange', 'red']
    xs,ys,er1,er2   =funcs().read_data('Fig1.xlsx', 'b_main', ['x', 'y', 'y_er_low', 'y_er_up'])

    for i in range(len(xs)):
        axes[1].errorbar([xs[i]],[ys[i]], yerr=([er1[i]],[er2[i]]), fmt='o',ms=30, color=colos[i],ecolor=colos[i],elinewidth=7,capsize=20, linestyle='-',linewidth=5, zorder=-i)
    axes[1].plot(xs,ys, linewidth=3.5, color='black', linestyle='-', zorder=-99)

    colos           =['lightblue', 'lightgreen', 'navajowhite', 'tomato']
    xs,ys,er1,er2   =funcs().read_data('Fig1.xlsx', 'b_massive', ['x', 'y', 'y_er_low', 'y_er_up'])

    for i in range(len(xs)):
        axes[1].errorbar([xs[i]],ys[i], yerr=([er1[i]],[er2[i]]), mfc='none', mec=colos[i],mew=5, fmt='o',ms=35, color=colos[i],ecolor=colos[i],elinewidth=7,capsize=20, linestyle='--',linewidth=5, zorder=-100+i)
    axes[1].plot(xs,ys, linewidth=3, color='black', linestyle='--', zorder=-100)

    axes[1].errorbar([-1],[-1], yerr=([0], [0]), fmt='o',ms=30, color='black',ecolor='black',elinewidth=7,capsize=20, linestyle='-',linewidth=5, label='Main sample')
    axes[1].errorbar([-1],[-1], yerr=([0], [0]), mfc='none', mec='black',mew=5, fmt='o',ms=35, color='black',ecolor='black',elinewidth=7,capsize=20,  linestyle='--',linewidth=5, label='Massive sample')

    style().sub_draw_bias(axes[1], 'Sig')

    #--------------------------------

    names   =[['samp_0', 'samp_1', 'samp_2', 'samp_3'], ['halo_1', 'halo_2', 'halo_3']]
    colos   =[['blue', 'green', 'orange', 'red'], ['dimgray', 'darkred', 'turquoise']]
    labs    =[['Diffuse dwarfs', '$7\leq\Sigma_*$<15', '$15\leq\Sigma_*$<25', 'Compact dwarfs'], [r'Groups, 11.0$\leq$log $M_{\rm h}$ (M$_{\odot}$)<12.0', r'Groups, 12.0$\leq$log $M_{\rm h}$ (M$_{\odot}$)<13.0', r'Groups, 13.0$\leq$log $M_{\rm h}$ (M$_{\odot}$)']]

    for o in range(len(names)):
        if o == 0:
            xs,xe1,xe2,ys,ye1,ye2   =funcs().read_data('Fig1.xlsx', 'c_main', ['x', 'x_er_low', 'x_er_up', 'y', 'y_er_low', 'y_er_up'])

            for i in range(len(xs)):
                ax.errorbar([xs[i]],[ys[i]], xerr=([xe1[i]], [xe2[i]]), yerr=([ye1[i]], [ye2[i]]), fmt='o', ms=e_size, color=colos[o][i],ecolor=colos[o][i],elinewidth=el_wid,capsize=c_size, zorder=i, label=labs[o][i])

        else:
            xs,ys,ye1,ye2   =funcs().read_data('Fig1.xlsx', 'c_groups', ['x', 'y', 'y_er_low', 'y_er_up'])
        
            for i in range(len(xs)):
                ax.errorbar([xs[i]],[ys[i]], yerr=([ye1[i]], [ye2[i]]), fmt='D', mew=3, ms=e_size, mfc='none', mec=colos[o][i], color=colos[o][i],ecolor=colos[o][i],elinewidth=el_wid,capsize=c_size, zorder=i)

                hms,bs          =funcs().read_data('Fig1.xlsx', 'c_hb_'+names[o][i], ['x', 'y'])
                ax.plot(hms,bs, color=colos[o][i],linewidth=l_wid,linestyle='--', zorder=i)

    style().sub_draw_bias(axes[2], 'hb')

    save= fig_rou+'Fig1.pdf'
    plt.subplots_adjust()
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

def efig1():
    fig,axes    =plt.subplots(2, 2)
    fig.set_size_inches(f_size*2, f_size*2)

    ax0,ax1,ax2,ax3 =axes[0,0],axes[0,1],axes[1,0],axes[1,1]

    ax0.text(0.09, 0.93, 'a', transform=ax0.transAxes,fontsize=50, ha='center', va='center', color='black', weight='bold')
    ax1.text(0.09, 0.93, 'b', transform=ax1.transAxes,fontsize=50, ha='center', va='center', color='black', weight='bold')
    ax2.text(0.09, 0.93, 'c', transform=ax2.transAxes,fontsize=50, ha='center', va='center', color='black', weight='bold')
    ax3.text(0.09, 0.93, 'd', transform=ax3.transAxes,fontsize=50, ha='center', va='center', color='black', weight='bold')

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'a_0', ['x', 'y'])
    ax0.scatter(xs,ys, color='green',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'a_1', ['x', 'y'])
    ax0.scatter(xs,ys, color='blue',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'a_2', ['x', 'y'])
    ax0.scatter(xs,ys, color='red',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'b_0', ['x', 'y'])
    ax1.scatter(xs,ys, color='green',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'b_1', ['x', 'y'])
    ax1.scatter(xs,ys, color='blue',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'b_2', ['x', 'y'])
    ax1.scatter(xs,ys, color='red',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'c_0', ['x', 'y'])
    ax2.scatter(xs,ys, color='green',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'c_1', ['x', 'y'])
    ax2.scatter(xs,ys, color='blue',  s=s_size)

    xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'c_2', ['x', 'y'])
    ax2.scatter(xs,ys, color='red',  s=s_size)

    style().sefig1_abc(ax0,0.6,  3, 250,      -0.05,  0.8, '$\Sigma_*$ (M$_{\odot}$pc$^{-2}$)', '$^{0.1}(g-r)$')
    style().sefig1_abc(ax1,1.6,  3, 250,      0.5,   3,   '$\Sigma_*$ (M$_{\odot}$pc$^{-2}$)', r'S$\acute{e}$rsic index')
    style().sefig1_abc(ax2,2.0,  -0.05,0.8,   0.5,   3,   '$^{0.1}(g-r)$', r'S$\acute{e}$rsic index')

    #--------------------------------

    labs    =['Total dwarfs', '0$\leq \Sigma_*<$7', '7$\leq \Sigma_*<$15', '15$\leq \Sigma_*<$25', '25$\leq \Sigma_*$']
    names   =['samp_4', 'samp_0', 'samp_1', 'samp_2', 'samp_3']
    colos   =['black', 'blue', 'green', 'orange', 'red']

    for i in range(len(names)):
        xs,ys   =funcs().read_data('Wang_EDfig1.xlsx', 'd_'+names[i], ['x', 'y'])
        ax3.plot(xs,ys, zorder=i,color=colos[i],linewidth=l_wid,label=labs[i])

    xs,low  =funcs().read_data('Wang_EDfig1.xlsx', 'd_low', ['x', 'y'])
    ax3.fill_between(xs,np.zeros(len(xs)), low, color='grey',alpha=0.2, label=r'$f_{\rm ctl}$($z$)')

    style().sefig1_d(ax3)

    save= fig_rou+'Wang_EDfig1.pdf'
    plt.subplots_adjust()
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

def efig2():
    labs0   =['samp_0', 'samp_1', 'samp_2']
    labs1   =['$0\leq\Sigma_*<7$', '$7\leq\Sigma_*<15$', '$15\leq\Sigma_*<25$']
    labs2   =['$z$-weighting', '$z$-matching']
    colos0  =['blue', 'green', 'orange']
    alabs   =['a', 'b', 'c', 'd', 'e', 'f']

    for k in range(2):
        fig,axes=plt.subplots(2, 3, sharex=True,gridspec_kw = {'hspace':0,  'height_ratios': [2, 1]})
        fig.set_size_inches(f_size*3, f_size)

        for oo in range(len(labs0)):
            ax1,ax2 =axes[0,oo],axes[1,oo]

            aid =k*3+oo
            ax1.text(0.3, 0.9, r''+alabs[aid]+'   '+labs2[k], transform=ax1.transAxes,fontsize=50, ha='center', va='center', color='black', weight='bold')

            if k == 0:
                names   =[labs0[oo], 'samp_3']
            else:
                names   =['m_'+labs0[oo], 'm_samp_3']

            labs    =[labs1[oo], '$25\leq\Sigma_*$']
            colos   =[colos0[oo],'red']

            for i in range(len(names)):
                xs,ys,er1,er2   =funcs().read_data('Wang_EDfig2.xlsx', alabs[aid]+'_'+names[i], ['x', 'y', 'y_er_low', 'y_er_up'])

                ax1.errorbar(xs,ys,yerr=(er1,er2), ms=e_size, zorder=-i,fmt='o', color=colos[i],ecolor=colos[i], linewidth=l_wid,elinewidth=el_wid, linestyle='-',capsize=c_size, label=labs[i])

            xs0,ys0,er10,er20   =funcs().read_data('Wang_EDfig2.xlsx', alabs[aid]+'_ratio', ['x', 'y', 'y_er_low', 'y_er_up'])

            ax2.errorbar(xs0,ys0,yerr=(er10,er20), ms=e_size, zorder=i,fmt='o', color=colos[0],ecolor=colos[0], linewidth=l_wid,elinewidth=el_wid, linestyle='-',capsize=c_size, label=labs[0]+'\n'+labs[1])

            ax2.plot([0.17, 0.9], [1.48, 1.48], color='black', linewidth=2, zorder=1000)

            b,b1,b2 =funcs().read_data('Wang_EDfig2.xlsx', alabs[aid]+'_fit', ['b', 'b1', 'b2'])

            xs  =xs0[np.where((xs0>2) & (xs0<10))[0]]
            ax2.fill_between(xs, [b[0]-b1[0]]*len(xs), [b[0]+b2[0]]*len(xs), zorder=100, color='grey', alpha=0.4, label='fit relative bias\n$'+str(round(b[0], 2))+'^{+'+str(round(b2[0], 2))+'}_{-'+str(round(b1[0], 2))+'}$')

            style().sefig2(ax1, ax2, oo)

        save= fig_rou+'fig_2pccf_dwarf_'+str(k)+'.pdf'
        plt.subplots_adjust(wspace=0, hspace=0)
        plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
        plt.close(fig)

    image1 = convert_from_path(frou+'Figs/fig_2pccf_dwarf_0.pdf')[0]
    image2 = convert_from_path(frou+'Figs/fig_2pccf_dwarf_1.pdf')[0]

    width1, height1 = image1.size
    width2, height2 = image2.size

    combined_width = max(width1, width2)
    combined_height = height1 + height2

    combined_image = Image.new('RGB', (combined_width, combined_height))
    combined_image.paste(image1, (0, 0))
    combined_image.paste(image2, (0, height1))
    combined_image.save(fig_rou+'Wang_EDfig2.pdf')

    return fig_rou+'Wang_EDfig2.pdf'

def efig3():
    colos   =['blue', 'green', 'orange', 'red']
    labs    =['Diffuse dwarfs', '$7\leq\Sigma_*$<15', '$15\leq\Sigma_*$<25', 'Compact dwarfs']
    lis     =['log $M_*$ (M$_{\odot}$)', '$^{0.1}$(g-r)']
    inters  =[8.5, 0.3]

    fig,axes=plt.subplots(1, 2)
    fig.set_size_inches(f_size*2, f_size*0.87)

    alabs   =['a', 'b']
    for o in range(len(lis)):
        ax      =axes[o]
        tlab    =str(inters[o])

        ax.text(0.07, 0.93, alabs[o], transform=ax.transAxes,fontsize=50, ha='center', va='center', color='black', weight='bold')

        l0s,l1s =[],[]
        for i in range(len(labs)):
            xs,ys   =funcs().read_data('Wang_EDfig3.xlsx', alabs[o]+'_l_'+str(i), ['x', 'y'])
            l0, =ax.plot(xs,ys, color=colos[i],linewidth=l_wid+i*0.5,linestyle='--', label=labs[i]+'\n'+lis[o]+'$\leq$'+tlab)

            xs,ys   =funcs().read_data('Wang_EDfig3.xlsx', alabs[o]+'_u_'+str(i), ['x', 'y'])
            l1, =ax.plot(xs,ys, color=colos[i],linewidth=l_wid,linestyle='-',  label=labs[i]+'\n'+lis[o]+'>'+tlab)

            l0s.append(l0)
            l1s.append(l1)

        style().sefig3(ax, l0s,l1s)

    save= fig_rou+'Wang_EDfig3.pdf'
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

def efig4():
    fig,axes    =plt.subplots(2, 3)
    fig.set_size_inches(f_size*3, f_size*2)

    alabs   =['a', 'b', 'c', 'd', 'e', 'f']
    for i in range(2):
        for j in range(3):
            tid =i*3 + j
            ax  =axes[i, j]
            ax.text(0.07, 0.93, alabs[tid], transform=ax.transAxes,fontsize=60, ha='center', va='center', color='black', weight='bold')

    labs    =[['z$\leq$0.025', 'z>0.025'], ['7.5$\leq$log $M_*$ (M$_{\odot}$)<8.5', '8.5$\leq$log $M_*$ (M$_{\odot}$)<9.0'], ['$^{0.1}$(g-r)$\leq$0.3', '0.3$\leq^{0.1}$(g-r)<0.6'], ['R.A.$\leq$181 degrees', 'R.A.>181 degrees']]
    colos   =['blue', 'green']

    for i in range(len(labs)):
        for j in range(len(labs[i])):
            xs,ys,er1,er2   =funcs().read_data('Wang_EDfig4.xlsx', alabs[i]+'_'+str(j), ['x', 'y', 'y_er_low', 'y_er_up'])

            axes[i//3, i%3].errorbar(xs,ys, yerr=(er1,er2), fmt='o',ms=e_size, color=colos[j],ecolor=colos[j],elinewidth=el_wid,capsize=c_size, linestyle='-',linewidth=l_wid, zorder=j, label=labs[i][j])

            style().sub_draw_bias(axes[i//3, i%3], 'ef4')

    axs     =[axes[1, 1], axes[1, 2]]
    labs    =[['No n-cut sample', 'n>1.6 sample'], ['No n-cut sample']]
    colos   =['blue', 'green']
    labs1   =['nn_cut', 'sn']

    for oo in range(len(labs)):
        ax      =axs[oo]

        for u in range(len(labs[oo])):
            xs,ys,er1,er2   =funcs().read_data('Wang_EDfig4.xlsx', alabs[oo+4]+'_'+str(u), ['x', 'y', 'y_er_low', 'y_er_up'])

            ax.errorbar(xs,ys, yerr=(er1,er2), fmt='o',ms=e_size, color=colos[u],ecolor=colos[u],elinewidth=el_wid,capsize=c_size, linestyle='-',linewidth=l_wid, zorder=10+u, label=labs[oo][u])

        style().sub_draw_bias(ax, labs1[oo])

    colos           =['lightblue']
    xs,ys,er1,er2   =funcs().read_data('Wang_EDfig4.xlsx', alabs[4]+'_'+str(2), ['x', 'y', 'y_er_low', 'y_er_up'])

    axs[0].errorbar(xs,ys, yerr=(er1,er2), fmt='o',ms=e_size, color='lightblue',ecolor='lightblue',elinewidth=el_wid,capsize=c_size, linestyle='-',linewidth=l_wid, label='Main sample')
    axs[0].legend(loc='upper right', fontsize=le_size, frameon=False)

    save= fig_rou+'Wang_EDfig4.pdf'
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

def efig5():
    fig,ax  = plt.subplots()
    fig.set_size_inches(f_size, f_size)

    names   =['samp_0', 'samp_3', 'sat_0.1', 'sat_0.298', 'sat_0.31']
    labs    =['Diffuse dwarfs', 'Compact dwarfs', 'Compact dwarfs\n+10% satellite contamination', 'Compact dwarfs\n+30% satellite contamination', 'Diffuse dwarfs\n+31% satellite contamination']
    colos   =['blue', 'red', 'dimgray', 'darkred', 'turquoise']
    lines   =['-', '-', ':', ':', ':']

    for i in range(len(names)):
        xs,ys,er1,er2   =funcs().read_data('Wang_EDfig5.xlsx', names[i], ['x', 'y', 'y_er_low', 'y_er_up'])

        if 'samp' in names[i]:
            ax.errorbar(xs,ys,yerr=(er1,er2), zorder=-i,fmt='o',ms=e_size,color=colos[i],ecolor=colos[i], elinewidth=el_wid, linestyle=lines[i], linewidth=l_wid, capsize=c_size, label=labs[i])

        elif 'sat' in names[i]:
            ax.errorbar(xs,ys,yerr=(er1,er2), zorder=-i,fmt='*',ms=e_size,mfc='none', mec=colos[i],color=colos[i],ecolor=colos[i], elinewidth=el_wid, linestyle=lines[i], linewidth=l_wid, capsize=c_size, label=labs[i])

    style().sub_draw_2pccf(ax, 1)

    save= fig_rou+'Wang_EDfig5.pdf'
    plt.subplots_adjust()
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

def efig6():
    names   =['samp_0', 'samp_1', 'samp_2', 'samp_3']
    labs    =['Diffuse dwarfs', '$7\leq\Sigma_*$<15', '$15\leq\Sigma_*$<25', 'Compact dwarfs']
    colos   =['blue', 'green', 'orange', 'red']
    alabs   =['a', 'b', 'c', 'd']

    fsm,fhm,scatter =funcs().read_data('Wang_EDfig6.xlsx', 'SHMR', ['x', 'y', 'y_er'])

    fig,axes    =plt.subplots(2, 2)
    fig.set_size_inches(f_size*2, f_size*2)

    for i in range(len(names)):
        ax          =axes[i//2, i%2]
        ax.text(0.1, 0.7, alabs[i], transform=ax.transAxes,fontsize=60, ha='center', va='center', color='black', weight='bold')

        if i == 0:
            sms     =[8.12, 7.45, 7.82, 8.35, 7.65, 8.11, 8.14, 8.14]
            hms     =[10.14, 10.08, 10.32, 9.92, 10.16, 10.00, 10.76, 10.39]
            er1     =[0.11, 0.19, 0.23, 0.18, 0.14, 0.17, 0.22, 0.35]
            er2     =[0.24, 0.38, 0.51, 0.32, 0.26, 0.30, 0.17, 0.44]

            ax.errorbar(sms,hms, yerr=(er1, er2), zorder=100, ms=e_size*0.9, fmt='o', color='cyan', ecolor='cyan', elinewidth=el_wid/2, capsize=c_size/5, label='Kong et al. 2022')

        xs,xer1,xer2,ys,yer1,yer2   =funcs().read_data('Wang_EDfig6.xlsx', names[i], ['x', 'x_er_low', 'x_er_up', 'y', 'y_er_low', 'y_er_up'])

        ax.errorbar(xs,ys, xerr=(xer1,xer2), yerr=(yer1,yer2), zorder=-100, ms=e_size/30, fmt='o', color=colos[i], ecolor='grey', elinewidth=el_wid/2, alpha=0.3, capsize=c_size/5)
        ax.scatter(xs,ys, marker='o',s=5.5*s_size,color=colos[i], zorder=1, label=labs[i])

        ax.plot(fsm,fhm, color='teal', zorder=-50)
        ax.fill_between(fsm, np.array(fhm)+np.array(scatter), np.array(fhm)-np.array(scatter), color='teal', alpha=0.5, zorder=-50, label='SHMR')

        style().sefig6(ax, r'log $M_{\rm h, HI}$ (M$_{\odot}$)')

    save= fig_rou+'Wang_EDfig6.pdf'
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

def efig7():
    names   =['samp_0', 'samp_1', 'samp_2', 'samp_3']
    labs    =['Diffuse dwarfs', '$7\leq\Sigma_*$<15', '$15\leq\Sigma_*$<25', 'Compact dwarfs']
    colos   =['blue', 'green', 'orange', 'red']

    fig,ax          =plt.subplots()
    fig.set_size_inches(f_size, f_size)

    for i in range(len(names)):
        xs,ys   =funcs().read_data('Wang_EDfig7.xlsx', names[i], ['x', 'y'])
        ax.scatter(xs,ys, color=colos[i], s=s_size, zorder=-i, alpha=0.4, label=labs[i])

        xs,ys   =funcs().read_data('Wang_EDfig7.xlsx', 'dis_'+names[i], ['x', 'y'])
        ax.plot(xs,ys, color=colos[i],linewidth=l_wid, zorder=i)

        style().sefig7(ax)

    save= fig_rou+'Wang_EDfig7.pdf'
    plt.savefig(save, bbox_inches='tight', pad_inches=0.2)
    plt.close(fig)

    return save

class style():
    def sub_draw_2pccf(self, ax, tlab):
        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        if tlab == 0:
            ux_1,ux_2,uy_1,uy_2 =0.07, 32, 2, 2000
            ax.legend(loc='upper right', fontsize=le_size*0.75, frameon=False)

        elif tlab == 1:
            ux_1,ux_2,uy_1,uy_2 =0.07, 32, 2, 400
            ax.legend(loc='upper right', fontsize=le_size/40*26, frameon=False)

        ax.axis([ux_1,ux_2,uy_1,uy_2])

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax1_1.set_xscale("log")
        ax1_2.set_yscale("log")

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax.set_xlabel(r'$r_{\rm p}$ ($h^{-1}\rm Mpc$)', fontsize=l_size)
        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax.set_ylabel(r'$w_{\rm p}(r_{\rm p})$ ($h^{-1}\rm Mpc$)', fontsize=l_size)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

    def sefig1_abc(self, ax,ty,ux_1,ux_2,uy_1,uy_2, xlab,ylab):
        if ty == 0.6 or ty == 1.6:
            ax.axhline(y=ty, color='grey', linewidth=5, zorder=100)
        else:
            ax.plot([-0.2, 0.6], [1.6, 1.6], color='grey', linewidth=l_wid, zorder=100)
            ax.axvline(x=0.6, color='grey', linewidth=5, zorder=100)

        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        ax.axis([ux_1,ux_2,uy_1,uy_2])

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        if "Sigma" in xlab:
            ax.set_xscale("log")
            ax1_1.set_xscale("log")
        else:
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax1_1.xaxis.set_major_locator(MultipleLocator(0.5))
            ax1_1.xaxis.set_minor_locator(MultipleLocator(0.1))

        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1_2.yaxis.set_major_locator(MultipleLocator(0.5))
        ax1_2.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax.set_xlabel(xlab,fontsize=l_size)
        ax.set_ylabel(ylab,fontsize=l_size)

        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

    def sefig1_d(self, ax):
        ax.legend(loc='upper right', fontsize=le_size, frameon=False)

        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        ux_1,ux_2,uy_1,uy_2 =0.001, 0.078, 0,70
        ax.axis([ux_1,ux_2,uy_1,uy_2])

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax.xaxis.set_major_locator(MultipleLocator(0.03))
        ax.xaxis.set_minor_locator(MultipleLocator(0.01))
        ax1_1.xaxis.set_major_locator(MultipleLocator(0.03))
        ax1_1.xaxis.set_minor_locator(MultipleLocator(0.01))

        ax.yaxis.set_major_locator(MultipleLocator(20))
        ax.yaxis.set_minor_locator(MultipleLocator(5))
        ax1_2.yaxis.set_major_locator(MultipleLocator(20))
        ax1_2.yaxis.set_minor_locator(MultipleLocator(5))

        ax.set_xlabel(r'$z$',fontsize=l_size)
        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax.set_ylabel(r'Probability Density Function',fontsize=l_size)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

    def sefig2(self, ax1, ax2, oo):
        ax1.legend(fontsize=37,loc='upper right', frameon=False)
        ax2.legend(fontsize=37,loc='upper left', frameon=False)
        ax2.axhline(y=1, color='grey', linewidth=3, linestyle=':', zorder=-10)

        ux_1,ux_2,uy_1,uy_2=0.07, 32, 2, 150
        dx_1,dx_2,dy_1,dy_2=0.07, 32, 0.5, 2.9

        ax1_1 = ax1.twiny()
        ax1_2 = ax1.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)

        ax2_1 = ax2.twiny()
        ax2_2 = ax2.twinx()
        ax2_1.set_xlim(dx_1,dx_2)
        ax2_2.set_ylim(dy_1,dy_2)

        ax1.axis([ux_1,ux_2,uy_1,uy_2])
        ax2.axis([dx_1,dx_2,dy_1,dy_2])

        ax1.spines['bottom'].set_linewidth(b_wid)
        ax1.spines['left'].set_linewidth(b_wid)
        ax1.spines['top'].set_linewidth(b_wid)
        ax1.spines['right'].set_linewidth(b_wid)

        ax2.spines['bottom'].set_linewidth(b_wid)
        ax2.spines['left'].set_linewidth(b_wid)
        ax2.spines['top'].set_linewidth(b_wid)
        ax2.spines['right'].set_linewidth(b_wid)

        ax1.set_xscale("log")
        ax1.set_yscale("log")
        ax1_1.set_xscale("log")
        ax1_2.set_yscale("log")
        ax2_1.set_xscale("log")

        if oo == 0:
            ax1.set_ylabel(r'$w_{\rm p}(r_{\rm p})$ ($h^{-1}\rm Mpc$)',fontsize=l_size)
            ax2.set_ylabel(r'$\rm Ratio$',fontsize=l_size)
        else:
            ax1.set_yticklabels([])
            ax2.set_yticklabels([])

        ax2.set_xlabel(r'$r_{\rm p}$ ($h^{-1}\rm Mpc$)', fontsize=l_size)

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])
        ax2_1.set_xticklabels([])
        ax2_2.set_yticklabels([])

        ax1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax2.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax2.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax2_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax2_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax2_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax2_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax2.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax2.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.yaxis.set_major_locator(MultipleLocator(0.5))
        ax2_2.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax2_2.yaxis.set_major_locator(MultipleLocator(0.5))

    def sefig3(self, ax, l0s,l1s):
        ux_1,ux_2   =0.008, 0.062
        uy_1,uy_2   =0.001,2.5
        ax.set_ylim(uy_1,uy_2)

        legend1 = ax.legend(handles=l0s, loc='lower left', fontsize=le_size*0.6, frameon=False)
        ax.add_artist(legend1)
        legend2 = ax.legend(handles=l1s, loc='upper right', fontsize=le_size*0.6, frameon=False)

        ax.set_yscale("log")

        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        ax.xaxis.set_major_locator(MultipleLocator(0.01))
        ax.xaxis.set_minor_locator(MultipleLocator(0.005))

        ax.set_xlabel(r'z',fontsize=l_size)
        ax.set_ylabel(r'$n(z)\ (n_0)$',fontsize=l_size)

        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)
        ax1_2.set_yscale("log")

        ax1_1.xaxis.set_major_locator(MultipleLocator(0.01))
        ax1_1.xaxis.set_minor_locator(MultipleLocator(0.005))

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

    def sefig6(self, ax, ylab):
        ux_1,ux_2,uy_1,uy_2 =7.4, 9.1, 9.5, 12.0
        ax.axis([ux_1,ux_2,uy_1,uy_2])

        ax.legend(loc='upper left', fontsize=le_size, frameon=False)

        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        ax.set_xlabel(r'log $M_*$ (M$_{\odot}$)', fontsize=l_size)
        ax.set_ylabel(ylab,fontsize=l_size)

        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_major_locator(MultipleLocator(0.5))
        ax1_1.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1_2.yaxis.set_major_locator(MultipleLocator(0.5))
        ax1_2.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

    def sefig7(self, ax):
        ux_1,ux_2,uy_1,uy_2 =7.5, 9.05, 8.1, 10.5
        ax.axis([ux_1,ux_2,uy_1,uy_2])
        ax.legend(loc='upper left', fontsize=le_size, frameon=False)

        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        ax.xaxis.set_major_locator(MultipleLocator(0.5))
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax.set_xlabel(r'log $M_*$ (M$_{\odot}$)',fontsize=l_size)
        ax.set_ylabel(r'log $M_{\rm HI}$ (M$_{\odot}$)',fontsize=l_size)

        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ax.get_xlim())
        ax1_2.set_ylim(ax.get_ylim())

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_major_locator(MultipleLocator(0.5))
        ax1_1.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1_2.yaxis.set_major_locator(MultipleLocator(0.5))
        ax1_2.yaxis.set_minor_locator(MultipleLocator(0.1))

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

    def sub_draw_bias(self, ax, tlab):
        ax.axhline(y=1, color='grey', linewidth=3, linestyle=':', zorder=-10)

        ax.spines['bottom'].set_linewidth(b_wid)
        ax.spines['left'].set_linewidth(b_wid)
        ax.spines['top'].set_linewidth(b_wid)
        ax.spines['right'].set_linewidth(b_wid)

        if 'Sig' in tlab:
            ax.set_xlabel(r'$\Sigma_*$ (M$_{\odot}$pc$^{-2}$)', fontsize=l_size)
            ux_1,ux_2,uy_1,uy_2 =4.3, 55, 0.77, 2.8
            ax.legend(loc='upper right', fontsize=le_size, frameon=False)

        elif 'hb' in tlab:
            ax.set_xlabel(r'log $M_{\rm h}$ (M$_{\odot}$)',fontsize=l_size)
            ux_1,ux_2,uy_1,uy_2 =9.7, 13.5, 0.77, 2.8
            ax.legend(loc='upper center', fontsize=le_size, frameon=False, bbox_to_anchor=(0.63, 1))

        elif 'nn_cut' in tlab:
            ax.set_xlabel(r'$\Sigma_*$ (M$_{\odot}$pc$^{-2}$)', fontsize=l_size)
            ux_1,ux_2,uy_1,uy_2 =4, 95, 0.7, 2.7
            ax.legend(loc='upper right', fontsize=le_size, frameon=False)

        elif 'sn' in tlab:
            ax.set_xlabel(r'S$\acute{e}$rsic index', fontsize=l_size)
            ux_1,ux_2,uy_1,uy_2 =0.7, 2.65, 0.6, 1.45
            ax.legend(loc='upper right', fontsize=le_size, frameon=False)

        elif 'ef4' in tlab:
            ax.set_xlabel(r'$\Sigma_*$ (M$_{\odot}$pc$^{-2}$)', fontsize=l_size)
            ux_1,ux_2,uy_1,uy_2 =4, 55, 0.6, 3.8
            ax.legend(loc='upper right', fontsize=le_size, frameon=False)

        ax.set_ylabel('Relative bias',fontsize=l_size)
        ax.axis([ux_1,ux_2,uy_1,uy_2])

        ax1_1 = ax.twiny()
        ax1_2 = ax.twinx()
        ax1_1.set_xlim(ux_1,ux_2)
        ax1_2.set_ylim(uy_1,uy_2)

        ax.yaxis.set_major_locator(MultipleLocator(0.5))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
        ax1_2.yaxis.set_major_locator(MultipleLocator(0.5))
        ax1_2.yaxis.set_minor_locator(MultipleLocator(0.1))

        if 'Sig' in tlab or 'nn_cut' in tlab or 'ef4' in tlab:
            ax.set_xscale("log")
            ax1_1.set_xscale("log")

            from matplotlib.ticker import NullFormatter

            ax.xaxis.set_major_formatter(NullFormatter())  # 主要刻度无文本
            ax.xaxis.set_minor_formatter(NullFormatter())  # 次要刻度无文本

            custom_ticks = [5, 10, 50]  # 自定义刻度位置
            custom_labels = ['5.0', '10.0', '50.0']  # 自定义刻度标签

            ax.set_xticks(custom_ticks)  # 设置主要刻度位置
            ax.set_xticklabels(custom_labels)
            ax1_1.set_xticks(custom_ticks)  # 设置主要刻度位置

        elif 'hb' in tlab:
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax1_1.xaxis.set_major_locator(MultipleLocator(0.5))
            ax1_1.xaxis.set_minor_locator(MultipleLocator(0.1))

        elif 'sn' in tlab:
            ax.xaxis.set_major_locator(MultipleLocator(0.5))
            ax.xaxis.set_minor_locator(MultipleLocator(0.1))
            ax1_1.xaxis.set_major_locator(MultipleLocator(0.5))
            ax1_1.xaxis.set_minor_locator(MultipleLocator(0.1))

            ax.yaxis.set_major_locator(MultipleLocator(0.2))
            ax.yaxis.set_minor_locator(MultipleLocator(0.05))
            ax1_2.yaxis.set_major_locator(MultipleLocator(0.2))
            ax1_2.yaxis.set_minor_locator(MultipleLocator(0.05))

        else:
            ax.set_xscale("log")
            ax1_1.set_xscale("log")

        ax.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax1_1.xaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=2,length=15)
        ax1_2.yaxis.set_tick_params(which='major',direction='in',width=4,length=25)
        ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=2,length=15)

        ax.xaxis.set_tick_params(labelsize=t_size, pad=l_pad)
        ax.yaxis.set_tick_params(labelsize=t_size, pad=l_pad)

        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

class funcs():
    def read_data(self, name, sheet_name, keys):
        df  =pd.read_excel(frou+name, sheet_name=sheet_name)
        lis =[]
        
        for i in keys:
            lis.append(df[i])
        
        return lis

    def make_dir(self, rou):
        if os.path.exists(rou)==False:
            os.makedirs(rou)

    def convert_to_jpg(self, path0):
        path            =path0[:len(path0)-3]
        image           =convert_from_path(path0)[0]
        width,height    =image.size

        combined_image = Image.new('RGB', (width,height))
        combined_image.paste(image, (0, 0))
        combined_image.save(path+'jpeg', dpi=(300, 300))

        f_size  =os.path.getsize(path+'jpeg') / (1024 ** 2)
        if f_size > 10:
            self.convert_size(path+'jpeg')

    def convert_size(self, path):
        img             =Image.open(path)
        width, height   =img.size

        new_width       =int(width * 0.3)
        new_height      =int(height * 0.3)
        resized_img     =img.resize((new_width, new_height), Image.ANTIALIAS)

        resized_img.save(path, dpi=(300, 300))

#End-----------------------------

frou        =os.getcwd()+'/data/'
fig_rou     =frou+'Figs/'
funcs().make_dir(fig_rou)

font_path   =frou+'arial.ttf'  # 替换为你的字体文件路径
custom_font =fm.FontProperties(fname=font_path)

warnings.filterwarnings("ignore", category=UserWarning, message="The PostScript backend does not support transparency")

Image.MAX_IMAGE_PIXELS =400000000

#--------------------------------

st  =time.time()

functions   =[fig1, efig1, efig2, efig3, efig4, efig5, efig6, efig7]

for i in range(len(functions)):
    tem     =functions[i]()

    if i != 0:
        funcs().convert_to_jpg(tem)    

print(time.time() - st)
