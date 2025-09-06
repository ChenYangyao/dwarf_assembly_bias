import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
import numpy as np
import os
import math
import random
import shutil
import sys
from astropy.io import fits
import time

pi=np.pi
c     =300000.0           
G     =6.754e-11             
ckm   =3.240779e-17      
ckg   =5.027854e-31

#==============================================================
#User manual===================================================

##Some useful functions======| check_increase(list) | find_gro(gid,gro_id) | ff1(line,Nmember,redshift,re,di,gid,gred,xra,xdec,ceorsa,xproperty) | maxmin(red1,red2) | ff(line,gid) | change1(rou,name) | rm_file(rou) | mkdir(rou,name) | mat(array,t) | theta(D1,D2,R1,R2) | to_dis(list,re,di) | ave_split(Nbin,para) | tcolos(num) | sort_list(1-11)(list(1-11))

## Reading and Output file===| read(1-20)(rou,name) | output(1-12)(list(1-12),rou,name) | output_shear(rou,Nmember,Nbin,redshift,re,di,gid,gred,hm,xra,xdec,ceorsa,abm) |  output_fit(1-12)(rou,name)

## Drawing pictures==========| class static_pic(rou,names,labs,xlab,title,colos,Nbin).main() | class len_pic(rou,name,colo,listyle,lab,title).paint() | class cross(rou,names).main() | class corre_pic(rou,name,colo,listyle,lab,title).paint() | region_pic(rou,ptitle,ra1,dec1,redl,lab1,colo1,x1,x2,y1,y2,sarea,dalpha)

##Calculating codes==========| class lensing().main |

#==============================================================
#Some small useful functions===================================

def cal_dis(lis, minx,maxx):
    Nbin    =30
    step    =(maxx-minx)/Nbin

    xs,ys,ys1=[],[],[]
    for i in range(Nbin):
        tem =[]
        for j in range(len(lis)):
            if (minx+i*step) <= lis[j] < (minx+(i+1)*step):
                tem.append(j)

        xs.append(minx+(i+0.5)*step)
        ys.append(tem)
        ys1.append(len(tem))

    tot =float(sum(ys1))
    ys1 =list(map(lambda i:i/tot, ys1))

    return xs,ys,ys1

def z_matching(rou,names):
    samps,dis   =[],[]
    for i in range(len(names)):
        ra,dec,red,nid,sm,ids   =read6(rou, names[i])
        xs,ys,ys1               =cal_dis(red, 0.0,0.08)

        samps.append([ra,dec,red,nid,sm,ids])
        dis.append([xs,ys,ys1])

    low         =[]
    for i in range(len(dis[0][0])):
        tem     =[]
        for j in range(len(dis)):
            tem.append(dis[j][2][i])
        low.append(min(tem))

    rou1    =rou+'z_matching/'
    make_dir(rou1)

    for i in range(len(dis)):
        tids    =[]
        for j in range(len(low)):
            tem     =int(len(dis[i][1][j])*low[j]/(dis[i][2][j]+0.0000001))
            tids    +=random.sample(dis[i][1][j], tem)

        if os.path.exists(rou1+'m_'+names[i]) == False:
            noutput6(tids, *samps[i], rou1, 'm_'+names[i])

        rou2    =rou1+'m_'+names[i]+'_small/'
        make_dir(rou2)
        shutil.copyfile(rou1+'m_'+names[i], rou2+'m_'+names[i])

        if os.path.exists(rou2+'m_'+names[i]+'.wrp') == False:
            os.system('/home/zwzhang/data/codes/cal_ab/2PCCF/crossdr72_small_scale '+rou2+'m_'+names[i]+' &')

        if os.path.exists(rou1+'m_'+names[i]+'.wrp') == False:
            os.system('/home/zwzhang/data/codes/cal_ab/2PCCF/crossdr72_large_scale '+rou1+'m_'+names[i]+' &')

def z_weighting(rou,names):
    samps,dis   =[],[]
    for i in range(len(names)):
        ra,dec,red,nid,sm,ids   =read6(rou, names[i])
        xs,ys,ys1               =cal_dis(red, 0.0,0.08)

        samps.append([ra,dec,red,sm,ids])
        dis.append([xs,ys,ys1])

    low     =dis[0][2]
    rou1    =rou+'z_weighting/'
    make_dir(rou1)

    for i in range(len(names)):
        wts     =np.zeros(len(samps[i][0]))
        for j in range(len(low)):
            tem     =low[j]/(dis[i][2][j]+0.0000001)
            for k in dis[i][1][j]:
                wts[k]  =tem

        if os.path.exists(rou1+'m_'+names[i]) == False:
            output6(samps[i][0],samps[i][1],samps[i][2],wts,samps[i][3],samps[i][4], rou1, 'm_'+names[i])

        rou2    =rou1+'m_'+names[i]+'_small/'
        make_dir(rou2)
        shutil.copyfile(rou1+'m_'+names[i], rou2+'m_'+names[i])

        if os.path.exists(rou2+'m_'+names[i]+'.wrp') == False:
            os.system('/home/zwzhang/data/codes/cal_ab/2PCCF_wts/crossdr72_small_scale '+rou2+'m_'+names[i]+' &')

        if os.path.exists(rou1+'m_'+names[i]+'.wrp') == False:
            os.system('/home/zwzhang/data/codes/cal_ab/2PCCF_wts/crossdr72_large_scale '+rou1+'m_'+names[i]+' &')

def split_proj_blocks(ras, decs):
    minx,maxx   =min(ras), max(ras)
    miny,maxy   =min(decs), max(decs)

    stepx       =1
    stepy       =1

    Nbinx       =int((maxx-minx)/float(stepx))
    Nbiny       =int((maxy-miny)/float(stepy))

    id1,id2     =[],[]

    blocks      =[[[] for i in range(Nbiny+1)] for j in range(Nbinx+1)]
   
    for i in range(len(ras)):
        tid1    =int((ras[i]-minx)/float(stepx))
        tid2    =int((decs[i]-miny)/float(stepy))

        id1.append(tid1)
        id2.append(tid2)

        blocks[tid1][tid2].append(i)

    return id1,id2,blocks

def find_sep(step):
    bases   =[0.01, 0.05, 0.1, 0.5, 1, 5, 10, 20, 30, 50, 80, 100, 150, 200, 300, 400, 600]
    bases   =list(reversed(bases))

    for i in range(len(bases)):
        if 3 <= (float(step)/bases[i]) <= 10:
            maj =bases[i]
            mio =bases[i]/2.0

            return maj,mio

    bases   =list(map(lambda i:i/2.0, bases))
    for i in range(len(bases)):
        if 3 <= (float(step)/bases[i]) <= 10:
            maj =bases[i]
            mio =bases[i]/2.0

            return maj,mio

    print(step, 'check it!!!')
    sys.exit(0)

def match_id(ids, i):
    aa=np.where(ids==i)
    tem=[]
    for j in range(len(aa)):
        for k in range(len(aa[j])):
            tem.append(aa[j][k])
    return tem

def cal_age(z, H0, WM, WV):
      WR = 0.        # Omega(radiation)
      WK = 0.        # Omega curvaturve = 1-Omega(total)
      c = 299792.458 # velocity of light in km/sec
      Tyr = 977.8    # coefficent for converting 1/H into Gyr
      DTT = 0.5      # time from z to now in units of 1/H0
      DTT_Gyr = 0.0  # value of DTT in Gyr
      age = 0.5      # age of Universe in units of 1/H0
      age_Gyr = 0.0  # value of age in Gyr
      zage = 0.1     # age of Universe at redshift z in units of 1/H0
      zage_Gyr = 0.0 # value of zage in Gyr
      DCMR = 0.0     # comoving radial distance in units of c/H0
      DCMR_Mpc = 0.0 
      DCMR_Gyr = 0.0
      DA = 0.0       # angular size distance
      DA_Mpc = 0.0
      DA_Gyr = 0.0
      kpc_DA = 0.0
      DL = 0.0       # luminosity distance
      DL_Mpc = 0.0
      DL_Gyr = 0.0   # DL in units of billions of light years
      V_Gpc = 0.0
      a = 1.0        # 1/(1+z), the scale factor of the Universe
      az = 0.5       # 1/(1+z(object))

      h = H0/100.
      WR = 4.165E-5/(h*h)   # includes 3 massless neutrino species, T0 = 2.72528
      WK = 1-WM-WR-WV
      az = 1.0/(1+1.0*z)
      age = 0.
      n=1000         # number of points in integrals
      for i in range(n):
        a = az*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        age = age + 1./adot

      zage = az*age/n
      zage_Gyr = (Tyr/H0)*zage
      DTT = 0.0
      DCMR = 0.0

      for i in range(n):
        a = az+(1-az)*(i+0.5)/n
        adot = np.sqrt(WK+(WM/a)+(WR/(a*a))+(WV*a*a))
        DTT = DTT + 1./adot
        DCMR = DCMR + 1./(a*adot)

      DTT = (1.-az)*DTT/n
      DCMR = (1.-az)*DCMR/n
      age = DTT+zage
      age_Gyr = age*(Tyr/H0)
      DTT_Gyr = (Tyr/H0)*DTT
      DCMR_Gyr = (Tyr/H0)*DCMR
      DCMR_Mpc = (c/H0)*DCMR
      return DTT_Gyr

def aoutputs1(list1,rou,tnum,num,name):
    output1(list1,rou,'list_'+str(num))
    count=[]
    for i in range(tnum):
        if os.path.exists(rou+'list_'+str(i))==True:
            count.append(i)
    if len(count)==tnum:
        dat=open(rou+name,'w+')
        for i in range(tnum):
            t1=read2(rou,'list_'+str(i))
            for j in range(len(t1)):
                print(t1[j],file=dat)
        dat.close()

def aoutputs2(list1,list2,rou,tnum,num,name):
    output2(list1,list2,rou,'list_'+str(num))
    count=[]
    for i in range(tnum):
        if os.path.exists(rou+'list_'+str(i))==True:
            count.append(i)
    if len(count)==tnum:
        dat=open(rou+name,'w+')
        for i in range(tnum):
            t1,t2=read2(rou,'list_'+str(i))
            for j in range(len(t1)):
                print(t1[j],t2[j],file=dat)
        dat.close()

def outs(count,list):
    nlist=[]
    for i in range(len(list)):
        tem=[]
        for j in count:
            tem.append(list[i][j])
        nlist.append(tem)
    return nlist

def noutput1(ids,list1,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],file=dat)
    dat.close()

def noutput2(ids,list1,list2,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],file=dat)
    dat.close()

def noutput3(ids,list1,list2,list3,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],file=dat)
    dat.close()

def noutput4(ids,list1,list2,list3,list4,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],file=dat)
    dat.close()

def noutput5(ids,list1,list2,list3,list4,list5,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],file=dat)
    dat.close()

def noutput6(ids,list1,list2,list3,list4,list5,list6,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],file=dat)
    dat.close()

def noutput7(ids,list1,list2,list3,list4,list5,list6,list7,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],file=dat)
    dat.close()

def noutput8(ids,list1,list2,list3,list4,list5,list6,list7,list8,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],file=dat)
    dat.close()

def noutput9(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],file=dat)
    dat.close()

def noutput10(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],file=dat)
    dat.close()

def noutput11(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],file=dat)
    dat.close()

def noutput12(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],file=dat)
    dat.close()

def noutput13(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],list13[i],file=dat)
    dat.close()

def noutput14(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,list14,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],list13[i],list14[i],file=dat)
    dat.close()

def noutput15(ids,list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,list14,list15,rou,name):
    dat=open(rou+name,'w+')
    for i in ids:
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],list13[i],list14[i],list15[i],file=dat)
    dat.close()

def check_increase(list):
    DD=0
    for i in range(len(list)-1):
        if list[i+1]>=list[i]:
            DD+=1
        else:
            return 0
    return 1  #increase

def find_gro(gid,gro_id):
    ids=[]
    y=mat(gid,gro_id)
    if y<(len(gid)-1):
        for i in range(y,len(gid)):
            if gid[i]==gro_id:
                ids.append(i)
            else:
                break
        for i in range(y,-1,-1):
            if gid[i]==gro_id:
                ids.append(i)
            else:
                break
    elif y==(len(gid)-1):
        for i in range(len(gid)-1,-1,-1):
            if gid[i]==gro_id:
                ids.append(i)
            else:
                break
    elif y==0:
        for i in range(len(gid)):
            if gid[i]==gro_id:
                ids.append(i)
            else:
                break
    return list(set(ids))

def ff1(line,Nmember,redshift,re,di,gid,gred,xra,xdec,ceorsa,xproperty):
    idx,property,IDS =[],[],[]
    lab,para,groupid,count1,count2=0,0,0,-1,0
    for i in range(len(gid)-line):
        if gid[line] == gid[line+i]:
            idx.append(line+i)
        else:
            break

    if len(idx) >= Nmember:
        count0=[]
        DD,cra,cdec = 0,0,0
        for i in idx:
            if ceorsa[i] == 1:
                if gred[i]<redshift:
                    count0.append(1)
                    cra,cdec  = xra[i],xdec[i]
                    DD = 1
        if len(count0)>1:
            print('NOooooo')
        if DD==1:
            for i in idx:
                if ceorsa[i]==2:
                    IDS.append(i)
                    property.append(abs(xproperty[i]))
            if len(property) >= (Nmember-1):
                property,ids=sort_list1(property)    

                y1=mat(re,gred[line])
                dis=di[y1]

                ras,decs,pdis=[],[],[]
                for i in range(Nmember-1):
                    ra1,dec1 = xra[IDS[ids[len(IDS)-i-1]]],xdec[IDS[ids[len(IDS)-i-1]]]
                    ras.append(ra1)
                    decs.append(dec1)

                for i in range(len(ras)):
                    pdis1=theta(cdec,decs[i],cra,ras[i])*dis
                    pdis.append(pdis1)

                para=sum(pdis)/len(pdis)
                groupid=gid[line]
                lab=1
                return len(idx)+line,lab,para,groupid,count1,count2
            else:
                return len(idx)+line,lab,para,groupid,count1,count2
        else:
            count1=gid[line]
            return len(idx)+line,lab,para,groupid,count1,count2
    else:
        count2=1
        return len(idx)+line,lab,para,groupid,count1,count2

def maxmin(red1,red2):
    maxr1 = max(red1)
    minr1 = min(red1)
    maxr2 = max(red2)
    minr2 = min(red2)
    rmax  = max(maxr1,maxr2)
    rmin  = min(minr1,minr2)
    return rmin,rmax

def ff(line,gid):
    idx,count = [],[]
    for j in range(len(gid)-line):
        if gid[line] == gid[line+j]:
            idx.append(line+j)
    return len(idx)+line,len(idx)

def change1(rou,name):
    if os.path.getsize(''+rou+name)!=0:
        dat = np.loadtxt(''+rou+name,unpack=True)
        ra  = dat[0][:]
        dec = dat[1][:]
        red = dat[2][:]
        sm  = dat[3][:]
        hm  = dat[4][:]
        dat1 = open(r''+rou+name+'_ra','w+')
        dat2 = open(r''+rou+name+'_dec','w+')
        dat3 = open(r''+rou+name+'_red','w+')
        dat4 = open(r''+rou+name+'_sm','w+')
        dat5 = open(r''+rou+name+'_hm','w+')
        for i in range(len(ra)):
            print(ra[i],file=dat1)
        dat1.close()
        for i in range(len(dec)):
            print(dec[i],file=dat2)
        dat2.close()
        for i in range(len(red)):
            print(red[i],file=dat3)
        dat3.close()
        for i in range(len(sm)):
            print(sm[i],file=dat4)
        dat4.close()
        for i in range(len(hm)):
            print(hm[i],file=dat5)
        dat5.close()
        return 1
    else:
        dat1 = open(r''+rou+name+'_ra','w+')
        dat2 = open(r''+rou+name+'_dec','w+')
        dat3 = open(r''+rou+name+'_red','w+')
        dat4 = open(r''+rou+name+'_sm','w+')
        dat5 = open(r''+rou+name+'_hm','w+')
        return 1

def rm_file(rou):
    if os.path.exists(rou):
        shutil.rmtree(rou,True)
        os.makedirs(rou)
    else:
        os.makedirs(rou)
    return 1

def make_dir(rou):
    if os.path.exists(rou)==False:
        os.makedirs(rou)
    else:
        return 'ok'

def mat(array,t):
    aa=check_increase(array)
    if aa!=1:
        print('Not increased!!!')
        sys.exit(0)

    low=0
    height=len(array)-1
    while low<=height:
        mid = int((low+height)/2)
        if array[mid]<t:
            low= mid+1
        elif array[mid]>t:
            height=mid-1
        else:
            return mid
    else:
        print('NOT RIGHT!!!',t)
        return -1

def theta(R1,D1, R2,D2):
    if D1==D2 and R1==R2:
        return 0
    else:
        xm0 = np.cos(pi/2.0+D1*pi/180.0)
        xm1 = np.cos(pi/2.0+D2*pi/180.0)
        xm2 = np.sin(pi/2.0+D1*pi/180.0)
        xm3 = np.sin(pi/2.0+D2*pi/180.0)
        xm4 = np.cos((R1-R2)*pi/180.0)
        aa1  = xm0*xm1+xm2*xm3*xm4
        if aa1>1:
            return 0.0
        else:
            return np.arccos(aa1)

def to_dis(list,re,di):
    dis = []
    for i in range(len(list)):
        y = mat(re,round(list[i],4))
        dis.append(di[y])
    return dis

def ave_split(Nbin,para):
    para1 =sorted(para)
    leng = len(para)
    step = leng//Nbin
    inter=[]
    """
    for i in range(Nbin+1):
        if 0<i<Nbin:
            inter.append(para1[i*step])
        elif i==Nbin:
            inter.append(para1[i*step-1])
        elif i==0:
            inter.append(para1[0])
    """

    for i in range(Nbin+1):
        if i!=0:
            inter.append(para1[i*step-1])
        else:
            inter.append(para1[0])

    idss = [[]]*Nbin
    for i in range(Nbin):
        tem = []
        for j in range(len(para)):
            if inter[i] <= para[j] <inter[i+1]:
                tem.append(j)
        idss[i]=tem
    return inter,idss

def tcolos(num):
    if num==1:
        return ['blue']
    elif num==2:
        return ['blue','red']
    elif num==3:
        return ['blue','red','darkgrey']
    elif num==4:
        return ['blue','green','red','darkgrey']
    elif num==5:
        return ['blue','green','violet','red','darkgrey']
    elif num==6:
        return ['blue','green','violet','orange','red','darkgrey']

def sort_list1(list1):
    ids=list(range(len(list1)))
    tt=list(zip(list1,ids))
    tx=sorted(tt)
    list1[:],ids[:]=zip(*tx)
    return list1,ids

def sort_list2(list1,list2):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,ids))
    tx=sorted(tt)
    list1[:],list2[:],ids[:]=zip(*tx)
    return list1,list2,ids

def sort_list3(list1,list2,list3):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],ids[:]=zip(*tx)
    return list1,list2,list3,ids

def sort_list4(list1,list2,list3,list4):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,ids

def sort_list5(list1,list2,list3,list4,list5):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,ids

def sort_list6(list1,list2,list3,list4,list5,list6):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,list6,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],list6[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,list6,ids

def sort_list7(list1,list2,list3,list4,list5,list6,list7):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,list6,list7,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],list6[:],list7[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,list6,list7,ids

def sort_list8(list1,list2,list3,list4,list5,list6,list7,list8):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,list6,list7,list8,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],list6[:],list7[:],list8[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,list6,list7,list8,ids

def sort_list9(list1,list2,list3,list4,list5,list6,list7,list8,list9):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,list6,list7,list8,list9,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],list6[:],list7[:],list8[:],list9[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,list6,list7,list8,list9,ids

def sort_list10(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],list6[:],list7[:],list8[:],list9[:],list10[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,ids

def sort_list11(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11):
    ids=list(range(len(list1)))
    tt=list(zip(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,ids))
    tx=sorted(tt)
    list1[:],list2[:],list3[:],list4[:],list5[:],list6[:],list7[:],list8[:],list9[:],list10[:],list11[:],ids[:]=zip(*tx)
    return list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,ids

#==============================================================
#Reading and Output file,the maximun columns of files is 20====

def read1(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    return l1

def read2(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    return l1,l2

def read3(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    return l1,l2,l3

def read4(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    return l1,l2,l3,l4

def read5(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    return l1,l2,l3,l4,l5

def read6(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    return l1,l2,l3,l4,l5,l6

def read7(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7

def read8(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8

def read9(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9

def read10(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10

def read11(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11

def read12(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12

def read13(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13

def read14(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14

def read15(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    l15=[float(item[14]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15

def read16(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    l15=[float(item[14]) for item in z]
    l16=[float(item[15]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16

def read17(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    l15=[float(item[14]) for item in z]
    l16=[float(item[15]) for item in z]
    l17=[float(item[16]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17

def read18(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    l15=[float(item[14]) for item in z]
    l16=[float(item[15]) for item in z]
    l17=[float(item[16]) for item in z]
    l18=[float(item[17]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18

def read19(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    l15=[float(item[14]) for item in z]
    l16=[float(item[15]) for item in z]
    l17=[float(item[16]) for item in z]
    l18=[float(item[17]) for item in z]
    l19=[float(item[18]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19

def read20(rou,name):
    with open(r''+rou+name,'r') as data:
        fp = data.readlines()

    z=[]
    for row in fp:
        i = row.split()
        z.append(i)
    l1=[float(item[0]) for item in z]
    l2=[float(item[1]) for item in z]
    l3=[float(item[2]) for item in z]
    l4=[float(item[3]) for item in z]
    l5=[float(item[4]) for item in z]
    l6=[float(item[5]) for item in z]
    l7=[float(item[6]) for item in z]
    l8=[float(item[7]) for item in z]
    l9=[float(item[8]) for item in z]
    l10=[float(item[9]) for item in z]
    l11=[float(item[10]) for item in z]
    l12=[float(item[11]) for item in z]
    l13=[float(item[12]) for item in z]
    l14=[float(item[13]) for item in z]
    l15=[float(item[14]) for item in z]
    l16=[float(item[15]) for item in z]
    l17=[float(item[16]) for item in z]
    l18=[float(item[17]) for item in z]
    l19=[float(item[18]) for item in z]
    l20=[float(item[19]) for item in z]
    return l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,l15,l16,l17,l18,l19,l20

def output1(list1,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],file = dat)
    dat.close()
    return 1

def output2(list1,list2,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],file=dat)
    dat.close()
    return 1

def output3(list1,list2,list3,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],file=dat)
    dat.close()
    return 1

def output4(list1,list2,list3,list4,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],file=dat)
    dat.close()
    return 1

def output5(list1,list2,list3,list4,list5,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],file=dat)
    dat.close()
    return 1

def output6(list1,list2,list3,list4,list5,list6,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],file=dat)
    dat.close()
    return 1

def output7(list1,list2,list3,list4,list5,list6,list7,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],file=dat)
    dat.close()
    return 1

def output8(list1,list2,list3,list4,list5,list6,list7,list8,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],file=dat)
    dat.close()
    return 1

def output9(list1,list2,list3,list4,list5,list6,list7,list8,list9,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],file=dat)
    dat.close()
    return 1

def output10(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],file=dat)
    dat.close()
    return 1

def output11(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],file=dat)
    dat.close()
    return 1

def output12(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],file=dat)
    dat.close()
    return 1

def output13(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],list13[i],file=dat)
    dat.close()
    return 1

def output14(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,list14,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],list13[i],list14[i],file=dat)
    dat.close()
    return 1

def output15(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,list13,list14,list15,rou,name):
    dat = open(r''+rou+name,'w+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],list13[i],list14[i],list15[i],file=dat)
    dat.close()
    return 1

def aoutput1(list1,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],file = dat)
    dat.close()
    return 1

def aoutput2(list1,list2,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],file=dat)
    dat.close()
    return 1

def aoutput3(list1,list2,list3,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],file=dat)
    dat.close()
    return 1

def aoutput4(list1,list2,list3,list4,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],file=dat)
    dat.close()
    return 1

def aoutput5(list1,list2,list3,list4,list5,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],file=dat)
    dat.close()
    return 1

def aoutput6(list1,list2,list3,list4,list5,list6,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],file=dat)
    dat.close()
    return 1

def aoutput7(list1,list2,list3,list4,list5,list6,list7,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],file=dat)
    dat.close()
    return 1

def aoutput8(list1,list2,list3,list4,list5,list6,list7,list8,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],file=dat)
    dat.close()
    return 1

def aoutput9(list1,list2,list3,list4,list5,list6,list7,list8,list9,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],file=dat)
    dat.close()
    return 1

def aoutput10(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],file=dat)
    dat.close()
    return 1

def aoutput11(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],file=dat)
    dat.close()
    return 1

def aoutput12(list1,list2,list3,list4,list5,list6,list7,list8,list9,list10,list11,list12,rou,name):
    dat = open(r''+rou+name,'a+')
    for i in range(len(list1)):
        print(list1[i],list2[i],list3[i],list4[i],list5[i],list6[i],list7[i],list8[i],list9[i],list10[i],list11[i],list12[i],file=dat)
    dat.close()
    return 1

def output_shear(rou,Nmember,Nbin,redshift,re,di,gid,gred,hm,xra,xdec,ceorsa,abm):
    line,count2=0,0
    para,groupid,check=[],[],[]
    while(line<len(gid)):
        line,lab,par,groid,count1,lab1=ff1(line,Nmember,redshift,re,di,gid,gred,xra,xdec,ceorsa,abm)
        if lab==1:
            para.append(par)
            groupid.append(groid)
        if count1!=-1:
            check.append(count1)
        if lab1==1:
            count2+=1

    print(check,"No central galaxy")
    print(len(para),len(groupid),'|',count2,"halo does not have 3 galaxies")
 
    lrac,ldecc,labmc,lhmc,lredc,lras,ldecs,labms,lhms,lreds = [],[],[],[],[],[],[],[],[],[] ## In here,abm can be changed into other parameters,such as stellar mass .etc
    leftgid=[i for i in list(set(gid)) if i not in list(set(groupid))] 

    for i in leftgid:
        nids=find_gro(gid,i)
        for j in nids:
            if ceorsa[j]==1:
                lrac.append(xra[j])
                ldecc.append(xdec[j])
                labmc.append(abm[j])
                lhmc.append(hm[j])
                lredc.append(gred[j])
            elif ceorsa[j]==2:
                lras.append(xra[j])
                ldecs.append(xdec[j])
                labms.append(abm[j])
                lhms.append(hm[j])
                lreds.append(gred[j]) 

    print('Left is done')

    output5(lrac,ldecc,lredc,labmc,lhmc,rou+'static/','lcens')
    output5(lras,ldecs,lreds,labms,lhms,rou+'static/','lsats')

    inter,idss=ave_split(Nbin,para)

    for i in range(Nbin):
        print(len(leftgid),len(idss[i]))
        output1(idss[i],rou,'ids_'+str(i))

        name0=['cen_ra0','cen_ra1','cen_ra2','cen_ra3','cen_ra4','cen_ra5','cen_ra6','cen_ra7']
        name1=['sat_ra0','sat_ra1','sat_ra2','sat_ra3','sat_ra4','sat_ra5','sat_ra6','sat_ra7']
        name2=['cen_dec0','cen_dec1','cen_dec2','cen_dec3','cen_dec4','cen_dec5','cen_dec6','cen_dec7']
        name3=['sat_dec0','sat_dec1','sat_dec2','sat_dec3','sat_dec4','sat_dec5','sat_dec6','sat_dec7']
        name4=['cen_red0','cen_red1','cen_red2','cen_red3','cen_red4','cen_red5','cen_red6','cen_red7']
        name5=['sat_red0','sat_red1','sat_red2','sat_red3','sat_red4','sat_red5','sat_red6','sat_red7']
        name6=['cen_hm0','cen_hm1','cen_hm2','cen_hm3','cen_hm4','cen_hm5','cen_hm6','cen_hm7']
        name7=['sat_hm0','sat_hm1','sat_hm2','sat_hm3','sat_hm4','sat_hm5','sat_hm6','sat_hm7']
        name8=['cen_abm0','cen_abm1','cen_abm2','cen_abm3','cen_abm4','cen_abm5','cen_abm6','cen_abm7']
        name9=['sat_abm0','sat_abm1','sat_abm2','sat_abm3','sat_abm4','sat_abm5','sat_abm6','sat_abm7']
        name12=['tes_c0','tes_c1','tes_c2','tes_c3','tes_c4','tes_c5','tes_c6','tes_c7']
        name13=['tes_s0','tes_s1','tes_s2','tes_s3','tes_s4','tes_s5','tes_s6','tes_s7']

        cen_ra,cen_dec,cen_red,cen_hm,cen_abm = [],[],[],[],[]
        sat_ra,sat_dec,sat_red,sat_hm,sat_abm = [],[],[],[],[]

        for j in idss[i]:
            nids=find_gro(gid,groupid[j])
            for k in nids:
                if ceorsa[k]==1:
                    cen_ra.append(xra[k])
                    cen_dec.append(xdec[k])
                    cen_red.append(gred[k])
                    cen_hm.append(hm[k])
                    cen_abm.append(abm[k])
                elif ceorsa[k]==2:
                    sat_ra.append(xra[k])
                    sat_dec.append(xdec[k])
                    sat_red.append(gred[k])
                    sat_hm.append(hm[k])
                    sat_abm.append(abm[k])

        dis_c = to_dis(cen_red,re,di)
        dis_s = to_dis(sat_red,re,di)

        oo = output1(cen_ra,rou+'static/',name0[i])
        oo = output1(cen_dec,rou+'static/',name2[i])
        oo = output1(cen_red,rou+'static/',name4[i])
        oo = output1(cen_hm,rou+'static/',name6[i])
        oo = output1(cen_abm,rou+'static/',name8[i])
        oo = output1(sat_ra,rou+'static/',name1[i])
        oo = output1(sat_dec,rou+'static/',name3[i])
        oo = output1(sat_red,rou+'static/',name5[i])
        oo = output1(sat_hm,rou+'static/',name7[i])
        oo = output1(sat_abm,rou+'static/',name9[i])

        oo = output4(cen_ra,cen_dec,cen_red,dis_c,rou,name12[i])
        oo = output4(sat_ra,sat_dec,sat_red,dis_s,rou,name13[i])
        
    return 1

def output_fit1(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],file=dat)
    dat.close()

def output_fit2(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],file=dat)
    dat.close()

def output_fit3(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],file=dat)
    dat.close()

def output_fit4(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],file=dat)
    dat.close()

def output_fit5(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],file=dat)
    dat.close()

def output_fit6(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],file=dat)
    dat.close()

def output_fit7(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],fp[1].data[i][6],file=dat)
    dat.close()

def output_fit8(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],fp[1].data[i][6],fp[1].data[i][7],file=dat)
    dat.close()

def output_fit9(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],fp[1].data[i][6],fp[1].data[i][7],fp[1].data[i][8],file=dat)
    dat.close()

def output_fit10(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],fp[1].data[i][6],fp[1].data[i][7],fp[1].data[i][8],fp[1].data[i][9],file=dat)
    dat.close()

def output_fit11(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],fp[1].data[i][6],fp[1].data[i][7],fp[1].data[i][8],fp[1].data[i][9],fp[1].data[i][10],file=dat)
    dat.close()

def output_fit12(rou,name):
    fp = fits.open(''+rou+name+'.fits')

    dat = open(r''+rou+name+'.txt','w+')
    for i in range(len(fp[1].data)):
        print(fp[1].data[i][0],fp[1].data[i][1],fp[1].data[i][2],fp[1].data[i][3],fp[1].data[i][4],fp[1].data[i][5],fp[1].data[i][6],fp[1].data[i][7],fp[1].data[i][8],fp[1].data[i][9],fp[1].data[i][10],fp[1].data[i][11],file=dat)
    dat.close()

#==============================================================
#Drawing pictures==============================================

class static_pic:
    def __init__(self,rou,names,labs,xlab,title,colos,Nbin):
        self.names  = names
        self.labs   = labs
        self.xlab   = xlab
        self.title  = title
        self.colos  = colos
        self.rou    = rou
        self.Nbin   = Nbin

    def re(self,rou,name):
        if os.path.getsize(''+self.rou+name)!=0:
            with open(r''+self.rou+name,'r') as data:
                fp =data.readlines()
            z=[]
            for row in fp:
                i=row.split()
                z.append(i)
            lis = [float(item[0]) for item in z]
            return lis
        else:
            return []

    def maxmin(self,list):
        tmax,tmin,tmean,tnum = [],[],[],[]
        for i in range(len(list)):
            if len(list[i])!=0:
                tmax.append(max(list[i]))
                tmin.append(min(list[i]))
                tmean.append(np.mean(list[i]))
                tnum.append(len(list[i]))
        return min(tmin),max(tmax),tnum,tmean

    def draw(self,pos,count,num,mean,minb,maxb):
        fig,ax = plt.subplots()
        fig.set_size_inches(18.5,10.5)

        for i in range(len(count)):
            if len(count[i])!=0:
                plt.plot(pos,count[i],color=self.colos[i],linewidth=3.5,label=self.labs[i]+',ave='+str(round(mean[i],2))+',num='+str(num[i]))

        plt.xlabel(r'$'+self.xlab+'$',fontsize=30,labelpad=0)
        plt.ylabel(r'$ PDF $',fontsize=30,labelpad=0)
        ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.025))       

        save=''+self.rou+self.title+'.png'

        plt.title(self.title,fontsize=30)
        plt.legend(fontsize=30,loc='upper right')
        plt.ylim([0,0.4])
        plt.xticks(fontsize=30)
        plt.yticks(fontsize=30)
        plt.savefig(save)
        plt.close(fig)
        return 1

    def main(self):
        tem,inter,pos,count=[],[],[],[[]]*self.Nbin
        for i in range(len(self.names)):
            tem.append(self.re(self.rou,self.names[i]))
        tmin,tmax,lnum,lmean = self.maxmin(tem)
        step = (tmax-tmin)/self.Nbin
        for i in range(self.Nbin+1):
            inter.append(i*step+tmin)

        for i in range(len(tem)):
            tem2=[]
            if len(tem[i])!=0:
                for j in range(self.Nbin):
                    tem1 = []
                    for k in range(len(tem[i])):
                        if inter[j]<=tem[i][k]<inter[j+1]:
                            tem1.append(k)
                    tem2.append(len(tem1)/len(tem[i]))
            print(sum(tem2))
            count[i]=tem2

        for i in range(self.Nbin):
            pos.append((inter[i]+inter[i+1])/2)
        aa = self.draw(pos,count,lnum,lmean,tmin,tmax)
        return 1

class len_pic:
    def __init__(self,rou,name,colo,listyle,lab,title):
        self.name     = name
        self.colo     = colo
        self.listyle  = listyle
        self.lab      = lab
        self.title    = title
        self.rou      = rou

    def mkdir(self,path):
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        else:
            return path

    def re(self,aa):
        lis_b=[]
        lis_c=[]
        lis_d=[]
        lis_e=[]
        lis_f=[]
        ftop =[]
        fbot =[]
        pos  =[]
        res  =[]

        with open(r''+self.rou+aa,'r') as data:
            fp = data.readlines()
        lis_a=[]
        for row in fp:
            i = row.split()
            lis_a.append(i)
        xaa = [float(item[0]) for item in lis_a]
        xaa = list(map(lambda i:math.log(i,10),xaa))
        xbb = [float(item[1]) for item in lis_a]
        xcc = [float(item[2]) for item in lis_a]

        aa,bb,cc=[],[],[]
        for i in range(len(xaa)):
            if xbb[i]!=np.nan and xcc[i]!=np.nan:
                aa.append(xaa[i])
                bb.append(xbb[i])
                cc.append(xcc[i])

        for i in range(len(bb)):
            if bb[i]>0:
                lis_b.append(i)
                if cc[i]>0:
                    lis_c.append(i)
        for i in lis_c:
            ttop = (bb[i]+cc[i])
            tbot = (bb[i]-cc[i])
            if ttop>0 and tbot>0:
                lis_d.append(i)
                lis_e.append(ttop)
                lis_f.append(tbot)
        tem = [i for i in range(len(bb)) if i not in lis_d]

        for i in lis_d:
            pos.append(aa[i])
            res.append(bb[i])

        lis_e = list(map(lambda i:math.log(i,10),lis_e))
        lis_f = list(map(lambda i:math.log(i,10),lis_f))
        res   = list(map(lambda i:math.log(i,10),res))

        for i in range(len(res)):
            ftop.append(lis_e[i]-res[i])
            fbot.append(res[i]-lis_f[i])

        lis_a1 , lis_b1 , tbot1 , ttop1 = self.lef(aa,bb,cc,tem)

        lx ,ly = self.sor(lis_a1 , lis_b1 , pos , res)

        return pos,res,ftop,fbot,tem,lis_a1,lis_b1,tbot1,ttop1,lx,ly

    def lef(self,a,b,c,rr):
        lis_a=[]
        lis_b=[]
        lis_c=[]
        for i in rr:
            lis_a.append(a[i])
            lis_b.append(b[i])
            lis_c.append(c[i])

        ttop=[]
        tbot=[]
        tval=[]
        for i in range(len(lis_a)):
            if lis_b[i]>0:
                ttop.append(math.log(lis_b[i]+lis_c[i],10)-math.log(lis_b[i],10))
                tbot.append(10)
                tval.append(math.log(lis_b[i],10))
            elif lis_b[i]<=0 and lis_c[i]+lis_b[i]>0:
                ttop.append(math.log(lis_b[i]+lis_c[i],10)+2)
                tbot.append(0)
                tval.append(-2)
            elif lis_b[i]<=0:
                ttop.append(0)
                tbot.append(0)
                tval.append(-2)
        return lis_a,tval,tbot,ttop

    def sor(self,aa,bb,a,b):
        id1 = []
        id2 = []
        c=(aa+a)
        d=sorted(c)
        x1 = np.zeros(len(c))
        y1 = np.zeros(len(c))
        for i in aa:
            for j in range(len(c)):
                if i==d[j]:
                    id1.append(j)
        for i in a:
            for j in range(len(c)):
                if i==d[j]:
                    id2.append(j)
        for i in range(len(id1)):
            x1[id1[i]] = aa[i]
            y1[id1[i]] = bb[i]
        for i in range(len(id2)):
            x1[id2[i]] = a[i]
            y1[id2[i]] = b[i]
        return x1,y1

    def plot(self,name,colo,listyle,lab):
        pos,res,ftop,fbot,tem,lis_a1,lis_b1,tbot1,ttop1,lx,ly = self.re(name)

        plt.errorbar(pos,res,yerr=(fbot,ftop),fmt='o',ecolor=colo,color=colo,elinewidth=4,capsize=4,label=listyle+','+lab)
        plt.errorbar(lis_a1,lis_b1,yerr=(tbot1,ttop1),fmt='o',ecolor=colo,color=colo,elinewidth=3.5,capsize=4)
        plt.plot(lx,ly,color=colo,linestyle=listyle,linewidth=4)

        return 1

    def paint(self):
        fig, ax=plt.subplots()
        fig.set_size_inches(18.5,10.5)

        for i in range(len(self.name)):
            cc = self.plot(self.name[i],self.colo[i],self.listyle[i],self.lab[i])

        plt.xlabel(r'$log(R/[Mpc/h])$',fontsize=30,labelpad=0)
        plt.ylabel(r'$log(ESD/[M_{sun}/pc^2])$',fontsize=30,labelpad=0)
        ax_d = ax.twinx()
        ax_c = ax.twiny()

        x1,x2=0.01,6
        y1,y2=-1.5,2.3
        ax_c.set_xlim(x1,x2)
        ax_d.set_ylim(y1,y2)

        ax.set_xlim(x1,x2)
        ax.set_ylim(y1,y2)

        ax.set_xscale("log")
        ax_c.set_xscale("log")

        #ax_c.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax_d.yaxis.set_minor_locator(MultipleLocator(0.1))
        #ax_c.xaxis.set_tick_params(which='major',direction='in')
        #ax_c.xaxis.set_tick_params(which='minor',direction='in')
        ax_d.yaxis.set_tick_params(which='major',direction='in')
        ax_d.yaxis.set_tick_params(which='minor',direction='in')
        #ax.xaxis.set_tick_params(which='major',direction='in')
        #ax.xaxis.set_tick_params(which='minor',direction='in')
        ax.yaxis.set_tick_params(which='major',direction='in')
        ax.yaxis.set_tick_params(which='minor',direction='in')

        ax_c.xaxis.set_tick_params(labelsize=25)
        ax_d.yaxis.set_tick_params(labelsize=25)

        ax_c.set_xticklabels([])
        ax_d.set_yticklabels([])

        ax.xaxis.set_tick_params(labelsize=25)
        ax.yaxis.set_tick_params(labelsize=25)

        #ax.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax.yaxis.set_minor_locator(MultipleLocator(0.1))
    
        print('ok')
        save= self.rou+self.name[len(self.name)-1]+'.png'
        plt.title(self.title,fontsize=30)

        ax.legend(fontsize=30)
        #ax.grid(which='minor')
        #ax.grid(which='major')
        plt.savefig(save)
        plt.savefig(save)
        plt.close(fig)
        return 1

class cross:
    def __init__(self,rou,names):
        self.rou   = rou
        self.names = names

    def mkdir(self,path):
        folder  =os.path.exists(path)
        if not folder:
            os.makedirs(path)
        return path

    def read(self,name):
        with open(r''+self.rou+name+'.wrp','r') as data:
            fp = data.readlines()
        z=[]
        for row in fp:
            i = row.split()
            z1=[]
            for j in i:
                z1.append(float(j))
            z.append(z1)
        return z

    def process(self,name1,name2):
        z1=self.read(name1)
        z2=self.read(name2)
        path = self.mkdir(self.rou+'out/')
        xs1,ys1,er1=[],[],[]
        xs2,ys2,er2=[],[],[]
        t1,t2=[],[]
        ave,err,xx = [],[],[]
        for i in range(len(z1)):
            tem1=[]
            xs1.append(z1[i][0])
            ys1.append(z1[i][1]*2)
            er1.append(z1[i][2])
            for j in range(len(z1[i])):
                if j>=3:
                    tem1.append(z1[i][j])
            t1.append(tem1)

        for i in range(len(z2)):
            tem2=[]
            xs2.append(z2[i][0])
            ys2.append(z2[i][1]*2)
            er2.append(z2[i][2])
            for j in range(len(z2[i])):
                if j>=3:
                    tem2.append(z2[i][j])
            t2.append(tem2)

        for i in range(len(t1)):
            tem=[]
            for j in range(len(t1[i])):
                if t1[i][j]>0 and t2[i][j]>0:
                    tem.append(math.log(t1[i][j]/t2[i][j],10))
            if tem!=[]:
                #print(math.log(ys1[i]/ys2[i],10),np.mean(tem))
                if ys1[i]>0 and ys2[i]>0:
                    ave.append(math.log(ys1[i]/ys2[i],10))
                    err.append(np.std(tem,ddof=0))
                    xx.append(xs1[i])

        dat1 = open(r''+self.rou+'out/'+name1,'w+')
        for i in range(len(xs1)):
            print(xs1[i],ys1[i],er1[i],file = dat1)
        dat1.close()

        dat2 = open(r''+self.rou+'out/'+name2,'w+')
        for i in range(len(xs2)):
            print(xs2[i],ys2[i],er2[i],file = dat2)
        dat2.close()

        dat3 = open(r''+self.rou+'out/err_'+name1+name2,'w+')
        for i in range(len(xx)):
            print(xx[i],ave[i],err[i],file = dat3)
        dat3.close()
        return 1

    def main(self):
        for i in range(1,len(self.names)):
            aa=self.process(self.names[i],self.names[0])

class corre_pic:
    def __init__(self,rou,name,colo,listyle,lab,title):
        self.name     = name
        self.colo     = colo
        self.listyle  = listyle
        self.lab      = lab
        self.title    = title
        self.rou      = rou

    def mkdir(self,path):
        isExist = os.path.exists(path)
        if not isExist:
            os.makedirs(path)
        else:
            return path

    def re(self,aa):
        lis_b=[]
        lis_c=[]
        lis_d=[]
        lis_e=[]
        lis_f=[]
        ftop =[]
        fbot =[]
        pos  =[]
        res  =[]

        with open(r''+self.rou+aa,'r') as data:
            fp = data.readlines()
        lis_a=[]
        for row in fp:
            i = row.split()
            lis_a.append(i)
        aa = [float(item[0]) for item in lis_a]
        aa = list(map(lambda i:math.log(i,10),aa))
        bb = [float(item[1]) for item in lis_a]
        cc = [float(item[2]) for item in lis_a]

        for i in range(len(bb)):
            if bb[i]>0:
                lis_b.append(i)
                if cc[i]>0:
                    lis_c.append(i)
        for i in lis_c:
            ttop = (bb[i]+cc[i])
            tbot = (bb[i]-cc[i])
            if ttop>0 and tbot>0:
                lis_d.append(i)
                lis_e.append(ttop)
                lis_f.append(tbot)
        tem = [i for i in lis_b if i not in lis_d]

        for i in lis_d:
            pos.append(aa[i])
            res.append(bb[i])

        lis_e = list(map(lambda i:math.log(i,10),lis_e))
        lis_f = list(map(lambda i:math.log(i,10),lis_f))
        res   = list(map(lambda i:math.log(i,10),res))

        for i in range(len(res)):
            ftop.append(lis_e[i]-res[i])
            fbot.append(res[i]-lis_f[i])

        lis_a1 , lis_b1 , tbot1 , ttop1 = self.lef(aa,bb,cc,tem)

        lx ,ly = self.sor(lis_a1 , lis_b1 , pos , res)

        return pos,res,ftop,fbot,tem,lis_a1,lis_b1,tbot1,ttop1,lx,ly,aa,bb,cc

    def lef(self,a,b,c,rr):
        lis_a=[]
        lis_b=[]
        lis_c=[]
        for i in rr:
            lis_a.append(a[i])
            lis_b.append(b[i])
            lis_c.append(c[i])

        ttop=[]
        tbot=[]
        for i in range(len(lis_a)):
            ttop.append(math.log(lis_b[i]+lis_c[i],10)-math.log(lis_b[i],10))
            tbot.append(10)
        lis_b = list(map(lambda i:math.log(i,10),lis_b))
        return lis_a,lis_b,tbot,ttop

    def sor(self,aa,bb,a,b):
        id1 = []
        id2 = []
        c=(aa+a)
        d=sorted(c)
        x1 = np.zeros(len(c))
        y1 = np.zeros(len(c))
        for i in aa:
            for j in range(len(c)):
                if i==d[j]:
                    id1.append(j)
        for i in a:
            for j in range(len(c)):
                if i==d[j]:
                    id2.append(j)
        for i in range(len(id1)):
            x1[id1[i]] = aa[i]
            y1[id1[i]] = bb[i]
        for i in range(len(id2)):
            x1[id2[i]] = a[i]
            y1[id2[i]] = b[i]
        return x1,y1

    def plot(self,name,colo,listyle,lab):
        pos,res,ftop,fbot,tem,lis_a1,lis_b1,tbot1,ttop1,lx,ly,xx1,yy1,zz1 = self.re(name)
        self.ax1.errorbar(pos,res,yerr=(fbot,ftop),fmt='o',ecolor=colo,color=colo,elinewidth=3.5,capsize=4,label=listyle+','+lab)
        self.ax1.errorbar(lis_a1,lis_b1,yerr=(tbot1,ttop1),fmt='o',ecolor=colo,color=colo,elinewidth=3.5,capsize=4)
        self.ax1.plot(lx,ly,color=colo,linewidth=3.5,linestyle=listyle)
        self.ax1.legend(fontsize=24,loc=3)
        return 1

    def plot1(self):
        for j in range(1,len(self.name)):
            dat = np.loadtxt(''+self.rou+'err_'+self.name[j]+self.name[0],unpack=True)
            xs = dat[0][:]
            ys = dat[1][:]
            er = dat[2][:]

            xs = list(map(lambda i:math.log(i,10)+0.009*j,xs))
            self.ax2.errorbar(xs,ys,yerr=er,fmt='o',ecolor=self.colo[j],color=self.colo[j],linewidth=3.5,elinewidth=3.5,linestyle='-',label=self.lab[j]+'/'+self.lab[0])
            self.ax2.legend(fontsize=24,loc=3)
        return 1

    def paint(self):
        fig, (ax1,ax2)=plt.subplots(2,1,sharex=True,gridspec_kw = {'hspace':0})
        self.ax1,self.ax2 = ax1,ax2
        fig.set_size_inches(18.5,10.5)

        for i in range(len(self.name)):
            cc = self.plot(self.name[i],self.colo[i],self.listyle[i],self.lab[i])

        pics = self.plot1()
        ax2.axhline(y=0,linewidth=3,color='black')

        ax1.axis([-2.0,1.8,0.2,4.8])
        ax2.axis([-2.0,1.8,-0.2,0.2])

        ax1_1 = ax1.twiny()
        ax1_2 = ax1.twinx()
        ax1_1.set_xlim(-2.0,1.8)
        ax1_2.set_ylim(0.2,4.8)

        ax2_1 = ax2.twinx()
        ax2_1.set_ylim(-0.2,0.2)

        ax1_1.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1_2.yaxis.set_minor_locator(MultipleLocator(0.2))
        ax1_1.xaxis.set_tick_params(which='major',direction='in')
        ax1_1.xaxis.set_tick_params(which='minor',direction='in')
        ax1_2.yaxis.set_tick_params(which='major',direction='in')
        ax1_2.yaxis.set_tick_params(which='minor',direction='in')
        #ax1_1.xaxis.set_tick_params(labelsize=25)
        #ax1_2.yaxis.set_tick_params(labelsize=25)
        ax1_1.set_xticklabels([])
        ax1_2.set_yticklabels([])

        ax2_1.yaxis.set_minor_locator(MultipleLocator(0.05))
        ax2_1.yaxis.set_tick_params(which='major',direction='in')
        ax2_1.yaxis.set_tick_params(which='minor',direction='in')
        #ax2_1.yaxis.set_tick_params(labelsize=25)
        ax2_1.set_yticklabels([])

        ax2.set_xlabel(r'$log(r_{p}/(h^{-1}Mpc))$',fontsize=30,labelpad=0)
        ax1.set_ylabel(r'$log(w_{p}(r_{p})/(h^{-1}Mpc))$',fontsize=25,labelpad=0)
        ax1.xaxis.set_tick_params(which='major',direction='in')
        ax1.xaxis.set_tick_params(which='minor',direction='in')
        ax1.yaxis.set_tick_params(which='major',direction='in')
        ax1.yaxis.set_tick_params(which='minor',direction='in')
        ax1.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax1.yaxis.set_minor_locator(MultipleLocator(0.2))
        ax2.xaxis.set_tick_params(which='major',direction='in')
        ax2.xaxis.set_tick_params(which='minor',direction='in')
        ax2.yaxis.set_tick_params(which='major',direction='in')
        ax2.yaxis.set_tick_params(which='minor',direction='in')
        ax2.xaxis.set_minor_locator(MultipleLocator(0.1))
        ax2.yaxis.set_minor_locator(MultipleLocator(0.05))

        save= self.rou+'pic_'+self.lab[len(self.lab)-1]+'.png'

        fig.suptitle(self.title,fontsize=30)
        ax1.xaxis.set_tick_params(labelsize=25)
        ax1.yaxis.set_tick_params(labelsize=25)
        ax2.xaxis.set_tick_params(labelsize=25)
        ax2.yaxis.set_tick_params(labelsize=25)

        plt.savefig(save)
        plt.close(fig)
        return 1

def region_pic(rou,ptitle,ra1,dec1,redl,lab1,colo1,x1,x2,y1,y2,sarea,dalpha):
    fig,ax = plt.subplots()
    fig.set_size_inches(18.5,10.5)

    for i in range(len(ra1)):
        dec1[i] = list(map(lambda i:np.sin(i/180.0*pi),dec1[i]))
        plt.scatter(ra1[i],dec1[i],color=colo1[i],s=sarea,alpha=dalpha,label=lab1[i]+',num='+str(len(ra1[i]))+','+str(round(min(redl[i]),3))+'<z<='+str(round(max(redl[i]),3)))

    ax.xaxis.set_minor_locator(MultipleLocator(10))
    ax.yaxis.set_minor_locator(MultipleLocator(0.1))

    #ax.axhline(y=0.95,linewidth=3,color='black')
    #ax.axhline(y=-0.3,linewidth=3,color='black')

    #ax.axvline(99,linewidth=3,color='black')
    #ax.axvline(283,linewidth=3,color='black')

    plt.title(ptitle,fontsize=30)

    plt.xlabel(r'$α_{J}$',fontsize=30)
    plt.ylabel(r'$sin(δ_{J})$',fontsize=30)
    plt.xticks(fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(loc = 'lower left',fontsize=30)
    plt.axis([x1,x2,y1,y2])
    save = rou+'region_'+lab1[0]+'.png'
    plt.savefig(save)
    plt.close(fig)
    return 1

#==============================================================

class lensing:
    def __init__(self,Rmin,Rmax,Nbin):
        self.H0=100.0
        self.dH0=c/self.H0
        self.fac=(c*c*ckg)/(4.*pi*G*ckm)     
        
        self.Nsamp=2500
        self.Rmin=Rmin
        self.Rmax=Rmax
        self.Nbin=Nbin

    def bin_step(self):
        ftemp=[]
        xtemp=(math.log(self.Rmax,10)-math.log(self.Rmin,10))/Nbin
        for i in range(self.Nbin+1):
            ytemp=math.log(self.Rmin,10)+i*xtemp
            ftemp.append(10**ytmp)
        return ftemp

    def region(self,sec,nra,ndec):
        up,low=[],[]
        for i in range(nra):
            up.append(1+i*ndec)
            low.append(nra+i*ndec)
    
        secs=[]
        if sec<=ndec and sec not in up and sec not in low:
            secs=[sec,sec-1,sec+1,sec+ndec,sec+ndec+1,sec+ndec-1]
        elif sec<=ndec and sec in up:
            secs=[sec,sec+1,sec+ndec,sec+ndec+1]
        elif sec<=ndec and sec in low:
            secs=[sec,sec-1,sec+ndec,sec+ndec-1]
        elif sec>=(nra-1)*ndec and sec not in up and sec not in low:
            secs=[sec,sec+1,sec-1,sec-ndec,sec-ndec+1,sec-ndec-1]
        elif sec>=(nra-1)*ndec and sec in up:
            secs=[sec,sec+1,sec-ndec,sec-ndec+1]
        elif sec>=(nra-1)*ndec and sec in low:
            secs=[sec,sec-1,sec-ndec,sec-ndec-1]
        elif sec in up and sec>ndec and sec<(nra-1)*ndec:
            secs=[sec,sec+1,sec-ndec,sec+ndec,sec-ndec+1,sec+ndec+1]
        elif sec in low and sec>ndec and sec<(nra-1)*ndec:
            secs=[sec,sec-1,sec-ndec,sec+ndec,sec-ndec-1,sec-ndec-1]
        elif sec>ndec and sec<(nra-1)*ndec and sec not in up and sec not in low:
                secs=[sec,sec+1,sec-1,sec-ndec,sec+ndec,sec-ndec+1,sec-ndec-1,sec+ndec+1,sec+ndec-1]
        return secs

    def main(self):
        start=time.clock()

        print('Start loading the files...')
        ral,decl,zl,dl=read4('/data9/zwzhang/','')
        print('Finish loading the lens file')
        print('============================')
        
        ras,decs,e1,e2,res,er1,er2,zs,ds,spa=read10('/data9/zwzhang/Dats/','sorted_sourcewmap3')
        print('Finish loading the source file')
        print('============================')

        secs,slines,elines,ramin,ramax,decmin,decmax=read7('/data9/zwzhang/Dats/','info_sorted_sourcewmap3')
        print('Finish loading the info_sorted_sourcewmap3 file')
        print('End of loading the files',time.clock()-start)
        print('=============================')

        print('Start to create the list...')
        
        lists=[]
        for i in range(len(ral)):
            tem=[]
            tsec=0
            for j in range(len(secs)-1):
                if ramin[j]<=ral[i]<ramax[j] and decmin[j]<=decl[i]<decmax[j]:
                    tsec=secs[j]
            nsecs=self.region(tsec,91,28)
            for j in range(len(nsecs)):
                sl,el=0,0
                for k in range(len(secs)):
                    if nsecs[j]==secs[k]:
                        sl=slines[k]
                        el=elines[k]
                        break
                for k in range(sl,el):
                    sep1=abs(ral[i]-ras[k])
                    sep2=abs(decl[i]-decs[k])
                    if sep1<=4.0 and sep2<=4.0 and zs[k]>(zl[i]+0.1) and res[k]>1.0/3.0:
                        tem.append(k)
                if len(tem)!=0:
                    lists.append(tem)
                else:
                    lists.append(tem)
                    print('There is no source galaxy...')

        print('Finish creating the list...',time.clock()-start)
        print('=============================')
        print('Start calculating the shear...')

        gm1,gm2,wgt=[],[],[]
        e1=list(map(lambda i:i>=0,e1))
        resp=2.0*(1.0-np.var(e1))
        bins=self.bin_step()

        for i in range(len(ral)):
            tgm1,tmg2,twgt=np.zeros(len(bins)-1),np.zeros(len(bins)-1),np.zeros(len(bins))
            for j in lists[i]:
                Sig=self.fac*ds[j]/(dl[i]*(ds[j]-dl[i]))/(1.0+zl[i])**2
                wt=1.0/(0.1648+er1[j]**2+er2[j]**2)/Sig**2

                xm0=np.cos(pi/2.0-decs[j]*pi/180.0)
                xm1=np.cos(pi/2.0-decl[i]*pi/180.0)
                xm2=np.sin(pi/2.0-decs[j]*pi/180.0)
                xm3=np.sin*(pi/2.0-decl[i]*pi/180.0)
                xm4=np.cos((ras[j]-ral[i])*pi/180.0)
                the=np.arccos(xm0*xm1+xm2*xm3*xm4)

                tpsin=np.sin*(((ras[j]-ral[i])*np.cos(decl[j]*pi/180.0))*pi/180.0)
                tpcin=np.sin*((decs[j]-decl[i])*pi/180.0)
                sph=(2.0*tpsin**2)/np.sin(the)**2-1.0
                cph=(2.0*tpsin*tpcin)/np.sin(the)**2

                for k in range(len(bins)-1):
                    tdis=dl[i]*the*(1.0+zl[i])
                    if bins[k]<=tdis<bins[k+1]:
                        ang=spa[j]*pi/90.0
                        ep=np.cos(ang)*e1[j]-np.sin(ang)*e2[j]
                        em=np.sin(ang)*e1[j]+np.cos(ang)*e2[j]
                        e45=cph*ep+sph*em
                        et=-sph*ep+cph*em

                        tgm1[k]+=et*wt*Sig/rs[j]
                        tgm2[k]+=45*wt*Sig/rs[j]
                        twgt[k]+=wt
            gm1.append(tgm1)
            gm2.append(tgm2)
            wgt.append(twgt)
            
        print('Finish calculating the shear...',time.clock()-start)        
        lists=[]
        ras,decs,e1,e2,res,er1,er2,zs,ds,spa=0,0,0,0,0,0,0,0,0,0
        ral,decl,zl,dl=0,0,0,0
        print('Finish releasing the memory...')
        print('=============================')
        print('Start bootstraping...')

        g1err,g2err,gamm1,gamm2,inter=[],[],[],[],[]
        for i in range(len(bins)-1):
            inter.append((bins[i]+bins[i+1])/2)
            tga1,tga2,twt=np.zeros(len(bins)-1),np.zeros(len(bins)-1),np.zeros(len(bins)-1)
            nga1,nga2,nwt=np.zeros(len(bins)-1),np.zeros(len(bins)-1),np.zeros(len(bins)-1)
        
            for j in range(len(gm1)):
                nga1[i]+=gm1[j][i]
                nga2[i]+=gm2[j][i]
                nwt[i]+=wgt[j][i]

            gamm1.append(nga1[i]/nwt[i]/resp)
            gamm2.append(nga2[i]/nwt[i]/resp)

            ga1,ga2=[],[]
            for j in range(self.Nsamp):
                ibt=np.random.randint(0,len(ral),size=len(ral))
                for k in ibt:
                    tga1[i]+=gm1[k][i]
                    tga2[i]+=gm2[k][i]
                    twt[i]+=wgt[k][i]
    
                for k in range(len(tga1)):
                    ga1.append(tga1[k]/twt[k]/resp)
                    ga2.append(tga2[k]/twt[k]/resp)

            g1err.append(np.var(ga1))
            g2err.append(np.var(ga2))

        print('Finish bootstraping...',time.clock()-start) 
        print('=============================')
        print('Start outputing the results')

        output5(inter,gamm1,g1err,gamm2,g2err)
        print('Everything is done!',time.clock()-start)
        print('=============================')

def squre(reg, Nbiny, minx, miny, step):
    num     =(reg % Nbiny)
    num1    =(reg // Nbiny)

    lowy    =(miny+num*step)
    upy     =(miny+(num+1)*step)

    lowx    =(minx+num1*step)
    upx     =(minx+(num1+1)*step)

    return [(lowx+upx)/2, (lowy+upy)/2]

def draw_blocks(oblocks, extra, runners, len_blocks, sub, control, Nbiny, minx, miny, step):
    fig, ax=plt.subplots()
    fig.set_size_inches(21.5,10.5)

    #==============================================================

    oxs, oys=[],[]
    for i in range(len(oblocks)):
        if len(oblocks[i])!=0:
            for j in oblocks[i]:
                tp=(squre(j, Nbiny, minx, miny, step))
                oxs.append(tp[0])
                oys.append(tp[1])

    exs, eys=[],[]
    if len(extra) != 0:
        for i in extra:
            tp=(squre(i, Nbiny, minx, miny, step))
            exs.append(tp[0])
            eys.append(tp[1])

    rxs, rys=[],[]
    if len(runners) != 0:
        for i in range(len(runners)):
                tp=(squre(runners[i], Nbiny, minx, miny, step))
                rxs.append(tp[0])
                rys.append(tp[1])

    lxs, lys=[],[]
    for i in len_blocks:
        tp=(squre(i, Nbiny, minx, miny, step))
        lxs.append(tp[0])
        lys.append(tp[1])

    #==============================================================

    ax.scatter(lxs,lys, color='red', s=25, label='total lenses')

    if len(oxs) != 0:
        ax.scatter(oxs,oys, marker='x', color='lightblue', s=17, label='total source region')

    labs    =['minus', 'plus']
    if control in labs:
        ax.scatter(exs, eys, marker='o', color='orange', s=5, label='used sources at a time')

        tp=squre(runners[0], Nbiny, minx, miny, step)
        ax.scatter(rxs, rys, color='indigo', s=19, label='used lenses at a time')

    if control == 'blines':
        if len(exs) != 0:
            ax.scatter(exs, eys, color='black', s=23, label='broundary')

    #==============================================================
    ux_1,ux_2,uy_1,uy_2=0, 360, -90, 90

    step=2
    for i in range(181):
        tem0=(ux_1+i*step)
        tem1=(uy_1+i*step)
        if tem0 <= ux_2:
            ax.axvline(x=tem0,linewidth=0.5,color='gray')
        if tem1 <= uy_2:
            ax.axhline(y=tem1,linewidth=0.5,color='gray')

    ax.axis([ux_1,ux_2,uy_1,uy_2])
    ax1_1 = ax.twiny()
    ax1_2 = ax.twinx()
    ax1_1.set_xlim(ux_1,ux_2)
    ax1_2.set_ylim(uy_1,uy_2)

    ax1_1.set_xticklabels([])
    ax1_2.set_yticklabels([])

    bwith=3.5
    ax.spines['bottom'].set_linewidth(bwith)
    ax.spines['left'].set_linewidth(bwith)
    ax.spines['top'].set_linewidth(bwith)
    ax.spines['right'].set_linewidth(bwith)

    ax.xaxis.set_tick_params(which='major',direction='in',width=3,length=12)
    ax.xaxis.set_tick_params(which='minor',direction='in',width=1.5,length=7)
    ax.yaxis.set_tick_params(which='major',direction='in',width=3,length=12)
    ax.yaxis.set_tick_params(which='minor',direction='in',width=1.5,length=7)

    ax1_1.xaxis.set_tick_params(which='major',direction='in',width=3,length=12)
    ax1_1.xaxis.set_tick_params(which='minor',direction='in',width=1.5,length=7)
    ax1_2.yaxis.set_tick_params(which='major',direction='in',width=3,length=12)
    ax1_2.yaxis.set_tick_params(which='minor',direction='in',width=1.5,length=7)

    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    ax1_1.xaxis.set_major_locator(MultipleLocator(20))
    ax1_1.xaxis.set_minor_locator(MultipleLocator(4))

    ax.yaxis.set_major_locator(MultipleLocator(20))
    ax.yaxis.set_minor_locator(MultipleLocator(4))
    ax1_2.yaxis.set_major_locator(MultipleLocator(20))
    ax1_2.yaxis.set_minor_locator(MultipleLocator(4))

    ax.set_xlabel(r'RA',fontsize=30)
    ax.set_ylabel(r'DEC',fontsize=30)

    ax.xaxis.set_tick_params(labelsize=30)
    ax.yaxis.set_tick_params(labelsize=30)

    ax.legend(fontsize=21, loc='upper left')

    if control not in labs:
        plt.title('Num sources='+str(len(oxs)), fontsize=30)
    else:
        plt.title('Num sources='+str(len(exs))+', Num lenses='+str(len(rxs)), fontsize=30)

    save= '/home/zwzhang/data/Dats/box/fig_'+control+'_'+str(sub)+'.png'

    plt.subplots_adjust(wspace=0,hspace=0)
    plt.savefig(save)
    plt.close(fig)


#==============================================================
#END===========================================================
