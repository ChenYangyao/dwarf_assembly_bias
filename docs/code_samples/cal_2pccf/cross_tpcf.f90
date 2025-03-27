subroutine cross_tpcf(omegam,omegal,h0,nq,qra,qdec,qred,qwei, &
                      n,ra,dec,red,wei,nr,rra,rdec,rred,rwei,ppwei, &
                      nbs,sepsmin,dseps,logseps,nbp,seppmin,dsepp,logsepp, &
                      nbv,sepvmin,dsepv,logsepv,xif)
!------------------------------------------------------------------------------
! Purpose:
!   Compute redshift-space cross two-point correlation functions.
!
! Description:
!   This routine is to compute the redshift-space cross 2PCF for a  sample of
!   particles (hereafter Sample Q) with respect to a reference sample (Sample
!   D for the real sample and Sample R for the random),  using the estimators
!            Nr     QD(s)
!   xi(s) = ---- * ------- - 1,
!            Nq     QR(s)
!   and
!                Nr     QD(rp,pi)
!   xi(rp,pi) = ---- * ----------- - 1.
!                Nq     QR(rp,pi)
!   Listed below is the decription on the varables involed:
!   xi(s),xi(rp,pi)--redshift-space 2PCFs
!   s----------------seperation in redshift space
!   rp---------------seperation perpendicular to the line of sight
!   pi---------------seperation parallel to the line of sight
!   Nq---------------number of particles in Q
!   Nr---------------number of particles in R
!   QD(s),QD(rp,pi)--number of Q-D pairs with seperations s or (rp,pi)
!   QR(s),QR(rp,pi)--number of Q-R pairs with seperations s or (rp,pi)
!
!   Here the errors are estimated using the bootstrap  resampling  technique.
!
! Inputs:
!   Variable----Type--------Dimension---Description----------------------------
!   omegam------real*4------1-----------Matter density at z=0
!   omegal------real*4------1-----------Valumm density at z=0
!   h0----------real*4------1-----------Hubble's constant at z=0
!   nq----------integer-----1-----------Number of particles in Q
!   qra,qdec----real*8------nq----------RA and Dec of particles in Q
!   qred--------real*4------nq----------Redshift of particles in Q
!   qwei--------real*4------nq----------Weight of particles in Q
!   n-----------integer-----1-----------Number of particles in D
!   ra,dec------real*8------n-----------RA and Dec of particles in D
!   red---------real*4------n-----------Redshift of particles in D
!   wei---------real*4------n-----------Weight of particles in D
!   nr----------integer-----1-----------Number of particles in R
!   rra,rdec----real*8------nr----------RA and Dec of particles in R
!   rred--------real*4------nr----------Redshift of particles in R
!   rwei--------real*4------nr----------Weight of particles in R
!   ppwei-------real*4------1-----------External function to compute the weight
!                                       for each pair.  This weight is  usually
!                                       used to correct for fiber collisions.
!   nbs---------integer-----1-----------Number of s bins
!   sepsmin-----real*4------1-----------Minimum of s
!   dseps-------real*4------1-----------Bin size for s
!   logseps-----integer-----1-----------Set to 1 to bin s in log space
!   nbp---------integer-----1-----------Number of rp bins
!   seppmin-----real*4------1-----------Minimum of rp
!   dsepp-------real*4------1-----------Bin size for rp
!   logsepp-----integer-----1-----------Set to 1 to bin rp in log space
!   nbv---------integer-----1-----------Number of pi bins
!   sepvmin-----real*4------1-----------Minimum of pi
!   dsepv-------real*4------1-----------Bin size for pi
!   logsepv-----integer-----1-----------Set to 1 to bin pi in log space
!   xif---------char*200----1-----------To specify the names of output files
!
! Outputs:
!   The outputs will be written into a number of ASCII/binary file,  with the
!   names of the files are specified by "xif"(see below). The detailed format
!   for each file can be found in file [xif].log
!   Name----------Type--------Description--------------------------------------
!   [xif].cnt-----Binary------Pair counts
!   [xif].xis-----ASCII-------xi(s)
!   [xif].wrp-----ASCII-------Projected 2PCF w(rp)
!   [xif].xipv----ASCII-------xi(rp,pi)
!   [xif].log-----ASCII-------A log file for the computation
!
! Warning:
!   If file [xif].log already exists, the routine will STOP !!!
!
! Subroutines/functions called:
!
! Internal subroutines/functions:
!
! Revision history:
!   03-Oct-2005  Written by Cheng Li, Garching
!------------------------------------------------------------------------------
implicit none
!============
! Definitions
  !-------
  ! Inputs
  real(kind=4) :: omegam,omegal,h0
  integer :: nq,n,nr
  real(kind=4) :: qra(nq),qdec(nq),ra(n),dec(n),rra(nr),rdec(nr)
  real(kind=4) :: qred(nq),qwei(nq),red(n),wei(n),rred(nr),rwei(nr)
  real(kind=4),external :: ppwei
  integer :: nbs,nbp,nbv,logseps,logsepp,logsepv
  real(kind=4) :: sepsmin,dseps,seppmin,dsepp,sepvmin,dsepv
  character(len=200) :: xif
  !----------
  ! Constants
  real(kind=4),parameter :: vc=299792.46 ! speed of light
  integer,parameter :: nbts=100 ! number of bootstrap samples
  integer,parameter :: mxh1=100,mxh2=100,mxh3=100 ! dimension of HOC
  !------------
  ! Seperations
  real(kind=4),allocatable :: seps(:),sepp(:),sepv(:)
  !------------------
  ! Bootstrap weights
  real(kind=4),allocatable :: qwbts(:,:),wbts(:,:)
  !------------
  ! Pair counts
  real(kind=8),allocatable :: qds(:),qdpv(:,:),bqds(:,:),bqdpv(:,:,:)
  real(kind=8),allocatable :: qrs(:),qrpv(:,:)
  !-----------
  ! Boundaries
  real(kind=4) :: ramin,ramax,decmin,decmax
  real(kind=4) :: dcommin,dcommax
  common /bound/ramin,ramax,decmin,decmax,dcommin,dcommax
  !-----------------
  ! Other variables
  character(len=40) :: sysdate
  real(kind=4) :: t0,cpu0
  real(kind=4),external :: cputime
  real(kind=4),allocatable :: qdcom(:),dcom(:),rdcom(:)
  real(kind=4),external :: comdis
  character(len=80) :: ppf
  integer :: nok,i,j,k,logunit,ounit1,ounit2,ounit3
  real(kind=4) :: sf,bxi(nbts),xi,xierr,bwrp(nbts),wrp,wrperr
  real(kind=4),external :: stddev

!============
! Now, begin!

  !------------------
  ! Open the log file
  call get_unit(logunit)
  open(logunit,file=trim(xif)//'.log',status='replace')
  call systime(sysdate)
  write(logunit,'(a)') sysdate
  t0=cputime()
  
  !--------------------------------
  ! Write down the input parameters
  write(logunit,'(a)') '----------------------'
  write(logunit,'(a)') 'Cosmology parameters: '
  write(logunit,'(a,f5.2)') 'Matter density at z=0   : ',omegam
  write(logunit,'(a,f5.2)') 'Valumm density at z=0   : ',omegal
  write(logunit,'(a,f5.1)') "Hubble's constant at z=0: ",h0
  write(logunit,'(a)') '----------'
  write(logunit,'(a)') 'Sample Q: '
  write(logunit,'(a,i7)')   'Number of particles: ',nq
  write(logunit,'(a,2(g14.8,1x))') 'Range of RA        : ',minval(qra), maxval(qra)
  write(logunit,'(a,2(g14.8,1x))') 'Range of Dec       : ',minval(qdec),maxval(qdec)
  write(logunit,'(a,2(g12.6,1x))') 'Range of  redshift : ',minval(qred),maxval(qred)
  write(logunit,'(a)') '----------'
  write(logunit,'(a)') 'Sample D: '
  write(logunit,'(a,i7)')   'Number of particles: ',n
  write(logunit,'(a,2(g14.8,1x))') 'Range of RA        : ',minval(ra), maxval(ra)
  write(logunit,'(a,2(g14.8,1x))') 'Range of Dec       : ',minval(dec),maxval(dec)
  write(logunit,'(a,2(g12.6,1x))') 'Range of  redshift : ',minval(red),maxval(red)
  write(logunit,'(a)') '----------'
  write(logunit,'(a)') 'Sample R: '
  write(logunit,'(a,i7)')   'Number of particles: ',nr
  write(logunit,'(a,2(g14.8,1x))') 'Range of RA        : ',minval(rra), maxval(rra)
  write(logunit,'(a,2(g14.8,1x))') 'Range of Dec       : ',minval(rdec),maxval(rdec)
  write(logunit,'(a,2(g12.6,1x))') 'Range of  redshift : ',minval(rred),maxval(rred)
  write(logunit,'(a)') '-------------'
  write(logunit,'(a)') 'Seperations: '
  write(logunit,'(a,i3,1x,2(f6.2,1x),i2)') 's : ',nbs,sepsmin,dseps,logseps
  write(logunit,'(a,i3,1x,2(f6.2,1x),i2)') 'rp: ',nbp,seppmin,dsepp,logsepp
  write(logunit,'(a,i3,1x,2(f6.2,1x),i2)') 'pi: ',nbv,sepvmin,dsepv,logsepv
  write(logunit,'(a)') '----------------------------------------------'

  !----------------
  ! Set separations
  write(logunit,'(a)') 'Setting separations...'
  allocate(seps(nbs+1)); allocate(sepp(nbp+1)); allocate(sepv(nbv+1))
  call binsize(nbs,seps,sepsmin,dseps,logseps)
  call binsize(nbp,sepp,seppmin,dsepp,logsepp)
  call binsize(nbv,sepv,sepvmin,dsepv,logsepv)

  !---------------------------
  ! Compute comoving distances
  write(logunit,'(a)') 'Computing comoving distances...'
  allocate(qdcom(nq)); allocate(dcom(n)); allocate(rdcom(nr))
  do i=1,nq
     qdcom(i)=vc/h0*comdis(qred(i),omegam,omegal)
  end do
  do i=1,n
     dcom(i)=vc/h0*comdis(red(i),omegam,omegal)
  end do
  do i=1,nr
     rdcom(i)=vc/h0*comdis(rred(i),omegam,omegal)
  end do

  !-------------------------
  ! Write out the boundaries
  write(logunit,'(a)') 'Writting out the boundaries...'
  ramin=0
  ramax=360
  decmin=min(minval(qdec),minval(dec),minval(rdec))-0.1
  decmax=max(maxval(qdec),maxval(dec),maxval(rdec))+0.1
  decmin=max(decmin,-90.0)
  decmax=min(decmax, 90.0)
  dcommin=min(minval(qdcom),minval(dcom),minval(rdcom))-1.
  dcommax=max(maxval(qdcom),maxval(dcom),maxval(rdcom))+1.
  dcommin=max(dcommin,0.0)
  write(logunit,'(a,2(g14.8,a))') ' RA  : (',ramin,',',ramax,')'
  write(logunit,'(a,2(g14.8,a))') ' Dec : (',decmin,',',decmax,')'
  write(logunit,'(a,2(g12.6,a))') ' Dcom: (',dcommin,',',dcommax,')'

  !------------
  ! Count pairs
  write(logunit,'(a)') 'Getting bootstrap weights for Q and D catalogs...'
  allocate(qwbts(nq,nbts)); allocate(wbts(n,nbts))
  call wboots90(nq,nbts,qwbts)
  call wboots90(n,nbts,wbts)
  write(logunit,'(a)') 'Counting cross pairs in Q and D catalogs...'
  allocate(qds(nbs)); allocate(bqds(nbs,nbts))
  allocate(qdpv(nbp,nbv)); allocate(bqdpv(nbp,nbv,nbts))
  ppf=''
  cpu0=cputime()
  call sph_drcountb(nq,qra,qdec,qdcom,qwei,n,ra,dec,dcom,wei, &
     nbs,seps,nbp,sepp,nbv,sepv,nbts,qwbts,wbts,qds,qdpv,bqds,bqdpv, &
     mxh1,mxh2,mxh3,ppwei,1,ppf,0.,0.,0.,0.,nok)
  write(logunit,'(a,f7.1,a)') 'Finished counting QD: ',cputime()-cpu0,' seconds'

  !-------------------------
  write(logunit,'(a)') 'Counting cross pairs in Q and R catalogs...'
  allocate(qrs(nbs)); allocate(qrpv(nbp,nbv))
  ppf=''
  cpu0=cputime()
  call sph_drcount2(nq,qra,qdec,qdcom,qwei,nr,rra,rdec,rdcom,rwei, &
     nbs,seps,nbp,sepp,nbv,sepv,qrs,qrpv,mxh1,mxh2,mxh3,ppwei,0, &
     ppf,0.,0.,0.,0.,nok)
  write(logunit,'(a,f7.1,a)') 'Finished counting QR: ',cputime()-cpu0,' seconds'

  !----------------------------------
  ! Save the pair counts
  write(logunit,'(a)') 'Saving the pair counts...'
  call get_unit(ounit1)
  open(ounit1,file=trim(xif)//'.cnt',status='replace',form='unformatted')
  write(ounit1) nq,n,nr,nbts,nbs,sepsmin,dseps,logseps, &
     nbp,seppmin,dsepp,logsepp,nbv,sepvmin,dsepv,logsepv
  write(ounit1) bqds,bqdpv,qds,qdpv,qrs,qrpv
  close(ounit1)
  write(logunit,'(a)') 'The pair counts saved in file '//trim(xif)//'.cnt'

  !----------------------------------------
  ! Correct for the case of qrs=0, qrpv=0
  write(logunit,'(a)') 'Correcting for the cases of QR=0...'
  call extrapolate_dr(nbs,seps,qrs,0.2,1.,3.)
  call extrapolate_drpv(nbp,nbv,sepp,sepv,qrpv,0.2,1.,3.)

  !---------------------
  write(logunit,'(a)') 'Computing 2PCFs...'
  call get_unit(ounit1)
  open(ounit1,file=trim(xif)//'.xis',status='replace')
  sf=float(nr)/float(n)
  write(logunit,'(a)') '--------xi(s)------------'
  do i=1,nbs
     if(qrs(i)>0.1) then
        xi=qds(i)/qrs(i)*sf-1
        do j=1,nbts
           bxi(j)=bqds(i,j)/qrs(i)*sf-1
        end do
        xierr=stddev(nbts,bxi)
     else
        xi=-1
        xierr=-1
     end if
     write(ounit1,'(4(g12.6,1x))') 0.5*(seps(i)+seps(i+1)),xi,xierr,(bxi(k),k=1,nbts)
     write(logunit,'(3(g12.6,1x))') 0.5*(seps(i)+seps(i+1)),xi,xierr
  end do
  close(ounit1)
  write(logunit,'(a)') 'xi(s) written in file '//trim(xif)//'.xis'
  write(logunit,'(a)') '-------w(rp)------------'
  call get_unit(ounit1)
  open(ounit1,file=trim(xif)//'.wrp',status='replace')
  call get_unit(ounit2)
  open(ounit2,file=trim(xif)//'.xipv',status='replace')
  do i=1,nbp
     wrp=0
     bwrp=0
     do j=1,nbv
        if(qrpv(i,j)>0.1) then
           xi=qdpv(i,j)/qrpv(i,j)*sf-1
           do k=1,nbts
              bxi(k)=bqdpv(i,j,k)/qrpv(i,j)*sf-1
           end do
           xierr=stddev(nbts,bxi)
        else
           xi=-1
           bxi=-1
           xierr=-1
        end if
        wrp=wrp+xi*dsepv
        bwrp=bwrp+bxi*dsepv
        write(ounit2,'(104(g12.6,1x))') 0.5*(sepp(i)+sepp(i+1)), &
             0.5*(sepv(j)+sepv(j+1)),xi,xierr,(bxi(k),k=1,nbts)
     end do
     if(sum(qrpv(i,:))>0.1)then
        wrp = (sum(qdpv(i,:))/sum(qrpv(i,:))*sf-1)*(sepv(nbv+1)-sepv(1))
        do k=1,nbts
           bwrp(k)=(sum(bqdpv(i,:,k))/sum(qrpv(i,:))*sf-1)*(sepv(nbv+1)-sepv(1))
        end do
     else
        wrp = -1.0
        bwrp = -1.0
     end if
     wrperr=stddev(nbts,bwrp)
     write(ounit1,'(1003(g12.6,1x))') 0.5*(sepp(i)+sepp(i+1)),wrp,wrperr,(bwrp(k),k=1,nbts)
     write(logunit,'(3(g12.6,1x))') 0.5*(sepp(i)+sepp(i+1)),wrp,wrperr
  end do
  close(ounit1)
  close(ounit2)
  write(logunit,'(a)') 'w(rp) written in file '//trim(xif)//'.wrp'
  write(logunit,'(a)') 'xi(rp,pi) written in file '//trim(xif)//'.xipv'

  write(logunit,'(a,f7.1,a)') 'Time totally: ',cputime()-t0,' seconds'

  close(logunit)
  return
end subroutine cross_tpcf


subroutine sph_drcountb(npt,ra,dec,dcom,wei,nrpt,rra,rdec,rdcom,rwei, &
     nbs,seps,nbp,sepp,nbv,sepv,nbts,wbts,rwbts,dds,ddpv,bdds,bddpv, &
     mxh1,mxh2,mxh3,wfib,iwfib,ppf,rplim1,rplim2,rvlim1,rvlim2,nok)
  real(kind=4) :: ra(npt),dec(npt),rra(nrpt),rdec(nrpt)
  real(kind=4) :: ral,rau,decl,decu,hc1,hc2
  real(kind=4) :: dcom(npt),rdcom(nrpt),wei(npt),rwei(nrpt)
  real(kind=4) :: seps(nbs+1),sepp(nbp+1),sepv(nbv+1)
  real(kind=4) :: wbts(npt,nbts),rwbts(nrpt,nbts)
  real(kind=8) :: dds(nbs),ddpv(nbp,nbv)
  real(kind=8) :: bdds(nbs,nbts),bddpv(nbp,nbv,nbts)

  real(kind=4),allocatable :: seps2(:),sepp2(:)
  integer,allocatable :: hoc(:,:,:),ll(:)
  real(kind=4),external :: wfib
  character(len=80) :: ppf

  allocate(seps2(nbs+1)); allocate(sepp2(nbp+1))
  seps2=seps*seps
  sepp2=sepp*sepp

  !------------------
  ! reset the counts
  dds=0
  ddpv=0
  bdds=0
  bddpv=0

  rpmax=max(sepp(nbp+1),seps(nbs+1))
  rvmax=max(sepv(nbv+1),seps(nbs+1))

  !-----------------------------------
  ! set the boundaries of the survey
  call sph_boundary(ral,rau,decl,decu,dcoml,dcomu)
  nc1=mxh1
  nc2=mxh2
  nc3=(dcomu-dcoml)/rvmax
  if(nc3>mxh3) nc3=mxh3
  hc1=(decu-decl)/float(nc1)
  hc2=(rau-ral)/float(nc2)
  hc3=(dcomu-dcoml)/float(nc3)

  !----------------------------------------------------
  ! set up the head-of-chain and the linked list tables
  !write(*,'(a)') 'Setting up the HOC and LL tables...'
  allocate(hoc(mxh1,mxh2,mxh3)); allocate(ll(nrpt))
  call sph_hocll(hoc,mxh1,mxh2,mxh3,nc1,nc2,nc3,ll,nrpt,rra,rdec,rdcom)

  !------------------------
  ! calculate the counts
  !write(*,'(a,3(i4,1x))') 'Counting pairs in HOC with dimensions: ',nc1,nc2,nc3
  nok=0
  if(trim(ppf)/='') then ! open a file to save indices of pairs
     open(23,file=trim(ppf)//'.drpp',status='replace',form='unformatted')
     write(23) npt,nrpt
  end if
  do i=1,npt
     iq1=(dec(i)-decl)/hc1+1.
     iq2=(ra(i)-ral)/hc2+1.
     iq3=(dcom(i)-dcoml)/hc3+1.
     lp_jq3: do jq3=iq3-1,iq3+1
        if(jq3<1.or.jq3>nc3) cycle lp_jq3
        rcl=(jq3-1)*hc3+dcoml
        stm2=sthmax2(rpmax,dcom(i),rcl)
        dltdec=2.*asind(stm2)
        jq1m=dltdec/hc1+1.
        lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
           if(jq1>nc1.or.jq1<1) cycle lp_jq1
           if(jq1==iq1) then
              dltra=dalp(stm2,dec(i),jq1,hc1,decl,1)
           else
              dltra=dalp(stm2,dec(i),jq1,hc1,decl,0)
           end if
           jq2m=dltra/hc2+1.
           jq2max=iq2+jq2m
           jq2min=iq2-jq2m
           if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
           lp_jq2: do jq2=jq2min,jq2max
              if(jq2>nc2) then
                 jq2t=jq2-nc2
              else if(jq2<1) then
                 jq2t=jq2+nc2
              else
                 jq2t=jq2
              end if
              j=hoc(jq1,jq2t,jq3)
              do while(j/=0)
                 call sph_countbb(i,j,ra(i),dec(i),dcom(i),wei(i), &
                      rra(j),rdec(j),rdcom(j),rwei(j),wfib,iwfib, &
                      nbs,seps2,nbp,sepp2,nbv,sepv, &
                      npt,nrpt,nbts,wbts,rwbts,dds,ddpv,bdds,bddpv, &
                      rvmax,ick,iok,rp,rv,wpp)
                 if(trim(ppf)/=''.and.iok==1) then
                    if(rp>rplim1.and.rp<rplim2) then
                       if(rv>rvlim1.and.rv<rvlim2) then
                          write(23) i,j,rp,rv,wpp
                          nok=nok+1
                       end if
                    end if
                 end if
                 if(ick==1.and.jq3>iq3) cycle lp_jq2
                 j=ll(j)
              end do
           end do lp_jq2
        end do lp_jq1
     end do lp_jq3
  end do
  if(trim(ppf)/='') then
     close(23)
     !write(*,'(i8,a)') nok,' cross pairs saved in file '//trim(ppf)//'.drpp'
  end if
end subroutine sph_drcountb
subroutine sph_drcount2(npt,ra,dec,dcom,wei,nrpt,rra,rdec,rdcom,rwei, &
     nbs,seps,nbp,sepp,nbv,sepv,dds,ddpv,mxh1,mxh2,mxh3,wfib,iwfib, &
     ppf,rplim1,rplim2,rvlim1,rvlim2,nok)
  real(kind=4) :: ra(npt),dec(npt),rra(nrpt),rdec(nrpt)
  real(kind=4) :: ral,rau,decl,decu,hc1,hc2
  real(kind=4) :: dcom(npt),rdcom(nrpt),wei(npt),rwei(nrpt)
  real(kind=4) :: seps(nbs+1),sepp(nbp+1),sepv(nbv+1)
  real(kind=8) :: dds(nbs),ddpv(nbp,nbv)

  real(kind=4),allocatable :: seps2(:),sepp2(:)
  integer,allocatable :: hoc(:,:,:),ll(:)
  real(kind=4),external :: wfib
  character(len=80) :: ppf

  allocate(seps2(nbs+1)); allocate(sepp2(nbp+1))
  seps2=seps*seps
  sepp2=sepp*sepp

  !------------------
  ! reset the counts
  dds=0
  ddpv=0

  rpmax=max(sepp(nbp+1),seps(nbs+1))
  rvmax=max(sepv(nbv+1),seps(nbs+1))

  !-----------------------------------
  ! set the boundaries of the survey
  call sph_boundary(ral,rau,decl,decu,dcoml,dcomu)
  nc1=mxh1
  nc2=mxh2
  nc3=(dcomu-dcoml)/rvmax
  if(nc3>mxh3) nc3=mxh3
  hc1=(decu-decl)/float(nc1)
  hc2=(rau-ral)/float(nc2)
  hc3=(dcomu-dcoml)/float(nc3)

  !----------------------------------------------------
  ! set up the head-of-chain and the linked list tables
  !write(*,'(a)') 'Setting up the HOC and LL tables...'
  allocate(hoc(mxh1,mxh2,mxh3)); allocate(ll(nrpt))
  call sph_hocll(hoc,mxh1,mxh2,mxh3,nc1,nc2,nc3,ll,nrpt,rra,rdec,rdcom)

  !------------------------
  ! calculate the counts
  !write(*,'(a,3(i4,1x))') 'Counting pairs in HOC with dimensions: ',nc1,nc2,nc3
  nok=0
  if(trim(ppf)/='') then ! open a file to save indices of pairs
     open(23,file=trim(ppf)//'.drpp',status='replace',form='unformatted')
     write(23) npt,nrpt
  end if
  do i=1,npt
     iq1=(dec(i)-decl)/hc1+1.
     iq2=(ra(i)-ral)/hc2+1.
     iq3=(dcom(i)-dcoml)/hc3+1.
     lp_jq3: do jq3=iq3-1,iq3+1
        if(jq3<1.or.jq3>nc3) cycle lp_jq3
        rcl=(jq3-1)*hc3+dcoml
        stm2=sthmax2(rpmax,dcom(i),rcl)
        dltdec=2.*asind(stm2)
        jq1m=dltdec/hc1+1.
        lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
           if(jq1>nc1.or.jq1<1) cycle lp_jq1
           if(jq1==iq1) then
              dltra=dalp(stm2,dec(i),jq1,hc1,decl,1)
           else
              dltra=dalp(stm2,dec(i),jq1,hc1,decl,0)
           end if
           jq2m=dltra/hc2+1.
           jq2max=iq2+jq2m
           jq2min=iq2-jq2m
           if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
           lp_jq2: do jq2=jq2min,jq2max
              if(jq2>nc2) then
                 jq2t=jq2-nc2
              else if(jq2<1) then
                 jq2t=jq2+nc2
              else
                 jq2t=jq2
              end if
              j=hoc(jq1,jq2t,jq3)
              do while(j/=0)
                 call sph_count2(ra(i),dec(i),dcom(i),wei(i), &
                      rra(j),rdec(j),rdcom(j),rwei(j),wfib,iwfib, &
                      nbs,seps2,nbp,sepp2,nbv,sepv, &
                      dds,ddpv,rvmax,ick,iok,rp,rv,wpp)
                 if(trim(ppf)/=''.and.iok==1) then
                    if(rp>rplim1.and.rp<rplim2) then
                       if(rv>rvlim1.and.rv<rvlim2) then
                          write(23) i,j,rp,rv,wpp
                          nok=nok+1
                       end if
                    end if
                 end if
                 if(ick==1.and.jq3>iq3) cycle lp_jq2
                 j=ll(j)
              end do
           end do lp_jq2
        end do lp_jq1
     end do lp_jq3
  end do
  if(trim(ppf)/='') then
     close(23)
     !write(*,'(i8,a)') nok,' cross pairs saved in file '//trim(ppf)//'.drpp'
  end if
end subroutine sph_drcount2
subroutine sph_countbb(ipt,jpt,rai,deci,ri,wi,raj,decj,rj,wj,wfib,iwfib, &
     nbs,seps2,nbp,sepp2,nbv,sepv, &
     npt,nrpt,nbts,wbts,rwbts,dds,ddpv,bdds,bddpv,rvmax,ick,iok,rp,rv,wpp)

  real(kind=4) :: rai,deci,raj,decj
  real(kind=4),external :: wfib
  real(kind=8) :: dds(nbs),ddpv(nbp,nbv)
  real(kind=8) :: bdds(nbs,nbts),bddpv(nbp,nbv,nbts)
  real(kind=4) :: seps2(nbs+1),sepp2(nbp+1),sepv(nbv+1)
  real(kind=4) :: wbts(npt,nbts),rwbts(nrpt,nbts)
  real(kind=4) :: minus1,pi,shth2,rp2,r2
  minus1=-1
  pi=acos(minus1)
  pi=pi/180

  shth2=(sind(0.5*(deci-decj)))**2+ &
       cosd(deci)*cosd(decj)*(sind(0.5*(rai-raj)))**2
  rv=abs(ri-rj)
  rp2=4.*ri*rj*shth2
  r2=rv**2+rp2
  rp=sqrt(rp2)

  !--------------------
  ! Weight of this pair
  wpp=wi*wj*(iwfib*(wfib(shth2)-1)+1)

  iok=0
  if(r2>seps2(nbs+1)) goto 2
  do i=nbs,1,-1
     if(r2>seps2(i)) then
        dds(i)=dds(i)+wpp
        do j=1,nbts
           bdds(i,j)=bdds(i,j)+wpp*wbts(ipt,j)*rwbts(jpt,j)
        end do
        goto 2
     end if
  end do

2 if(rp2>sepp2(nbp+1).or.rv>sepv(nbv+1)) goto 3
  do i=nbp,1,-1
     if(rp2>sepp2(i)) then
        do j=nbv,1,-1
           if(rv>sepv(j)) then
              iok=1
              ddpv(i,j)=ddpv(i,j)+wpp
              do k=1,nbts
                 bddpv(i,j,k)=bddpv(i,j,k) &
                      +wpp*wbts(ipt,k)*rwbts(jpt,k)
              enddo
              goto 3
           endif
        enddo
        goto 3
     endif
  enddo

3 ick=0
  if(rv>rvmax) ick=1
end subroutine sph_countbb
subroutine sph_count2(rai,deci,ri,wi,raj,decj,rj,wj,wfib,iwfib, &
     nbs,seps2,nbp,sepp2,nbv,sepv,dds,ddpv,rvmax,ick,iok,rp,rv,wpp)
  real(kind=4) :: rai,deci,raj,decj
  real(kind=4),external :: wfib
  real(kind=8) :: dds(nbs),ddpv(nbp,nbv)
  real(kind=4) :: seps2(nbs+1),sepp2(nbp+1),sepv(nbv+1)
  real(kind=4) :: minus1,pi,shth2,rp2,r2
  minus1=-1
  pi=acos(minus1)
  pi=pi/180

  shth2=(sind(0.5*(deci-decj)))**2+ &
       cosd(deci)*cosd(decj)*(sind(0.5*(rai-raj)))**2
  rv=abs(ri-rj)
  rp2=4.*ri*rj*shth2
  r2=rv**2+rp2
  rp=sqrt(rp2)

  !--------------------
  ! Weight of this pair
  wpp=wi*wj*(iwfib*(wfib(shth2)-1)+1)

  iok=0
  if(r2>seps2(nbs+1)) goto 2
  do i=nbs,1,-1
     if(r2>seps2(i)) then
        dds(i)=dds(i)+wpp
        goto 2
     end if
  end do

2 if(rp2>sepp2(nbp+1).or.rv>sepv(nbv+1)) goto 3
  do i=nbp,1,-1
     if(rp2>sepp2(i)) then
        do j=nbv,1,-1
           if(rv>sepv(j)) then
              iok=1
              ddpv(i,j)=ddpv(i,j)+wpp
              goto 3
           endif
        enddo
        goto 3
     endif
  enddo

3 ick=0
  if(rv>rvmax) ick=1
end subroutine sph_count2
subroutine extrapolate_dr(nbin,sep,dr,seplim,sep1,sep2)
  real(kind=4) :: sep(nbin+1),s(nbin)
  real(kind=8) :: dr(nbin)

  k1=1
  k2=1
  do i=1,nbin
     s(i)=0.5*(sep(i)+sep(i+1))
     if(s(i)<sep1) k1=i+1
     if(s(i)<sep2) k2=i+1
  end do
  k1=min(k1,nbin)
  k2=max(k2,k1)

  x1=s(k1)
  x2=s(k2)
  y1=dr(k1)
  y2=dr(k2)
  yk=alog10(y2/y1)/alog10(x2/x1)
  do i=1,nbin
     if(s(i)>seplim) exit
     dr(i)=dr(k1)*(s(i)/x1)**yk
  end do

end subroutine extrapolate_dr
subroutine extrapolate_drpv(nbp,nbv,sepp,sepv,crpv,rplim,rp1,rp2)
  real(kind=4) :: sepp(nbp+1),sepv(nbv+1),rp(nbp)
  real(kind=8) :: crpv(nbp,nbv)

  k1=1
  k2=1
  do i=1,nbp
     rp(i)=0.5*(sepp(i)+sepp(i+1))
     if(rp(i)<rp1) k1=i+1
     if(rp(i)<rp2) k2=i+1
  end do
  k1=min(k1,nbp)
  k2=max(k2,k1)
  k=(k1+k2)/2

  do i=1,nbp
     if(rp(i)>rplim) exit
     do j=1,nbv
        crpv(i,j)=crpv(k,j)*(rp(i)/rp(k))**2
     end do
  end do

end subroutine extrapolate_drpv
