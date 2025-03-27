
function sthmax2(rpmax,rg,rkl)
  if(rkl==0)then
     sthmax2=1.
  else
     sthmax2=rpmax/sqrt(2.*rg*rkl)
     if(sthmax2>1) sthmax2=1.
  endif
end function sthmax2

real(kind=4) function wfib_dr72(shth2)
  real(kind=4) :: shth2,minus1,pi
  minus1=-1
  pi=acos(minus1)

  theta=asin(sqrt(shth2))*2*3600.*180./pi

  p9 = 89.697894

  if(theta<55) then
     p0 =  0.41748608
     p1 = -0.13445859
     p2 = -0.028616099
  else if(theta<p9) then
     p0 =  1.0618320
     p1 = -0.024856426
     p2 =  0.016569177
  else
     p0 =  1.0000000
     p1 = -0.065890517
     p2 = -0.0035219729
  end if

  wfib_dr72 = 1.0 / (p0 + p1 * exp(p2*theta))

end function wfib_dr72

real(kind=4) function wfib_dr4(shth2)
  real(kind=4) :: shth2,minus1,pi
  minus1=-1
  pi=acos(minus1)

  theta=sqrt(shth2)*2*3600.*180./pi

  wfib_sdss=1.0
  if(theta<55) then
     wfib_dr4= &
     (1+52.014168*theta**(-0.92056626))/(6.1709135*theta**(-0.48324151))
  else
     wfib_dr4= &
     (1+26.109340*theta**(-0.72920537))/(1+24.121585*theta**(-0.72637058))
  end if

end function wfib_dr4

real(kind=4) function wfib_tdf(shth2)
  real(kind=4) :: shth2

  wfib_tdf=1.0

end function wfib_tdf

subroutine sph_hocll(hoc,mxh1,mxh2,mxh3,nc1,nc2,nc3,ll,np,phi,th,r)
  real(kind=4) :: phi(np),th(np),phl,phu,thtl,thtu,hc1,hc2
  real(kind=4) :: r(np)
  integer hoc(mxh1,mxh2,mxh3),ll(np),q(3)

  if(nc1>mxh1.or.nc2>mxh2.or.nc3>mxh3) &
       stop 'nc1 or nc2 or nc3 too large in HOCLL_SPH !'

  !------------------------------------
  ! set the head of chain table to zero
  hoc=0

  !-------------------------------------
  ! set up the boundaries of the survey
  call sph_boundary(phl,phu,thtl,thtu,rl,ru)
  hc1=(thtu-thtl)/float(nc1)
  hc2=(phu-phl)/float(nc2)
  hc3=(ru-rl)/float(nc3)

  !------------------------------
  ! begin to fill the two tables 
  do i=1, np
     ! find the cell coordinate for the particle i
     q(1)=int((th(i)-thtl)/hc1)+1.
     q(2)=int((phi(i)-phl)/hc2)+1.
     q(3)=int((r(i)-rl)/hc3)+1.

     ! check the coordinates
     if(q(1)>nc1.or.q(1)<1) cycle
     if(q(2)>nc2) then
        q(2)=q(2)-nc2
     else if(q(2)<1) then
        q(2)=q(2)+nc2
     end if
     if(q(3)>nc3.or.q(3)<1) cycle

     ! fill the tables
     ll(i)=hoc(q(1),q(2),q(3))
     hoc(q(1),q(2),q(3))=i

     ! change the orders in increasing r(i) values;
     ! i.e. the head of the chain has the smallest r!
     if(ll(i)==0) cycle
     if(r(i)<=r(ll(i))) cycle
     hoc(q(1),q(2),q(3))=ll(i)
     i0=ll(i)
10   if(ll(i0)/=0)then
        if(r(i)<=r(ll(i0)))then
           ll(i)=ll(i0)
           ll(i0)=i
           cycle
        else
           i0=ll(i0)
           goto 10
        endif
     else
        ll(i0)=i
        ll(i)=0
     endif
  end do
end subroutine sph_hocll
function dalp(sthmax2,decg,i,hc1,decl,isele)
  real(kind=4) :: decg,hc1,decl
  real(kind=4) :: cdeccl,cdeccu,cdeccmin,atmp,btmp,sdalp2,sdalp22

  cdeccl=cosd(decl+(i-1)*hc1)
  cdeccu=cosd(decl+i*hc1)
  cdeccmin=min(cdeccu,cdeccl)
  
  if(cdeccmin==0) then
     dalp=180.
     return
  endif

  if(isele==1) then
     sdalp2=sthmax2/sqrt(cosd(decg)*cdeccmin)
  else
     atmp=sind((decg-(decl+(i-1)*hc1))*0.5)
     btmp=sind((decg-(decl+i*hc1))*0.5)
     atmp=atmp**2
     btmp=btmp**2
     atmp=min(atmp,btmp)
     sdalp22=(sthmax2**2-atmp)/(cosd(decg)*cdeccmin)
     if(sdalp22<0) sdalp22=0.
     sdalp2=sqrt(sdalp22)
  endif

  if(sdalp2>1)then
     dalp=180.
  else
     dalp=asind(sdalp2)*2.
  endif
end function dalp

subroutine systime(str)
  implicit none
  character(len=20) :: date,time
  character(len=40) :: str

  call date_and_time(date,time)
  str = 'Date:'//date(1:4)//'-'//date(5:6)//'-'//date(7:8) &
    //' Time:'//time(1:2)//':'//time(3:4)//':'//time(5:6)

  return
end subroutine systime

subroutine sph_boundary(ral,rau,decl,decu,dcoml,dcomu)
  real(kind=4) :: ral,rau,decl,decu,ramin,ramax,decmin,decmax
  common /bound/ramin,ramax,decmin,decmax,dcommin,dcommax
  ral=ramin
  rau=ramax
  decl=decmin
  decu=decmax
  dcoml=dcommin
  dcomu=dcommax
end subroutine sph_boundary

subroutine get_unit ( iunit )

!*******************************************************************************
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Modified:
!
!    18 September 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer IUNIT, the free unit number.
!
  implicit none

  integer i
  integer ios
  integer iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end
