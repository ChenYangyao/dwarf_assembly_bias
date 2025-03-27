      function cputime()
      real tarray(2)

      cputime=etime(tarray)
      return
      end

      subroutine binsize(nbin,sep,sepmin,dsep,logsep)
      real sep(nbin+1)
      do i=1,nbin+1
         if(logsep == 1) then
            sep(i)=sepmin*10**((i-1)*dsep)
         else
            sep(i)=sepmin+(i-1)*dsep
         endif
      enddo
      end

! Compute comoving line-of-sight distances (for c/H_0=1).
! z---------redshit
! omegam----Omega_M at z=0
! omegal----Omega_Lambda at z=0
      function comdis(z,omegam,omegal)
      external cdfunc
      common /fluc/omegam0,omegal0

      omegam0=omegam
      omegal0=omegal

      ul=1.
      dl=1./sqrt(1+z)
      call qromb(cdfunc,dl,ul,comdis)

      comdis=2*comdis
      return
      end

      !-----------------
      function cdfunc(y)
      common /fluc/omegam,omegal

      cdfunc=1./sqrt(omegam+omegal*y**6+(1.-omegam-omegal)*y**2)
      return
      end

      subroutine wboots90(npt,nsamp,w)
      real(kind=4) :: w(npt,nsamp),sd

      do isamp=1,nsamp
        do ipt=1,npt
          w(ipt,isamp)=0
        enddo
      enddo

      do isamp=1,nsamp
         do ipt=1,npt
            call random_number(sd)
            isel=floor(sd*(npt-1))+1
            w(isel,isamp)=w(isel,isamp)+1
         enddo
      enddo

      end

!------------------------------------------------------------------------
!     Compute the sum of an array
!     npt-------number of data points
!     x---------1D array
      real function total(npt,x)
      integer npt,i
      real x(npt),y

      y = 0
      do i=1,npt
         y = y+x(i)
      enddo

      total = y
      return
      end

!------------------------------------------------------------------------
!     Compute the mean of an array
!     npt-------number of data points
!     x---------1D array
      real function avg(npt,x)
      implicit none
      integer npt
      real x(npt),y,total
      external total

      y = total(npt,x)/npt

      avg = y
      return
      end

!------------------------------------------------------------------------
!     Compute the standard deviation of an array
!     npt--------number of data points
!     x----------1D array
      real function stddev(npt,x)
      implicit none
      integer npt,i
      real x(npt),y,xmean,avg,x2(npt),total,sum
      external avg,total
      
      if (npt.le.1) then
         y = 0
      else
         xmean = avg(npt,x)
         do i=1,npt
            x2(i) = (x(i)-xmean)**2
         end do
         sum = total(npt,x2)
         if(sum.le.0)then
            y = 0
         else
            y = sqrt(sum/(npt-1))
         end if
      end if

      stddev = y
      return
      end
!-----
!   PI
      function pi()
      !pi = acos(-1.)
      pi = 3.1415926535
      end

!-----------------------------------
!   Same as the intrinsic function cos() but the input angle should be 
!   in unit of degree
      function cosd(deg)
      cosd = cos(deg * pi() / 180.)
      end

!-----------------------------------
!   Same as the intrinsic function sin() but the input angle should be 
!   in unit of degree
      function sind(deg)
      sind = sin(deg * pi() / 180.)
      end

!-----------------------------------
!   Same as the intrinsic function acos() but the output angle is 
!   in unit of degree
      function acosd(x)
      acosd = acos(x) * 180. / pi()
      end

!-----------------------------------
!   Same as the intrinsic function asin() but the output angle is
!   in unit of degree
      function asind(x)
      asind = asin(x) * 180. / pi()
      end
