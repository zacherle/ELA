module output_list

use iso_c_binding,      only: C_INT, C_FLOAT
use, intrinsic :: iso_fortran_env, only : real64

implicit none

contains

subroutine write_hy3(lu, hypo, startpoint, marr, arrs, hyp, modfn, co, info)

use layers,             only: model_error, reading_error
use const_raddeg,       only: RAD2DEG, DEG2RAD
use calc_covar,         only: get_errellipse
use version,            only: pversion
use arrivals,           only: cFileArr, cRecordArr
use hypfile,            only: cFileHyp, cRecordHyp
!use datetime_module     only: datetime
use gather_module,      only: collect_obs_simul

implicit none

INTERFACE
  real(C_FLOAT) function maxgap(naz, az) BIND(C)
    USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_INT,C_FLOAT
    IMPLICIT NONE
    integer(C_INT),value :: naz
    real(kind=C_FLOAT),dimension(naz),intent(IN)  :: az
  end function maxgap

  function localmagnitude(hypo,r)
    use arrivals,             only    : cRecordArr
    implicit none
    real, dimension(3), intent(IN)    :: hypo
    type(cRecordArr), intent(IN)      :: r
    real                              :: localmagnitude
  end function localmagnitude
END INTERFACE

integer lu
real,dimension(4),intent(IN) :: hypo
real,dimension(4),intent(IN) :: startpoint
integer, intent(IN)          :: marr
type(cFileArr),intent(IN)    :: arrs
type(cFileHyp),intent(IN)    :: hyp
character(LEN=*),intent(IN)  :: modfn
real(real64),dimension(4,4),intent(IN) :: co
integer, intent(IN)          :: info

   type(cRecordArr)   :: a
   type(cRecordHyp)   :: r

   character*22  whole_date
   integer ita(9)
   integer stime
   integer i,j
   integer isec,msec

   real    dxer, dyer, dzer, dter, l1, l2, theta
   logical, dimension(4)  :: ee_nan=.false.
   real    az_theta
   real    rmsres

   real(kind=C_FLOAT),dimension(marr)  :: az
   real    gap

   real(real64) :: fi, rla
   real    meridian_con

   integer year_ref,month_ref,day_ref,hour_ref,minute_ref
   real t_ref
   integer year_orig, month_orig, day_orig, hour_orig, minute_orig

real(real64),dimension(7,marr) :: g
integer                        :: ngather
real :: aobs, acal, ares, depi, dhypo, aaz, amag, dx,dy,dz, atoa

integer ::  nmag
real    ::  smag, smag2

   call arrs%todatetime(hypo(4),year_orig,month_orig,day_orig,hour_orig,minute_orig,isec,msec)
   year_orig = mod(year_orig,100)
   write (whole_date,&
   &'(2(i2.2,"-"),i2.2,2x,2(i2.2,":"),i2.2,".",i3.3)')&
   &year_orig,month_orig,day_orig,hour_orig,minute_orig,isec,msec

   t_ref = 0.0
   call arrs%todatetime(t_ref,year_ref,month_ref,day_ref,hour_ref,minute_ref,isec,msec)
   year_ref = mod(year_ref,100)
! --------------------------------------------------------------------

   write (lu,'("program       :ELA, ",a)') pversion
   write (lu,'("model         :",a)') modfn(1:lnblnk(modfn))       ! model name
   write (lu,'("model error   :",f5.3," s")') model_error
   write (lu,'("reading error :",f5.3," s")') reading_error

   stime = time ()  !system time
   call ltime (stime,ita)
   write (lu,'("create time   :",i2.2,"-",i2.2,"-",i2.2,1x,&
              &i2.2,":",i2.2,":",i2.2)')&
              & mod(ita(6),100),ita(5)+1,ita(4),ita(3),ita(2),ita(1)
   write (lu,'("event",9x,":",a)') hyp%filename(1:lnblnk(hyp%filename))

   write (lu,'("start(x,y,z,t):(",&
   &f7.2,",",f7.2,",",f7.2,",",f7.2,")")') startpoint

   write (lu,'("reference time:",&
   &i2.2,"-",i2.2,"-",i2.2,1x,i2.2,":",i2.2)')&
   &year_ref,month_ref,day_ref,hour_ref,minute_ref

! header for station data
   write (lu,&
   &'("----------------------------------------------------------",&
   &  "-------------")')
   write (lu,'(" sta   |obs. t.|cal. t.|res. |amplitude"&
   &,"|freq|w| epi |hypo |azm|ain|xmag")')
   write(lu, '("       |  [s]  |  [s]  | [s] |  [m/s]  "&
   &,"|[Hz]| |[km] |[km] |[o]|[o]|    ")')
   write (lu,&
   &'("----------------------------------------------------------",&
   &  "-------------")')

   meridian_con = mconvergence(hypo(1),hypo(2))
! station data
   call collect_obs_simul(hypo,arrs%n_arr,g,ngather)
   nmag=0
   smag=0.0
   smag2=0.0
   j=0
   do i=1,arrs%n_arr
     a = arrs%arr(i)
     r = hyp%rec(a%id)
     aobs=a%trec   !g(5,i)
     acal=g(1,i)+hypo(4)
     ares=aobs-acal
     dx=hypo(1)-a%X
     dy=hypo(2)-a%Y
     dz=hypo(3)-a%Z
     dhypo = sqrt(dx**2+dy**2+dz**2)
     depi  = sqrt(dx**2+dy**2)
     aaz = real(mod(720.0+atan2(dy,dx)*RAD2DEG-180.0-meridian_con     ,360.0))
     if (a%wt > 0.0) then
         j=j+1
         az(j)=aaz
     end if
     amag=localmagnitude(hypo(1:3),a)
     if (amag > -9.9) then
         nmag=nmag+1
         smag=smag+amag
         smag2=smag2+amag**2
     end if
     atoa=g(7,i)

     write(lu,100, advance='no') r%sta, r%phase, aobs, acal, ares,&
             &a%amp,a%freq,nint(r%wt),depi,dhypo,nint(aaz),nint(atoa)
     if (a%phase(1:1).eq.'S') then
         write(lu,'(f4.1)') amag
     else
         write(lu,*) 
     end if
   end do
   gap=maxgap(j,az)
   !sum((res*w)**2)/sum(w)
   rmsres=sqrt(sum(((g(5,:)-g(1,:)-hypo(4))*g(6,:))**2)/sum(g(6,:)))

100   format (a5,' ',a1,'|',&
          &f7.2,   '|',f7.2,   '|',f5.3,'|',1pe9.2,'|',0pf4.1,&
          &'|',i1,'|',0pf5.1,'|',0pf5.1,'|',i3,'|',i3,'|')

! error ellipse
  call get_errellipse(co,dxer,dyer,dzer,dter,l1,l2,theta)
  ee_nan=.true.
  if (dxer>0 .and. dxer<100) ee_nan(1)=.false. 
  if (dyer>0 .and. dyer<100) ee_nan(2)=.false. 
  if (dzer>0 .and. dzer<50) ee_nan(3)=.false. 
  if (dter>0 .and. dter<5) ee_nan(4)=.false. 

! coordinates local --> Krovak
!   theta = mod(360.0+theta+p_fi,360.0)
! coordinates Krovak --> geofraphic
!   az_theta = theta-nangle
   az_theta = theta-180.0-meridian_con
   az_theta = mod(360.0+az_theta,360.0)

   call XY2FL (hypo(2)*1000, hypo(1)*1000, fi, rla)

   write(lu,'(//,"hypocenter data:",/,"----------------")')

   write(lu,'("origin time          t:",2x,a22,$)') whole_date
   if (ee_nan(4)) then
      write(lu,'(1x,"+-    NaN")')
   else
      write(lu,'(1x,"+-",1x,f6.3)') dter
   endif

   write(lu,'("x-coordinate         x:",2x,f7.2,$)') hypo(1)
   if (ee_nan(1)) then
      write(lu,'(1x,"+-    NaN   km",$)')
   else
      write(lu,'(1x,"+-",1x,f6.2,3x,"km",$)') dxer
   endif
   write(lu,'(7x,"(fi:",f10.6," deg)")') fi

   write(lu,'("y-coordinate         y:",2x,f7.2,$)') hypo(2)
   if (ee_nan(2)) then
      write(lu,'(1x,"+-    NaN   km",$)')
   else
      write(lu,'(1x,"+-",1x,f6.2,3x,"km",$)') dyer
   endif
   write(lu,'(3x,"(lambda:",f10.6," deg)")') rla

   write(lu,'("depth                z:",2x,f7.2,$)') hypo(3)
   if (ee_nan(3)) then
      write(lu,'(1x,"+-    NaN   km")')
   else
      write(lu,'(1x,"+-",1x,f6.2,3x,"km")') dzer
   endif

   write(lu,'("magnitude           ml:",2x,f7.2,1x,"+-",1x,f6.2)') smag/nmag, sqrt((smag2-smag**2/nmag)/nmag)
   write(lu,'("rms of time residuals :",8x,f6.2,8x,"s")') rmsres
   write(lu,'("angular gap           :",11x,i3,8x,"deg")') nint(gap)
   write(lu,'("info                  :",11x,i3)') info

   if(ee_nan(1) .or. ee_nan(2)) then
      write(lu,'("error ellipse axis l1 :",10x," NaN",8x,"km")')
      write(lu,'("              axis l2 :",10x," NaN",8x,"km")')
   else
      write(lu,'("error ellipse axis l1 :",8x,f6.2,8x,"km")') l1
      write(lu,'("              axis l2 :",8x,f6.2,8x,"km")') l2
   endif

   if(ee_nan(1) .or. ee_nan(2)) then
      write(lu,'("              theta   :",2x,"NaN deg (to grid)",$)')
      write(lu,'(4x,"(azimuth:",2x,"NaN deg)")')
   else
      write(lu,'("              theta   :",2x,f6.1," deg (to grid)",$)') theta
      write(lu,'(4x,"(azimuth:",f6.1," deg)")') az_theta
   endif

   return

contains

real function mconvergence(X, Y)

    real X,Y   ! X,Y [km]

    mconvergence = 0.008257*Y+2.373*Y/X

end function mconvergence


end subroutine write_hy3

end module output_list
