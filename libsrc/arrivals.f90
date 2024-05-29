module Arrivals

use datetime_module,   only: datetime, timedelta
use hypfile

implicit none
private

type, public     :: cRecordArr
   integer          :: id = 0
   character(LEN=5) :: sta = ''
   character(LEN=5) :: phase = ''
   real             :: trec = 0.0
   real             :: wt = 0.0
   real             :: amp = 0.0
   real             :: freq = 0.0
   integer          :: key = 0
   real             :: X
   real             :: Y
   real             :: Z
contains
!procedure ::
end type cRecordArr


type, public      :: cFileArr
   integer           :: n_arr = 0
   type (datetime)   :: reftime
   real              :: sumwt = 0.0
   real              :: sumwt2 = 0.0
   real              :: maxwt = 0.0
   type(cRecordArr), allocatable,  dimension(:)  :: arr
contains
procedure :: fill => cfilearr_fill
procedure :: toDateTime => cfilearr_todatetime
end type cFileArr

type(cFileArr),pointer,public,save   :: parrs

contains

subroutine cfilearr_todatetime(this,rtime,year,month,day,hour,minute,isec,msec)
implicit none
class(cFileArr)      :: this

real, intent(IN)     :: rtime
integer, intent(OUT) :: year,month,day,hour,minute,isec,msec

type(datetime)       :: pcas
integer              :: rtime_m

    rtime_m = nint(rtime*1000)
    write(*,*) rtime_m
    pcas = this%reftime + timedelta(0,0,0,0,rtime_m)
    write(*,*) pcas%isoformat()
    year =  pcas%getYear()
    month = pcas%getMonth()
    day   = pcas%getDay()
    hour  = pcas%getHour()
    minute= pcas%getMinute()
    isec  = pcas%getSecond()
    msec  = pcas%getMillisecond()

end subroutine cfilearr_todatetime

!
subroutine cfilearr_fill(this, hyp)
!
use stations,   only: nstat, stat_name, xstat, ystat, zstat
!
   implicit none

   class(cFileArr)      :: this
   type (cFileHyp)      :: hyp

!
   integer i
   integer j
!
   real    pwt
   real    psum, psum2, mwt
   type(datetime)             :: pcas, mcas
   type(timedelta)            :: dcas
   type(cRecordHyp)           :: r
   type(cRecordArr)           :: a
! reftime
   mcas = datetime(3000)
   do i=1,hyp%n_rec
     r = hyp%rec(i)
     pcas = datetime(r%getYear(),r%month,r%day,&
        &r%hour,r%minute,r%isec,r%msec)
     if (pcas .lt. mcas) mcas=pcas
   end do
!   if(mcas%getSecond() .lt. 30) then
!      mcas = mcas - timedelta(0,0,1)
!   endif
   this%reftime = datetime(mcas%getYear(),mcas%getMonth(),mcas%getDay(),mcas%getHour(),mcas%getMinute())
   write(*,*) this%reftime%isoformat()

if (allocated(this%arr)) deallocate(this%arr)
allocate(this%arr(hyp%n_rec))

   this%n_arr=0
   do i=1,hyp%n_rec
     r = hyp%rec(i)
!     if (r%wt .eq. 4) cycle
     a%id = r%id
     a%sta = r%sta
     a%phase = r%phase
     a%amp = r%amp_v
     if (r%period > 0) then
       a%freq = 1.0/r%period
     else
       a%freq = 99.9
     end if

! compute weights
     pwt = r%wt
! weights refer (are inverse proportional) to the standard deviation
     if (hyp%hyr) then
        if (pwt .gt. 0) then
             pwt = pwt/1000.0
             pwt = 1. / pwt
        else if (pwt .lt. 0) then
             pwt = 0.0
        else
             pwt = 0.0
        endif
    else
             pwt=(4.-pwt)/4.
    endif
    a%wt = pwt

     pcas = datetime(r%getYear(),r%month,r%day,&
        &r%hour,r%minute,r%isec,r%msec)
     dcas = pcas - this%reftime
     a%trec = real(dcas%total_seconds())
! key
     a%key=0
     do j=1,nstat
         if (a%sta .eq. stat_name(j)) then
             a%key=j
             a%X=xstat(j)
             a%Y=ystat(j)
             a%Z=zstat(j)
         endif
     end do
     if (a%key .eq. 0) then
         write (*,'(1x,"Station ",a5," not found.")')&
               &   this%arr(i)%sta
         cycle
     end if
     this%arr(i)=a
     this%n_arr=this%n_arr+1
   end do   !i n_rec

!
!   sum weights
!
   psum=0.0
   psum2=0.0
   mwt = 0.0
   do  i=1,this%n_arr
      pwt=this%arr(i)%wt
      psum=psum+pwt
      psum2=psum2+pwt*pwt
      if (pwt > mwt) mwt = pwt
   end do  !i n_arr
   this%sumwt=psum
   this%sumwt2=psum2
   this%maxwt=mwt
!
return

end subroutine cfilearr_fill

!      if (type(i) .eq. 'S') then
!      else if (type(i) .eq. 'P') then
!      else
!         write (*,'(1x," : Wrong type of arrival in record ",i2)') i
!      endif
end module Arrivals
