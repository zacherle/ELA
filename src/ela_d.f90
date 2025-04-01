
subroutine ela_d()

use calc_covar,   only: cov_matrix1
use arrivals,     only: parrs
use hypfile,      only: hyp
use fw_tdx,       only: get_res_w, get_G_w, get_topt
use output_list,  only: name_output_hy3, write_hy3

use, intrinsic :: iso_fortran_env, only : real64

implicit none

interface
   subroutine loc_hypo_lm(hypo, m, fvec, fjac, fix_depth, info)
use, intrinsic :: iso_fortran_env, only : real64
      real(real64), dimension(4), intent(INOUT) :: hypo
      integer, intent(IN)                       :: m
      real(real64), dimension(m), intent(OUT)   :: fvec(m)
      real(real64), dimension(m,4), intent(OUT) :: fjac(m, 4)
      logical, intent(IN)  :: fix_depth
      integer, intent(OUT) :: info
   end subroutine loc_hypo_lm
end interface

integer i,j

integer n0

character*255 line
character*4 answer4
logical init_nea

logical  :: fix_x, fix_y, fix_depth 

real(real64), dimension(4)  :: hypo
real, dimension(4)  :: startpt
real                :: x0,y0,z0,t0

integer                               :: marr, narr
real(real64), dimension(:), allocatable   :: fvec
real(real64), dimension(:), allocatable   :: w
real(real64), dimension(:,:), allocatable :: fjac
integer                               :: info
real(real64), dimension(4,4) :: covM

integer :: lu
integer :: ios
character (256) :: iom


   fix_x=.false.
   fix_y=.false.
   fix_depth=.false.
   startpt=0.0

40 continue
! start point
   write(*,'(1x,"Start point [X,Y,Z]",/&
   &     ,7x,"   orig. time given by minimizing procedure",/&
   &     ,7x,"   space coord. = 0 ... value of the nearest station",/&
   &     ,7x,"   [0,0,0,0]:_")')
   read (*,'(a)',end=40) line
   if (line.eq.' ') then
      init_nea = .true.
   else
      read (line,*,err=45,end=45) startpt
   endif
45 continue

   write(*,*) startpt,init_nea
50 continue
   write (*,'(1x,"Enter fixed coordinates of the epicenter:",/,9x,&
   &"Fixed coordinates XYZ   -  ''X or Y or XYZ''",/,9x,&
   &"Fixed depth             -  ''Z''   [ ]:_")')
!
   read (*,'(a)',err=50,end=50) answer4
   if (answer4.eq.' ') then
   else
!
      if ( index(answer4,'X').ne.0.or.index(answer4,'x').ne.0) then
         fix_x=.true.
         fix_y=.true.
         fix_depth=.true.
      end if
      if ( index(answer4,'Y').ne.0.or.index(answer4,'y').ne.0) then
         fix_x=.true.
         fix_y=.true.
         fix_depth=.true.
      end if
      if ( index(answer4,'Z').ne.0.or.index(answer4,'z').ne.0) then
         fix_depth=.true.
      end if
   end if
   write(*,*) fix_x,fix_y,fix_depth

marr=parrs%n_arr
allocate(fvec(marr))
allocate(w(marr))
allocate(fjac(marr,4))

if (fix_x.or.fix_y) then
! solve the forward problem
    hypo=startpt
    call get_topt(marr,4,hypo)
else
! not xyz fix
narr=0
do i=1,marr
   if (parrs%arr(i)%wt .gt. 0) then
           narr=narr+1
   end if
end do

if (narr .lt. 3) then
      write (*,'(1x,"# of arrivals in hyp_file  <  3")')
      return
end if

!  search the nearest station
n0=1
t0=1e20
do i=1,marr
    if (parrs%arr(i)%phase(1:1) .eq. 'S') cycle
    if (parrs%arr(i)%trec .gt. t0) cycle
    t0=parrs%arr(i)%trec
    n0=i
end do
!  coord of the nearest station
    x0=parrs%arr(n0)%X
    y0=parrs%arr(n0)%Y
    z0=parrs%arr(n0)%Z

   hypo=[x0+0.1, y0+0.1, 7.0, t0-0.5]
   hypo(3) = startpt(3)
   call loc_hypo_lm(hypo, marr, fvec, fjac, fix_depth, info)
end if ! not xyz fix

! covariance
call get_res_w(marr,4,hypo,fvec,w)
call get_G_w(marr,4,hypo,fjac,w)
fvec=fvec*w
do j=1,4
    fjac(:,j)=fjac(:,j)*w
end do

covM=cov_matrix1(marr,4,fvec,fjac,w)
do i=1,4
write (*,'(4(f16.6," "))') covM(i,:)
end do

if(len(trim(name_output_hy3)) > 3) then
   open (file=trim(name_output_hy3), newunit=lu, iostat=ios, iomsg=iom, status='UNKNOWN')
   if (ios/=0) then
       write(*,*) 'iostat = ', ios
       write(*,*) 'iomsg: '//trim(iom)
       stop
   end if
   call write_hy3(lu,real(hypo),startpt,marr,parrs,hyp,'model',covM,info)

   close (UNIT = lu,status='KEEP')
end if

end subroutine ela_d


