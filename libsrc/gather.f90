module gather_module

use, intrinsic :: iso_fortran_env, only : real64
implicit none

contains

subroutine collect_obs_simul(srcxyz,m,gather,ngather)
!observed arrivals, weights and
!calculated travel times, dt/dx for all stations

use arrivals,           only: parrs  !, cFileArr, cRecordArr
use twopoint

implicit none

real,dimension(3), intent(IN)       :: srcxyz
integer, intent(IN)                 :: m
!gather(1,:) calc tt
!gather(2,:) calc dt/dx
!gather(3,:) calc dt/dy
!gather(4,:) calc dt/dz
!gather(5,:) obs arrival
!gather(6,:) obs weight
real(real64),dimension(7,m),intent(OUT) :: gather
integer, intent(OUT)                    :: ngather

integer           :: i,j,narr
real,dimension(3) :: obsxyz
character(LEN=5)  :: ph
real,dimension(4) :: td
real              :: toa

narr=parrs%n_arr
if (narr .gt. m) then
    narr=m
end if
j=0
do i=1,narr
    obsxyz(1)=parrs%arr(i)%X
    obsxyz(2)=parrs%arr(i)%Y
    obsxyz(3)=parrs%arr(i)%Z
    ph=trim(parrs%arr(i)%phase)
    call two_point_td(ph,srcxyz,obsxyz,td,toa)
    !if (td(1) .le. 0.0) cycle
    j=j+1
    !write(*,*) ph,srcxyz,obsxyz,td,toa
    !tcal(i)=td(1)+dly(key(i))  !  travel time
    gather(1,j)=td(1)              !  travel time
    gather(2,j)=td(2)              !  dt/dx
    gather(3,j)=td(3)              !  dt/dy
    gather(4,j)=td(4)              !  dt/dz
    !  dt/dt_o = 1
    if (td(1) .gt. 0.0) then
        gather(5,j)=parrs%arr(i)%trec
        gather(6,j)=parrs%arr(i)%wt
        gather(7,j)=toa
    else
        gather(5,j)=parrs%arr(i)%trec
        gather(6,j)=0.0
        gather(7,j)=toa
    end if
end do
ngather=j
end subroutine collect_obs_simul

end module gather_module
