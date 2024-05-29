module twopoint

implicit none
private get1Dvel
public  two_point_td

contains

function get1Dvel(depth,nl,d,v)

implicit none
real                           :: get1Dvel
real,               intent(IN) :: depth
integer,            intent(IN):: nl
real, dimension(:), intent(IN):: d
real, dimension(:), intent(IN):: v

integer    i
real       vel
real       vair                        !velocity of sound wave in air
parameter (vair = 0.300)

if (depth .lt. d(1)) then
    get1Dvel=vair
    return
end if

vel=vair
do i=nl,1,-1
    if (depth .ge. d(i)) exit  !break cycle
end do
if (i .gt. 0) then
    vel=v(i)
else
    vel=vair
end if

get1Dvel=vel

end function get1Dvel

subroutine two_point_td(ph,srcxyz,obsxyz,td,toa)
! Compute travel time, derivations and take off angle
! z-axes downward
use const_raddeg
use layers
implicit none

real,dimension(3), intent(IN)  :: srcxyz    !coordinates of the hypocenter
real,dimension(3), intent(IN)  :: obsxyz    !coordinates of the station
character(*), intent(IN)   :: ph        !seismic phase
!      td(1) travel time
!      td(2:4) derivation on the x,y,z in the hypocenter
real,dimension(4), intent(OUT) :: td
! toa - ray take off angle of the hypocenter relative to the z-axis
real, intent(OUT)              :: toa

interface
include 'rt1Dlayers.fh'
end interface

   type(cLayers), pointer  :: lptr=>null()
   integer :: nl        !number of layers
   real(8) toas
   real(8) toac
   real(8) raypar
   integer type_of_wave
   real    v_hypo       !velocity in hypoc.
   real    delta        !distance
   real    tt           !travel time 1D

   lptr=>get_layers(ph)
   if (associated(lptr) .eqv. .false.) then
      td(1)=-1.0
      td(2:4)=[0.0,0.0,0.0]
      return
   end if
   nl=lptr%n_layers
   delta = sqrt((obsxyz(1)-srcxyz(1))**2+(obsxyz(2)-srcxyz(2))**2)
   if (delta .lt. 0.0001) then
      delta = 0.0001
   endif

!  ray traceing in 1D
   call raytr1D (delta,srcxyz(3), nl,lptr%z,lptr%v, tt, raypar,type_of_wave)

!  Get the velocity at the hypocenter to calculate the derivative dt/dz
   v_hypo = get1Dvel(srcxyz(3),nl,lptr%z,lptr%v)

!  travel time
   td(1) = tt

!  derivation by x  of travel time
   td(2)=real(-raypar*(obsxyz(1)-srcxyz(1))/delta)

!  derivation by y
   td(3)=real(-raypar*(obsxyz(2)-srcxyz(2))/delta)

!  derivation by z
   td(4)=0.0
   toas=raypar*v_hypo
   toac = dsqrt(1D0-toas*toas)
   if (isnan(toac)) then
      toas=1.D0
      toac=0.D0
   end if
   if (type_of_wave .eq. 1) then
      toac=-toac
   end if
   td(4)=real(toac/v_hypo)

!  the directional cosine of the ray to the z-axis
   toa=real(acos(-toac)*RAD2DEG)
   return
end subroutine two_point_td

end module twopoint
