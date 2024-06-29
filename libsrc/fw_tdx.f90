module fw_tdx

!use povrch,             only: relief8
use gather_module,      only: collect_obs_simul

implicit none

private
public Xp, Xo, tdx3, tdxder, get_res_w, get_G_w

integer, parameter :: dp=kind(0d0)

contains

function Xo(x,zmin)
! transform from bounded parameter space Xp to optimization procedure space Xo
    real(dp), dimension(4)             :: Xo
    real(dp), dimension(4), intent(IN) :: x
    real(dp), intent(IN) :: zmin
    real(dp) o
! reflective transformation zp -> zo
    o=sqrt((x(3)-zmin+1)**2-1)
    !o=abs(x(3))
    !o=o+10.0
    Xo=[x(1)-1100.0,x(2)-700.0,o,x(4)]
end function Xo

function Xp(x,zmin)
! transform from optimization procedure space Xo to bounded parameter space Xp
    real(dp), dimension(4)             :: Xp
    real(dp), dimension(4), intent(IN) :: x
    real(dp), intent(IN) :: zmin
    real(dp) p
! reflective transformation zo -> zp
    p=x(3)
    !p=p-10.0
    p=zmin-1+sqrt(p**2+1)
    !p=abs(x(3))
    Xp=[x(1)+1100,x(2)+700.0,p,x(4)]
end function Xp

function dpdo(z_o)
! the derivative of the z-coordinate transformation 
! df(y)/dx = df(y)/dy * dy(x)/dx
! x = z_o, y = z_p(z_o) = sqrt(z_o^2+1)-1
    real(dp) dpdo
    real(dp) z_o
    real(dp) r
    r=z_o/sqrt(1+z_o**2)
    ! prevention of zero
    r=r+sign(0.001_dp,z_o)
    dpdo=r
end function dpdo

      
subroutine tdxder(m, n, x, fvec, fjac, ldfjac, iflag)

implicit none
integer, intent(in) :: m, n, ldfjac, iflag
real(dp), dimension(n), intent(in) :: x
real(dp), dimension(m), intent(out) :: fvec
real(dp), dimension(ldfjac,n), intent(out) :: fjac

real(dp), dimension(4) :: xp_
real(dp), dimension(m) :: w

   !transform from optimization x to bounded parametric xp space
!   xp_=Xp(x,relief8(x(1),x(2)))
   xp_=Xp(x,0d0)
! write(*,*) "A: ", xp_
!   call get_topt(m,n,xp_)
! write(*,*) "B: ", xp_

   if (iflag .eq. 1) then
      call get_res_w(m,n,xp_,fvec,w)
      fvec=fvec*w

   else if (iflag .eq. 2) then
      call get_G_w(m,n,xp_,fjac,w)

      ! transform dt/dz from bounded xp -> optimization xo space
      ! dpdo = dz_p(z_o)/dz_o 
      fjac(:,3)=fjac(:,3)*dpdo(x(3))

      fjac(:,1)=fjac(:,1)*w
      fjac(:,2)=fjac(:,2)*w
      fjac(:,3)=fjac(:,3)*w
      fjac(:,4)=fjac(:,4)*w
!      fjac(:,4)=0.0_dp

   else
      write(*,*) xp_

   end if
end subroutine tdxder

subroutine tdx3(m, n, x, fvec, fjac, ldfjac, iflag)

implicit none
integer, intent(in) :: m, n, ldfjac, iflag
real(dp), dimension(n), intent(in) :: x
real(dp), dimension(m), intent(out) :: fvec
real(dp), dimension(ldfjac,n), intent(out) :: fjac

real(dp), dimension(m) :: w

   if (iflag .eq. 1) then
      call get_res_w(m,n,x,fvec,w)
      fvec=fvec*w

   else if (iflag .eq. 2) then
      call get_G_w(m,n,x,fjac,w)

      fjac(:,1)=fjac(:,1)*w
      fjac(:,2)=fjac(:,2)*w
      fjac(:,3)= 0d0
      fjac(:,4)=fjac(:,4)*w

   else
      write(*,*) x

   end if

end subroutine tdx3

subroutine get_topt(m, n, x)

implicit none
integer, intent(in) :: m, n
real(dp), dimension(n), intent(inout) :: x

real(dp),dimension(7,m) :: g
integer                 :: ngather

   call collect_obs_simul(real(x),m,g,ngather)
   if (ngather < m) then
       write(*,*) 'ngather<m ', ngather,m
       ngather=m
   end if

   x(4)=sum((g(5,:)-g(1,:))*g(6,:))/sum(g(6,:))
end subroutine get_topt

subroutine get_res_w(m, n, x, fvec, w)

implicit none
integer, intent(in) :: m, n
real(dp), dimension(n), intent(in) :: x
real(dp), dimension(m), intent(out) :: fvec
real(dp), dimension(m), intent(out) :: w

real(dp),dimension(7,m) :: gather
integer                 :: ngather

   call collect_obs_simul(real(x),m,gather,ngather)
   if (ngather < m) then
       write(*,*) 'ngather<m ', ngather,m
       ngather=m
   end if

   fvec=0d0
   w=0d0
   fvec(1:ngather)= gather(5,1:ngather)-x(4)-gather(1,1:ngather)
   w(1:ngather)= gather(6,1:ngather)
end subroutine get_res_w


subroutine get_G_w(m, n, x, fjac,  w)

implicit none
integer, intent(in) :: m, n
real(dp), dimension(n), intent(in) :: x
real(dp), dimension(m,n), intent(out) :: fjac
real(dp), dimension(m), intent(out) :: w

real(dp),dimension(7,m) :: gather
integer                 :: ngather

   call collect_obs_simul(real(x),m,gather,ngather)
   if (ngather < m) then
       write(*,*) 'ngather<m ', ngather,m
       ngather=m
   end if

   fjac=0d0
   w=0d0
   fjac(1:ngather,1)= -gather(2,1:ngather)
   fjac(1:ngather,2)= -gather(3,1:ngather)
   fjac(1:ngather,3)= -gather(4,1:ngather)
   fjac(1:ngather,4)= -1d0
   w(1:ngather)= gather(6,1:ngather)
end subroutine get_G_w

end module fw_tdx
