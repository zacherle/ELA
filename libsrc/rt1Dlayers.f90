
!**************************************************************************
! This subroutine is a modification of TRVDRV
! writen by J.P. Eaton (HYPOLAYR, 1969),
! which is part of  HYPO71, LEE AND LAHR (USGS OPEN-FILE REPORT 75-311, 1975).
!**************************************************************************

!subroutine raytr1D (s_x,s_z, nl,d,v, nray,ray_x,ray_z, tt,raypar,type_of_wave)
subroutine raytr1D (s_x,s_z, nl,d,v, tt,raypar,type_of_wave)

implicit none

real,               intent(IN)  :: s_x
real,               intent(IN)  :: s_z

! 1D model
integer           , intent(IN)  :: nl
real, dimension(:), intent(IN)  :: d
real, dimension(:), intent(IN)  :: v

!!ray
!integer,            intent(OUT) :: nray
!real, dimension(:), intent(OUT) :: ray_x    !2*m_layer
!real, dimension(:), intent(OUT) :: ray_z

real,               intent(OUT) :: tt
real(8),            intent(OUT) :: raypar
integer,            intent(OUT) :: type_of_wave

!all coordinates are aligned with the depth-layer system
!z-axis is downward
!-----------------------------------------------------------------
!     type_of_wave  =  1 ... refracted wave
!                      2 ... direct wave in first layer
!                      3 ... head wave along top of layer JL
!                      4 ... direct wave below first layer
!-----------------------------------------------------------------


!integer, parameter :: m_layer=100      ! m_layer=nl
!ray
integer               :: nray
real, dimension(2*nl) :: ray_x    !2*m_layer
real, dimension(2*nl) :: ray_z


   real xovmax
   real delxtr
   real dellit
   real delbig
   real xtest
   real tdir
   real tmin
   real tdc
   real tdj1
   real ub
   real ul
   real xlit
   real u
!
   integer ii,jj,l,k,ll,m
   integer j1
   integer jl
   real   tinj(nl+1),didj(nl+1),tr(nl+1)
   real, dimension(nl)   :: thk
   real, dimension(nl)   :: vsq
   integer i

   type_of_wave=0

   do i=1,nl-1
      thk(i) = d(i+1)-d(i)
   end do

   do i=1,nl
      vsq(i) = v(i)**2
   end do

   tmin=999.99

   if (nl.eq.1) then
      jl=nl
      go to 90
   endif

   jl=nl
   do l=1,nl
      if (d(l) .gt. s_z) then
         jj=l
         jl=l-1    !nearest higher interface
         exit      !break cycle
      endif
   end do
   ! s_z-d(jl) distance from focal to the nearest higher interface

   if (jl .eq. nl) then
!  only direct wave exist
      go to 100
   endif

!
!  test for refracted wave existence
!
   do l=jj,nl
      tinj(l)=hytid(jl,l)-(s_z-d(jl))*sqrt(vsq(l)-vsq(jl))/(v(l)*v(jl))
      didj(l)=hydid(jl,l)-(s_z-d(jl))*v(jl)/sqrt(vsq(l)-vsq(jl))
   end do
   xovmax=v(jj)*v(jl)*(tinj(jj)-hytid(jl,jl))/(v(jj)-v(jl))
   do m=jj,nl
      tr(m)=tinj(m)+s_x/v(m)
   end do
!
!  searching for the minimum arrival time for refracted wave
!
   tmin=999.99
   do m=jj,nl
      if (tr(m) .gt. tmin) cycle
      if (didj(m) .gt. s_x) cycle
      k=m
      tmin=tr(m)
   end do
!
!  on distance less then XOVMAX may be a direct wave arrival
!
   if (s_x .lt. xovmax) go to 90
!
!  travel time & derivatives for refracted wave
!
80 continue

!
!  set type of wave
!
   type_of_wave = 1
   !write(*,*) type_of_wave
!
!  no. of ray_xnts of intersection (of linear sections)
!
   nray = 2 * k - jj
!
!  set coordinates of linear sections; ray_x  ... profil coordinate
!                                      ray_z  ... z coordinate
!
   ray_x(1) = (thk(jl) - (s_z-d(jl)) ) * v(jl) / sqrt(vsq(k) - vsq(jl))
   ray_z(1) = d(jl+1)
   tt = sqrt(ray_x(1)**2+(ray_z(1)-s_z)**2)/v(jl)
!
!  travel from focal to the top of layer k
!
   do ii = 2,k-jj+1
      ray_x(ii) = ray_x(ii-1) + thk(jl+ii-1) * v(jl+ii-1) / sqrt(vsq(k) - vsq(jl+ii-1))
      ray_z(ii) = d(jl+ii)
      tt = tt + sqrt((ray_x(ii)-ray_x(ii-1))**2+(ray_z(ii)-ray_z(ii-1))**2)/v(jl+ii-1)
   end do
!
!  travel along the top of layer k
!
   ray_x(k-jj+2) = ray_x(k-jj+1) + s_x - didj(k)
   ray_z(k-jj+2) = d(k)
   tt = tt + (ray_x(k-jj+2)-ray_x(k-jj+1))/v(k)
!
!  from the top of layer k to the station
!
   do ii = k-jj+3,2*k-jj
      ray_x(ii) = ray_x(ii-1) + thk(2*k-jj+2-ii) * v(2*k-jj+2-ii) / sqrt(vsq(k)-vsq(2*k-jj+2-ii))
      ray_z(ii) = d(2*k-jj+2-ii)
      tt = tt + sqrt((ray_x(ii)-ray_x(ii-1))**2+(ray_z(ii)-ray_z(ii-1))**2)/v(2*k-jj+2-ii)
   end do

   tt = tt + sqrt((s_x-ray_x(2*k-jj))**2+(d(1)-ray_z(2*k-jj))**2)/v(1)
!write(*,*) k,tmin,tt
!   if (commut) then
!      toas=-v(1)/v(k)
!   else
!      toas=-v(jl)/v(k)
!   endif
   raypar=1/v(k)
   go to 260
!
!  calculation for direct wave
!
90 if (jl .ne. 1) go to 100
!
!  direct wave in first layer
!
   tdj1=sqrt((s_z-d(1))**2+s_x**2)/v(1)
   if (tdj1 .ge. tmin) then
!
!  refracted wave is faster
!
      go to 80
   endif
!
!  travel time & derivatives for direct wave in first layer
!

   type_of_wave = 2
   !write(*,*) type_of_wave
   nray=0

!   if (commut) then
!      toas=-s_x/sqrt((s_z-d(1))**2+s_x**2)
!   else
!      toas=s_x/sqrt((s_z-d(1))**2+s_x**2)
!   endif
   raypar=s_x/sqrt((s_z-d(1))**2+s_x**2)/v(1)

   tt = tdj1

   go to 260
!
!  find a direct wave that will emerge at the station
!
100 continue
   !sin alfa_max ... maximal angle of incidence
   ub=s_x/sqrt(s_x**2+(s_z-d(jl))**2)
   !sin alfa_min ... minimal angle
   xlit=s_x*(s_z-d(jl))/(s_z-d(1))
   ul=xlit/sqrt(xlit**2+(s_z-d(jl))**2)
   if (isnan(ul)) then
   ! top of layer
           ! ub = 1
           ! xlit = 0
           ! ul = NaN
           ul=0.0
           delbig = 0.0
           dellit = 0.0
           goto 105
   end if

!  distance of arrival of direct wave with angle of incidence
!  arcsin(ub) [resp. arcsin(ul)]

   delbig=(s_z-d(jl))*ub/sqrt(1-ub**2)
   dellit=(s_z-d(jl))*ul/sqrt(1-ul**2)

105 continue
   j1=jl-1
   do l=1,j1
      delbig=delbig+(thk(l)*ub)/sqrt(vsq(jl)/vsq(l)-ub**2)
      dellit=dellit+(thk(l)*ul)/sqrt(vsq(jl)/vsq(l)-ul**2)
   end do
!
!  delbig < s_x ... only direct wave along top of layer jl
!
   if (delbig.lt.s_x) then
      u=1.
      go to 190
   endif
!
!  iteration cyclus for searching appropriate angle of incidence
!  (wave hit station with selected error {+- 0.001 km})
!
   do ll=1,30
      if (delbig-dellit .lt. 0.001) go to 180
      u=(ub+ul)/2.
      delxtr=(s_z-d(jl))*u/sqrt(1-u**2)
      do  l=1,j1
         delxtr=delxtr+(thk(l)*u)/sqrt(vsq(jl)/vsq(l)-u**2)
      end do
      xtest=s_x-delxtr
      if (abs(xtest) .le. 0.001) go to 190
      if ( xtest < 0.0) then
!  s_x < delxtr
         ub=u
         delbig=delxtr
      else
!  s_x > delxtr
         ul=u
         dellit=delxtr
      endif
      if (ll.lt.10 .and. 1.0-u .lt. 0.00001) go to 190
   end do

180 u=(ub+ul)/2.
190 if (1.0-u .gt. 0.00001) go to 220
!
!  if u is too near 1, compute tdir as wave along the top of layer jl
!
    tdc=hytid(jl,jl)+s_x/v(jl)
    if (jl .eq. nl) go to 210       !focal in the lowest layer ...
    !only direct wave exist
    if (tdc .ge. tmin) go to 80     !refracted wave is faster
210 continue

!
!  along top of layer jl
!
   type_of_wave = 3
   !write(*,*) type_of_wave
   nray = jl - 1
   ray_x(1) = s_x - hydid(jl,jl)
!
   ray_z(1) = d(jl)

   tt = ray_x(1)/v(jl)
   do ii = jl-1,2,-1
      ray_x(jl+1-ii) = ray_x(jl-ii) + thk(ii) * v(ii) / sqrt(vsq(jl) - vsq(ii))
      ray_z(jl+1-ii) = d(ii)
      tt = tt + sqrt((ray_x(jl+1-ii)-ray_x(jl-ii))**2+(ray_z(jl+1-ii)-ray_z(jl-ii))**2)/v(ii)
   end do
   tt = tt + sqrt((s_x-ray_x(jl-1))**2+(d(1)-ray_z(jl-1))**2)/v(1)

!   if (commut) then
!      toas=-v(1)/v(jl)
!   else
!      toas=1D0
!   endif
   raypar=1/v(jl)

   go to 260
!
!  travel time & derivatives for direct wave below first layer
!

220 continue
   !travel time in layer concerning focal
   tdir=(s_z-d(jl))/(v(jl)*sqrt(1.0-u**2))
   do l=1,j1
      tdir=tdir+(thk(l)*v(jl))/(vsq(l)*sqrt(vsq(jl)/vsq(l)-u**2))
   end do
   if (jl .eq. nl) go to 245       !only direct wave exist
   if (tdir .ge. tmin) go to 80    !refracted wave is faster
245 continue

   type_of_wave = 4
   !write(*,*) type_of_wave
   nray = jl - 1
   ray_x(1) = (s_z-d(jl))*u / sqrt(1 - u**2)
   if (isnan(ray_x(1))) ray_x(1)=0.0
   ray_z(1) = d(jl)

   do ii = jl-1,2,-1
      ray_x(jl+1-ii) = ray_x(jl-ii) + (thk(ii) * u) / sqrt(vsq(jl) / vsq(ii) - u**2)
      ray_z(jl+1-ii) = d(ii)
   end do

!   if (commut) then
!      toas=-v(1)/v(jl)*u
!   else
!      toas=u
!   endif
   raypar=u/v(jl)

   tt = tdir

260 continue

return

contains

function hydid (j,ml)
   implicit none

   real hydid
   integer j,ml

   real flj
   real dm
   integer l

   hydid=0.
   if (ml.eq.1) return
   do l=1,ml-1
      dm=thk(l)*v(l)/sqrt(vsq(ml)-vsq(l))
      flj=1.
      if (l.ge.j) flj=2.
      hydid=hydid+flj*dm
   end do
end function hydid

function hytid (j,ml)
   implicit none

   real hytid
   integer j,ml

   real flj
   real tim
   integer l

   hytid=0.
   if (ml.eq.1) return
   do l=1,ml-1
      tim=thk(l)*sqrt(vsq(ml)-vsq(l))/(v(l)*v(ml))
      flj=1.
      if (l.ge.j) flj=2.
      hytid=hytid+flj*tim
   end do
end function hytid

end subroutine raytr1D

