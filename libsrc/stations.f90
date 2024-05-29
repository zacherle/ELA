module stations

implicit none

   integer, parameter :: nStation=480

   integer            :: nstat = 0
   character(LEN=5), dimension(nStation) :: stat_name
   real, dimension(nStation) :: xstat
   real, dimension(nStation) :: ystat
   real, dimension(nStation) :: zstat
   real, dimension(nStation) :: dly

save

contains

subroutine read_stations(lu)

implicit none
integer         lu

integer        :: ios
character(256) :: iom
integer         i

!  number of stations
!   read (lu, * ,iostat=ios,err=50,end=50,iomsg=iom) nstat

!  do a test on allowed number of stations:
!  nstat requested; nStation dimensioned
!   if(nstat.gt.nStation)then
!      write (*,'(1x,a,": Too many stations. Max:",i3)') nStation
!      stop
!   endif

!  stations
   i=0
   do while (i .lt. 480)
      i = i+1
      read (lu, * ,iostat=ios,err=50,end=60,iomsg=iom)&
      &stat_name(i),xstat(i),ystat(i),zstat(i),dly(i)
!  z-coordinate in model file is upward --> set to downward
      zstat(i)=-zstat(i)
      !write(*,*) i, stat_name(i)
   end do
return

50 continue

if (ios/=0) then
    write(*,*) 'iostat = ', ios
    write(*,*) 'iomsg: '//trim(iom)
    stop
end if

60 continue
nstat=i
end subroutine read_stations

end module stations
