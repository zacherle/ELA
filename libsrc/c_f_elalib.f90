module c_f_elalib
use :: iso_c_binding
use :: HypFile
use :: Arrivals

implicit none

public filehyp_read
public filesta_read
public filemod_read

type(cFileHyp),save     :: hyp
type(cFileArr),target,save     :: arrs

public c_f_str
contains

subroutine C_F_str(c_str,f_str, s_len)
   type(C_PTR), value,intent(in) :: c_str
!character(kind=c_char),dimension(*), intent(in) :: c_str
   character(*),intent(out) :: f_str
   integer,intent(in) :: s_len

   integer i
   character,dimension(:),pointer :: p_str
   call c_f_pointer(c_str,p_str,[1])
   do i=1,s_len
      if (p_str(i) .eq. c_null_char) exit
      f_str(i:i) = p_str(i)
   end do
end subroutine C_F_str


function filehyp_read(hypname) result(rmsg) bind(c,name='filehyp_read')
   type(C_PTR),value,intent(in) :: hypname
   integer(c_int) :: rmsg

   character(LEN=255) :: name_s=' '

   call C_F_str(hypname, name_s, 255)
   rmsg=hyp%load(name_s)
   call arrs%fill(hyp)
   parrs=>arrs

end function filehyp_read

function filesta_read(staname) result(rmsg) bind(c,name='filesta_read')
use :: stations
   type(C_PTR), value,intent(in) :: staname
   integer(c_int) :: rmsg     

   character(LEN=255) :: name_s=' '
   integer,parameter :: lu=12
   integer :: ios
   character (256) :: iom
   rmsg=1
   call C_F_str(staname, name_s, 255)
   open (file=name_s, unit=lu, iostat=ios, iomsg=iom, status='OLD')
   if (ios/=0) then
       write(*,*) 'iostat = ', ios
       write(*,*) 'iomsg: '//trim(iom)
       stop
   end if
   rmsg=0
   call read_stations(lu)

end function filesta_read

function filemod_read(modname) result(rmsg) bind(c,name='filemod_read')
use :: layers
   type(C_PTR), value,intent(in) :: modname
   integer(c_int) :: rmsg     

   character(LEN=255) :: name_s=' '
   integer :: lu
   integer :: ios
   character (256) :: iom

   type(cLayers),pointer :: lptr
   integer          :: iphase=0
   integer          :: i,j, nj

   rmsg=1
   call C_F_str(modname, name_s, 255)
   open (file=name_s, newunit=lu, iostat=ios, iomsg=iom, status='OLD')
   if (ios/=0) then
       write(*,*) 'iostat = ', ios
       write(*,*) 'iomsg: '//trim(iom)
       stop
   end if
   rmsg=0
   iphase=0
   do
      iphase=iphase+1
      lptr=>velmod(iphase)
      rmsg=lptr%load(lu)
      !lptr=cLayers(phase='')
      !rmsg=layers_load(lptr,lu)
      !velmod(iphase)=lptr
      write(*,*) rmsg
      if (rmsg .ne. 0) exit
   end do
   do i=1,iphase
      if (velmod(i)%n_layers .gt. 0) then
         nj = velmod(i)%n_layers     
         write(*,*) velmod(i)%phase
         !write(*,'(i3,1x,f6.3,1x,f6.3)') (j,velmod(i)%z(j), velmod(i)%v(j), j=1,nj)
      end if
   end do
   close (UNIT = lu)

end function filemod_read

function filehy3_write(hy3name) result(rmsg) bind(c,name='filehy3_write')
use output_list,        only: write_hy3
   type(C_PTR), value,intent(in) :: hy3name
   integer(c_int) :: rmsg     

   character(LEN=255) :: name_s=' '
   integer :: lu
   integer :: ios
   character (256) :: iom

   rmsg=1
   call C_F_str(hy3name, name_s, 255)
   open (file=name_s, newunit=lu, iostat=ios, iomsg=iom, status='UNKNOWN')
   if (ios/=0) then
       write(*,*) 'iostat = ', ios
       write(*,*) 'iomsg: '//trim(iom)
       stop
   end if
   !call write_hy3(lu,hypo,startpt,marr,arrs,hyp,modfn,covM,info)

   close (UNIT = lu,status='KEEP')

end function filehy3_write


end module c_f_elalib

!program testhyp
!
!use HypFile
!use Arrivals
!
!integer              :: rmsg
!type(cFileHyp)     :: hyp
!type(cFileArr)     :: arrs
!
!rmsg = hyp%load("test.hyp")
!rmsg = hyp%dump("testo.hyp")
!call arrs%fill(hyp)
!end program testhyp

