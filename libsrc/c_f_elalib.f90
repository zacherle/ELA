module c_f_elalib
use :: iso_c_binding

implicit none

public filehyp_read
public filesta_read
public filemod_read

!type(cFileHyp),save     :: hyp
!type(cFileArr),target,save     :: arrs

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
use :: HypFile,  only: cFileHyp, hyp
use :: Arrivals, only: cFileArr, parrs
   type(C_PTR),value,intent(in) :: hypname
   integer(c_int) :: rmsg

   character(LEN=255) :: name_s=' '

   call C_F_str(hypname, name_s, 255)
   rmsg=hyp%load(name_s)
   call parrs%fill(hyp)
   !parrs=>arrs

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

subroutine setpar_name_o(hy3name) bind(c,name='setpar_name_o')
use output_list,        only: name_output_hy3
   type(C_PTR), value,intent(in) :: hy3name

   call C_F_str(hy3name, name_output_hy3, 255)

end subroutine setpar_name_o

subroutine setpar_model(model_err, reading_err, modname) bind(c,name='setpar_model')
use layers,  only: model_error, reading_error, name_model
   real(c_double),value,intent(in) :: model_err, reading_err
   type(C_PTR), value,intent(in) :: modname

   model_error = real(model_err)
   reading_error = real(reading_err)
   call C_F_str(modname, name_model, 255)

end subroutine setpar_model

end module c_f_elalib
