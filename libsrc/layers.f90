module layers

implicit none
public

real         :: model_error=0.0
real         :: reading_error=0.016

type         :: cLayers
   character(LEN=5) :: phase=' '
   integer    :: n_layers=0
   real, allocatable, dimension(:) :: z
   real, allocatable, dimension(:) :: v
contains
procedure :: load => layers_load
end type cLayers

type(cLayers),dimension(5),target  :: velmod

save

public get_layers

contains

function layers_load(lay,iu) result(rmsg)

class(cLayers), intent(inout) :: lay
integer, intent(in)           :: iu        
integer                       :: rmsg

integer              :: ios
character(len=255)   :: iom
character(LEN=40)    :: cbuf
!character,dimension(40) :: cbuf
integer              :: icbuf
integer              :: i
real, dimension(500) :: z_tmp
real, dimension(500) :: v_tmp
integer          :: iphase
character(LEN=5) :: phase

iphase=0
i=0
do
  read (UNIT=iu, FMT='(A)', IOSTAT=ios, iomsg=iom) cbuf
  if  (ios > 0) then  !ios > 0 chyba
    close (UNIT = iu)
    print *, "Encountered Error:"
    print *, iom
    rmsg=ios
    return
  end if
  if (ios < 0) then  !ios = -1 konec souboru
    if (i > 0) then
       lay%n_layers=i
       allocate(lay%z(i))
       allocate(lay%v(i))
       lay%z=z_tmp(1:i)
       lay%v=v_tmp(1:i)
    end if
    rmsg=ios
    !exit             !break cycle
    return
  end if
  if(ios == 0) then  !ios == 0  ok
    icbuf=index(cbuf,'#')
    if (icbuf > 0 .and. icbuf < 5) cycle
    icbuf=index(cbuf,'wave')
    if (icbuf .gt. 0) then  ! keyword wave
       read(cbuf(5:),'(A)') phase
       iphase=iphase+1
       if (iphase .eq. 1) then
          lay%phase=adjustl(phase)
       else if (iphase > 1) then
          backspace iu
          lay%n_layers=i
          allocate(lay%z(i))
          allocate(lay%v(i))
          lay%z=z_tmp(1:i)
          lay%v=v_tmp(1:i)
          exit 
       end if
    else   
      i = i + 1
      read(cbuf,*) z_tmp(i),v_tmp(i)      
    end if
  end if
end do
rmsg=ios
end function layers_load

function get_layers(ph) result(lptr)

   character(*),intent(IN) :: ph
   type(cLayers),pointer :: lptr
   integer               :: iphase=0
   integer               :: i

   lptr=>null()
   do i=1,5
      if (velmod(i)%n_layers .gt. 0) then
         if (index(velmod(i)%phase,ph).eq.1) then
            lptr=>velmod(i)
            return
         end if
      end if
   end do
end function get_layers

end module layers
