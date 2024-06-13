module HypFile

implicit none
private

integer :: iu
integer :: ou


type, public        :: cRecordHyp
      integer          :: id     = 0  
      character(LEN=5) :: sta    = ''
      character(LEN=4) :: phase  = ''
      integer          :: ichan  = 0
      integer          :: year   = 0
      integer          :: month  = 0
      integer          :: day    = 0
      integer          :: hour   = 0
      integer          :: minute = 0
      integer          :: isec   = 0
      integer          :: msec   = 0
      integer          :: micros = 0
      real             :: wt     = 0
      real             :: amp_d  = 0
      integer          :: pol    = 0
      real             :: amp_v  = 0
      real             :: period = 0
contains
procedure :: rrec => recordhyp_read
procedure :: wrec => recordhyp_write
procedure :: wrec_hyr => recordhyp_write_hyr
procedure :: getYear => recordhyp_getyear
end type cRecordHyp

type, public       :: cFileHyp
      integer      :: n_rec
      logical      :: hyr
      type(cRecordHyp), allocatable, dimension(:) :: rec
      character(LEN=:), allocatable :: filename
contains

procedure :: load  => filehyp_load
procedure :: dump  => filehyp_dump
end type cFileHyp

type(cFileHyp),public     :: hyp
save

contains

function recordhyp_read(r) result(res)

class(cRecordHyp), intent(out) ::  r
integer                        ::  res

character(LEN=100) :: cbuf
integer            :: iline

   read (iu,'(a)',iostat=res) cbuf
   if(res /= 0) then
     return
   endif

   read (cbuf,'(a5)') r%sta

   r%phase = ''
      iline=index(cbuf(6:),'P');
      if (iline .gt. 0) then
         r%phase = 'P'
      endif
   if (iline .le. 0) then
      iline=index(cbuf(6:),'S');
      if (iline .gt. 0) then
         r%phase = 'S'
      endif
   endif
   if (iline .le. 0) then
      res = -1
      return
   endif

   iline = iline + 6

   read (cbuf(iline:),*)&
      &r%ichan, r%year, r%month, r%day,&
      &r%hour, r%minute, r%isec, r%msec,&
      &r%micros, r%wt, r%amp_d, r%pol,&
      &r%amp_v, r%period

end function recordhyp_read

subroutine recordhyp_write(r)

class(cRecordHyp), intent(in) ::  r

    write(ou,101) r%sta, r%phase, r%ichan,&
     &mod(r%year,100), r%month, r%day,&
     &r%hour, r%minute, r%isec, r%msec,&
     &r%micros, int(r%wt), r%amp_d, r%pol,&
     &r%amp_v, r%period

101   format(a5,1x,a1,1x,i2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,&
      &1x,i3.3,1x,i3.3,1x,i1,1x,1pe12.6,2x,i1,1x,1pe12.6,1x,1pe12.6)

end subroutine recordhyp_write

subroutine recordhyp_write_hyr(r)

class(cRecordHyp), intent(in) ::  r

    write(ou,101) r%sta, r%phase, r%ichan,&
     &mod(r%year,100), r%month, r%day,&
     &r%hour, r%minute, r%isec, r%msec,&
     &r%micros, (r%wt), r%amp_d, r%pol,&
     &r%amp_v, r%period

101   format(a5,1x,a1,1x,i2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,1x,i2.2,&
      &1x,i3.3,1x,i3.3,1x,f5.0,1x,1pe12.6,2x,i1,1x,1pe12.6,1x,1pe12.6)

end subroutine recordhyp_write_hyr

function recordhyp_getyear(r) result(year)

class(cRecordHyp), intent(in) ::  r
integer           :: year

    year = mod(r%year,100)
    if (year .gt. 70) then
      year = year + 1900
    else
      year= year + 2000
    endif

end function recordhyp_getyear


function filehyp_load(hyp,hypname) result(rmsg)

class(cFileHyp), intent(out)  :: hyp
integer                       :: rmsg
character(LEN=*)              :: hypname

integer                       :: res
integer                       :: ios
character (256)               :: iom
character(LEN=120)            :: cbuffer
type(cRecordHyp)              :: r
integer                       :: i

hyp%filename = trim(hypname)

open (NEWUNIT=iu,FILE=hyp%filename,FORM='FORMATTED',STATUS='OLD',IOSTAT=ios,IOMSG=iom)
if (ios/=0) then
    write(*,*) 'iostat = ', ios
    write(*,*) 'iomsg: '//trim(iom)
    rmsg=ios
    stop
end if

hyp%hyr = .false.
i=len_trim(hyp%filename)
if ( hyp%filename(i-3:i) .eq. '.hyr' .or. &
    &hyp%filename(i-3:i) .eq. '.HYR'     ) then
  hyp%hyr = .true.
endif

hyp%n_rec = 0
do
  read (UNIT =iu, FMT ='(A)', IOSTAT =res) cbuffer
  if(res == 0) then  !res == 0  ok
    hyp%n_rec = hyp%n_rec + 1
    cycle
  endif
  if(res < 0) then   !res = -1 konec souboru
    exit             !break cycle
  endif
  close (UNIT = iu)  !res > 0 chyba
  rmsg=res
  stop
end do

allocate(hyp%rec(hyp%n_rec))

rewind(iu)

do i=1,hyp%n_rec
  res=r%rrec()
  if(res == 0) then  !res == 0  ok
    r%id=i
    hyp%rec(i)=r
    cycle            
  endif
  if(res < 0) then   !res = -1 konec souboru
    exit             !break cycle
  endif
  close (UNIT = iu)  !res > 0 chyba
  rmsg=res
  stop
end do

close (UNIT = iu)
rmsg=0

end function filehyp_load

function filehyp_dump(hyp,hypname) result(rmsg)

class(cFileHyp), intent(in)   :: hyp
integer                       :: rmsg
character(LEN=*)              :: hypname

integer                       :: res
integer                       :: ios
character (256)               :: iom
type(cRecordHyp)              :: r
integer                       :: i

open (NEWUNIT =ou, FILE = hypname,  ACTION='write', STATUS='replace', IOSTAT=ios,IOMSG=iom)
if (ios/=0) then
    write(*,*) 'iostat = ', ios
    write(*,*) 'iomsg: '//trim(iom)
    rmsg=ios
    stop
end if

do i=1,hyp%n_rec
  r=hyp%rec(i)
  if(hyp%hyr) then
    call r%wrec_hyr()
  else
    call r%wrec()
  endif
end do

close (UNIT = ou)
rmsg=0

end function filehyp_dump

end module HypFile


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

