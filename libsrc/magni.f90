
function localmagnitude(hypo,r)
! formula Scherbaum , Stoll ,1983
use arrivals,             only    : cRecordArr
implicit none

real, dimension(3), intent(IN)    :: hypo
type(cRecordArr), intent(IN)      :: r
real                              :: localmagnitude

real, parameter  :: min_dist = 5.0
real             :: hypodist, mag

mag=-9.9
if (r%amp > 0.0 .and.&
    &(r%phase(1:1) .eq. 'S' .or.&
    & r%phase(1:1) .eq. 's' .or. r%phase(1:1) .eq. 'L')) then

    hypodist=sqrt((r%X-hypo(1))**2+(r%Y-hypo(2))**2+(r%Z-hypo(3))**2)
    mag=log10(r%amp/6.283/r%freq*2.8*1e6/0.6325)+0.1+1.4*log10(hypodist)

end if
localmagnitude=mag

end function localmagnitude
