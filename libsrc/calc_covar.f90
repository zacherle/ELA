module calc_covar
use, intrinsic :: iso_fortran_env, only : real64

use layers,       only: model_error, reading_error
use const_raddeg

implicit none

contains

subroutine eigv_sym4(a,w)

integer, parameter :: n = 4
integer, parameter :: lda = n
real(real64), dimension(lda,n) :: a
real(real64), dimension(n)     :: w

integer, parameter :: lwmax = 200
integer            :: info, lwork
real(real64), dimension(lwmax) :: work

!  evaluate eigenvalues
!
![in]	JOBZ	CHARACTER*1
!          = 'N':  Compute eigenvalues only;
!          = 'V':  Compute eigenvalues and eigenvectors.
!
![in]	UPLO	CHARACTER*1
!          = 'U':  Upper triangle of A is stored;
!          = 'L':  Lower triangle of A is stored.
!
![in]	N	INTEGER
!          The order of the matrix A.  N >= 0.
!
![in,out]	A
!          A is DOUBLE PRECISION array, dimension (LDA, N)
!          On entry, the symmetric matrix A.  If UPLO = 'U', the
!          leading N-by-N upper triangular part of A contains the
!          upper triangular part of the matrix A.  If UPLO = 'L',
!          the leading N-by-N lower triangular part of A contains
!          the lower triangular part of the matrix A.
!          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
!          orthonormal eigenvectors of the matrix A.
!          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
!          or the upper triangle (if UPLO='U') of A, including the
!          diagonal, is destroyed.
!
![in]	LDA	INTEGER
!          The leading dimension of the array A.  LDA >= max(1,N).
!
![out]	W
!          W is DOUBLE PRECISION array, dimension (N)
!          If INFO = 0, the eigenvalues in ascending order.
!
![out]	WORK
!          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
!          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
!
![in]	LWORK	INTEGER
!          The length of the array WORK.  LWORK >= max(1,3*N-1).
!          For optimal efficiency, LWORK >= (NB+2)*N,
!          where NB is the blocksize for DSYTRD returned by ILAENV.
!
!          If LWORK = -1, then a workspace query is assumed; the routine
!          only calculates the optimal size of the WORK array, returns
!          this value as the first entry of the WORK array, and no error
!          message related to LWORK is issued by XERBLA.
!
![out]	INFO	INTEGER
!          = 0:  successful exit
!          < 0:  if INFO = -i, the i-th argument had an illegal value
!          > 0:  if INFO = i, the algorithm failed to converge; i
!                off-diagonal elements of an intermediate tridiagonal
!                form did not converge to zero.

      lwork = -1
!     call dsyev( 'Vectors', 'Upper', n, a, lda, w, work, lwork, info )
      call dsyev( 'Vectors', 'Lower', n, a, lda, w, work, lwork, info )
      lwork = int( work( 1 ) )
!write(*,*) lwork
!      call dsyev( 'Vectors', 'Upper', n, a, lda, w, work, lwork, info )
      call dsyev( 'Vectors', 'Lower', n, a, lda, w, work, lwork, info )
      if( info.gt.0 ) then
         write(*,*)'The algorithm failed to compute eigenvalues.'
         stop
      end if

end subroutine eigv_sym4

function cov_matrix1(m,n,resw,Gw,w)
!covariance matrix
integer, intent(IN)          :: m,n
real(real64), dimension(m), intent(IN)   :: resw
real(real64), dimension(m), intent(IN)   :: w
real(real64), dimension(m,n), intent(IN) :: Gw
real(real64), dimension(n,n) :: cov_matrix1

real(real64), parameter :: kappa=1.0d-08

real(real64), dimension(n,n) :: GG,eigV,covM,pmat !,eigL
real(real64), dimension(n) :: eig_l
real(real64)               :: varD
integer                    :: marrv
real(real64)               :: sumw
integer :: j

!GG=G*G'
GG=MATMUL(TRANSPOSE(Gw),Gw)

!eigenvalues (ortogonal) decomposition
!V*diag(l)*V'=GG
eigV=GG
call eigv_sym4(eigV,eig_l)
!test
!   eigL=0.0
!   do i=1,4
!         eigL(i,i)=eig_l(i)
!   end do
!   GG=MATMUL(eigV,MATMUL(eigL,TRANSPOSE(eigV)))

if (minval(eig_l)/maxval(eig_l)<kappa) then
        eig_l=eig_l+kappa
end if

!covariance form eigenvalues
!C=V*diag(1/l)*V'
do j=1,n
    pmat(j,:)=eigV(:,j)/eig_l(j)
end do
covM=MATMUL(eigV,pmat)
!test
!   eigL=0.0
!   do i=1,4
!         eigL(i,i)=1/eig_l(i)
!   end do
!   covM=MATMUL((eigV),MATMUL(eigL,TRANSPOSE(eigV)))


sumw=sum(w)
marrv=count(w>0)
varD=dot_product(resw,resw)/sumw
if (marrv>4) then
        varD=varD*marrv/(marrv-4)
else
        varD=varD*marrv
end if

! reading error
   if (varD .lt. reading_error**2) then
      varD = reading_error**2
   endif

! model error
   varD=varD+model_error**2


cov_matrix1=covM*varD

end function cov_matrix1


subroutine get_errellipse(co,dxer,dyer,dzer, dter,l1,l2,theta)

real(real64), dimension(4,4),intent(IN) :: co
real,intent(OUT) ::   dxer,dyer,dzer,dter
real,intent(OUT) ::   theta
real,intent(OUT) ::   l1,l2

real(real64)    deter
real(real64)    d11,d21,d22
real(real64)    al,bl
real(real64)    tl

! error ellipse for epicenter
   deter=co(1,1)*co(2,2)-co(1,2)*co(2,1)
   d11=co(2,2)/deter
   d22=co(1,1)/deter
   d21=-co(2,1)/deter
   theta=atan(2*d21/(d11-d22))/2
   al=d11*cos(theta)**2 + 2*d21*cos(theta)*sin(theta) + d22*sin(theta)**2
   bl=d11*sin(theta)**2 - 2*d21*cos(theta)*sin(theta) + d22*cos(theta)**2

   l1=sqrt(1/al)
   l2=sqrt(1/bl)

! major axis of the error ellipse to the x-axis
   theta = real(theta*RAD2DEG)

! l1 is major axis
   if (l2 .gt. l1 ) then
      tl=l1
      l1=l2
      l2=tl
      theta=theta+90.0
   endif

   dxer=sqrt(abs(co(1,1)))
   dyer=sqrt(abs(co(2,2)))

   dzer=sqrt(abs(co(3,3)))
   dter=sqrt(abs(co(4,4)))

end subroutine get_errellipse

end module calc_covar
