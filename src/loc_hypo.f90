
subroutine loc_hypo_lm(hypo, m, fvec, fjac, fix_depth, info)

!use povrch,             only: relief8
use fw_tdx, only: Xp, Xo, tdx3, tdxder  !, tdxdif

implicit none

integer, parameter :: dp=kind(1.0d0)

interface

   subroutine fdjac2(fcn, m, n, x, Fvec, Fjac, Ldfjac, Iflag, Epsfcn, Wa)
use iso_fortran_env, only: wp => real64
        implicit none
        procedure() :: fcn !! the user-supplied subroutine which
                            !! calculates the functions.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                            !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                            !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                 !! which specifies the leading dimension of the array fjac.
        integer, intent(inout) :: Iflag !! an integer variable which can be used to terminate
                                 !! the execution of fdjac2. see description of [[func2]].
        real(wp), intent(in) :: Epsfcn !! an input variable used in determining a suitable
                                  !! step length for the forward-difference approximation. this
                                  !! approximation assumes that the relative errors in the
                                  !! functions are of the order of epsfcn. if epsfcn is less
                                  !! than the machine precision, it is assumed that the relative
                                  !! errors in the functions are of the order of the machine
                                  !! precision.
        real(wp), intent(inout) :: x(n) !! an input array of length n.
        real(wp), intent(in) :: Fvec(m) !! an input array of length m which must contain the
                                   !! functions evaluated at x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array which contains the
                                            !! approximation to the jacobian matrix evaluated at x.
        real(wp), intent(inout) :: Wa(m) !! a work array of length m.
   end subroutine fdjac2

   subroutine lmder1(fcn, m, n, x, Fvec, Fjac, Ldfjac, Tol, Info, Ipvt, Wa, Lwa)
use iso_fortran_env, only: wp => real64
        implicit none

        procedure() :: fcn
        !procedure(fcn_lmder) :: fcn !! user-supplied subroutine which
                                    !! calculates the functions and the jacobian.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(in) :: Ldfjac !! a positive integer input variable not less than m
                                     !! which specifies the leading dimension of the array fjac.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows.
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  algorithm estimates that the relative error
                                    !!    in the sum of squares is at most tol.
                                    !!  * ***info = 2***  algorithm estimates that the relative error
                                    !!    between x and the solution is at most tol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
                                    !!    jacobian to machine precision.
                                    !!  * ***info = 5***  number of calls to fcn with iflag = 1 has
                                    !!    reached 100*(n+1).
                                    !!  * ***info = 6***  tol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  tol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than 5*n+m.
        integer, intent(out) :: Ipvt(n) !! an integer output array of length n. ipvt
                                       !! defines a permutation matrix p such that jac*p = q*r,
                                       !! where jac is the final calculated jacobian, q is
                                       !! orthogonal (not stored), and r is upper triangular
                                       !! with diagonal elements of nonincreasing magnitude.
                                       !! column j of p is column ipvt(j) of the identity matrix.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                   !! when the algorithm estimates either that the relative
                                   !! error in the sum of squares is at most tol or that
                                   !! the relative error between x and the solution is at
                                   !! most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(out) :: Fjac(Ldfjac, n) !! an output m by n array. the upper n by n submatrix
                                                !! of fjac contains an upper triangular matrix r with
                                                !! diagonal elements of nonincreasing magnitude such that
                                                !!```
                                                !!        t     t           t
                                                !!       p *(jac *jac)*p = r *r,
                                                !!```
                                                !! where p is a permutation matrix and jac is the final
                                                !! calculated jacobian. column j of p is column ipvt(j)
                                                !! (see below) of the identity matrix. the lower trapezoidal
                                                !! part of fjac contains information generated during
                                                !! the computation of r.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

    end subroutine lmder1


    subroutine lmdif1(fcn, m, n, x, Fvec, Tol, Info, Iwa, Wa, Lwa)
use iso_fortran_env, only: wp => real64
        implicit none

        procedure() :: fcn !! the user-supplied subroutine which
                                !! calculates the functions.
        integer, intent(in) :: m !! a positive integer input variable set to the number
                                !! of functions.
        integer, intent(in) :: n !! a positive integer input variable set to the number
                                !! of variables. n must not exceed m.
        integer, intent(out) :: Info !! an integer output variable. if the user has
                                    !! terminated execution, info is set to the (negative)
                                    !! value of iflag. see description of fcn. otherwise,
                                    !! info is set as follows:
                                    !!
                                    !!  * ***info = 0***  improper input parameters.
                                    !!  * ***info = 1***  algorithm estimates that the relative error
                                    !!    in the sum of squares is at most tol.
                                    !!  * ***info = 2***  algorithm estimates that the relative error
                                    !!    between x and the solution is at most tol.
                                    !!  * ***info = 3***  conditions for info = 1 and info = 2 both hold.
                                    !!  * ***info = 4***  fvec is orthogonal to the columns of the
                                    !!    jacobian to machine precision.
                                    !!  * ***info = 5***  number of calls to fcn has reached or
                                    !!    exceeded 200*(n+1).
                                    !!  * ***info = 6***  tol is too small. no further reduction in
                                    !!    the sum of squares is possible.
                                    !!  * ***info = 7***  tol is too small. no further improvement in
                                    !!    the approximate solution x is possible.
        integer, intent(in) :: Lwa !! a positive integer input variable not less than
                                  !! m*n+5*n+m.
        integer, intent(inout) :: Iwa(n) !! an integer work array of length n.
        real(wp), intent(in) :: Tol !! a nonnegative input variable. termination occurs
                                   !! when the algorithm estimates either that the relative
                                   !! error in the sum of squares is at most tol or that
                                   !! the relative error between x and the solution is at
                                   !! most tol.
        real(wp), intent(inout) :: x(n) !! an array of length n. on input x must contain
                                       !! an initial estimate of the solution vector. on output x
                                       !! contains the final estimate of the solution vector.
        real(wp), intent(out) :: Fvec(m) !! an output array of length m which contains
                                        !! the functions evaluated at the output x.
        real(wp), intent(inout) :: Wa(Lwa) !! a work array of length lwa.

    end subroutine lmdif1

    pure real(wp) function enorm(n, x)
use iso_fortran_env, only: wp => real64

        implicit none

        integer, intent(in) :: n !! a positive integer input variable.
        real(wp), intent(in) :: x(n) !! an input array of length n.
    end function enorm

end interface

real(dp), dimension(4), intent(INOUT) :: hypo
integer, intent(IN)                   :: m
real(dp), dimension(m), intent(OUT)   :: fvec(m)
real(dp), dimension(m,4), intent(OUT) :: fjac(m, 4)
logical, intent(IN)  :: fix_depth
integer, intent(OUT) :: info

real(dp),dimension(4) :: x
real(dp) :: tol
integer :: ipvt(size(x))
real(dp), allocatable :: wa(:)
!integer i, Iflag
!real(dp) eps

tol = sqrt(epsilon(1._dp))

if (fix_depth) then
        ! 5*n+m
        allocate(wa(5*size(x) + size(fvec)))
        x = hypo
        call lmder1(tdx3, m, size(x), x, fvec, fjac, size(fjac, 1), 100*tol, &
                info, ipvt, wa, size(wa))
        deallocate(wa)
else
        ! m*n+5*n+m.
        allocate(wa(size(fvec)*size(x)+5*size(x) + size(fvec)))
!        x = Xo(hypo,relief8(hypo(1),hypo(2)))
        x = Xo(hypo,0d0)
        call lmder1(tdxder, m, size(x), x, fvec, fjac, size(fjac, 1), tol, &
                info, ipvt, wa, size(wa))
        deallocate(wa)
!        x=Xp(x,relief8(x(1),x(2)))
        x = Xp(x,0d0)
end if

print 1000, enorm(m, fvec), info, x
1000 format(5x, 'FINAL L2 NORM OF THE RESIDUALS', d15.7 // &
            5x, 'EXIT PARAMETER', 16x, i10              // &
            5x, 'FINAL APPROXIMATE SOLUTION'            // &
            5x, 4d15.7)
hypo=x

end subroutine loc_hypo_lm
