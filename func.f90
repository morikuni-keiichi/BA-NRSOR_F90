module func
  use floating_point_kinds
  implicit none
contains
  function nrm2(x, k) result(norm)

!*  Purpose
!*  =======
!*
!*  nrm2 returns the euclidean norm of a vector via the function
!*  name, so that
!*
!*     nrm2 := sqrt( x'*x )
!*
!*  Further Details
!*  ===============
!*
!*  -- This version written on 25-October-1982.
!*     Modified on 14-October-1993 to inline the call to DLASSQ.
!*     Sven Hammarling, Nag Ltd.
!      translated into Fortran 90/95 and changed by Keiichi Morikun
!
! References:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, and Pete Stewart,
!    LINPACK User's Guide,
!    Society for Industrial and Applied Mathematics (SIAM),
!		 Philadelphia, 1979,
!	   ISBN-10: 089871172X,
!    ISBN-13: 978-0-898711-72-1.
!
!    Charles Lawson, Richard Hanson, David Kincaid, and Fred Krogh,
!    Algorithm 539,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, Pages 308-323.
!*
!*  =====================================================================

  real(real64), intent(in) :: x(:)
  real(real64) absxi, norm, scale, ssq
  integer i, k

  scale = zero
  ssq = one

  do i = 1, k
    if (x(i) /= zero) then
      absxi = abs(x(i))
      if (scale <= absxi) then
        ssq = one + ssq*(scale/absxi)**2
        scale = absxi
      else
        ssq = ssq + (absxi/scale)**2
      endif
    endif
  enddo

  norm = scale*sqrt(ssq)

  return

  end function nrm2


  subroutine rotg(da, db, c, s)

!*  Purpose
!*  =======
!*
!*     rotg construct givens plane rotation.
!*
!*  Further Details
!*  ===============
!*
!*     jack dongarra, linpack, 3/11/78.
!      translated into Fortran 90/95 and changed by Keiichi Morikuni
!
!  References:
!
!    Jack Dongarra, Jim Bunch, Cleve Moler, and Pete Stewart,
!    LINPACK User's Guide,
!    Society for Industrial and Applied Mathematics (SIAM),
!		 Philadelphia, 1979,
!	   ISBN-10: 089871172X,
!    ISBN-13: 978-0-898711-72-1.
!
!    Charles Lawson, Richard Hanson, David Kincaid, and Fred Krogh,
!    Algorithm 539,
!    Basic Linear Algebra Subprograms for Fortran Usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, September 1979, Pages 308-323.
!*
!*  =====================================================================

  real(real64), intent(inout) :: da, db
  real(real64), intent(out) :: c, s
  real(real64) r,roe,scale,z

  roe = db
  if (abs(da) > abs(db)) roe = da

     scale = abs(da) + abs(db)

     if (scale /= zero) go to 10

     c = one
     s = zero
     r = zero
     z = zero
     go to 20

  10 r = scale*sqrt((da/scale)**2 + (db/scale)**2)

      if (roe<zero) r = -r
      c = da / r
      s = db / r
      z = one

      if (abs(da) >  abs(db)) z = s

      if (abs(db) >= abs(da) .and. c /= zero) z = one / c

  20 da = r

  db = z
  return

  end subroutine rotg

end module func