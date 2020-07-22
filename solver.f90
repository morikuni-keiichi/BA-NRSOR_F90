module solver
  use func
  use,intrinsic :: iso_fortran_env
  use floating_point_kinds
  implicit none
  private :: NRSOR, atNRSOR
  public :: BAGMRES
contains
  subroutine BAGMRES(AC, ia, jp, m, n, &
                      x0, b, &
                      at, nin, omg, &
                      tol, omax, rmax, &
                      x, &
                      iter, riter, relres, conv, verbose)

  real(real64) H(omax+1, omax), V(n, omax+1)
  real(real64), intent(in) :: b(:), AC(:)
  real(real64), intent(out) :: relres(:), x(:)  
  real(real64) c(omax), g(omax+1), r(m), s(omax), & 
            tmp_x(n), w(n), x0(n), y(omax), Aei(n)
  real(real64), intent(in) :: tol
  real(real64), intent(inout) :: omg
  real(real64) :: beta, beta0, eps, inprod, nrmATr, nrmATr0, tmp

  integer, intent(in) :: ia(:), jp(:)
  integer, intent(in) :: at, m, n, omax, rmax, verbose
  integer, intent(inout) :: nin
  integer, intent(out) :: iter, riter, conv
  integer :: bd = 0, i, iter_tot = 1, j, k, k1, k2, l, p
  
  eps = epsilon(eps)
  conv = 0

!	The initial approximate solution is equal to zero
  x0(1:n) = zero

! The 2-norm of a_j, j=1,...,n, required for NR-SOR
  do j = 1, n
    k1 = jp(j)
    k2 = jp(j+1)-1
    inprod = sum(AC(k1:k2)*AC(k1:k2))
    if (inprod /= zero) then
      Aei(j) = one / inprod
    else
      write(*, *) 'warning: ||aj|| = 0.0, j =', j      
      stop
    endif
  enddo

! nrmATr0 = ||A^T b||_2
  w(1:n) = zero
  do j = 1, n
    k1 = jp(j)
    k2 = jp(j+1)-1
    w(j) = sum(AC(k1:k2)*b(ia(k1:k2)))
  enddo
  nrmATr0 = nrm2(w(1:n), n)

  ! Main loop
  do p = 0, rmax

    r(1:m) = zero
    do j = 1, n
      tmp = x0(j)
      do l = jp(j), jp(j+1)-1
        i = ia(l)
        r(i) = r(i) + tmp*AC(l)
      enddo
    enddo

  ! r0 = b - A*x0
    r(1:m) = b(1:m) - r(1:m)

    if (at > 0) then
      if (verbose == 1) then
        write(*, '(a)') '  Automatic parameter tuning started'
      endif
      call atNRSOR(AC, ia, jp, m, n, r, Aei, nin, omg, w)
      if (verbose == 1) then
        write(*, '(a)') '  Automatic parameter tuning finished'
      endif
    else if (at == 0) then
    ! NR-SOR inner-iteration preconditioning without automatic parameter tuning
      call NRSOR(AC, ia, jp, n, r, Aei, nin, omg, w)
    endif

    if (p == 0) then
      beta = nrm2(w(1:n), n)
      beta0 = beta
    else
      beta = nrm2(w(1:n), n)
    endif

  ! Normalization of the vector v_1
    tmp = one / beta
    V(1:n, 1) = tmp * w(1:n)

    g(1) = beta

  ! Restart loop
    do k = 1, omax

      r(1:m) = zero
      do j = 1, n
        tmp = V(j, k)
        do l = jp(j), jp(j+1)-1
          i = ia(l)
          r(i) = r(i) + tmp*AC(l)
        enddo
      enddo

    ! NR-SOR inner-iteration preconditioning
      call NRSOR(AC, ia, jp, n, r, Aei, nin, omg, w)

    ! Modified Gram-Schmidt orthogonalization
      do i = 1, k
        tmp = sum(w(1:n) * V(1:n, i))
        w(1:n) = w(1:n) - tmp*V(1:n, i)
        H(i, k) = tmp
      enddo

      tmp = nrm2(w(1:n), n)

      if (tmp > eps) then
      ! Normalization of the vector v_{k+1}
        H(k+1, k) = tmp
        tmp = one / tmp
        V(1:n, k+1) = tmp * w(1:n)
      else
        write(*, *) 'BREAKDOW at step', k
        write(*, *) 'h_{k+1, k} =', H(k+1, k)
        bd = 1
      endif

    ! Application of Givens rotations
      do i = 1, k-1
        tmp =  c(i)*H(i, k) + s(i)*H(i+1, k)
        H(i+1, k) = -s(i)*H(i, k) + c(i)*H(i+1, k)
        H(i, k) = tmp
      enddo

    ! Construction of Givens rotations
      call rotg(H(k, k), H(k+1, k), c(k), s(k))

      H(k, k) = one / H(k, k)

      g(k+1) = -s(k)*g(k)
      g(k) = c(k)*g(k)

      relres(iter_tot) = abs(g(k+1)) / beta0

    !	Convergence check
      if (relres(iter_tot) < tol .or. & 
         k == omax .or. bd == 1) then

        !	Derivation of the approximate solution x_k
        !	Backward substitution
        y(k) = g(k) * H(k, k)
        do i = k-1, 1, -1
          y(i) = (g(i) - sum(H(i, i+1:k) * y(i+1:k))) * H(i, i)
        enddo

        x(1:n) = x0(1:n) + matmul(V(1:n, 1:k), y(1:k))

        if (bd == 1) then
          iter = iter_tot
          Riter = p
          return
        endif

        r(1:m) = zero
        do j = 1, n
          tmp = x(j)
          do l = jp(j), jp(j+1)-1
            i = ia(l)
            r(i) = r(i) + tmp*AC(l)
          enddo
        enddo

        ! r_k = b - A*x_k
        r(1:m) = b(1:m) - r(1:m)

      ! w = A^T r
        do j = 1, n
          k1 = jp(j)
          k2 = jp(j+1)-1
          w(j) = sum(AC(k1:k2) * r(ia(k1:k2)))
        enddo

        if (nrm2(w(1:n), n) < tol*nrmATr0) then
          
          conv = 1
          iter = iter_tot
          Riter = p          

          return

        endif

      endif

      iter_tot = iter_tot + 1

    enddo

    x0(1:n) = x(1:n)

  enddo

  iter = iter_tot - 1
  Riter = p - 1

  end subroutine BAGMRES

  subroutine NRSOR(AC, ia, jp, n, rhs, Aei, nin, omg, x)

  real(real64), intent(in) :: AC(:), Aei(:)
  real(real64), intent(inout) :: rhs(:)
  real(real64), intent(out) :: x(:)
  integer, intent(in) :: ia(:), jp(:)
  integer, intent(in) :: n, nin
  real(real64), intent(in) :: omg
  real(real64) d
  integer i, j, k, k1, k2, l

  x(1:n) = zero

  do k = 1, nin
    do j = 1, n
      k1 = jp(j)
      k2 = jp(j+1)-1
      d = omg * sum(rhs(ia(k1:k2))*AC(k1:k2)) * Aei(j)
      x(j) = x(j) + d
      if (k == nin .and. j == n) return
      do l = k1, k2
        i = ia(l)
        rhs(i) = rhs(i) - d*AC(l)
      enddo
    enddo
  enddo

  end subroutine NRSOR

  subroutine atNRSOR(AC, ia, jp, m, n, rhs, Aei, nin, omg, x)

  real(real64), intent(in) :: AC(:), Aei(:)
  real(real64), intent(inout) :: rhs(:)
  real(real64), intent(out) :: x(:)
  integer, intent(in) :: ia(:), jp(:)
  integer, intent(in) :: m, n
  integer, intent(out) :: nin
  real(real64) :: d, e, omg, res1, res2 = zero, y(n), tmprhs(m)
  integer i, ii, j, k, k1, k2, l

  tmprhs(1:m) = rhs(1:m)

  x(1:n) = zero
  y(1:n) = zero

! 	Tune the number of inner iterations
  do i = 1, 50

    do j = 1, n
      k1 = jp(j)
      k2 = jp(j+1)-1
      d = sum(rhs(ia(k1:k2))*AC(k1:k2)) * Aei(j)
      x(j) = x(j) + d
      do l = k1, k2
        ii = ia(l)
        rhs(ii) = rhs(ii) - d*AC(l)
    enddo
  enddo
  
  d = maxval(abs(x(1:n)))
  e = maxval(abs(x(1:n) - y(1:n)))

  if (e < 1.0d-1 * d .or. i == 50) then
    nin = i
    exit
  endif

  y(1:n) = x(1:n)

  enddo

! 	Tune the relaxation parameter
  do k = 19, 1, -1

    omg = 1.0d-1 * dble(k)

    rhs(1:m) = tmprhs(1:m)

    x(1:n) = zero

    do i = 1, nin
      do j = 1, n
        k1 = jp(j)
        k2 = jp(j+1)-1
        d = omg * sum(rhs(ia(k1:k2))*AC(k1:k2)) * Aei(j)
        x(j) = x(j) + d
        do l = k1, k2
          ii = ia(l)
          rhs(ii) = rhs(ii) - d*AC(l)
        enddo
      enddo
    enddo

    res1 = nrm2(rhs(1:m), m)

    if (k < 19) then
      if (res1 > res2) then
        omg = omg + 1.0d-1
        x(1:n) = y(1:n)
        return
      else if (k == 1) then
        omg = 0.1d0
        return
      endif
    endif

    res2 = res1

    y(1:n) = x(1:n)

  enddo

  end subroutine atNRSOR
end module solver
