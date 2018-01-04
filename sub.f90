module sub
  use func
  use globvar
  implicit none
  public :: read_prm, read_mat, output, discard
contains
subroutine read_prm(nin, omg, omax, output_mode, rmax, tol, at)
  use globvar
  implicit none

  real(dp), intent(out) :: omg, tol
  integer, intent(out) :: at, nin, omax, output_mode, rmax
  integer :: fi = 11, i, is, narg
  character*100 argv

  narg = command_argument_count()
  if (narg > 8) then
    stop "Error: Too many arguments."
  else if (narg < 8) then
    stop "Error: Too few arguments."
  endif

  do i = 1, narg
      call getarg(i, argv)
      if (index(argv(1:5), "--at=") == 1) then
         read(argv(6:len_trim(argv)), *) at
         write(*, *) "at=", at
      else if (index(argv(1:6), "--nin=") == 1) then
         read(argv(7:len_trim(argv)), *) nin
         write(*, *) "nin=", nin
      else if (index(argv(1:6), "--omg=") == 1) then
         read(argv(7:len_trim(argv)), *) omg
         write(*, *) "omax(1)=", omg
      else if (index(argv(1:6), "--tol=") == 1) then
         read(argv(7:len_trim(argv)), *) tol
         write(*, *) "tol=", tol
      else if (index(argv(1:7), "--omax=") == 1) then
         read(argv(8:len_trim(argv)), *) omax
         write(*, *) "omax=", omax
      else if (index(argv(1:7), "--rmax=") == 1) then
         read(argv(8:len_trim(argv)), *) rmax
         write(*, *) "rmax=", rmax
      else if (index(argv(1:14), "--output_mode=") == 1) then
         read(argv(15:len_trim(argv)), *) output_mode
         write(*, *) "output_mode=", output_mode
      else if (index(argv(1:5), "--fi=") /= 1) then
         write(*, *) argv
         stop "Unknown argument"
      endif
   enddo

end subroutine read_prm
!---------------------------------------------------------------------------

subroutine read_mat(AC, ia, jp, m, n, x0, b)
  use globvar
  implicit none

  real(dp), allocatable, intent(out) :: AC(:), b(:), x0(:)

  integer, allocatable, intent(out) :: ia(:), jp(:)
  integer, intent(inout) :: m, n
  integer is, j, l, nnz
  integer :: fi_AC = 14, fi_jp = 15, fi_ia = 16

  open(fi_AC, file = 'RANDL7/AC.ccs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open AC file'
  read(fi_AC, *) m, n, nnz

  allocate(AC(nnz))
  allocate(jp(n+1))
  allocate(ia(nnz))
  allocate(b(m))
  allocate(x0(n))

! Load AC
  read(fi_AC, *) (AC(l), l = 1, nnz)
  close(fi_AC)

! Load jp
  open(fi_jp, file ='RANDL7/jp.ccs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open jp file'
  read(fi_jp, *) (jp(j), j = 1, n+1)
  close(fi_jp)

! Load ia
  open(fi_ia, file = 'RANDL7/ia.ccs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open ia file'
  read(fi_ia, *) (ia(l), l = 1, nnz)
  close(fi_ia)

!	b : = random
!	call random_seed
  call random_number(b(1:m))
  b(1:m) = 2.0d0*b(1:m) - one

end subroutine read_mat


subroutine output(AC, ia, jp, m, n, &
                  b, &
                  nin, omg, &
                  x, &
                  iter, omax, riter, relres, t_tot, &
                  output_mode)
  use func
  use globvar
  implicit none

  real(dp), intent(in) :: AC(:), b(:)  
  real(dp), intent(in) :: t_tot, x(n)
  real(dp), intent(inout) :: relres(iter)
  real(dp) r(m), ATr(n)
  real(dp), intent(in) :: omg
  real(dp) tmp, nrm_b, nrm_r, nrmATb, nrmATr

  integer, intent(in) :: ia(:), jp(:)
  integer, intent(in) :: iter, m, n, nin, omax, output_mode, riter
  integer :: info = 10, reshis = 11, sol = 12
  integer i, is, j, l, k, k1, k2

  nrm_b = nrm2(b(1:m), m)

  do j = 1, n
    k1 = jp(j)
    k2 = jp(j+1)-1
    ATr(j) = sum(AC(k1:k2) * b(ia(k1:k2)))
  enddo

  nrmATb = nrm2(ATr(1:n), n)

  r(1:m) = zero
  do j = 1, n
    tmp = x(j)
    do l = jp(j), jp(j+1)-1
      i = ia(l)
      r(i) = r(i) + tmp*AC(l)
    enddo
  enddo

  r(1:m) = b(1:m) - r(1:m)

  nrm_r = nrm2(r(1:m), m)

  do j = 1, n
    k1 = jp(j)
    k2 = jp(j+1)-1
    ATr(j) = sum(AC(k1:k2) * r(ia(k1:k2)))
  enddo

  nrmATr = nrm2(ATr(1:n), n)

  if (output_mode == 0) then

    write(*, '(a, f6.3)') ' omega: ', omg
    write(*, *) '# of outer iterations: ', iter
    write(*, *) '# of inner iterations: ', nin
    write(*, *) '# of restarts: ', riter
    write(*, *) 'Restart cycle: ', omax
    write(*, '(a, f16.5)') ' CPU time: ', t_tot
    write(*, *) 'Relative residual: ', relres(iter)
    write(*, *) 'Actual relative residual (ATr): ', nrmATr / nrmATb
    write(*, *) 'Actual relative residual (r): ', nrm_r / nrm_b

  else

    write(*, '(f6.3, f10.5, i4)') omg, t_tot, iter

  endif

  open(info, file='info.dat', action='write', iostat=is, status='replace')
  if (is /= 0) stop 'cannot open info.dat file'
  write(info, '(a, f6.3)') ' omega: ', omg
  write(info, *) '# of outer iterations: ', iter
  write(info, *) '# of inner iterations: ', nin
  write(info, *) '# of restarts: ', riter
  write(info, *) 'Restart cycle: ', omax
  write(info, '(a, f16.5)') ' CPU time: ', t_tot
  write(info, *) 'Relative residual: ', relres(iter)
  write(info, *) 'Actual relative residual (ATr): ', nrmATr / nrmATb
  write(info, *) 'Actual relative residual (r): ', nrm_r / nrm_b
  close(info)

  open(reshis, file='reshis.dat', action='write', iostat=is, status='replace')
  if (is /= 0) stop 'cannot open reshis.dat file'
  do k = 1, iter
    write(reshis, *) k, relres(k)
  enddo
  close(reshis)

  open(sol, file='solution.dat', action='write', iostat=is, status='replace')
  if (is /= 0) stop 'cannot open solution.dat file'
  do j = 1, n
    write(sol, *) x(j)
  enddo
  close(sol)

end subroutine output

subroutine discard(AC, ia, jp, b)
  use globvar
  implicit none
  real(dp), allocatable, intent(inout) :: AC(:), b(:)
  integer,  allocatable, intent(inout) :: ia(:), jp(:)

  deallocate(AC, jp, ia)
  deallocate(b)

end subroutine discard
end module sub