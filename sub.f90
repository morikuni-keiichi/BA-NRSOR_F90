module sub
  use func
  use,intrinsic :: iso_fortran_env
  use floating_point_kinds, only : zero, one
  implicit none
  public :: read_prm, read_mat, output, discard
contains
subroutine read_prm(nin, omg, omax, verbose, rmax, tol, at, directory)
  implicit none

  real(real64), intent(inout) :: omg, tol
  integer, intent(inout) :: at, nin, omax, rmax, verbose
  character(*), intent(out) :: directory
  integer :: i, ierror, length, narg, satus
  integer :: set_at = 0, set_fi = 0, set_nin = 0, set_omg = 0, set_tol = 0, set_omax = 0, set_rmax = 0
  character*100 argv

  narg = command_argument_count()

  do i = 1, narg
    call get_command_argument(i, argv, length, satus)
    if (satus > 0) then
      error stop "in get_command_argument"
    endif
    if (index(argv(1:5), "--at=") == 1) then        ! Flag for automatic parameter tuning
      read(argv(6:len_trim(argv)), '(i1)', iostat=ierror) at          
      if (ierror /= 0) then          
        error stop "option --at="
      endif
      set_at = 1
      if (verbose == 1) then
        if (at == 0) then
          write(*, '(a)') "  Disabled automatic paramete tuning"
        else
          write(*, '(a)') "  Enabled automatic paramete tuning"
        endif
      endif
    endif
  enddo

  do i = 1, narg
    call get_command_argument(i, argv, length, satus)
    if (satus > 0) then
      error stop "in get_command_argument"
    endif
    if (index(argv(1:5), "--at=") == 1) then        ! Flag for automatic parameter tuning
    elseif (index(argv(1:6), "--nin=") == 1) then  ! Maximum number of inner iterations
      read(argv(7:len_trim(argv)), '(i10)', iostat=ierror) nin
      if (ierror /= 0) then          
        error stop "option --nin="
      endif
      set_nin = 1
      if (verbose == 1) then
        if (at == 0) then
          write(*, '(a, i0)') "  # of inner iterations: ", nin
        else
          write(*, '(a, i0)') "  Maximum # of inner iterations: ", nin
        endif
      endif
    else if (index(argv(1:6), "--omg=") == 1) then  ! Acceleration (relaxation) parameter 
      read(argv(7:len_trim(argv)), '(f16.10)', iostat=ierror) omg
      if (ierror /= 0) then          
        error stop "option --omg="
      endif
      set_omg = 1
      if (verbose == 1) then
        if (at == 0) then
          write(*, '(a, f0.2)') "  Acceleration parameter (omega): ", omg
        else
          write(*, '(a, f0.2)') "  Initial acceleration parameter (omega): ", omg
          endif
      endif
    else if (index(argv(1:6), "--tol=") == 1) then  ! Stopping criterion
      read(argv(7:len_trim(argv)), '(e16.10e2)', iostat=ierror) tol
      if (ierror /= 0) then          
        error stop "option --tol="
      endif
      set_tol = 1
      if (verbose == 1) then
        write(*, '(a, es8.2e2)') "  Stopping criterion: ", tol
      endif
    else if (index(argv(1:7), "--omax=") == 1) then ! Maximum # of outer iterations
       read(argv(8:len_trim(argv)), '(i10)', iostat=ierror) omax
       if (ierror /= 0) then          
          error stop "option --omax="
       endif
       set_omax = 1
       if (verbose == 1) then
        write(*, '(a, i0)') "  Maximum # of outer iterations: ", omax
      endif
    else if (index(argv(1:7), "--rmax=") == 1) then ! Maximum # of restart cycles
      read(argv(8:len_trim(argv)), '(i10)', iostat=ierror) rmax
      if (ierror /= 0 ) then          
        error stop "option --rmax="
      endif
      set_rmax = 1
      if (verbose == 1) then
        write(*, '(a, i0)') "  Maximum # of restarts: ", rmax
      endif
    else if (index(argv(1:2), "-v") == 1) then          
    else if (index(argv(1:12), "--directory=") == 1) then      
      read(argv(13:len_trim(argv)), '(a)', iostat=ierror) directory
      if (ierror /= 0 ) then          
        error stop "option --directory="
      endif
      set_fi = 1
      if (verbose == 1) then
        write(*, '(a, a)') "  Directory: ", trim(directory)
      endif
    else
      error stop "Unknown argument"
    endif
  enddo  

! Set the default values of parameters (if not set above)
  if (set_at == 0) then
    at = 1
    if (verbose == 1) then
      write(*, *) "  Enebled automatic paramete tuning (default)"
    endif
  endif

  if (set_nin == 0) then
    nin = 50
    if (verbose == 1) then
      write(*, *) "  Maximum # of inner iterations (default): ", nin
    endif
  endif

  if (set_omg == 0) then
    omg = 1.0d0
    if (verbose == 1) then
      write(*, *) "  Initial Value of omega (default): ", omg
    endif
  endif    

  if (set_tol == 0) then
    tol = 1.0d-8
    if (verbose == 1) then
      write(*, *) "  Stoppin  g criterion (default)", tol
    endif
  endif

  if (set_omax == 0) then
    omax = 800
    if (verbose == 1) then
      write(*, *) "  Maximum # of outer iterations (default): ", omax
    endif
  endif

  if (set_rmax == 0) then
    rmax = 0
    if (verbose == 1) then
      write(*, *) "  Maximum # of restarts (default): ", rmax
    endif
  endif  

  if (set_fi == 0) then
    directory = "RANDL7"
    if (verbose == 1) then
      write(*, *) "  fi: ", directory
    endif
  endif    

end subroutine read_prm
!---------------------------------------------------------------------------

subroutine read_mat(directory, AC, ia, jp, m, n, x0, b)  
  implicit none

  character(*), intent(in) :: directory
  real(real64), allocatable, intent(out) :: AC(:), b(:), x0(:)

  integer, allocatable, intent(out) :: ia(:), jp(:)
  integer, intent(inout) :: m, n
  integer is, j, l, nnz
  integer :: fi_AC = 14, fi_jp = 15, fi_ia = 16

  open(fi_AC, file = trim(directory)//'/AC.ccs', action = 'read', iostat = is, status="old")
  if (is /= 0) error stop 'cannot open file AC.ccs'
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
  open(fi_jp, file = trim(directory)//'/jp.ccs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open jp file'
  read(fi_jp, *) (jp(j), j = 1, n+1)
  close(fi_jp)

! Load ia
  open(fi_ia, file = trim(directory)//'/ia.ccs', action = 'read', iostat = is, status="old")
  if (is /= 0) stop 'cannot open ia file'
  read(fi_ia, *) (ia(l), l = 1, nnz)
  close(fi_ia)

!	b : = random
  call random_seed
  call random_number(b(1:m))
  b(1:m) = 2.0d0*b(1:m) - one

end subroutine read_mat


subroutine output(AC, ia, jp, m, n, &
                  b, &
                  nin, omg, &
                  x, &
                  iter, omax, riter, relres, t_tot, &
                  verbose)
  use func  
  implicit none

  real(real64), intent(in) :: AC(:), b(:)  
  real(real64), intent(in) :: t_tot, x(n)
  real(real64), intent(inout) :: relres(iter)
  real(real64) r(m), ATr(n)
  real(real64), intent(in) :: omg
  real(real64) tmp, nrm_b, nrm_r, nrmATb, nrmATr

  integer, intent(in) :: ia(:), jp(:)
  integer, intent(in) :: iter, m, n, nin, omax, verbose, riter
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

  if (verbose == 1) then

    write(*, '(a)') 'Results:'
    write(*, '(a, f0.2)') '  omega: ', omg
    write(*, '(a, i0)') '  # of outer iterations: ', iter
    write(*, '(a, i0)') '  # of inner iterations: ', nin
    write(*, '(a, i0)') '  # of restarts: ', riter
    write(*, '(a, i0)') '  Restart cycle: ', omax
    write(*, '(a, f0.2)') '   CPU time: ', t_tot
    write(*, '(a, es8.2e2)') '  Relative residual: ', relres(iter)
    write(*, '(a, es8.2e2)') '  Actual relative residual (ATr): ', nrmATr / nrmATb
    write(*, '(a, es8.2e2)') '  Actual relative residual (r): ', nrm_r / nrm_b

  else

    write(*, '(a)') 'Results (omega time iter):'
    write(*, '(f0.2, a, f0.2, a, i0)') omg, ' ', t_tot, ' ', iter
   ! w rite(logcsv, '(i0, a, i0, a, f0.2, a, f0.2)') iter, ',', nin, ',', omg, ',', t_tot

  endif

  open(info, file='info.dat', action='write', iostat=is, status='replace')
  if (is /= 0) stop 'cannot open info.dat file'
  write(info, '(a, f0.2)') '  omega: ', omg
  write(info, '(a, i0)') '  # of outer iterations: ', iter
  write(info, '(a, i0)') '  # of inner iterations: ', nin
  write(info, '(a, i0)') '  # of restarts: ', riter
  write(info, '(a, i0)') '  Restart cycle: ', omax
  write(info, '(a, f0.2)') '   CPU time: ', t_tot
  write(info, '(a, es8.2e2)') '  Relative residual: ', relres(iter)
  write(info, '(a, es8.2e2)') '  Actual relative residual (ATr): ', nrmATr / nrmATb
  write(info, '(a, es8.2e2)') '  Actual relative residual (r): ', nrm_r / nrm_b
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
  implicit none
  real(real64), allocatable, intent(inout) :: AC(:), b(:)
  integer,  allocatable, intent(inout) :: ia(:), jp(:)

  deallocate(AC, jp, ia)
  deallocate(b)

end subroutine discard
end module sub