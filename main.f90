program main
  use,intrinsic :: iso_fortran_env
  use floating_point_kinds, only : zero, one
  use sub
  use solver
  implicit none

  real(real64), allocatable :: b(:), relres(:), x(:), x0(:), AC(:)
  integer,  allocatable :: jp(:), ia(:)
  real(real64) omg, t1, t2, t_cri, t_tot, tol
  integer :: at, err, iter, logcsv = 20, m, n, nin, omax, output_mode, riter, rmax

! Read parameters  
  write(*, *) "Started setting the values of parameters"
  call read_prm(nin, omg, omax, output_mode, rmax, tol, at)
  write(*, *) "Finished setting the values of parameters"
  write(*, *)

! Read matrix data
  call read_mat(AC, ia, jp, m, n, x0, b)

  allocate(x(n), stat=err)
  if (err /= 0) then
    error stop "allocate x"
  endif
  allocate(relres(omax*(rmax+1)))

  if (err /= 0) then
    error stop "allocate relres"
  endif
  
  open(logcsv, file='log.csv', position='append')

  t_tot = zero
  t_cri = zero

  write(*, *) 'Executing BA-GMRES'
  call cpu_time(t1)
  
  call BAGMRES(AC, ia, jp, m, n, &
                x0, b, &
                at, nin, omg, &
                tol, omax, rmax, &
                x, &
                iter, riter, relres)

  call cpu_time(t2)
  write(*, *) 'BA-GMRES terminated'
  write(*, *) 

  t_tot = t_tot + t2 - t1 ! total CPU Time

  write(*, *) 'Outputting results'
  call output(AC, ia, jp, m, n, &
              b, &
              nin, omg, & 
              x, &
              iter, omax, riter, relres, t_tot, &
              output_mode)
  write(logcsv, '(i9, a1, i2, a1, f4.1, a1, f16.4)') iter, ',', nin, ',', omg, ',', t_tot
  close(logcsv)

  call discard(AC, ia, jp, b)

end program main