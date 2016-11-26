program main
	use globvar
	use solver
	implicit none

	real(dp), allocatable :: RelRes(:), x(:)
	real(dp) t1, t2, t_cri, t_tot
	integer :: logcsv = 20, Iter, Riter

! Read parameters
	call read_prm()

! Read matrix data
	call read_mat()

	allocate(x(n))
	allocate(RelRes(omax*(rmax+1)))

	open(logcsv, file='log.csv', position='append')

	t_tot = zero
	t_cri = zero
	write(*, *) 'Executing'
	call cpu_time(t1)
	call BAGMRES(Iter, Riter, RelRes, x)
	call cpu_time(t2)
	write(*, *) 'Terminated'
	t_tot = t_tot + t2 - t1
	write(*, *) 'Outputting results'
	call output(Iter, Riter, t_tot, RelRes, x)
	write(logcsv, '(i9, a1, i2, a1, f4.1, a1, f16.4)') Iter, ',', nin, ',', omg, ',', t_tot

	close(logcsv)

	call discard()

end program main
