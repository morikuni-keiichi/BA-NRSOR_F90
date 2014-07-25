subroutine read_prm()
	use globvar
	implicit none

	integer :: fi = 11, is
	
	open(fi, file='prm.dat', action='read', iostat = is, status="old")
	if (is /= 0) stop 'cannot open a parameter file'
	read(fi, *) cri_res 
	read(fi, *) omg
	read(fi, *) omax 
	read(fi, *) nin 
	read(fi, *) rmax
	read(fi, *) at
	read(fi, *) omode
	close(fi)

	end subroutine read_prm
!---------------------------------------------------------------------------

subroutine read_mat()
	use globvar
	implicit none

	real(dp) :: two = 2.0d0
	integer is, j, l
	integer :: fi_AC = 14, fi_jp = 15, fi_ia = 16

	open(fi_AC, file = 'RANDL7/AC.ccs', action = 'read', iostat = is, status="old")
	if (is /= 0) stop 'cannot open AC file'
	read(fi_AC, *) m, n, nnz

	allocate(AC(nnz))
	allocate(jp(n+1))
	allocate(ia(nnz))
	allocate(b(m))
	allocate(Aei(n))

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
	b(1:m) = two*b(1:m) - one

end subroutine read_mat 


subroutine output(Iter, Riter, t_tot, RelRes, x)
	use func
	use globvar
	implicit none

	real(dp), intent(in) :: t_tot, x(n)
	real(dp), intent(inout) :: RelRes(Iter)
	real(dp) r(m), ATr(n)
	real(dp) tmp, norm_b, norm_r, norm_ATb, norm_ATr
	integer, intent(in) :: Iter, Riter	
	integer :: info = 10, reshis=11, sol = 12
	integer i, is, j, l, k, k1, k2

	norm_b = nrm2(b(1:m), m)

	do j = 1, n
		k1 = jp(j)
		k2 = jp(j+1)-1
		ATr(j) = sum(AC(k1:k2) * b(ia(k1:k2)))
	enddo
	
	norm_ATb = nrm2(ATr(1:n), n)

	r(1:m) = zero
	do j = 1, n
		tmp = x(j)
		do l = jp(j), jp(j+1)-1
			i = ia(l)
			r(i) = r(i) + tmp*AC(l)
		enddo
	enddo

	r(1:m) = b(1:m) - r(1:m)	

	norm_r = nrm2(r(1:m), m) 
	
	do j = 1, n
		k1 = jp(j)
		k2 = jp(j+1)-1
		ATr(j) = sum(AC(k1:k2) * r(ia(k1:k2)))
	enddo

	norm_ATr = nrm2(ATr(1:n), n)

	tmp = one / norm_ATb
	RelRes(1:Iter) = tmp * RelRes(1:Iter)

	if (omode == 0) then

		write(*, '(a, f6.3)') ' omega: ', omg
		write(*, *) '# of outer iterations: ', Iter
		write(*, *) '# of inner iterations: ', nin	
		write(*, *) '# of restarts: ', Riter	
		write(*, *) 'Restart cycle: ', omax
		write(*, '(a, f16.5)') ' CPU time: ', t_tot
		write(*, *) 'Relative residual: ', RelRes(Iter)
		write(*, *) 'Actual relative residual (ATr): ', norm_ATr / norm_ATb
		write(*, *) 'Actual relative residual (r): ', norm_r / norm_b

	else

		write(*, '(f6.3, f10.5, i4)') omg, t_tot, Iter
	
	endif

	open(info, file='info.dat', action='write', iostat=is, status='replace')
	if (is /= 0) stop 'cannot open info.dat file'
	write(info, '(a, f6.3)') ' omega: ', omg
	write(info, *) '# of outer iterations: ', Iter
	write(info, *) '# of inner iterations: ', nin
	write(info, *) '# of restarts: ', Riter	
	write(info, *) 'Restart cycle: ', omax
	write(info, '(a, f16.5)') ' CPU time: ', t_tot
	write(info, *) 'Relative residual: ', RelRes(Iter)
	write(info, *) 'Actual relative residual (ATr): ', norm_ATr / norm_ATb	
	write(info, *) 'Actual relative residual (r): ', norm_r / norm_b	
	close(info)
	
	open(reshis, file='reshis.dat', action='write', iostat=is, status='replace')
	if (is /= 0) stop 'cannot open reshis.dat file'
	do k = 1, Iter
		write(reshis, *) k, RelRes(k)
	enddo
	close(reshis)
	
	open(sol, file='solution.dat', action='write', iostat=is, status='replace')
	if (is /= 0) stop 'cannot open solution.dat file'
	do j = 1, n
		write(sol, *) x(j)
	enddo
	close(sol)

end subroutine output

subroutine discard()
	use globvar
	implicit none

	deallocate(AC, jp, ia)
	deallocate(b)
	deallocate(Aei)

end subroutine discard
