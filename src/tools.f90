subroutine trSTe(nSv,Sv,nTe,Te,nTe_ptr,Te_ptr,tr_STe)
	double precision	:: tr_STe
	integer	:: nSv,nTe,Te(nTe),nTe_ptr,Te_ptr(nTe_ptr)
	double precision	:: Sv(nSv)
	! internal variables
	integer	:: j,h,k
	tr_STe = 0.
	do j = 1, nTe_ptr - 1
		k = j * (j - 1) / 2
		do h = Te_ptr(j), Te_ptr(j + 1) - 1
			if(Te(h).le.j)	tr_STe = tr_STe + Sv(Te(h) + k)
		end do
	end do
	tr_STe = 2. * tr_STe
end subroutine trSTe

subroutine trSTeSTe(nSvh,Svh,nTe,Te,nTe_ptr,Te_ptr,ncl,cl,tr_SvhTeSvhTe)
	integer	:: nSvh,nTe,Te(nTe),nTe_ptr,Te_ptr(nTe_ptr),ncl,cl(ncl)
	double precision	:: Svh(nSvh),tr_SvhTeSvhTe
	!internal variables
	integer	:: i,j,h,k
	double precision	:: B(ncl**2)
	B = 0.
	k = 0
	do j = 1, ncl
		do i = 1, ncl
			k = k + 1
			do h = Te_ptr(cl(j)), Te_ptr(cl(j) + 1) - 1
				if(cl(i).le.Te(h)) then
					B(k) = B(k) + Svh(cl(i) + Te(h) * (Te(h) - 1) / 2)
				else 
					B(k) = B(k) + Svh(Te(h) + cl(i) * (cl(i) - 1) / 2)
				end if
			end do
		end do
	end do
	tr_SvhTeSvhTe = 0.
	do j = 1, ncl
		do i = j + 1, ncl
			tr_SvhTeSvhTe = tr_SvhTeSvhTe + B(i + ncl * (j - 1)) * B(j + ncl * (i - 1))
		end do
	end do
	tr_SvhTeSvhTe = 2. * tr_SvhTeSvhTe
	do j = 1, ncl
		tr_SvhTeSvhTe = tr_SvhTeSvhTe + B(j + ncl * (j - 1))**2
	end do
end subroutine trSTeSTe

subroutine updateSvh_v(p,nSvh,Svh,dth,nTv,Tv_pkg,Tv_rw)
	integer	:: p,nSvh,nTv,Tv_pkg(nTv),Tv_rw(nTv)
	double precision	:: Svh(nSvh),dth
	!internal variables
	integer	:: i0,i1,k,i,j,n,jj,ii
	double precision	:: vk, B(nSvh)	
	do k = 1, nTv
		vk = dth / (1. + dth * Svh(Tv_pkg(k)))
		i0 = 1 + Tv_rw(k) * (Tv_rw(k) - 1) / 2
		i1 = Tv_rw(k) * (Tv_rw(k) + 1) / 2
		n = 0
		do j = i0, i1
			do i = i0, j
				n = n + 1
				B(n) =  Svh(i) * Svh(j)
			end do
		end do
		j = i1
		do jj = Tv_rw(k), p - 1
			j = j + jj
			do i = i0, i1
				n = n + 1
				B(n) =  Svh(i) * Svh(j)
			end do
			i = i1
			do ii = Tv_rw(k), jj
				i = i + ii
				n = n + 1
				B(n) =  Svh(i) * Svh(j)				
			end do
		end do
		Svh = Svh - vk * B
	end do
end subroutine updateSvh_v

subroutine updateSvh_e(p,nSv,Svh,dth,nTe,Te,nTe_ptr,Te_ptr)
	integer	:: p,nSv,nTe,Te(nTe),nTe_ptr,Te_ptr(nTe_ptr)
	double precision	:: Svh(nSv),dth
	!internal variables
	integer	:: i,j,k
	double precision	:: Sh(p,p),v1(p),v2(p),vk
	k = 0
	do j = 1, p
		do i = 1, j - 1
			k = k + 1
			Sh(i,j) = Svh(k)
			Sh(j,i) = Svh(k)
		end do
		k = k + 1
		Sh(j,j) = Svh(k)
	end do
	do k = 1, p
		if(Te_ptr(k).ne.Te_ptr(k + 1)) then
			v1 = 0.
			do j = Te_ptr(k), Te_ptr(k + 1) - 1
				v1 = v1 + Sh(:,Te(j))
			end do
			vk = dth / (1. + dth * v1(k))
			v2 = Sh(k,:)
			do j = 1, p
				do i = 1, p
					Sh(i,j) = Sh(i,j) - vk * v1(i) * v2(j)
				end do
			end do
		end if
	end do
	k = 0
	do j = 1, p
		do i = 1, j
			k = k + 1
			Svh(k) = Sh(i,j)
		end do 
	end do
end subroutine updateSvh_e