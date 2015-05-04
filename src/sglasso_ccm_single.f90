!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						SUBROUTINE USED TO COMPUTE THE SGLASSO ESTIMATOR BY THE CCM ALGORITHM (SINGLE VERSION)
!
!	AUTHOR: Luigi Augugliaro
!	VERSION: 1.0-0
!	DATA: 23/04/2014
!
!	DESCRIPTION: the subroutine sglasso_ccm_single implements the cyclic coordinate minimization algorithm proposed
!						Augugliaro et al., (2014) to fit a weighted l1-norm penalized RCON(V, E) model. 
!
!	ARGUMENTS: see arguments in sglasso_ccm_path.f90
subroutine sglasso_ccm_single(nSv,Sv,nTv,Tv_pkg,Tv_rw,nv,Tv_ptr,nTe,Te,nTe_ptr,Te_ptr,ne,Te_scn,Te_ptr_scn, &
	nstep,trnc,tol,rho,grd,th,w,n,conv)
	integer	:: nSv,nTv,Tv_pkg(nTv),Tv_rw(nTv),nv,Tv_ptr(nv+1),nTe, Te(nTe),nTe_ptr,Te_ptr(nTe_ptr),ne, &
		Te_scn(ne+1),Te_ptr_scn(ne+1),nstep,n,conv
	double precision	:: grd(nv+ne),Sv(nSv),trnc,tol,rho,th(nv+ne),w(ne)
	!internal variables
	logical	:: A(ne),flg
	integer	:: p,i,j,h,k,i0,i1,j0,j1,h0,h1,Te_scn_sz(ne), &
		Te_cl_ptr(ne+1),Te_cl_sz(ne),Te_cl(ne*(Tv_ptr(nv+1)-1)),df,max_df
	double precision	:: Svh(nSv),tr_SvTv(nv),tr_SvTe(ne),max_rho,Du,oth(nv+ne), &
		dth,tr_SvhTeSvhTe,tr_SvhTvSvhTv,th_diff,Du_sgn
	
	A = .false.
	Svh = 0.
	oth = 0.
	df = nv
	max_df = nv + ne
	p = Tv_ptr(nv + 1) - 1
	
	k = 0
	Te_cl_ptr(1) = 1
	do i = 1, ne
		Te_scn_sz(i) = Te_scn(i + 1) - Te_scn(i)
		j = 0
		do h = Te_ptr_scn(i), Te_ptr_scn(i + 1) - 2
			j = j + 1
			if(Te_ptr(h).ne.Te_ptr(h + 1)) then
				k = k + 1
				Te_cl(k) = j
			end if
		end do
		Te_cl_ptr(i + 1) = k + 1
		Te_cl_sz(i) = Te_cl_ptr(i + 1) - Te_cl_ptr(i)
	end do
	
	do i = 1, nv
		i0 = Tv_ptr(i)
		i1 = Tv_ptr(i + 1) - 1
		tr_SvTv(i) = sum(Sv(Tv_pkg(i0:i1)))
		th(i) = (i1 - i0 + 1)/tr_SvTv(i)
		Svh(Tv_pkg(i0:i1)) = 1/th(i)
	end do
	
	max_rho = 0.
	do i = 1, ne
		i0 = Te_scn(i)
		i1 = Te_scn(i + 1) - 1
		j0 = Te_ptr_scn(i)
		j1 = Te_ptr_scn(i + 1) - 1
		call trSTe(nSv,Sv,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1),tr_SvTe(i))
		Du = - tr_SvTe(i)
		grd(nv + i) = Du
		max_rho = max(max_rho,abs(Du) / w(i))
		if(abs(Du) / w(i).ge.rho) then
			A(i) = .true.
			df = df + 1
		end if
	end do
	
	if(rho.ge.max_rho) return
	
	n = 0
	do
		oth(1:nv) = th(1:nv)
		do i = 1, ne
			if(A(i))	oth(nv + i) = th(nv + i)
		end do
		
		if(df.gt.nv) then
			do i = 1, ne
				if(A(i)) then
					i0 = Te_scn(i)
					i1 = Te_scn(i + 1) - 1
					j0 = Te_ptr_scn(i)
					j1 = Te_ptr_scn(i + 1) - 1
					h0 = Te_cl_ptr(i)
					h1 = Te_cl_ptr(i + 1) - 1
					do
						n = n + 1
						if(n.ge.nstep) then
							conv = 1
							return
						end if					
						call trSTe(nSv,Svh,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1),Du)
						Du = Du - tr_SvTe(i)
						grd(nv + i) = Du
						Du_sgn = dsign(dble(1),Du)
						Du = Du - rho * Du_sgn * w(i)
                        if(abs(Du).le.tol) exit
						call trSTeSTe(nSv,Svh,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1), &
							Te_cl_sz(i),Te_cl(h0:h1),tr_SvhTeSvhTe)
						dth = Du / tr_SvhTeSvhTe
                        if(isnan(dth)) then
                            conv = 2
                            return
                        end if
						th(nv + i) = th(nv + i) + dth
						call updateSvh_e(p,nSv,Svh,dth,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1))
					end do
				end if
			end do
		end if
		
		do i = 1, nv
			i0 = Tv_ptr(i)
			i1 = Tv_ptr(i + 1) - 1
			do
				n = n + 1
				if(n.ge.nstep) then
					conv = 1
					return
				end if
				Du = sum(Svh(Tv_pkg(i0:i1))) - tr_SvTv(i)	
				grd(i) = Du
				if(abs(Du).le.tol) exit		
				tr_SvhTvSvhTv = 0.
				do h = i0 + 1, i1
					do k = i0, h - 1
						tr_SvhTvSvhTv = tr_SvhTvSvhTv + Svh(int(Tv_rw(k) + Tv_rw(h) * (Tv_rw(h) - 1) / 2))**2
					end do
				end do
				tr_SvhTvSvhTv = 2 * tr_SvhTvSvhTv + sum(Svh(Tv_pkg(i0:i1))**2)
				dth = Du / tr_SvhTvSvhTv
                if(isnan(dth)) then
                    conv = 2
                    return
                end if
				th(i) = th(i) + dth
				call updateSvh_v(p,nSv,Svh,dth,Tv_ptr(i+1)-Tv_ptr(i),Tv_pkg(i0:i1),Tv_rw(i0:i1))
			end do
		end do
		
		th_diff = sum(abs(th(1:nv) - oth(1:nv)))
		if(th_diff/df.lt.tol) then
			do i = 1, ne
				if(A(i)) then
					th_diff = th_diff + abs(th(nv + i) - oth(nv + i))
				end if
			end do
		end if
		
		if(th_diff/df.lt.tol) then
			flg = .true.
            do i = 1, ne
                if(A(i).and.abs(th(nv + i)).gt.trnc) then
                    if((th(nv + i) * grd(nv + i)).lt.0.) then
                        flg = .false.
                        call updateSvh_e(p,nSv,Svh,-th(nv + i),Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1))
                        th(nv + i) = 0.
                        A(i) = .false.
                        df = df - 1
                    end if
                end if
            end do
            if(flg) then
                do i = 1, ne
                    if(.not.A(i)) then
                        call trSTe(nSv,Svh,Te_scn_sz(i),Te(Te_scn(i):(Te_scn(i + 1) - 1)), Tv_ptr(nv+1), &
                            Te_ptr(Te_ptr_scn(i):(Te_ptr_scn(i + 1) - 1)),Du)
                        Du = Du - tr_SvTe(i)
                        grd(nv + i) = Du
                        if(abs(Du)/w(i) - rho .gt. 0.) then
                            flg = .false.
                            A(i) = .true.
                            df = df + 1
                        end if
                    end if
                end do
            end if
			if(flg) exit
		end if
	end do
	
end subroutine sglasso_ccm_single
