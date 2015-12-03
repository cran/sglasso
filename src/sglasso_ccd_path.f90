!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!						SUBROUTINE USED TO COMPUTE THE SGLASSO ESTIMATOR BY THE CCD ALGORITHM (PATH VERSION)
!
!	AUTHOR: Luigi Augugliaro
!	VERSION: 1.2-2
!	DATA: 01/12/2015
!
!	DESCRIPTION: the subroutine sglasso_ccd_path implements the cyclic coordinate descent algorithm proposed
!						Augugliaro et al., (2014) to fit a weighted l1-norm penalized RCON(V, E) model. 
!
!	ARGUMENTS: see arguments in sglasso_ccm_path.f90
subroutine sglasso_ccd_path(nSv,Sv,nTv,Tv_pkg,Tv_rw,nv,Tv_ptr,nTe,Te,nTe_ptr,Te_ptr,ne,Te_scn,Te_ptr_scn,nstep,trnc,tol, &
	rho,nrho,min_rho,grd,th,w,pnl_flg,n,conv)
	integer	:: nSv,nTv,Tv_pkg(nTv),Tv_rw(nTv),nv,Tv_ptr(nv+1),nTe,Te(nTe),nTe_ptr,Te_ptr(nTe_ptr),ne, &
		Te_scn(ne+1),Te_ptr_scn(ne+1),nstep,nrho,pnl_flg(ne),n,conv
	double precision	:: Sv(nSv),trnc,tol,rho(nrho),min_rho,grd(nv+ne,nrho),th(nv+ne,nrho),w(ne)
	!internal variables
	logical	:: A(ne),flg
	integer	:: p,i,j,h,k,i0,i1,j0,j1,h0,h1,Te_scn_sz(ne), &
		Te_cl_ptr(ne+1),Te_cl_sz(ne),Te_cl(ne*(Tv_ptr(nv+1)-1)),nr,df,max_df
	double precision	:: Svh(nSv),tr_SvTv(nv),tr_SvTe(ne),max_rho,Du,oth(nv+ne),ath(nv+ne), &
		dth,tr_SvhTeSvhTe,tr_SvhTvSvhTv,th_diff,Du_sgn,cf,arho,agrd(nv+ne)
		
	do i = 1, ne
   	if (pnl_flg(i).eq.0) then
			A(i) = .true.
      else 
			A(i) = .false.
		end if
	end do
	Svh = 0.
	ath = 0.
	oth = 0.
	agrd = 0.
	p = Tv_ptr(nv + 1) - 1
	df = nv + ne - sum(pnl_flg)
	max_df = nv + ne
	
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
		ath(i) = (i1 - i0 + 1)/tr_SvTv(i)
		Svh(Tv_pkg(i0:i1)) = 1/ath(i)
	end do
			
	max_rho = 0.
	do i = 1, ne
		i0 = Te_scn(i)
		i1 = Te_scn(i + 1) - 1
		j0 = Te_ptr_scn(i)
		j1 = Te_ptr_scn(i + 1) - 1
		call trSTe(nSv,Sv,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1),tr_SvTe(i))
		Du = - tr_SvTe(i)
		agrd(nv + i) = Du
		if(pnl_flg(i).eq.1) max_rho = max(max_rho,abs(Du) / w(i))
	end do

	if(min_rho.gt.max_rho) return

	rho(1) = max_rho
	cf = dexp((dlog(min_rho) - dlog(max_rho)) / (nrho - 1))
	do i = 2, nrho
		rho(i) = rho(i - 1) * cf
	end do
	
	n = 0
	do nr = 1, nrho

		arho = rho(nr)
		do i = 1, ne
			if(.not.A(i).and.pnl_flg(i).eq.1) then
				if(abs(agrd(nv + i))/w(i).gt.arho) then
					A(i) = .true.
					df = df + 1
				end if
			end if
		end do

		do
			oth(1:nv) = ath(1:nv)
			do i = 1, ne
				if(A(i))	oth(nv + i) = ath(nv + i)
			end do

			if(df.gt.nv) then
				do i = 1, ne
					if(A(i)) then
               	n = n + 1
                  if(n.ge.nstep) then
                  	if(nr.eq.1) then
                     	nrho = nr
                     else
                     	nrho = nr - 1
                     end if
                     conv = 1
                     return
						end if
						i0 = Te_scn(i)
						i1 = Te_scn(i + 1) - 1
						j0 = Te_ptr_scn(i)
						j1 = Te_ptr_scn(i + 1) - 1
						h0 = Te_cl_ptr(i)
						h1 = Te_cl_ptr(i + 1) - 1
						call trSTe(nSv,Svh,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1),Du)
						Du = Du - tr_SvTe(i)
						agrd(nv + i) = Du
   	         	if(pnl_flg(i).eq.1) then
   	            	Du_sgn = dsign(dble(1),Du)
							Du = Du - arho * Du_sgn * w(i)
						end if
						call trSTeSTe(nSv,Svh,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1), &
						Te_cl_sz(i),Te_cl(h0:h1),tr_SvhTeSvhTe)
						dth = Du / tr_SvhTeSvhTe
						if(dth.ne.dth) then
							conv = 2
							return
						end if
						ath(nv + i) = ath(nv + i) + dth
						call updateSvh_e(p,nSv,Svh,dth,Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1))
					end if
				end do
			end if

			do i = 1, nv
				n = n + 1
				if(n.ge.nstep) then
					if(nr.eq.1) then
						nrho = nr
					else
						nrho = nr - 1
					end if
					conv = 1
					return
				end if
				i0 = Tv_ptr(i)
				i1 = Tv_ptr(i + 1) - 1
				Du = sum(Svh(Tv_pkg(i0:i1))) - tr_SvTv(i)
				agrd(i) = Du
				tr_SvhTvSvhTv = 0.
				do h = i0 + 1, i1
					do k = i0, h - 1
						tr_SvhTvSvhTv = tr_SvhTvSvhTv + Svh(int(Tv_rw(k) + Tv_rw(h) * (Tv_rw(h) - 1) / 2))**2
					end do
				end do
				tr_SvhTvSvhTv = 2 * tr_SvhTvSvhTv + sum(Svh(Tv_pkg(i0:i1))**2)
				dth = Du / tr_SvhTvSvhTv
				if(dth.ne.dth) then
					conv = 2
					return
				end if
				ath(i) = ath(i) + dth
				call updateSvh_v(p,nSv,Svh,dth,Tv_ptr(i+1)-Tv_ptr(i),Tv_pkg(i0:i1),Tv_rw(i0:i1))
			end do

			th_diff = sum(abs(ath(1:nv) - oth(1:nv)))
			if(th_diff/df.lt.tol) then
				do i = 1, ne
					if(A(i)) then
						th_diff = th_diff + abs(ath(nv + i) - oth(nv + i))
					end if
				end do
			end if

			if(th_diff/df.lt.tol) then
				flg = .true.
				do i = 1, ne
					if(A(i).and.abs(ath(nv + i)).gt.trnc.and.pnl_flg(i).eq.1) then
						if((ath(nv + i) * agrd(nv + i)).lt.0.) then
							flg = .false.
							call updateSvh_e(p,nSv,Svh,-ath(nv + i),Te_scn_sz(i),Te(i0:i1),Tv_ptr(nv+1),Te_ptr(j0:j1))
							ath(nv + i) = 0.
							A(i) = .false.
							df = df - 1
						end if
					end if
				end do
				if(flg) then
					do i = 1, ne
						if(.not.A(i)) then
							call trSTe(nSv,Svh,Te_scn_sz(i),Te(Te_scn(i):(Te_scn(i + 1) - 1)), &
								Tv_ptr(nv+1),Te_ptr(Te_ptr_scn(i):(Te_ptr_scn(i + 1) - 1)),Du)
							Du = Du - tr_SvTe(i)
							agrd(nv + i) = Du
							if(abs(Du)/w(i) - arho .gt. 0.) then
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

		th(:, nr) = ath
		grd(:, nr) = agrd
		
	end do
end subroutine sglasso_ccd_path
