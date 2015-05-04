subroutine gdf_fun(n,p,X,S,nrho,Kh,gdf,info)
	integer	::n,p,nrho,info
	double precision	:: X(n,p),S(p,p),Kh(p,p,nrho),gdf(nrho)
	!internal variables
	integer	:: i,j,ii,jj,k,nr
	logical	:: A(p,p,nrho)
	double precision	::	Sk(p,p,n),S_Sk(p,p,n),gdf_partial,Sh(p,p)
	do nr = 1, nrho
		do j = 1, p
			do i = 1, j
				if(abs(Kh(i,j,nr)).ne.0.) then
					A(i,j,nr) = .true.
					A(j,i,nr) = .true.
					else 
						A(i,j,nr) = .false.
						A(j,i,nr) = .false.
				end if
			end do	
		end do
	end do
	do k = 1, n
		do j = 1, p
			do i = 1, j
				if(any(A(i,j,:))) then
					Sk(i,j,k) = X(k,i) * X(k,j)
					S_Sk(i,j,k) = S(i,j) - Sk(i,j,k)
					Sk(j,i,k) = Sk(i,j,k)
					S_Sk(j,i,k) = S_Sk(i,j,k)
					else
						Sk(i,j,k) = 0.
						S_Sk(i,j,k) = 0.
						Sk(j,i,k) = 0.
						S_Sk(j,i,k) = 0.
				end if
			end do
		end do
	end do
	do nr = 1, nrho
		gdf(nr) = 0.
		Sh = Kh(:,:,nr)
		call dpotrf('U',p,Sh,p,info)
		if(info.ne.0) return
		call dpotri('U',p,Sh,p,info)
		if(info.ne.0) return
		do j = 1, p
			do i = 1, j - 1
				Sh(j,i) = Sh(i,j)
			end do
		end do
		do k = 1, n
			do j = 1, p
				do i = 1, p
					if(A(i,j,nr)) then
						gdf_partial = 0.
						do jj = 1, p
							do ii = 1, p
								if(A(ii,jj,nr)) then
									gdf_partial = gdf_partial + Kh(i,jj,nr) * S_Sk(jj,ii,k) * Kh(ii,j,nr)
								end if
							end do
						end do
						gdf(nr) = gdf(nr) + (Sh(i,j) - Sk(i,j,k)) * gdf_partial
					end if
				end do
			end do
		end do		
		gdf(nr) = gdf(nr) / (n - 1)
	end do
end subroutine gdf_fun
