!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine skubo_w(nR, norb, norb_ex, nv, nc, nRvec, R, B, hhop, shop, npointstotal, rkx, &
rky, rkz, fk_ex, e_ex)

implicit real*8 (a-h,o-z)

dimension nRvec(nR,3)
dimension R(3,3)
dimension B(norb,3)
dimension hhop(norb,norb,nR)
dimension shop(norb,norb,nR)
dimension rkx(npointstotal)
dimension rky(npointstotal)
dimension rkz(npointstotal)
dimension fk_ex(norb_ex,norb_ex)
dimension e_ex(norb_ex)

dimension rhop(3,nR,norb,norb)
dimension sderhop(3,nR,norb,norb)
dimension hderhop(3,nR,norb,norb)

dimension hkernel(norb,norb)
dimension skernel(norb,norb)
dimension hderkernel(3,norb,norb)
dimension sderkernel(3,norb,norb)
dimension akernel(3,norb,norb)
dimension pgaugekernel(3,norb,norb)

dimension hk_ev(norb,norb)
dimension e(norb)
dimension pgauge(3,norb,norb)
dimension vjseudoa(3,norb,norb) 
dimension vjseudob(3,norb,norb)
dimension vme(3,norb,norb)

allocatable vme_ex(:,:,:)

allocatable wp(:)
allocatable wn_sp(:)
allocatable wn_ex(:)
allocatable skubo_sp_int(:,:,:)
allocatable skubo_ex_int(:,:,:)
allocatable sigma_w_sp(:,:,:)
allocatable sigma_w_ex(:,:,:)

complex*16 sderhop
complex*16 hderhop
complex*16 fk_ex
complex*16 hkernel
complex*16 skernel
complex*16 hderkernel
complex*16 sderkernel
complex*16 akernel
complex*16 pgaugekernel
complex*16 hk_ev
complex*16 pgauge
complex*16 vjseudoa
complex*16 vjseudob
complex*16 vme

complex*16 vme_ex 

complex*16 skubo_sp_int
complex*16 skubo_ex_int 
complex*16 sigma_w_sp
complex*16 sigma_w_ex

character(100) type_broad
character(100) material_name
character(100) file_name_sp
character(100) file_name_ex
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_kubo_parameters(w0,wrange,nw,eta,type_broad,material_name, &
file_name_sp,file_name_ex)
eta=eta/27.211385d0
call crossproduct(R(1,1),R(1,2),0.0d0,R(2,1), &
R(2,2),0.0d0,cx,cy,cz)         
vcell=sqrt(cx**2+cy**2+cz**2)

nkubo_sp=norb**2*npointstotal
nkubo_ex=norb_ex

allocate (wp(nw))
allocate (wn_sp(nkubo_sp))
allocate (wn_ex(nkubo_ex))
allocate (vme_ex(3,nkubo_ex,2))
allocate (skubo_sp_int(3,3,nkubo_sp))
allocate (skubo_ex_int(3,3,nkubo_ex))
allocate (sigma_w_sp(3,3,nw))
allocate (sigma_w_ex(3,3,nw))

wp=0.0d0
wn_sp=0.0d0
wn_ex=0.0d0
skubo_sp_int=0.0d0
skubo_ex_int=0.0d0
sigma_w_sp=0.0d0
sigma_w_ex=0.0d0
vme_ex=0.0d0
!Brillouin zone sampling
kcount=1


call hoppings_observables_TB(norb,nR,nRvec,R,B,shop,hhop,rhop,sderhop,hderhop)
do ibz=1,npointstotal
	write(*,*) 'point:',ibz,npointstotal         
	rkxp=rkx(ibz)
	rkyp=rky(ibz)
	rkzp=rkz(ibz)
			
	!get matrices in the \alpha, \alpha' basis (orbitals,k)    		
	call get_vme_kernels(rkxp,rkyp,rkzp,nR,nRvec,norb,R, &
	hkernel,skernel,shop,hhop,rhop,sderhop,hderhop,sderkernel,hderkernel,akernel, &
	pgaugekernel)


	!velocity matrix elements
	call get_eigen_vme(rkxp,rkyp,rkzp,nR,nRvec,norb,R, &
	hkernel,skernel,akernel,sderkernel,hderkernel, &
	pgaugekernel,hk_ev,e,pgauge,vjseudoa,vjseudob,vme) 

	write(*,*) ibz,vme(1,1,1)
	
	!pause
	!compute exciton velocity matrix element
	do iex=1, nkubo_ex
		wn_ex(iex)=e_ex(iex)
		do nj=1, 3
			vme_ex(nj,iex,1)=vme_ex(nj,iex,1)+fk_ex(ibz,iex)*vme(nj,1,2)
			vme_ex(nj,iex,2)=vme_ex(nj,iex,2)+fk_ex(ibz,iex)*vme(nj,2,1)
		end do
	end do

	write(*,*) ibz,vme_ex(1,1,1)

	!get strength for kubo SP
	call get_kubo_intens(nv,npointstotal,vcell,norb,nkubo_sp,e,kcount,vme,wn_sp,skubo_sp_int)
end do


!fill kubo oscillators of EXCITONS
do nn=1,nkubo_ex  
	!save oscillator stregths
	do nj=1,3
		do njp=1,3		  
			skubo_ex_int(nj,njp,nn)=1.0d0/(dble(npointstotal)*vcell) &
			*conjg(vme_ex(nj,nn,1))*vme_ex(njp,nn,1)/wn_ex(nn)   !pick the correct order of operators
		end do	
	end do
end do


!broad the delta points of sp and exciton
do ialpha=1,3
	do ialphap=1,3
		call broad_vector(nkubo_sp,wn_sp,skubo_sp_int(ialpha,ialphap,:), &
		w0/27.211385d0,wrange/27.211385d0,nw,wp,sigma_w_sp(ialpha,ialphap,:),eta)
	end do
end do
do ialpha=1,3
	do ialphap=1,3
		call broad_vector(nkubo_ex,wn_ex,skubo_ex_int(ialpha,ialphap,:), &
		w0/27.211385d0,wrange/27.211385d0,nw,wp,sigma_w_ex(ialpha,ialphap,:),eta)
	end do
end do


!write frequency dependent conductivity
open(50,file=file_name_sp)
open(60,file=file_name_ex)
do iw=1,nw
	feps=1.0d0
	write(50,*) wp(iw)*27.211385d0,-realpart(feps*sigma_w_sp(1,1,iw)), &
				-realpart(feps*sigma_w_sp(1,2,iw)), &
				-realpart(feps*sigma_w_sp(1,3,iw)), &
				-realpart(feps*sigma_w_sp(2,1,iw)), &
				-realpart(feps*sigma_w_sp(2,2,iw)), &
				-realpart(feps*sigma_w_sp(2,3,iw)), &
				-realpart(feps*sigma_w_sp(3,1,iw)), &
				-realpart(feps*sigma_w_sp(3,2,iw)), &
				-realpart(feps*sigma_w_sp(3,3,iw))                
	write(60,*) wp(iw)*27.211385d0,realpart(feps*sigma_w_ex(1,1,iw)), &
				realpart(feps*sigma_w_ex(1,2,iw)), &
				realpart(feps*sigma_w_ex(1,3,iw)), &
				realpart(feps*sigma_w_ex(2,1,iw)), &
				realpart(feps*sigma_w_ex(2,2,iw)), &
				realpart(feps*sigma_w_ex(2,3,iw)), &
				realpart(feps*sigma_w_ex(3,1,iw)), &
				realpart(feps*sigma_w_ex(3,2,iw)), &	
				realpart(feps*sigma_w_ex(3,3,iw))						  
end do
close(50)
close(60)


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_kubo_parameters(w0,wrange,nw,eta,type_broad,material_name, &
file_name_sp,file_name_ex)
implicit real*8 (a-h,o-z)

character(100) type_broad
character(100) material_name
character(100) file_name_sp
character(100) file_name_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
open(10,file='kubo_w.in')
read(10,*)
read(10,*) w0
read(10,*)
read(10,*) wrange 
read(10,*)
read(10,*) nw
read(10,*) 
read(10,*) eta
read(10,*) 
read(10,*) type_broad
read(10,*) 
read(10,*) material_name
read(10,*) 
read(10,*) file_name_sp
read(10,*) file_name_ex
close(10)
	
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    	  
subroutine hoppings_observables_TB(norb,nR,nRvec,R,B,shop,hhop, &
rhop,sderhop,hderhop)
implicit real*8 (a-h,o-z)
dimension R(3,3),B(norb,3),nRvec(nR,3),phop(3,nR,norb,norb),rhop(3,nR,norb,norb)
dimension shop(norb,norb,nR),hhop(norb,norb,nR)
dimension sderhop(3,nR,norb,norb),hderhop(3,nR,norb,norb)

complex*16 hderhop,sderhop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
rhop=0.0d0
hderhop=0.0d0
sderhop=0.0d0
do iR=1,nR
	do ialpha=1,norb
		do ialphap=1,ialpha
			Rx=nRvec(iR,1)*R(1,1)+nRvec(iR,2)*R(2,1)
			Ry=nRvec(iR,1)*R(1,2)+nRvec(iR,2)*R(2,2)
			Rz=0.0d0

			hderhop(1,iR,ialpha,ialphap)=complex(0.0d0,Rx)*hhop(ialpha,ialphap,iR)
			hderhop(2,iR,ialpha,ialphap)=complex(0.0d0,Ry)*hhop(ialpha,ialphap,iR)
			hderhop(3,iR,ialpha,ialphap)=0.0d0
		
			sderhop(1,iR,ialpha,ialphap)=complex(0.0d0,Rx)*shop(ialpha,ialphap,iR)
			sderhop(2,iR,ialpha,ialphap)=complex(0.0d0,Ry)*shop(ialpha,ialphap,iR)
			sderhop(3,iR,ialpha,ialphap)=0.0d0

		end do
	end do
end do

!point-like interactions
do nn=1,norb
	rhop(1,1,nn,nn)=B(nn,1)
	rhop(2,1,nn,nn)=B(nn,2)
	rhop(3,1,nn,nn)=B(nn,3)
end do
	

return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
subroutine get_vme_kernels(rkx,rky,rkz,nR,nRvec,norb,R, &
hkernel,skernel,shop,hhop,rhop,sderhop,hderhop,sderkernel,hderkernel,akernel, &
pgaugekernel)
implicit real*8 (a-h,o-z)

dimension R(3,3)
dimension nRvec(nR,3)

dimension shop(norb,norb,nR)
dimension hhop(norb,norb,nR)
dimension rhop(3,nR,norb,norb)
dimension sderhop(3,nR,norb,norb)
dimension hderhop(3,nR,norb,norb)

dimension hkernel(norb,norb)
dimension skernel(norb,norb)	  
dimension sderkernel(3,norb,norb)
dimension hderkernel(3,norb,norb)
dimension akernel(3,norb,norb)
dimension pgaugekernel(3,norb,norb)

complex*16 hderhop,sderhop
complex*16 phase,factor,hkernel,skernel
complex*16 sderkernel,hderkernel,akernel
complex*16 pgaugekernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  	  
!evaluate \alpha \alpha' matrices for a given k

hkernel=0.0d0
hderkernel=0.0d0
skernel=0.0d0
sderkernel=0.0d0
akernel=0.0d0
pgaugekernel=0.0d0
do ialpha=1,norb
do ialphap=1,ialpha   

	do iRp=1,nR
		Rx=dble(nRvec(iRp,1))*R(1,1)+dble(nRvec(iRp,2))*R(2,1)
		Ry=dble(nRvec(iRp,1))*R(1,2)+dble(nRvec(iRp,2))*R(2,2)
		Rz=0.0d0
		phase=complex(0.0d0,rkx*Rx+rky*Ry+rkz*Rz)
		factor=exp(phase)     
					
		hkernel(ialpha,ialphap)=hkernel(ialpha,ialphap)+ &
		factor*hhop(ialpha,ialphap,iRp)                
		skernel(ialpha,ialphap)=skernel(ialpha,ialphap)+ &
		factor*shop(ialpha,ialphap,iRp)   
					
		do nj=1,3 
			sderkernel(nj,ialpha,ialphap)=sderkernel(nj,ialpha,ialphap)+ &
			factor*sderhop(nj,iRp,ialpha,ialphap)
		
			hderkernel(nj,ialpha,ialphap)=hderkernel(nj,ialpha,ialphap)+ &
			factor*hderhop(nj,iRp,ialpha,ialphap) 

			akernel(nj,ialpha,ialphap)=akernel(nj,ialpha,ialphap)+ &
			factor*(rhop(nj,iRp,ialpha,ialphap)+ &
			complex(0.0d0,1.0d0)*sderhop(nj,iRp,ialpha,ialphap))            
		end do  
	end do
	
	do nj=1,3
		pgaugekernel(nj,ialpha,ialphap)=hkernel(ialpha,ialphap)* &
		complex(0.0d0,1.0d0)* &
		(rhop(nj,1,ialphap,ialphap)-rhop(nj,1,ialpha,ialpha))
	end do
	
	!complex conjugate
	do nj=1,3
		hkernel(ialphap,ialpha)=conjg(hkernel(ialpha,ialphap))
		skernel(ialphap,ialpha)=conjg(skernel(ialpha,ialphap))
		sderkernel(nj,ialphap,ialpha)=conjg(sderkernel(nj,ialpha,ialphap))
		hderkernel(nj,ialphap,ialpha)=conjg(hderkernel(nj,ialpha,ialphap))
		akernel(nj,ialphap,ialpha)=conjg(akernel(nj,ialpha,ialphap))+ &
		complex(0.0d0,1.0d0)*conjg(sderkernel(nj,ialpha,ialphap))     
		pgaugekernel(nj,ialphap,ialpha)=conjg(pgaugekernel(nj,ialpha,ialphap))        
	end do
	
end do
end do  
	
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
subroutine get_eigen_vme(rkx,rky,rkz,nR,nRvec,norb,R, &
hkernel,skernel,akernel,sderkernel,hderkernel,pgaugekernel, &
hk_ev,e,pgauge,vjseudoa,vjseudob,vme) 
implicit real*8 (a-h,o-z)

dimension R(3,3),nRvec(nR,3),hkernel(norb,norb),skernel(norb,norb)
dimension sderkernel(3,norb,norb),hderkernel(3,norb,norb)
dimension pgaugekernel(3,norb,norb),pgauge(3,norb,norb)
dimension hk_alpha(norb,norb),hk_ev(norb,norb),e(norb)
dimension akernel(3,norb,norb)

dimension vjseudoa(3,norb,norb),vjseudob(3,norb,norb),vme(3,norb,norb)

complex*16 hkernel,akernel,skernel,sderkernel,hderkernel
complex*16 hk_alpha,hk_ev,vjseudoa,vjseudob,vme,pgaugekernel,pgauge

complex*16 amu,amup,aux1,factor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  

e=0.0d0
call diagoz(norb,e,hkernel) 
hk_ev(:,:)=hkernel(:,:)

!phase election
do j=1,norb
aux1=0.0d0
do i=1,norb
	aux1=aux1+hk_ev(i,j)
end do
!argument of the sym
arg=atan2(aimag(aux1),realpart(aux1))
factor=exp(complex(0.0d0,-arg))
!write(*,*) 'sum is now:',aux1*factor
hk_ev(:,j)=hk_ev(:,j)*factor
end do   

write(*,*) 'computing velocity matrix elements'      
vme=0.0d0
vjseudoa=0.0d0
vjseudob=0.0d0
pgauge=0.0d0
do nn=1,norb
do nnp=1,nn        
	!momentums and A and B term
	do ialpha=1,norb
	do ialphap=1,norb
		amu=hk_ev(ialpha,nn)
		amup=hk_ev(ialphap,nnp)
		do nj=1,3      
		pgauge(nj,nn,nnp)=pgauge(nj,nn,nnp)+ &
		conjg(amu)*amup*pgaugekernel(nj,ialpha,ialphap)
		
		vjseudoa(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+ &
		conjg(amu)*amup*hderkernel(nj,ialpha,ialphap)
		
		vjseudob(nj,nn,nnp)=vjseudob(nj,nn,nnp)+conjg(amu)*amup* &
		(e(nn)*akernel(nj,ialpha,ialphap)-e(nnp)*conjg(akernel(nj,ialphap,ialpha)))* &
		complex(0.0d0,1.0d0)         
		end do                      
	end do
	
	end do
	do nj=1,3
	vme(nj,nn,nnp)=vjseudoa(nj,nn,nnp)+vjseudob(nj,nn,nnp)
	
	pgauge(nj,nnp,nn)=conjg(pgauge(nj,nn,nnp))
	vme(nj,nnp,nn)=conjg(vme(nj,nn,nnp))
	vjseudoa(nj,nnp,nn)=conjg(vjseudoa(nj,nn,nnp))
	vjseudob(nj,nnp,nn)=conjg(vjseudob(nj,nn,nnp))
	end do	
end do
end do

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_kubo_intens(nv,npointstotal,vcell,norb,nkubo_sp,e,kcount,vme,wn_sp, &
skubo_sp_int)
implicit real*8 (a-h,o-z)

dimension vme(3,norb,norb)
dimension skubo_sp_int(3,3,nkubo_sp)
dimension wn_sp(nkubo_sp)
dimension e(norb)

complex*16 vme,skubo_sp_int,factor1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
do nn=1,norb	
!fermi disrtibution
if (nn.le.nv) then
	fnn=1.0d0
else
	fnn=0.0d0
end if		    
do nnp=1,norb
	wn_sp(kcount)=e(nn)-e(nnp)
	!fermi distribution
	if (nnp.le.nv) then
	fnnp=1.0d0
	else
	fnnp=0.0d0
	end if                    
	!DEDICE PREFACTOR WITH OCCUPATION
	if ((fnn-fnnp).eq.0.0d0) then !UPDATE: energies removed from this factor
	factor1=0.0d0         
	else
	factor1=(fnn-fnnp)/(wn_sp(kcount)) 
	end if           	  
	!save oscillator stregths
	do nj=1,3
	do njp=1,3		  
		skubo_sp_int(nj,njp,kcount)=2.0d0*factor1/(dble(npointstotal)*vcell) &
		*vme(nj,nn,nnp)*vme(njp,nnp,nn)		        		  
	end do	
	end do	
	kcount=kcount+1			
end do
end do	  	  

return
end	  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
subroutine broad_vector(n,wn,fn,w0,wrange,nw,wp,fw_c,eta)
implicit real*8 (a-h,o-z)

dimension wn(n),fn(n),fn_real(n),wp(nw),wp_aux(nw),fw(nw),fw_c(nw)
complex*16 fw_c,fn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
wp=0.0d0
fw_c=0.0d0
do i=1,n
fn_real(i)=realpart(fn(i))
end do

do i=1,nw
wp(i)=(w0+wrange/dble(nw)*dble(i-1))
do inn=1,n
	fw_c(i)=fw_c(i)+fn_real(inn)*(aimag(1.0d0/(wp(i)-wn(inn)-complex(0.0d0,1.0d0)*eta)))
	!write(*,*) sc**(-1.0d0)*wn(inn),fn(inn)
end do
!pause
!write(*,*) sc**(-1.0d0)*wp(i),fw_c(i)
end do


return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine crossproduct(ax,ay,az,bx,by,bz,cx,cy,cz)
implicit real*8 (a-h,o-z)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
cx=ay*bz-az*by
cy=az*bx-ax*bz
cz=ax*by-ay*bx     
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         diagoz
!   INPUTS:       h matrix to diagonalize
!                 n dimension of h
!   OUTPUTS:      w; eigenvalues of h
!                 h;  gives eigenvectors by columns as output
!   DESCRIPTION:  this subroutine uses Lapack libraries to diagonalize
!                 an hermitian complex matrix.
!
!     Juan Jose Esteve-Paredes                28.11.2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine diagoz(n,w,h)
!finding the eigenvalues of a complex matrix using LAPACK
implicit real*8 (a-h,o-z)
!declarations, notice double precision
dimension w(n)
complex*16 h(n,n),WORK(2*n)
dimension RWORK(3*n-2)
character*1 JOBZ,UPLO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!find the solution using the LAPACK routine ZGEEV
JOBZ='V'
UPLO='U'
LWORK=2*n
	
call zheev(JOBZ, UPLO, n, h, n, w, WORK, LWORK, RWORK, INFO)
end
