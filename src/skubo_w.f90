!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine skubo_w(nR,norb,norb_ex,nv_ex,nc_ex,nv,Rvec,R,B,hhop,shop,npointstotal,rkx, &
rky,rkz,fk_ex,e_ex,eigval_stack,eigvec_stack)
implicit real*8 (a-h,o-z)

!out of subroutine arrays 
dimension Rvec(nR,3)
dimension R(3,3)
dimension hhop(norb,norb,nR)
dimension shop(norb,norb,nR)
dimension rkx(npointstotal)
dimension rky(npointstotal)
dimension rkz(npointstotal)
dimension fk_ex(norb_ex,norb_ex)
dimension e_ex(norb_ex)
dimension B(norb,3)
dimension rhop(3,nR,norb,norb)
dimension eigval_stack(nv_ex + nc_ex,npointstotal)
dimension eigvec_stack(norb,nv_ex + nc_ex,npointstotal)

!sp and exciton variables
allocatable hk_ev(:,:)
allocatable e(:)
allocatable vme(:,:,:,:)
allocatable sigma_w_sp(:,:,:)

!exciton arrays
allocatable vme_ex(:,:,:)  
allocatable wp(:)
allocatable skubo_ex_int(:,:,:)
allocatable sigma_w_ex(:,:,:)


complex*16 hhop
complex*16 fk_ex
complex*16 eigvec_stack
complex*16 hk_ev
complex*16 vme
complex*16 vme_ex 
complex*16 skubo_ex_int 
complex*16 sigma_w_sp
complex*16 sigma_w_ex

character(100) type_broad
character(100) file_name_sp
character(100) file_name_ex


pi=acos(-1.0d0)

nbands = nv_ex + nc_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Here we work in atomic units
Rvec=Rvec/0.52917721067121d0
B=B/0.52917721067121d0 
R=R/0.52917721067121d0 
hhop=hhop/27.211385d0
rkx=rkx*0.52917721067121d0
rky=rky*0.52917721067121d0
rkz=rkz*0.52917721067121d0
e_ex=e_ex/27.211385d0
eigval_stack=eigval_stack/27.211385d0

call crossproduct(R(1,1),R(1,2),R(1,3),R(2,1), &
R(2,2),R(2,3),cx,cy,cz)         
vcell=sqrt(cx**2+cy**2+cz**2)

!call fill_nRvec(nR,R,Rvec,nRvec)
!get printing parameters
call get_kubo_parameters(w0,wrange,nw,eta,type_broad, &
file_name_sp,file_name_ex)
eta=eta/27.211385d0
 
!SP arrays
allocate (wp(nw))
allocate (sigma_w_sp(3,3,nw))
allocate (hk_ev(norb,nbands))
allocate (e(nbands))
allocate (vme(npointstotal,3,nbands,nbands))
wp=0.0d0
wn_sp=0.0d0
sigma_w_sp=0.0d0
do i=1,nw
  wp(i)=(w0+wrange/dble(nw)*dble(i-1))/27.211385d0
end do  

!exciton arrays
norb_ex_band=nv_ex*nc_ex !number of electron-hole pairs per k-point 
norb_ex_cut=norb_ex  !total number of optical transitions 
allocate (vme_ex(3,norb_ex_cut,2))
allocate (sigma_w_ex(3,3,nw))
allocate (skubo_ex_int(3,3,norb_ex_cut))
vme_ex=0.0d0
sigma_w_ex=0.0d0
skubo_ex_int=0.0d0

call exciton_oscillator_strength(nR,norb,norb_ex,nv_ex,nc_ex,nv,Rvec,R,B,hhop,shop,npointstotal,rkx, &
                                 rky,rkz,fk_ex,e_ex,eigval_stack,eigvec_stack,vme,vme_ex,.false.)


! Obtain SP Kubo
do ibz=1,npointstotal
  e(:) = eigval_stack(:, ibz)
          
  !get strength for kubo SP
  call get_kubo_intens(nv_ex,npointstotal,vcell,nbands,e,vme(ibz, :, :, :),nw,wp,sigma_w_sp,eta)
  
end do


!fill kubo oscillators of EXCITONS
do nn=1,norb_ex_cut
  !save oscillator stregths
  do nj=1,3
    do njp=1,3
      skubo_ex_int(nj,njp,nn)=pi/(dble(npointstotal)*vcell) &
      *conjg(vme_ex(nj,nn,1))*vme_ex(njp,nn,1)/e_ex(nn)   !pick the correct order of operators
    end do
  end do 
end do


!excitons
do ialpha=1,2
  do ialphap=1,2
    call broad_vector(type_broad,norb_ex_cut,e_ex,skubo_ex_int(ialpha,ialphap,:), &
    nw,wp,sigma_w_ex(ialpha,ialphap,:),eta)
  end do
end do 
 
!write frequency dependent conductivity	  
open(50,file=file_name_sp)
open(60,file=file_name_ex)
do iw=1,nw
  feps=1.0d0
  !feps=4.0d0*pi*1.0d0/137.035999084d0*100.0d0   !absorbance units
  !eps=4.0d0   !\sigma_0 units
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

!write exciton oscillator strengths
norb_ex_cut=nv_ex*nc_ex*npointstotal
write(60,*) '' ! Empty line
do iex=1,norb_ex_cut
  write(60,*) e_ex(iex)*27.211385d0,realpart(vme_ex(1,iex,1)), imagpart(vme_ex(1,iex,1)), &
              realpart(vme_ex(2,iex,1)), imagpart(vme_ex(2,iex,1)), &
              realpart(vme_ex(3,iex,1)), imagpart(vme_ex(2,iex,1))
end do

close(50)
close(60)

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exciton_oscillator_strength(nR,norb,norb_ex,nv_ex,nc_ex,nv,Rvec,R,B,hhop,shop,npointstotal,rkx, &
  rky,rkz,fk_ex,e_ex,eigval_stack,eigvec_stack,vme,vme_ex,convert_to_au)
  implicit real*8 (a-h,o-z)
  
  !out of subroutine arrays 
  dimension Rvec(nR,3)
  dimension R(3,3)
  dimension hhop(norb,norb,nR)
  dimension shop(norb,norb,nR)
  dimension rkx(npointstotal)
  dimension rky(npointstotal)
  dimension rkz(npointstotal)
  dimension fk_ex(norb_ex,norb_ex)
  dimension e_ex(norb_ex)
  dimension B(norb,3)
  dimension rhop(3,nR,norb,norb)
  dimension eigval_stack(nv_ex + nc_ex,npointstotal)
  dimension eigvec_stack(norb,nv_ex + nc_ex,npointstotal)
  dimension vme_ex(3,norb_ex,2)
  dimension vme(npointstotal,3,nv_ex + nc_ex,nv_ex + nc_ex)
  
  !sp and exciton variables
  allocatable sderhop(:,:,:,:)
  allocatable hderhop(:,:,:,:)
  allocatable hkernel(:,:)
  allocatable skernel(:,:)
  allocatable hderkernel(:,:,:)
  allocatable sderkernel(:,:,:)
  allocatable akernel(:,:,:)
  allocatable pgaugekernel(:,:,:)  
  allocatable hk_ev(:,:)
  allocatable e(:)
  allocatable pgauge(:,:,:)
  allocatable vjseudoa(:,:,:) 
  allocatable vjseudob(:,:,:)  
  
  complex*16 hhop
  complex*16 sderhop
  complex*16 hderhop
  complex*16 fk_ex
  complex*16 eigvec_stack
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

  logical convert_to_au


  if (convert_to_au) then
    Rvec=Rvec/0.52917721067121d0
    B=B/0.52917721067121d0 
    R=R/0.52917721067121d0 
    hhop=hhop/27.211385d0
    rkx=rkx*0.52917721067121d0
    rky=rky*0.52917721067121d0
    rkz=rkz*0.52917721067121d0
    e_ex=e_ex/27.211385d0
    eigval_stack=eigval_stack/27.211385d0
  end if

  nbands = nv_ex + nc_ex
  
  pi=acos(-1.0d0)
      
  !get unit cell volume
  call crossproduct(R(1,1),R(1,2),R(1,3),R(2,1), &
  R(2,2),R(2,3),cx,cy,cz)         
  vcell=sqrt(cx**2+cy**2+cz**2)
   
  !SP arrays
  allocate (sderhop(3,nR,norb,norb))
  allocate (hderhop(3,nR,norb,norb))
  allocate (hkernel(norb,norb))
  allocate (skernel(norb,norb))
  allocate (hderkernel(3,norb,norb))
  allocate (sderkernel(3,norb,norb))
  allocate (akernel(3,norb,norb))
  allocate (pgaugekernel(3,norb,norb))  
  allocate (hk_ev(norb,nbands))
  allocate (e(nbands))
  allocate (pgauge(3,nbands,nbands))
  allocate (vjseudoa(3,nbands,nbands)) 
  allocate (vjseudob(3,nbands,nbands))
  
  !exciton arrays
  norb_ex_band=nv_ex*nc_ex !number of electron-hole pairs per k-point 
  norb_ex_cut=norb_ex  !total number of optical transitions 

  vme_ex=0.0d0
  
  !getting some SP variables
  call hoppings_observables_TB(norb,nR,Rvec,shop,hhop,rhop,sderhop,hderhop)
  !11/05/2023 JJEP: fill rhop here. Easier to extend to DFT later 
  rhop=0.0d0
  do nn=1,norb
    do nj=1,3
      rhop(nj,1,nn,nn)=B(nn,nj)
    end do
  end do
  
  !Brillouin zone sampling	  
  !!$OMP PARALLEL DO PRIVATE(rkxp,rkyp,rkzp), &
  !!$OMP PRIVATE(hkernel,skernel,sderkernel,hderkernel,akernel), &  
  !!$OMP PRIVATE(hk_ev,e,pgaugekernel,pgauge,vjseudoa,vjseudob,vme), &
  !!$OMP PRIVATE(nj,iex,ib,nv_ip,nc_ip,j)
  do ibz=1,npointstotal
    !write(*,*) 'point:',ibz,npointstotal         
    rkxp=rkx(ibz)
    rkyp=rky(ibz)
    rkzp=rkz(ibz)
  
    hk_ev(:, :) = eigvec_stack(:, :, ibz)
    e(:) = eigval_stack(:, ibz)
            
    !get matrices in the \alpha, \alpha' basis (orbitals,k)    		
    call get_vme_kernels(rkxp,rkyp,rkzp,nR,Rvec,norb,hkernel,skernel,shop, &
    hhop,rhop,sderhop,hderhop,sderkernel,hderkernel,akernel,pgaugekernel)
    !velocity matrix elements
    call get_eigen_vme(norb,nbands,skernel,hkernel,akernel,hderkernel, &
    pgaugekernel,hk_ev,e,pgauge,vjseudoa,vjseudob,vme(ibz, :, :, :)) 

    !fill V_N
    !!$OMP critical	
    do nj=1,3
      do iex=1,norb_ex_cut
        !!! Iterate over band index explicitly (AJU 03-06-23)
        do ic=1,nc_ex
        do iv=1,nv_ex
        
          j = nc_ex * nv_ex * (ibz - 1) + nv_ex * (ic - 1) + iv
          !get valence/conduction indices      
          !nv_ip=(nv-nv_ex)+ib-int((ib-1)/nv_ex)*nv_ex
          !nc_ip=(nv+1)+int((ib-1)/nv_ex)
          !j=(ibz-1)*norb_ex_band+ib
      
          vme_ex(nj,iex,1)=vme_ex(nj,iex,1)+fk_ex(j,iex)*vme(ibz, nj,iv,ic + nv_ex)
          vme_ex(nj,iex,2)=vme_ex(nj,iex,2)+fk_ex(j,iex)*vme(ibz, nj,ic + nv_ex,iv)
  
        end do
        end do
      end do
    end do       
    
  end do
  
  return
end
  
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine fill_nRvec(nR,R,Rvec,nRvec)
implicit real*8 (a-h,o-z)
dimension R(3,3),Rvec(nR,3),nRvec(nR,3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Deprecated; Bravais lattice vectors are passed directly in cartesian coordinates (AJU 03-06-23)
nRvec=0
do inR=1,nR
  if (abs(Rvec(inR,3)).gt.0.1d0) then
    nz=10
  else
    nz=0
  end if
end do

kcount=0
do n3=-nz,nz !run n3 first (better for 2D)
  do n1=-30,30
    do n2=-30,30
      do inR=1,nR
        Rx=n1*R(1,1)+n2*R(2,1)+n3*R(3,1)
        Ry=n1*R(1,2)+n2*R(2,2)+n3*R(3,2)
        Rz=n1*R(1,3)+n2*R(2,3)+n3*R(3,3)
        d=sqrt((Rx-Rvec(inR,1))**2+(Ry-Rvec(inR,2))**2 &
        +(Rz-Rvec(inR,3))**2)
        !write(*,*) inR,n1,n2,n3,d
        if (d.lt.0.01d0) then
          !write(*,*) inR,n1,n2,n3,kcount
          nRvec(inR,1)=n1
          nRvec(inR,2)=n2
          nRvec(inR,3)=n3
          kcount=kcount+1
          if (kcount.eq.nR) goto 88
        end if
      end do
    end do
  end do
end do
88 continue
  
return
end
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_kubo_parameters(w0,wrange,nw,eta,type_broad, &
file_name_sp,file_name_ex)
implicit real*8 (a-h,o-z)

character(100) type_broad
character(100) file_name_sp
character(100) file_name_ex

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
read(10,*) file_name_sp
read(10,*) file_name_ex
close(10)
    
return
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    	  
subroutine hoppings_observables_TB(norb,nR,Rvec,shop,hhop, &
rhop,sderhop,hderhop)

implicit real*8 (a-h,o-z)
dimension Rvec(nR,3),rhop(3,nR,norb,norb)
dimension shop(norb,norb,nR),hhop(norb,norb,nR)
dimension sderhop(3,nR,norb,norb),hderhop(3,nR,norb,norb)

complex*16 hhop
complex*16 hderhop,sderhop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         
hderhop=0.0d0
sderhop=0.0d0
do iR=1,nR
  do ialpha=1,norb
  do ialphap=1,ialpha
    Rx=Rvec(iR,1)
    Ry=Rvec(iR,2)
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


return
end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   	        
subroutine get_vme_kernels(rkx,rky,rkz,nR,Rvec,norb, &
hkernel,skernel,shop,hhop,rhop,sderhop,hderhop,sderkernel,hderkernel,akernel, &
pgaugekernel)
implicit real*8 (a-h,o-z)

dimension Rvec(nR,3)

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

complex*16 hhop
complex*16 hderhop,sderhop
complex*16 phase,factor,hkernel,skernel
complex*16 sderkernel,hderkernel,akernel
complex*16 pgaugekernel
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  	  	  
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
    Rx=Rvec(iRp,1)
    Ry=Rvec(iRp,2)
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
subroutine get_eigen_vme(norb,nbands,skernel, &
hkernel,akernel,hderkernel,pgaugekernel, &
hk_ev,e,pgauge,vjseudoa,vjseudob,vme) 

implicit real*8 (a-h,o-z)

dimension skernel(norb,norb),hkernel(norb,norb)
dimension hderkernel(3,norb,norb)
dimension pgaugekernel(3,norb,norb),pgauge(3,nbands,nbands)
dimension hk_ev(norb,nbands),e(nbands)
dimension akernel(3,norb,norb)
dimension ecomplex(norb)
dimension vjseudoa(3,nbands,nbands),vjseudob(3,nbands,nbands),vme(3,nbands,nbands)

complex*16 ecomplex
complex*16 skernel,hkernel,akernel,hderkernel
complex*16 hk_ev,vjseudoa,vjseudob,vme,pgaugekernel,pgauge

complex*16 amu,amup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
  
!!! Remove diagonalization; instead use stored eigvecs and eigvals (AJU 29-05-23)
! e=0.0d0
! call diagoz(norb,e,hkernel) 
! !call diagoz_arbg(norb,ecomplex,skernel,hkernel)             
! !call order(norb,ecomplex,e,hkernel) 
! do ii=1,norb
! do jj=1,norb
!   hk_ev(ii,jj)=hkernel(ii,jj)
! end do
! end do

! !phase election
! do j=1,norb
!   aux1=0.0d0
!   do i=1,norb
!     aux1=aux1+hk_ev(i,j)
!   end do
  
!   !argument of the sym
!   arg=atan2(aimag(aux1),realpart(aux1))
!   factor=exp(complex(0.0d0,-arg))
!   !write(*,*) 'sum is now:',aux1*factor
!   do ii=1,norb
!     hk_ev(ii,j)=hk_ev(ii,j)*factor
! 	!hk_ev(ii,j)=hk_ev(ii,j)*1.0d0
!   end do      
! end do   

!write(*,*) 'computing velocity matrix elements'      
vme=0.0d0
vjseudoa=0.0d0
vjseudob=0.0d0
pgauge=0.0d0

!$OMP PARALLEL DO
do nn=1,nbands
  
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
    !if (nn.eq.2 .and. nnp.eq.1) then
    !write(*,*) nn,nnp,ialpha,ialphap,akernel(2,ialpha,ialphap),conjg(akernel(2,ialphap,ialpha))
    !pause
    !end if 
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
!$OMP END PARALLEL DO

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_kubo_intens(nv,npointstotal,vcell,nbands,e,vme,nw,wp,sigma_w_sp,eta)
implicit real*8 (a-h,o-z)

dimension vme(3,nbands,nbands)
dimension wp(nw),sigma_w_sp(3,3,nw)
dimension e(nbands)

complex*16 vme,skubo,sigma_w_sp
pi=acos(-1.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
do iw=1,nw  
  do nn=1,nbands
  !fermi disrtibution
    if (nn.le.nv) then
      fnn=1.0d0
    else
      fnn=0.0d0
    end if
    do nnp=1,nbands
      !fermi distribution
      if (nnp.le.nv) then
        fnnp=1.0d0
      else
        fnnp=0.0d0
      end if                    
      !DEDICE PREFACTOR WITH OCCUPATION
      if (abs(fnn-fnnp).lt.0.1d0) then !UPDATE: energies removed from this factor
        factor1=0.0d0         
      else
        factor1=(fnn-fnnp)/(e(nn)-e(nnp))
      end if
	  !lorentzian
      !delta_nnp=1.0d0/pi*aimag(1.0d0/(wp(iw)-e(nn)+e(nnp)-complex(0.0d0,eta)))
	  !exponential
    delta_nnp=1.0d0/eta*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta**2)*(wp(iw)-e(nn)+e(nnp))**2)
	  
          !save oscillator stregths
      do nj=1,3
        do njp=1,3
          skubo=vme(nj,nn,nnp)*vme(njp,nnp,nn)
          sigma_w_sp(nj,njp,iw)=sigma_w_sp(nj,njp,iw)+ &
          pi/(dble(npointstotal)*vcell)*factor1*skubo*delta_nnp
        end do
      end do
    		  
    end do
  end do
end do  
!pause

return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
subroutine broad_vector(type_broad,n,wn,fn,nw,wp,fw_c,eta)

implicit real*8 (a-h,o-z)

dimension wn(n),fn(n),fn_real(n),wp(nw),fw_c(nw)
complex*16 fw_c,fn
character(100) type_broad
pi=acos(-1.0d0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!wp=0.0d0
fw_c=0.0d0
do i=1,n
fn_real(i)=realpart(fn(i))
end do

do i=1,nw
!wp(i)=(w0+wrange/dble(nw)*dble(i-1))
do inn=1,n
  rbroad=1.0d0
  if (type_broad.eq.'lorentzian') then
    rbroad=1.0d0/pi*aimag(1.0d0/(wp(i)-wn(inn)-complex(0.0d0,1.0d0)*eta))
  end if
  if (type_broad.eq.'gaussian') then
    rbroad=1.0d0/eta*1.0d0/sqrt(2.0d0*pi)*exp(-0.5d0/(eta**2)*(wp(i)-wn(inn))**2)
  end if
  if (type_broad.eq.'exponential') then
    rbroad=0.5d0/eta*exp(-abs(wp(i)-wn(inn))/eta)
  end if
  fw_c(i)=fw_c(i)+fn_real(inn)*rbroad
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
return
end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   NAME:         diagoz_arbg
!   INPUTS:       h matrix to diagonalize
!                 s overlap matrix
!                 n dimension of h
!   OUTPUTS:      w; eigenvalues of h
!                 h;  gives eigenvectors by columns as output
!   DESCRIPTION:  this subroutine uses Lapack libraries to diagonalize
!                 an arbitrary complex matrix generalized problem
!
!     Juan Jose Esteve-Paredes                24.04.2017
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine diagoz_arbg(n,w,s,h)
!finding the eigenvalues of a complex matrix using LAPACK
      implicit real*8 (a-h,o-z)
!declarations, notice double precision
      complex*16 h(n,n),s(n,n),w(n),VL(n,n),VR(n,n),WORK(200*n)
      complex*16 ALPHA(n),BETA(n)      
      dimension RWORK(8*n),wreal(n)
      dimension ereal(n)
      integer INFO
      
      complex*16 haux(n,n),saux(n,n),ssuma
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      saux=s
      haux=h
      !find the solution using the LAPACK routine ZGGEEV     
      !e=0.0d0
      !do i=1,3
        !write(*,*) (s(i,j),j=1,3)
      !end do
      !pause
      !call zggev('V','V',n,s,n,h,n,ALPHA,BETA,VL,n,VR,n,WORK,200*n, RWORK,INFO)
      call zhegv( 1, 'V', 'U', n, h, n, s, n, wreal, WORK, 200*n, RWORK, INFO )
      !do i=1,n
        !write(*,*) wreal(i),1.0d0/wreal(i)
      !end do
      !do i=1,n
        !w(i)=ALPHA(i)/BETA(i)
        !write(*,*) ALPHA(i),BETA(i),w(i)
        !do j=1,n
          !h(i,j)=VR(i,j)
        !end do 
      !end do
      !pause


      !ANOTHER SIMILAR TO THE PREVIOUS ONE
      !call zgegv( 'N', 'V', n, h, n, s, n, ALPHA, BETA, &
      !QVL, n, VR, n, WORK, 200*n, RWORK, INFO )
      !write(*,*) 'Diagonalization=',INFO
      !write(*,*) INFO
      do i=1,n
        !w(i)=ALPHA(i)/BETA(i)
        w(i)=complex(wreal(i),0.0d0)
        !write(*,*) ALPHA(i),BETA(i),w(i)
        !do j=1,n
          !h(i,j)=VR(i,j)
        !end do 
        !write(*,*) wreal(i),h(i,1)
        !write(*,*) ALPHA(i),BETA(i),w(i)
      end do
      !pause
      do k=1,n
        ssuma=0.0d0
        do i=1,n
          do j=1,n
            ssuma=ssuma+conjg(h(i,k))*h(j,k)*saux(i,j)
          end do
        end do
        snorma=sqrt(realpart(ssuma)) 
        h(:,k)=1.0d0/snorma*h(:,k)       
        !if (k.eq.1) write(*,*) ssuma
      end do
      
      !write(*,*) 'Diagonalization=',INFO
      if (INFO.ne.0) then
        write(*,*) 'Bad diagonalization'
        write(*,*) 'Eigenvalues/ground state coefficients'
        do i=1,n
          write(*,*) wreal(i),h(i,1)
        end do        
        open(91,file='h_matlab_real.dat')
        open(92,file='s_matlab_real.dat')
        open(93,file='h_matlab_imag.dat')
        open(94,file='s_matlab_imag.dat')
        do i=1,n
          write(91,*) (realpart(haux(i,j)), j=1,n)
          write(92,*) (realpart(saux(i,j)), j=1,n)
          write(93,*) (aimag(haux(i,j)), j=1,n)
          write(94,*) (aimag(saux(i,j)), j=1,n)
        end do 
        close(91)
        close(92)
        close(93)
        close(94)  
        pause
      end if
      s=saux
      !do i=1,n
        !write(*,*) i,h(i,1),w(i)
      !end do
      !pause
    
      return
      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     !subroutine to order the complez eigenvalues and making them real
     !it also orders the eigenvectors
     subroutine order(n,ecomplex,e,hk)
     implicit real*8 (a-h,o-z)
     complex*16 ecomplex,hk,hkaux
     dimension ecomplex(n),e(n),eaux(n),hk(n,n),hkaux(n,n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     hkaux=hk
     do i=1,n
       eaux(i)=realpart(ecomplex(i))
     end do
     
     do i=1,n
       imin=minloc(eaux,1)      
       e(i)=eaux(imin)
       hk(:,i)=hkaux(:,imin)
       eaux(imin)=1.0d8
     end do
    
     return
     end