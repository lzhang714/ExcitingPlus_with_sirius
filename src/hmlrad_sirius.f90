
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine hmlrad_sirius
  
  use modmain
#ifdef _SIRIUS_
  use mod_sirius
#endif
!  
! 2019 L.Zhang@UFL created to comply to SIRIUS library, following the spirit of Exciting-Sirius interface. 
!
  implicit none
  
  ! local variables
  integer is,ia,ias,nr,ir,if1,if3,inonz,ireset1,ireset3
  integer l1,l2,l3,m1,m2,m3,lm1,lm2,lm3
  integer ilo,ilo1,ilo2,io,io1,io2,nalo1,maxnlo,i
  integer mid,info,lwork
  real(8) t1,t2,angular,rm,aa
  real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
  real(8) r2inv(nrmtmax),rmtable(nrmtmax)
  
  integer                    haaijSize          ! no hloloijSize here
  integer,    allocatable :: haloijSize(:)
  integer,    allocatable :: lfromlm(:),mfromlm(:)
  complex(8), allocatable :: haaij(:,:,:),haloij(:,:,:),hloloij(:,:,:)
  real(8),    allocatable :: haaintegrals(:,:,:,:,:),halointegrals(:,:,:,:),hlolointegrals(:,:,:)  
  
  logical valence_relativity
  
  logical :: usesplines  
  real(8) alpha
  parameter (alpha=1d0 / 137.03599911d0)
  integer polyord
  parameter (polyord=3)
  integer ipiv(polyord+1)
  real(8), allocatable :: poly(:,:),weight(:),ints(:),abscissa(:),work(:)
  
  ! here some variables for the sake of convenience.
  ! we have lm2l,lm2m in src/addons/init3.f90, same thing. 
  allocate (lfromlm(lmmaxvr))
  allocate (mfromlm(lmmaxvr))
  do l1 = 0,lmaxvr
    do m1 = -l1, l1
          lm1 = idxlm(l1, m1)
          lfromlm(lm1)=l1
          mfromlm(lm1)=m1
    enddo
  enddo
  
  ! if consider valence relativity: aa != 0
  valence_relativity = .false.
  if (valence_relativity) then
    aa = 0.5d0*alpha**2
  else
    aa = 0d0
  endif
  
  ! -------- poly spline stuff
  allocate(poly(0:polyord,0:polyord)) !  poly
  allocate(ints(0:polyord))           !  ints
  allocate(abscissa(0:polyord))       !  abscissa

  allocate(work(1))
  call dgetri(polyord+1,poly,polyord+1,ipiv,work,-1,info)
  lwork=int(work(1))
  deallocate(work)
  allocate(work(lwork))
  
  ! -------- determine size & allocate haaij/haloij/hloloij
  haaijSize=0
  Do is = 1, nspecies
    if1=0
    Do l1 = 0, lmaxapw         ! EXCITING used lmaxmat
      Do m1 = - l1, l1
            lm1 = idxlm (l1, m1)
            Do io1 = 1, apword(l1, is)
              if1=if1+1
            End Do
      End Do
    End Do
    if (if1.gt.haaijSize) haaijSize=if1
  Enddo
  ! APW-APW
  if (allocated(haaij)) deallocate(haaij)
  allocate(haaij(haaijSize,haaijSize,natmtot))
  haaij=dcmplx(0d0,0d0)
  ! EXCITING used lmaxmat in the 5th dimension of haaintegrals, 
  ! I changed it to lmaxapw b/c the 3rd dimension is lmaxapw.
  allocate(haaintegrals(lmmaxvr, apwordmax, 0:lmaxapw, apwordmax, 0:lmaxapw)) 
  haaintegrals(:,:,:,:,:)=1d100  
  
  if (allocated(haloij)) deallocate(haloij)
  if (allocated(haloijSize)) deallocate(haloijSize)
  allocate(haloijSize(nspecies))
  haloijSize=0
  maxnlo=0
  Do is = 1, nspecies
        ias=idxas (1, is)
        ilo=nlorb (is)
        if (ilo.gt.0) then
          l1 = lorbl (ilo, is)
          lm1 = idxlm (l1, l1)
          l3 = lorbl (1, is)
          lm3 = idxlm (l3, -l3)
          haloijSize(is)=idxlo (lm1, ilo, ias)- idxlo (lm3, 1, ias)+1
          if (maxnlo.lt.haloijSize(is)) maxnlo=haloijSize(is)
        endif
  Enddo
  if (maxnlo.gt.0) then 
        ! APW-LO
        allocate(haloij(maxnlo,haaijSize,natmtot))
        haloij=dcmplx(0d0,0d0)
        ! EXCITING used lmaxmat in the 3rd dimension of halointegrals, 
        ! I changed it to lmaxapw.
        Allocate (halointegrals(lmmaxvr, apwordmax, 0:lmaxapw, nlomax)) 
        ! LO-LO
        if (allocated(hloloij)) deallocate(hloloij)
        allocate(hlolointegrals(lmmaxvr,nlomax,nlomax))
        allocate(hloloij(maxnlo,maxnlo,natmtot))
        hloloij=dcmplx(0d0,0d0)
  else 
        write(*,'(" ERROR in hamiltonian_radial_integrals, maxnlo : ", i2.2)') maxnlo
  endif
  
  
  ! memo for EP
  !  lmaxvr: angular momentum cut-off for the muffin-tin density and potential
  ! lmmaxvr: (lmaxvr+1)^2
  ! lmaxmat: angular momentum cut-off for the outer-most loop in the hamiltonian 
  !          and overlap matrix setup (EX: lo in APW+lo method)
  
  
  ! ====================================== ! 
  !              loop species              !
  ! ====================================== !
  do is=1,nspecies
    
    ! nr, r2, r2inv
    nr=nrmt(is)
    do ir=1,nr
      r2(ir)=spr(ir,is)**2
      r2inv(ir)=1d0/r2(ir)
    enddo
    
    ! ---------------------- get weight
    allocate(weight(nr))
    weight=0d0

    mid=polyord/2
    do ir=2,nr-polyord-1
      abscissa(0:polyord)=spr(ir:ir+polyord,is)-spr(ir,is)  ! spr(spnrmax,nspecies) in genrmesh.f90
      poly(:,0)=1d0
      do io2=1,polyord
        do io1=0,polyord
              poly(io1,io2)=poly(io1,io2-1)*abscissa(io1)
        enddo
      enddo
      do io1=0,polyord
        ints(io1)=(poly(mid+1,io1)*abscissa(mid+1)-poly(mid,io1)*abscissa(mid))/dble(io1+1)
      enddo
      ! dgetrf computes an LU factorization of a general M-by-N matrix
      call dgetrf(polyord+1,polyord+1,poly,polyord+1,ipiv,info) 
      ! dgetri computes the inverse of a matrix using the LU factorization computed by DGETRF
      call dgetri(polyord+1,poly,polyord+1,ipiv,work,lwork,info) 
      call dgemm('N','N',1,polyord+1,polyord+1,1d0,ints(0),1,poly(0,0),polyord+1,1d0,weight(ir),1)
    enddo

    abscissa(0:polyord)=spr(1:1+polyord,is)-spr(1,is)
    poly(:,0)=1d0
    do io2=1,polyord
      do io1=0,polyord
        poly(io1,io2)=poly(io1,io2-1)*abscissa(io1)
      enddo
    enddo
    do io1=0,polyord
        ints(io1)=(poly(mid+1,io1)*abscissa(mid+1)-poly(0,io1)*abscissa(0))/dble(io1+1)
    enddo
    call dgetrf(polyord+1,polyord+1,poly,polyord+1,ipiv,info)
    call dgetri(polyord+1,poly,polyord+1,ipiv,work,lwork,info)
    call dgemm('N','N',1,polyord+1,polyord+1,1d0,ints(0),1,poly(0,0),polyord+1,1d0,weight(1),1)

    abscissa(0:polyord)=spr(nr-polyord:nr,is)-spr(nr-polyord,is)
    poly(:,0)=1d0
    do io2=1,polyord
      do io1=0,polyord
        poly(io1,io2)=poly(io1,io2-1)*abscissa(io1)
      enddo
    enddo
    do io1=0,polyord
      ints(io1)=(poly(polyord,io1)*abscissa(polyord)-poly(mid,io1)*abscissa(mid))/dble(io1+1)
    enddo
    call dgetrf(polyord+1,polyord+1,poly,polyord+1,ipiv,info)
    call dgetri(polyord+1,poly,polyord+1,ipiv,work,lwork,info)
    call dgemm('N','N',1,polyord+1,polyord+1,1d0,ints(0),1,poly(0,0),polyord+1,1d0,weight(nr-polyord),1)    
    
    
    ! ====================== ! 
    !       loop atoms       !
    ! ====================== !
    do ia=1,natoms(is)
      
      ias=idxas(ia,is)
      
      do ir = 1, nr
        rmtable(ir) = 1d0/(1d0-aa*veffmt (1, ir, ias)*y00)
      end do
        
      !---------------------------!
      !     APW-APW integrals     !
      !---------------------------!
      
      do l1 = 0, lmaxapw            ! EXCITING used lmaxmat here, I changed it to lmaxapw
        do io1 = 1, apword (l1,is)
          do l3 = 0, lmaxapw           ! EXCITING used lmaxmat here
            do io2 = 1, apword (l3, is)
                 
                if (l1.eq.l3) then
                  angular = dble(l1*(l1+1))
                  do ir = 1,nr
                     t1 = apwfr(ir, 1, io1, l1, ias) * apwfr(ir, 1, io2, l3, ias)
                     t2 = apwfr(ir, 2, io1, l1, ias) * apwfr(ir, 2, io2, l3, ias)
                     fr(ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir) + t1*veffmt(1,ir,ias)*y00) * r2(ir)
                  end do
                  call fderiv(-1, nr, spr(:,is), fr, gr, cf)
                  haaintegrals(1,io2,l3,io1,l1) = gr(nr) / y00
                else
                  haaintegrals(1,io2,l3,io1,l1) = 0.d0
                end if

                ! In Exciting, the double loop is under the condition usesplines=.false.
                do lm2 = 2, lmmaxvr
                  m2 = mfromlm(lm2)
                  l2 = lfromlm(lm2)
                  haaintegrals(lm2, io2, l3, io1, l1) = 0d0
                  do ir = 1, nr
                    haaintegrals(lm2, io2, l3, io1, l1) = haaintegrals(lm2, io2, l3, io1, l1) + &
                    & apwfr(ir, 1, io1, l1, ias) * apwfr(ir, 1, io2, l3, ias) * r2(ir) * veffmt(lm2, ir, ias) * weight(ir)
                  end do
                end do

            enddo
          enddo
        enddo
      enddo

      !---------------------------!
      !     lo-APW integtrals     !
      !------------------------- -!      
          
      do ilo = 1, nlorb(is)
        l1 = lorbl(ilo,is)
        do l3 = 0, lmaxapw           ! EXCITING used lmaxmat here
          do io = 1, apword(l3,is)
                     
                if (l1 .Eq. l3) then
                  angular=dble(l1*(l1+1))
                  do ir = 1, nr
                    !rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                    t1=apwfr(ir, 1, io, l1, ias)*lofr(ir, 1, ilo, ias)
                    t2=apwfr(ir, 2, io, l1, ias)*lofr(ir, 2, ilo, ias)
                    fr(ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir) + t1*veffmt(1,ir,ias)*y00) * r2(ir)
                  enddo
                  call fderiv(-1, nr, spr(:,is), fr, gr, cf)
                  halointegrals(1, io, l3, ilo) = gr (nr) / y00
                else
                  halointegrals(1, io, l3, ilo) = 0.d0
                endif

                ! In Exciting, the double loop is under the condition usesplines=.false.
                do lm2 = 2, lmmaxvr
                  m2 = mfromlm(lm2)
                  l2 = lfromlm(lm2)
                  halointegrals(lm2, io, l3, ilo)=0d0
                  do ir = 1, nr
                    halointegrals(lm2, io, l3, ilo) = halointegrals(lm2, io, l3, ilo) + &
                    & lofr(ir, 1, ilo, ias) * apwfr(ir, 1, io, l3, ias) * r2(ir) * veffmt(lm2, ir, ias) * weight(ir)
                  enddo
                enddo

          enddo
        enddo
      enddo

      !---------------------------!
      !      lo-lo integrals      !
      !---------------------------!
      
      do ilo1 = 1, nlorb(is)
        l1 = lorbl(ilo1,is)
        do ilo2 = 1, nlorb(is)
          l3 = lorbl(ilo2,is)
              
                if (l1 .Eq. l3) then
                  angular=dble(l1*(l1+1))
                  do ir = 1, nr
                    !rm=1d0/(1d0-a*veffmt (1, ir, ias)*y00)
                    t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                    t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                    fr(ir) = (0.5d0*t2*rmtable(ir) + 0.5d0*angular*t1*rmtable(ir)*r2inv(ir) + t1*veffmt(1,ir,ias)*y00) * r2(ir)
                  enddo
                  call fderiv(-1, nr, spr(:,is), fr, gr, cf)
                  hlolointegrals(1, ilo1, ilo2) = gr(nr) / y00
                else
                  hlolointegrals(1, ilo1, ilo2) = 0.d0
                endif
                    
                ! In Exciting, the double loop is under the condition usesplines=.false.
                do lm2 = 2, lmmaxvr
                  m2 = mfromlm(lm2)
                  l2 = lfromlm(lm2)
                  hlolointegrals(lm2, ilo1, ilo2)=0d0
                  do ir = 1, nr
                    hlolointegrals(lm2, ilo1, ilo2) = hlolointegrals(lm2, ilo1, ilo2) + &
                    &lofr(ir, 1, ilo1, ias) * lofr(ir, 1, ilo2, ias) * r2(ir) * veffmt(lm2, ir, ias) * weight(ir)
                  end do
                end do

        enddo
      enddo


      ! now pass haaintegrals/halointegrals/hlolointegrals to SIRIUS
      if (use_sirius_library.and..not.use_sirius_hmlrad) then
#ifdef _SIRIUS_
               ! pass apw-apw
               do l1 = 0, lmaxapw
                  do io1 = 1, apword(l1,is)
                     do l2 = 0, lmaxapw
                        do io2 = 1, apword(l2,is)
                          call sirius_set_h_radial_integrals( sctx, ias, lmmaxvr, &
                               & haaintegrals(1, io1, l1, io2, l2),  l1=l1,o1=io1,  l2=l2,o2=io2 )
                        end do
                     end do
                  end do
               end do
               ! pass apw-lo
               do ilo = 1, nlorb(is)
                  do l2 = 0, lmaxapw
                     do io2 = 1, apword(l2,is)
                          call sirius_set_h_radial_integrals( sctx, ias, lmmaxvr, &
                               & halointegrals(1, io2, l2, ilo),  ilo1=ilo,  l2=l2,o2=io2 )
                          call sirius_set_h_radial_integrals( sctx, ias, lmmaxvr, &
                               & halointegrals(1, io2, l2, ilo),  l1=l2,o1=io2,  ilo2=ilo )
                     end do
                  end do
               end do
               ! pass lo-lo
               do ilo1 = 1, nlorb(is)
                  do ilo2 = 1, nlorb(is)
                          call sirius_set_h_radial_integrals( sctx, ias, lmmaxvr, &
                               & hlolointegrals(1, ilo1, ilo2),  ilo1=ilo1,  ilo2=ilo2 )
                  end do
               end do
#else
               stop sirius_error
#endif
      endif
      
      ! ---------------------------------------------------------------- ! 
      ! EXCITING then does the angular integral calculation here, using  !
      ! haaintegrals/halointegrals/hlolointegrals and Gaunt coefficients !
      ! to get haaij/haloij/hloloij (which have major usage everywhere). !
      ! At the moment, we only produce what sirius needs.                !
      ! ---------------------------------------------------------------- !

    enddo ! ia
    
    ! ====================== ! 
    !       loop atoms       !
    ! ====================== !
    
    deallocate(weight)

  enddo ! is
  
  ! ====================================== ! 
  !              loop species              !
  ! ====================================== !  
  
  deallocate(lfromlm,mfromlm)  
  deallocate(poly,abscissa,ints,work)
  
  deallocate(haaintegrals)
  deallocate(halointegrals)
  deallocate(hlolointegrals)
  
  deallocate(haaij)
  deallocate(haloij)
  deallocate(hloloij)
      
  return
  
end subroutine





