
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: olprad
! !INTERFACE:

subroutine olprad

  use modmain
!#ifdef _SIRIUS_
!  use mod_sirius
!#endif
  
! !DESCRIPTION:
!   Calculates the radial overlap integrals of the APW and local-orbital basis
!   functions. In other words, for spin $\sigma$ and atom $j$ of species $i$, it
!   computes integrals of the form
!   $$ o^{\sigma;ij}_{qp}=\int_0^{R_i}u^{\sigma;ij}_{q;l_p}(r)v^{\sigma;ij}_p(r) r^2dr $$
!   and
!   $$ o^{\sigma;ij}_{pp'}=\int_0^{R_i}v^{\sigma;ij}_p(r)v^{\sigma;ij}_{p'}(r) r^2dr,\quad l_p=l_{p'} $$
!   where $u^{\sigma;ij}_{q;l}$ is the $q$th APW radial function for angular
!   momentum $l$; and $v^{\sigma;ij}_p$ is the $p$th local-orbital radial
!   function and has angular momentum $l_p$.
!
! !REVISION HISTORY:
!   Created November 2003 (JKD)
!EOP
!BOC

  implicit none

  ! local variables
  integer is,ia,ias,ir,nr
  integer l,ilo,ilo1,ilo2,io,io1,io2
  real(8) aa

  ! automatic arrays
  real(8) r2(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
  real(8) :: angular,t1,t2,rm,a,alpha
  parameter (alpha=1d0/137.03599911d0)
  
  ! for the linearized Koelling-Harmon (valence relativity) correction.
  ! only used when passing to sirius
  real(8), allocatable :: h1aa (:,:,:,:)
  real(8), allocatable :: h1loa (:,:,:)
  real(8), allocatable :: h1lolo (:,:,:)


  logical valence_relativity

  ! ========================================================
  ! one important boolean set ONLY in this subroutine.
  valence_relativity = .false.
  ! ========================================================


  if (allocated(h1aa)) deallocate(h1aa)
  allocate (h1aa(apwordmax,apwordmax,lmaxapw,natmtot))   ! LZ changed "lmaxapw" to "lmaxmat"
  if (allocated(h1loa)) deallocate(h1loa)
  allocate (h1loa(apwordmax,nlomax,natmtot))
  if (allocated(h1lolo)) deallocate(h1lolo)
  allocate (h1lolo(nlomax, nlomax, natmtot))
  
  !h1aa=0.0d0
  !h1loa=0.0d0
  !h1lolo=0.0d0

  ! -------- if consider valence relativity then aa!=0
  if (valence_relativity) then
    aa = 0.5d0*alpha**2
  else
    aa = 0d0
  endif

  ! ---------------------------- double loop is+ia
  do is=1,nspecies
    nr=nrmt(is)
    do ir=1,nr
      r2(ir)=spr(ir,is)**2
    end do
    do ia=1,natoms(is)
      ias=idxas(ia,is)

      !---------------------------!
      !     APW-APW integtrals    !   ONLY calculate val-rel correction term h1aa
      !---------------------------!

      if (valence_relativity) then
          do l=0, lmaxapw                   ! LZ changed "lmaxmat" to "lmaxapw"
            angular=dble(l*(l+1))
            do io1=1, apword(l,is)
              do io2=1, apword(l,is)
                  do ir=1, nr
                      rm=1d0/(1d0-a*veffmt(1, ir, ias)*y00)
                      t1=apwfr(ir, 1, io1, l, ias)*apwfr(ir, 1, io2, l, ias)
                      t2=apwfr(ir, 2, io1, l, ias)*apwfr(ir, 2, io2, l, ias)
                      fr(ir) = aa*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2)*r2(ir)
                  enddo
                  call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                  !h1aa(io1,io2,l,ias) = gr(nr)
              enddo ! io2
            enddo ! io1
          enddo ! l
      endif

      !---------------------------!
      !     APW-lo integtrals     !
      !---------------------------!
      
      do ilo=1,nlorb(is)
        l=lorbl(ilo,is)
        do io=1,apword(l,is)
            do ir=1,nr
              fr(ir)=apwfr(ir,1,io,l,ias)*lofr(ir,1,ilo,ias)*r2(ir)
            end do
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            oalo(io,ilo,ias)=gr(nr)
                    if (valence_relativity) then
                      angular=dble(l*(l+1))
                      do ir=1, nr
                        rm=1d0/(1d0-a*veffmt(1, ir, ias)*y00)
                        t1=apwfr(ir, 1, io, l, ias)*lofr(ir, 1, ilo, ias)
                        t2=apwfr(ir, 2, io, l, ias)*lofr(ir, 2, ilo, ias)
                        fr(ir) = (aa*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2))*r2(ir)
                      enddo
                      call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                      !h1loa(io,ilo,ias) = gr(nr)
                    endif
        enddo
      enddo

      !---------------------------!
      !      lo-lo integrals      !
      !---------------------------!
      
      do ilo1=1,nlorb(is)
        l=lorbl(ilo1,is)
        do ilo2=1,nlorb(is)
            if (lorbl(ilo2,is).eq.l) then
              do ir=1,nr
                fr(ir)=lofr(ir,1,ilo1,ias)*lofr(ir,1,ilo2,ias)*r2(ir)
              end do
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              ololo(ilo1,ilo2,ias)=gr(nr)
                    if (valence_relativity) then
                      angular=dble(l*(l+1))
                      do ir = 1, nr
                        rm=1d0/(1d0-a*veffmt(1, ir, ias)*y00)
                        t1=lofr(ir, 1, ilo1, ias)*lofr(ir, 1, ilo2, ias)
                        t2=lofr(ir, 2, ilo1, ias)*lofr(ir, 2, ilo2, ias)
                        fr(ir) = (aa*(0.5d0*t2*rm**2 + 0.5d0*angular*t1*rm**2/spr(ir,is)**2))*r2(ir)
                      enddo
                      call fderiv(-1,nr,spr(:,is),fr,gr,cf)
                      !h1lolo (ilo1,ilo2,ias) = gr(nr)
                    endif            
            endif
        enddo
      enddo
    
    enddo
  enddo
  ! ---------------------------- double loop is+ia






  ! ---------------------------- pass overlap integrals to sirius
  !   
  ! We pass only the bare "oalo" and "ololo" to Sirius. 
  ! 
  ! NOTE: 
  ! In Exciting-Sirius interface, the valence relativity corrections (h1aa,h1loa,h1lolo) are 
  ! also passed to Sirius. 
  ! Here I assumed when we change "valence_rel" from "none" to "iora" in call sirius_set_parameters, 
  ! sirius will calculate and apply the correction. So no corrections are passed to sirius. 
  
  ! About Koelling and Harmon, from https://euler.phys.cmu.edu/cluster/WIEN2k/2Basic_concepts.html
  ! The linearized augmented plane wave (LAPW) method is among the most accurate methods for performing electronic structure calculations for crystals. It is based on the density functional theory for the treatment of exchange and correlation and uses e.g. the local spin density approximation (LSDA). Several forms of LSDA potentials exist in the literature , but recent improvements using the generalized gradient approximation (GGA) are available too (see sec. 2.1). For valence states relativistic effects can be included either in a scalar relativistic treatment (Koelling and Harmon 77) or with the second variational method including spin-orbit coupling (Macdonald 80, Novák 97). Core states are treated fully relativistically (Desclaux 69). 
  
  
  if (use_sirius_library.and.pass_olprad_to_sirius) then
  
#ifdef _SIRIUS_
         do is = 1, nspecies
           do ia = 1, natoms(is)
             ias = idxas(ia, is)
               !---------------------------!
               !     APW-lo integtrals     !
               !---------------------------!
               do ilo = 1, nlorb(is)
                  l = lorbl(ilo, is)
                  do io = 1, apword(l, is)
                     call sirius_set_o_radial_integral(sctx, ias, oalo(io, ilo, ias), l, o1=io, ilo2=ilo)
                     call sirius_set_o_radial_integral(sctx, ias, oalo(io, ilo, ias), l, ilo1=ilo, o2=io)
                  end do
               end do
               !---------------------------!
               !      lo-lo integrals      !
               !---------------------------!
               do ilo1 = 1, nlorb(is)
                  l = lorbl(ilo1, is)
                  do ilo2 = 1, nlorb(is)
                     if (lorbl(ilo2, is) .eq. l) then
                        call sirius_set_o_radial_integral(sctx, ias, ololo(ilo1, ilo2, ias), l, ilo1=ilo1, ilo2=ilo2)
                     end if
                  end do
               end do
               !----------------------------------------!
               !     valence relativity corrections     !
               !----------------------------------------!
               if (valence_relativity) then
                  
                     do l = 0, lmaxapw                       ! LZ changed "lmaxmat" to "lmaxapw"
                       do io1 = 1, apword (l, is)
                         do io2 = 1, apword (l, is)
                           !call sirius_set_o1_radial_integral(sctx, ias, h1aa(io1, io2, l, ias), l1=l, o1=io1, l2=l, o2=io2)
                         enddo
                       enddo
                     enddo
                
                     do ilo = 1, nlorb (is)
                       l = lorbl (ilo, is)
                       do io = 1, apword (l, is)
                         !call sirius_set_o1_radial_integral(sctx, ias, h1loa(io, ilo, ias), l1=l, o1=io, ilo2=ilo)
                         !call sirius_set_o1_radial_integral(sctx, ias, h1loa(io, ilo, ias), ilo1=ilo, l2=l, o2=io)
                       enddo
                     enddo

                     do ilo1 = 1, nlorb (is)
                       do ilo2 = 1, nlorb (is)
                         !call sirius_set_o1_radial_integral(sctx, ias, h1lolo(ilo1, ilo2, ias), ilo1=ilo1, ilo2=ilo2)
                       enddo
                     enddo
               endif
            end do ! ia
         end do ! is
#else
         stop sirius_error
#endif
  
  endif
  ! ---------------------------- pass overlap integrals to sirius
  
  deallocate(h1aa)
  deallocate(h1loa)
  deallocate(h1lolo)
  
  return
  
end subroutine
!EOC
