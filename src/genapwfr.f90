
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: genapwfr
! !INTERFACE:
subroutine genapwfr
! !USES:
  use modmain
!#ifdef _SIRIUS_
!  use mod_sirius
!#endif
! !DESCRIPTION:
!   Generates the APW radial functions. This is done by integrating the scalar
!   relativistic Schr\"{o}dinger equation (or its energy deriatives) at the
!   current linearisation energies using the spherical part of the effective
!   potential. The number of radial functions at each $l$-value is given by the
!   variable {\tt apword} (at the muffin-tin boundary, the APW functions have
!   continuous derivatives up to order ${\tt apword}-1$). Within each $l$, these
!   functions are orthonormalised with the Gram-Schmidt method. The radial
!   Hamiltonian is applied to the orthonormalised functions and the results are
!   stored in the global array {\tt apwfr}.
!
! !REVISION HISTORY:
!   Created March 2003 (JKD)
!   Copied to equivalent atoms, February 2010 (A. Kozhevnikov and JKD)
!EOP
!BOC

  implicit none

  ! local variables
  integer is,ia,ja,ias,jas
  integer nr,ir,nn,l,io,jo
  real(8) t1
  
  ! automatic arrays
  logical done(natmmax)
  real(8) vr(nrmtmax),fr(nrmtmax),gr(nrmtmax),cf(4,nrmtmax)
  real(8) p0(nrmtmax,apwordmax),p1(nrmtmax),p1s(apwordmax)    ! LZ: Exciting has: p1 (nrmtmax, apwordmax)
  real(8) q0(nrmtmax,apwordmax),q1(nrmtmax,apwordmax)
  real(8) hp0(nrmtmax)

  do is=1,nspecies
    nr=nrmt(is)
    done(:)=.false.
    do ia=1,natoms(is)
      if (.not.done(ia)) then
        ias=idxas(ia,is)

        ! use spherical part of potential
        vr(1:nr)=veffmt(1,1:nr,ias)*y00
        do l=0,lmaxapw
          do io=1,apword(l,is)           ! note io1,io2 vs io,jo
            
            ! Integrate the radial Schrodinger equation. "sirius_radial_solver" no longer available....
            
!            if (use_sirius_library.and.use_sirius_radial_solver) then
!#ifdef _SIRIUS_
!              call sirius_radial_solver("none", nint(-spzn(is)), apwdm(io1, l, is), l, 0,&
!                                            &apwe(io,l,ias), nr, spr(1, is), vr(1), nn,&
!                                            &p0(1, io1), p1(1, io1), q0(1, io1), q1(1, io1))
!#endif
!            else

              call rschroddme(solsc,apwdm(io,l,is), l, 0, apwe(io,l,ias), nprad, nr, spr(:,is), &
                              &vr, nn, p0(:,io), p1, q0(:,io), q1(:,io))

              ! original EP and Exciting calls, put here for reference.
              !
              ! EP:       call rschroddme(solsc,apwdm( io,l,is), l, 0, apwe( io, l, ias), nprad,nr, &
              ! EP:                       spr(:,is),vr,nn,p0(:,io),p1,q0(:,io),q1(:,io))
              !
              ! Exciting: Call rschroddme(      apwdm(io1,l,is), l, 0, apwe(io1, l, ias), nr, spr(:, is), &
              ! Exciting:                    &vr, nn, p0(:, io1), p1(:, io1), q0(:, io1), q1(:, io1))
              !
              ! some differences: EP uses "io" while Exciting uses "io1";
              !                   EP has "solsc" and "nprad" in the call while Exciting doesn't have them; 
              !                   "p1" is 1D array in EP while it's 2D array in Exciting;  

!            end if    
            
            ! normalise radial functions
            fr(1:nr)=p0(1:nr,io)**2
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            t1=1.d0/sqrt(abs(gr(nr)))
            p0(1:nr,io)=t1*p0(1:nr,io)
            p1s(io)=t1*p1(nr)
            q0(1:nr,io)=t1*q0(1:nr,io)
            q1(1:nr,io)=t1*q1(1:nr,io)
            ! subtract linear combination of previous vectors
            do jo=1,io-1
              fr(1:nr)=p0(1:nr,io)*p0(1:nr,jo)
              call fderiv(-1,nr,spr(:,is),fr,gr,cf)
              t1=gr(nr)
              p0(1:nr,io)=p0(1:nr,io)-t1*p0(1:nr,jo)
              p1s(io)=p1s(io)-t1*p1s(jo)
              q0(1:nr,io)=q0(1:nr,io)-t1*q0(1:nr,jo)
              q1(1:nr,io)=q1(1:nr,io)-t1*q1(1:nr,jo)
            end do ! jo
            ! normalise radial functions
            fr(1:nr)=p0(1:nr,io)**2
            call fderiv(-1,nr,spr(:,is),fr,gr,cf)
            t1=abs(gr(nr))
            if (t1.lt.1.d-20) then
              write(*,*)
              write(*,'("Error(genapwfr): degenerate APW radial functions")')
              write(*,'(" for species ",I4)') is
              write(*,'(" atom ",I4)') ia
              write(*,'(" angular momentum ",I4)') l
              write(*,'(" and order ",I4)') io
              write(*,*)
              stop
            end if
            t1=1.d0/sqrt(t1)
            p0(1:nr,io)=t1*p0(1:nr,io)
            p1s(io)=t1*p1s(io)
            q0(1:nr,io)=t1*q0(1:nr,io)
            q1(1:nr,io)=t1*q1(1:nr,io)
            ! apply the Hamiltonian    
            ! LZ, here we calculate hp0 used below, EXCITING uses p1 & p0 instead
            call rschrodapp(solsc,l,nr,spr(:,is),vr,p0(:,io),q0(:,io),q1(:,io),hp0)
            ! divide by r and store in global array
            do ir=1,nr
              t1=1.d0/spr(ir,is)
              apwfr(ir,1,io,l,ias)=t1*p0(ir,io)
              apwfr(ir,2,io,l,ias)=t1*hp0(ir)
            end do
            ! derivative at the muffin-tin surface
            apwdfr(io,l,ias)=(p1s(io)-p0(nr,io)*t1)*t1
            
            ! pass apwfr to sirius.
            ! should I put this after copying to equivalent atoms ?             
            if (use_sirius_library.and.pass_apwfr_to_sirius) then
#ifdef _SIRIUS_
              call sirius_set_radial_function(sctx, ias, 0, apwfr(1, 1, io, l, ias), l=l, o=io)
              call sirius_set_radial_function(sctx, ias, 1, apwfr(1, 2, io, l, ias), l=l, o=io)
#else
              stop sirius_error
#endif
            endif
            
          end do  ! io
        end do  ! l
        
        
        done(ia)=.true. 
           
        ! copy to equivalent atoms
        do ja=1,natoms(is)
          if ((.not.done(ja)).and.(eqatoms(ia,ja,is))) then
            jas=idxas(ja,is)
            do l=0,lmaxapw
              do io=1,apword(l,is)
                apwfr(1:nr,1,io,l,jas)=apwfr(1:nr,1,io,l,ias)
                apwfr(1:nr,2,io,l,jas)=apwfr(1:nr,2,io,l,ias)
                apwdfr(io,l,jas)=apwdfr(io,l,ias)
              end do
            end do
            done(ja)=.true.
          end if
        end do      
      
      end if  ! if (.not.done(ia))
    end do  ! ia 
  end do  ! is

  return

end subroutine
!EOC
