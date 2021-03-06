
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: poteff
! !INTERFACE:

subroutine poteff

  use modmain
!#ifdef _SIRIUS_
!  use mod_sirius
!#endif

! !DESCRIPTION:
!   Computes the effective potential by adding together the Coulomb and
!   exchange-correlation potentials. See routines {\tt potcoul} and {\tt potxc}.
!
! !REVISION HISTORY:
!   Created April 2003 (JKD)
!EOP
!BOC

  implicit none

  ! local variables
  integer is,ia,ias,ir,offs,i,j
  real(8) ts0,ts1
  call timesec(ts0)
  call timer_start(t_pot)

  ! -------------------------------------------- 
  ! compute the Coulomb potential
  ! --------------------------------------------     
  
  if (use_sirius_library.and.use_sirius_vha) then
#ifdef _SIRIUS_
      ! offsets, I don't understand this line ......
      !offs = sum(ngr_loc_all(0:sirius_fft_comm_rank)) - ngr_loc + 1   
      offs = 1
      ! generate Coulomb potential by solving Poisson equation with SIRIUS. 
      call sirius_generate_coulomb_potential(gs_handler, vclmt(1,1,1), vclir(offs))
      !call gatherir(vclir(1))
      ! Exciting gets "vmad" from SIRIUS, do I need this? 
      call sirius_get_vha_el(gs_handler, vmad(1))
#else
      stop sirius_error
#endif
  else
      call potcoul
      ! Exciting only: potcoul calls zpotcoul to solve Poisson equation; vmad gets assigned value in zpotcoul; 
  endif
  
  ! --------------------------------------------
  ! compute the exchange-correlation potential
  ! --------------------------------------------
  
  if (use_sirius_library.and.use_sirius_vxc) then
#ifdef _SIRIUS_
      ! offs 
      !offs = sum(ngr_loc_all(0:sirius_fft_comm_rank)) - ngr_loc + 1
      offs = 1
      ! generate XC potential with SIRIUS.
      if (ndmag.eq.0) then
        call sirius_generate_xc_potential(gs_handler, bool(.false.), vxcmt(1,1,1), vxcir(1))
      else if (ndmag.eq.1) then
        call sirius_generate_xc_potential(gs_handler, bool(.false.), vxcmt(1,1,1), vxcir(1),&
                                         &bxcmt_z=bxcmt(1,1,1,1), bxcrg_z=bxcir(1,1))
      else 
        call sirius_generate_xc_potential(gs_handler, bool(.false.), vxcmt(1,1,1), vxcir(1),&
                                         &bxcmt_x=bxcmt(1,1,1,1), bxcrg_x=bxcir(1,1),&
                                         &bxcmt_y=bxcmt(1,1,1,2), bxcrg_y=bxcir(1,2),&
                                         &bxcmt_z=bxcmt(1,1,1,3), bxcrg_z=bxcir(1,3))

      endif
      !call gatherir(vxcir(1))
      !do j = 1, ndmag
      !  call gatherir(bxcir(1, j))
      !enddo
#else
      stop sirius_error
#endif
  else
      call potxc
  endif

  ! EXCITING only: shift interstitial part to reference value
  ! vclir(:) = vclir(:) + shift  
  
  ! --------------------------------------------
  ! compute vxc_val
  ! --------------------------------------------
  ! rho_val : flag for vxcval,               default false in readinput.f90
  ! pt_core : flag for partial core density, default false in readinput.f90
  if (rho_val.or.pt_core) then 
    call potxc_val
  endif
  
  ! ------------------------------------------------------------
  ! add Coulomb and exchange-correlation potentials together: 
  ! ------------------------------------------------------------

  ! muffin-tin part
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      do ir=1,nrmt(is)
        veffmt(:,ir,ias)=vclmt(:,ir,ias)+vxcmt(:,ir,ias)
      end do
    end do
  end do

  ! interstitial part
  veffir(:)=vclir(:)+vxcir(:)
  call timer_stop(t_pot)
  call timesec(ts1)
  timepot=timepot+ts1-ts0
  return

end subroutine
!EOC
