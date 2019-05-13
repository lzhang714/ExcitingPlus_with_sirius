
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: gencfun
! !INTERFACE:
subroutine gencfun
! !USES:
  use modmain
!#ifdef _SIRIUS_
!  use mod_sirius
!#endif
! !DESCRIPTION:
!   Generates the smooth characteristic function. This is the function which is
!   0 within the muffin-tins and 1 in the intersitial region and is constructed
!   from radial step function form-factors with $G<G_{\rm max}$. The form
!   factors are given by
!   $$ \tilde{\Theta}_i(G)=\begin{cases}
!    \frac{4\pi R_i^3}{3 \Omega} & G=0 \\
!    \frac{4\pi R_i^3}{\Omega}\frac{j_1(GR_i)}{GR_i} & 0<G\le G_{\rm max} \\
!    0 & G>G_{\rm max}\end{cases}, $$
!   where $R_i$ is the muffin-tin radius of the $i$th species and $\Omega$ is
!   the unit cell volume. Therefore the characteristic function in $G$-space is
!   $$ \tilde{\Theta}({\bf G})=\delta_{G,0}-\sum_{ij}\exp(-i{\bf G}\cdot
!    {\bf r}_{ij})\tilde{\Theta}_i(G), $$
!   where ${\bf r}_{ij}$ is the position of the $j$th atom of the $i$th species.
!
! !REVISION HISTORY:
!   Created January 2003 (JKD)
!EOP
!BOC

  implicit none

  ! local variables
  integer is,ia,ig,ifg,m
  real(8) t1
  complex(8) zt1
  ! allocatable arrays
  real(8), allocatable :: ffacg(:)
  complex(8), allocatable :: zfft(:)
  allocate(ffacg(ngrtot))
  allocate(zfft(ngrtot))

  ! allocate global characteristic function arrays
  if (allocated(cfunig)) deallocate(cfunig)
  allocate(cfunig(ngrtot))
  if (allocated(cfunir)) deallocate(cfunir)
  allocate(cfunir(ngrtot))
  
  
  ! ----------------------------------------------------------------------------    
  if (use_sirius_library.and.use_sirius_cfun) then
    if (allocated(cfunig)) deallocate(cfunig)
    allocate(cfunig(ngvec))

#ifdef _SIRIUS_
    !m = sum(ngr_loc_all(0:sirius_fft_comm_rank)) - ngr_loc + 1
    m = 1
    call sirius_get_step_function(sctx, cfunig(1), cfunir(m))
    !call gatherir(cfunir(1))
#else
    stop sirius_error
#endif

  else

    cfunig(:)=0.d0
    cfunig(1)=1.d0
    ! begin loop over species
    do is=1,nspecies
      ! generate the smooth step function form factors
          write(*,*)' debug flag, gencfun, 1 '
      call genffacg(is,ngrtot,ffacg)
          write(*,*)' debug flag, gencfun, 2 '
      ! loop over atoms
      do ia=1,natoms(is)
        do ig=1,ngrtot
          ! structure factor
          t1=-dot_product(vgc(:,ig),atposc(:,ia,is))
          zt1=cmplx(cos(t1),sin(t1),8)
          ! add to characteristic function in G-space
          cfunig(ig)=cfunig(ig)-zt1*ffacg(ig)
        end do
      end do
    end do

          write(*,*)' debug flag, gencfun, 3 '

    do ig=1,ngrtot
      ifg=igfft(ig)
      zfft(ifg)=cfunig(ig)
    end do

          write(*,*)' debug flag, gencfun, 4 '

    ! Fourier transform to real-space
    call zfftifc(3,ngrid,1,zfft)
    cfunir(:)=dble(zfft(:))
    deallocate(ffacg,zfft)

          write(*,*)' debug flag, gencfun, 5 '

  end if
  ! ----------------------------------------------------------------------------    

  return

end subroutine
!EOC

