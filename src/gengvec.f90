! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.
!BOP
! !ROUTINE: gengvec
! !INTERFACE:
subroutine gengvec
! !USES:
  use modmain
  use modtest
!#ifdef _SIRIUS_
!  use mod_sirius
!#endif
! !DESCRIPTION:
!   Generates a set of ${\bf G}$-vectors used for the Fourier transform of the
!   charge density and potential and sorts them according to length. Integers
!   corresponding to the vectors in lattice coordinates are stored, as well as
!   the map from these integer coordinates to the ${\bf G}$-vector index. A map
!   from the ${\bf G}$-vector set to the standard FFT array structure is also
!   generated. Finally, the number of ${\bf G}$-vectors with magnitude less than
!   {\tt gmaxvr} is determined.
!
! !REVISION HISTORY:
!   Created October 2002 (JKD)
!   Increased number of G-vectors to ngrtot, July 2007 (JKD)
!EOP
!BOC

  implicit none

  ! local variables
  integer ig,i1,i2,i3,j1,j2,j3,k
  real(8) v(3),t1

  ! allocatable arrays
  integer, allocatable :: idx(:)
  integer, allocatable :: iar(:)
  real(8), allocatable :: rar(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (use_sirius_library.and.use_sirius_gvec) then
  
#ifdef _SIRIUS_

    ! ---------------------------------------------------------------
    ! remind myself some EP definitions: 
    !   ivg:   G-vec in integer(lattice) coordi
    !   intgv: integer grid intervals for each direction, integer(3,2)
    !   ivgig: map from integer grid to G-vector array
    !   igfft: map from G-vector array to FFT array
    !   vgc:   G-vec in Cartesian coordi
    !   gc:    length of G-vectors, size=ngvec
    ! ---------------------------------------------------------------
    ! 
    ! 1st, get from sirius: ngvec
    ! 
    ! Note: sirius_get_num_gvec(sctx) will return int sim_ctx.gvec().num_gvec(), 
    !       Sirius documentation says it is "total number of G vectors", but don't confuse it with "ngrtot".
    !       This gives us ngvec rather than ngrtot. 
    ! 
    ngvec = sirius_get_num_gvec(sctx)
    ! 
    ! 2nd, allocate several arries in EP side, using the "ngvec"
    ! 
    If (allocated(ivg)) deallocate (ivg)
    Allocate(ivg(3,ngvec))
    If (allocated(ivgig)) deallocate (ivgig)
    Allocate(ivgig(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),intgv(3,1):intgv(3, 2)))
    If (allocated(igfft)) deallocate (igfft)
    Allocate(igfft(ngvec))
    If (allocated(vgc)) deallocate (vgc)
    Allocate(vgc(3,ngvec))
    If (allocated(gc)) deallocate (gc)
    Allocate(gc(ngvec))
    ! 
    ! 3rd, get from sirius: several G-vec arraies;
    !      get from sirius: the mapping between G-vec index and FFT index;
    ! 
    call sirius_get_gvec_arrays(sctx, ivg(1,1), vgc(1,1), gc(1), ivgig(intgv(1,1), intgv(2,1), intgv(3,1)))
    call sirius_get_fft_index(sctx, igfft(1))

              !write(*,*)'  '    
              !write(*,*)' debug flag, gengvec, 1 '
              !write(*,'(" debug flag, gengvec, ngrtot(from EP)     = ", I10    )') ngrtot
              !write(*,'(" debug flag, gengvec, ngvec(from SIRIUS)  = ", I10    )') ngvec
              !write(*,*)'  ' 
    
#else
    stop sirius_error
#endif

  else

    ! allocate local arrays
    allocate(idx(ngrtot))
    allocate(iar(ngrtot))
    allocate(rar(ngrtot))
    ! allocate global G-vector arrays
    if (allocated(ivg)) deallocate(ivg)
    allocate(ivg(3,ngrtot))
    if (allocated(ivgig)) deallocate(ivgig)
    allocate(ivgig(intgv(1,1):intgv(1,2),intgv(2,1):intgv(2,2),intgv(3,1):intgv(3,2)))
    if (allocated(igfft)) deallocate(igfft)
    allocate(igfft(ngrtot))
    if (allocated(vgc)) deallocate(vgc)
    allocate(vgc(3,ngrtot))
    if (allocated(gc)) deallocate(gc)
    allocate(gc(ngrtot))
    ! find length of each G-vector
    ig=0
    do i1=intgv(1,1),intgv(1,2)
      do i2=intgv(2,1),intgv(2,2)
        do i3=intgv(3,1),intgv(3,2)
          v(:)=dble(i1)*bvec(:,1)+dble(i2)*bvec(:,2)+dble(i3)*bvec(:,3)
          t1=v(1)**2+v(2)**2+v(3)**2
          ig=ig+1
          ! map from G-vector to (i1,i2,i3) index
          ivg(1,ig)=i1
          ivg(2,ig)=i2
          ivg(3,ig)=i3
          ! length of each G-vector
          gc(ig)=sqrt(t1)
        end do
      end do
    end do
    ! sort by vector length
    call sortidx(ngrtot,gc,idx)
    ! re-order arrays
    do ig=1,ngrtot
      rar(ig)=gc(ig)            ! copy gc to rar
    end do
    do ig=1,ngrtot 
      gc(ig)=rar(idx(ig))       ! re-assign rar to gc
    end do
    do k=1,3
      do ig=1,ngrtot
        iar(ig) = ivg(k,ig)       ! copy ivg to iar
      end do
      do ig=1,ngrtot
        ivg(k,ig) = iar(idx(ig))  ! re-assign iar to ivg
      end do
    end do
    ivgig(:,:,:) = 0
    do ig=1,ngrtot         ! loop ig, this gives ivgig,vgc,igfft, from sorted ivg
      i1 = ivg(1,ig)
      i2 = ivg(2,ig)
      i3 = ivg(3,ig)
      ! map from (i1,i2,i3) index to G-vector
      ivgig(i1,i2,i3) = ig
      ! assign G-vectors to global array
      vgc(:,ig) = dble(i1) * bvec(:,1) + dble(i2) * bvec(:,2) + dble(i3) * bvec(:,3)
      ! Fourier transform index
      if (i1.ge.0) then
        j1=i1
      else
        j1=ngrid(1)+i1
      end if
      if (i2.ge.0) then
        j2=i2
      else
        j2=ngrid(2)+i2
      end if
      if (i3.ge.0) then
        j3=i3
      else
        j3=ngrid(3)+i3
      end if
      igfft(ig)=j3*ngrid(2)*ngrid(1)+j2*ngrid(1)+j1+1
    end do ! loop ig
    ! find the number of vectors with G < gmaxvr
    ngvec=1
    do ig=ngrtot,1,-1
      if (gc(ig).lt.gmaxvr) then
        ngvec = ig
        goto 10
      end if
    end do
10  continue
    ! write number of G-vectors to test file
    call writetest(900,'number of G-vectors',iv=ngvec)
    ! 
    deallocate(idx,iar,rar)


              !write(*,*)'  '    
              !write(*,*)' debug flag, gengvec, 2 '
              !write(*,'(" debug flag, gengvec, ngrtot = ", I10    )') ngrtot
              !write(*,'(" debug flag, gengvec, ngvec  = ", I10    )') ngvec
              !write(*,*)'  ' 

  end if
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

              write(*,*)' -------------------------- '    
              write(*,*)' debug flag, gengvec done. '
              write(*,*)' -------------------------- ' 


  return

end subroutine
!EOC

