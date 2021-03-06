subroutine gw_mapping(iqrmap,qqnrmap,agqmap,akmap)
use modmain
use mod_addons_q

implicit none
integer,intent(inout) :: iqrmap(2,nkpt+nvq0-1)
integer,intent(inout) :: qqnrmap(nvq,2,nkpt+nvq0-1)
integer,intent(inout) :: agqmap(ngqmax,nvq)
integer,intent(inout) :: akmap(nsymcrys,nkptnr)
integer :: ist1, ist2, ig, ig1,igq,igqtmp,iv(3)
integer :: iqtmp, iktmp, iq1,iqibz
integer :: ik,isym,ik1,lspl,iq
real(8) :: s(3,3),v(3),v1(3),v2(3),t1,vak(3)
real(8) :: vibz(3),vbz(3),avbz(3),vg(3),vag(3),vqg0ag(3),vgq1(3)
integer :: vg0(3)

! At this moment, assume no parallelization for q-mesh, i.e. dim_q = 1
!

! iqrmap(1,iqibz): index of q_{IBZ} in the non-reduced q-mesh;
! iqrmap(2,iqibz): number of the set {q_{BZ}} related to this q_{IBZ}.
iqrmap(:,:)=0

! qqnrmap: set of q_{BZ} points generated by q_{IBZ} through 
!         q_{IBZ}+G0 = a*q_{BZ} where {a} are the symmetry operations. 
!         Note that q_{BZ}s are in non-reduced q-mesh. 
! qqnrmap(iqbz,1,iqibz): index of each q_{BZ} point in the
!                        non-reduced q-mesh, regarding q_{IBZ};
! qqnrmap(iqbz,2,iqibz): index of the symmetry operation "a" that yields 
!                        a*q_{BZ} = q_{IBZ}+G0
qqnrmap(:,:,:)=0

! agqmap(igq,iqbz): for a given q_{BZ}, maps G+q_{BZ} to (q_{IBZ}+G0+a*G)
!                   = G'+q_{IBZ}, where a is the symmetry operation.
agqmap(:,:)=0

! akmap(isym,ik): mapping non-reduced k to a*k in non-reduced k-mesh
akmap(:,:)=0

!print out the mapping between vklnr and vkl
if (mpi_grid_root()) then
 open(200,file='symmetryk.dat',form='formatted',status='replace')
 do ik=1,nkptnr
  call findkpt(vklnr(:,ik),isym,ik1) ! ik1, index of k_{IBZ}
  write(200,*) "ik:",ik
  write(200,*) "vklnr(:,ik):"
  write(200,*) vklnr(:,ik)
  write(200,*) "ik1:",ik1
  write(200,*) "vkl(:,ik1):"
  write(200,*) vkl(:,ik1)
  lspl=lsplsymc(isym)
  write(200,*) "lspl:",lspl
  write(200,*) "symlat(:,:,lspl):"
  write(200,*) dble(symlat(1,:,lspl))
  write(200,*) dble(symlat(2,:,lspl))
  write(200,*) dble(symlat(3,:,lspl))
  write(200,*)
 enddo
 call flushifc(200)
 close(200)

 open(201,file='symmetryq.dat',form='formatted',status='replace')
 open(202,file='agqmap.dat',form='formatted',status='replace')
 open(203,file='akmap.dat',form='formatted',status='replace')
endif

! mapping between q_{BZ} and q_{IBZ}
do iq=1,nvq
 if (vq_gamma(iq)) then
  ik1=1
  iqrmap(1,iq)=iq               ! index of qibz in the non-reduced q-mesh
  iqrmap(2,iq)=1                ! Nqbz for each q_{IBZ}
  ! 1:index of non-reduced q, 1:nvq0 are gamma points
  qqnrmap(iqrmap(2,iq),1,iq)=iq
  ! 2:index of identity operation
  qqnrmap(iqrmap(2,iq),2,iq)=1
 else
  call findkpt(vql(:,iq),isym,ik1)   ! ik1 > 1, a*q_{BZ} = q_{IBZ} + G0
  iqrmap(2,ik1+nvq0-1)=iqrmap(2,ik1+nvq0-1)+1
  qqnrmap(iqrmap(2,ik1+nvq0-1),1,ik1+nvq0-1)=iq
  qqnrmap(iqrmap(2,ik1+nvq0-1),2,ik1+nvq0-1)=lsplsymc(isym)
  if (lsplsymc(isym).eq.1) iqrmap(1,ik1+nvq0-1)=iq
 endif

 if (mpi_grid_root()) then
  write(201,*) "iq:",iq
  write(201,*) "vql(:,iq):"
  write(201,*) vql(:,iq)
  if (.not. vq_gamma(iq)) then
   write(201,*) "ik1+nvq0-1:",ik1+nvq0-1
  else
   write(201,*) "ik1+nvq0-1:",iq
  endif
  write(201,*) "vkl(:,ik1):"
  write(201,*) vkl(:,ik1)
  if (.not. vq_gamma(iq)) then
   write(201,*) "qqnrmap(nq,2,iqibz):",&
              &qqnrmap(iqrmap(2,ik1+nvq0-1),2,ik1+nvq0-1)  ! isym
   lspl=qqnrmap(iqrmap(2,ik1+nvq0-1),2,ik1+nvq0-1)
  else
   write(201,*) "qqnrmap(nq,2,iq):",qqnrmap(iqrmap(2,iq),2,iq)
   lspl=1
  endif
  write(201,*) "symlat(:,:,lspl):"
  write(201,*) dble(symlat(1,:,lspl))
  write(201,*) dble(symlat(2,:,lspl))
  write(201,*) dble(symlat(3,:,lspl))
  write(201,*)
  call flushifc(201)
 endif
enddo

!! find out qBZ+G -> qIBZ+G' = qIBZ + G0 + a*G
do iq=1,nkpt+nvq0-1
 iqibz=iqrmap(1,iq)  ! index of qIBZ in non-reduced q-mesh
 vibz(:)=vql(:,iqibz) ! qIBZ vector
 do iq1=1,iqrmap(2,iq)  ! n(qIBZ)
  iqtmp=qqnrmap(iq1,1,iq)  !index of qBZ
  vbz(:)=vql(:,iqtmp)   ! qBZ vector
  lspl=qqnrmap(iq1,2,iq)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,vbz,avbz)
  v1(:)=avbz(:)-vibz(:)
  call r3frac(epslat,v1,vg0)  ! v2: G0

  do ig=1,ngq(iqtmp)
   vg(:)=dble(ivg(:,igqig(ig,iqtmp)))   ! G
   call r3mtv(s,vg,vag)    ! aG
   vqg0ag(:)=vibz(:)+vag(:)+dble(vg0(:))   ! qIBZ + aG + G0

   do ig1=1,ngq(iqibz)
    vgq1(:)=vibz(:)+dble(ivg(:,igqig(ig1,iqibz)))  ! G'+qIBZ
    t1=abs(vgq1(1)-vqg0ag(1))+abs(vgq1(2)-vqg0ag(2))+abs(vgq1(3)-vqg0ag(3))
    if (t1.lt.epslat) then
     agqmap(ig,iqtmp)=ig1   ! index for G'+qIBZ
     exit
    endif
   enddo

   if (mpi_grid_root()) then
    if ((agqmap(ig,iqtmp).le.0).or.(agqmap(ig,iqtmp).gt.ngq(iqibz))) then
     write(202,*) "cannot find ig1!"
     write(202,*) "aG+G0+qIBZ:",vqg0ag(:)
     write(202,*) "iqbz:",iqtmp
     write(202,*) "ngq(qibz),ngq(qbz):",ngq(iqibz),ngq(iqtmp)
    else
     write(202,*) "iqbz:",iqtmp
     write(202,*) "ig:",ig  ! qBZ + G
     write(202,*) "aG+G0+qIBZ:",vqg0ag(:)
     write(202,*) "G:"
     write(202,*) ivg(:,igqig(ig,iqtmp))   ! G
     write(202,*) "lspl:",lspl
     write(202,*) "agqmap(ig,iqtmp):",agqmap(ig,iqtmp)  ! G'+qIBZ
     write(202,*) "aG+G0:"
     write(202,*) ivg(:,igqig(agqmap(ig,iqtmp),iqibz))  ! G'
     write(202,*)
    endif
   endif
  enddo !ig
 enddo !iq1
enddo !iq

!! find out ak -> k'
do ik=1,nkptnr
 v1(:)=vklnr(:,ik)    ! kBZ
 do isym=1,nsymcrys
  lspl=lsplsymc(isym)
  s(:,:)=dble(symlat(:,:,lspl))
  call r3mtv(s,v1,vak)   ! ak
  call r3frac(epslat,vak,iv)   ! (a*k)_{BZ}
  do ik1=1,nkptnr
   v1(:)=vklnr(:,ik1)  ! k'
   t1=abs(vak(1)-v1(1))+abs(vak(2)-v1(2))+abs(vak(3)-v1(3))
   if (t1.lt.epslat) then
    akmap(lspl,ik)=ik1   ! index for ak=k'
    exit
   endif
  enddo

  if (mpi_grid_root()) then
   if ((akmap(lspl,ik).eq.0).or.(akmap(lspl,ik).gt.nkptnr)) &
      write(203,*) "cannot find ik1!"
   write(203,*) "ik:",ik
   write(203,*) "vklnr:"
   write(203,*) vklnr(:,ik)
   write(203,*) "lspl:",lspl
   write(203,*) "symlat(:,:,lspl):"
   write(203,*) symlat(1,:,lspl)
   write(203,*) symlat(2,:,lspl)
   write(203,*) symlat(3,:,lspl)
   write(203,*) "akmap(lspl,ik):",akmap(lspl,ik)
   write(203,*) vklnr(:,akmap(lspl,ik))
   write(203,*) "ak:"
   write(203,*) vak(:)
   write(203,*)
  endif
 enddo !isym
enddo !ik

if (mpi_grid_root()) then
 close(201)
 close(202)
 close(203)
endif

return
end subroutine
