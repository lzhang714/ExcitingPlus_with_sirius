subroutine sic_genvme_dotp(t0only)
use modmain
use mod_sic
implicit none
logical, intent(in) :: t0only
integer nwtloc,iloc,i,n,j,n1,j1,vl(3),jp,ispn
real(8) pos1(3),pos2(3)
complex(8), allocatable :: wvtp(:,:,:)

! compute matrix elements of SIC potential
!  sic_vme(i) = <(W*V)_n|W_{n1,T}>
!  i = {n,n1,T}
allocate(wvtp(s_ntp,s_nr_min,nspinor))
sic_vme=zzero
if (t0only) then
  nwtloc=mpi_grid_map(sic_wantran%nwt0,dim_k)
else
  nwtloc=mpi_grid_map(sic_wantran%nwt,dim_k)
endif
jp=-1
do iloc=1,nwtloc
  if (t0only) then
    j=mpi_grid_map(sic_wantran%nwt0,dim_k,loc=iloc)  
    i=sic_wantran%iwt0(j)
  else
    i=mpi_grid_map(sic_wantran%nwt,dim_k,loc=iloc)
  endif
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  pos1(:)=wanpos(:,n)
  pos2(:)=wanpos(:,n1)+vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  if (j.ne.jp) then
    jp=j
    do ispn=1,nspinor
! convert to spherical coordinates
      call zgemm('T','N',s_ntp,s_nr_min,lmmaxwan,zone,s_ylmf,lmmaxwan,&
        &s_wvlm(1,1,ispn,j),lmmaxwan,zzero,wvtp(1,1,ispn),s_ntp)
    enddo
  endif
  sic_vme(i)=s_spinor_dotp(pos1,pos2,s_wvlm(1,1,1,j),wvtp,s_wlm(1,1,1,j1))
! filter small matrix elements
  if (abs(sic_vme(i)).lt.1d-9) sic_vme(i)=zzero
enddo
call mpi_grid_reduce(sic_vme(1),sic_wantran%nwt,dims=(/dim_k/),all=.true.)
deallocate(wvtp)
return
end subroutine
