subroutine sic_wan_pot(fout)
use modmain
use mod_sic
use mod_nrkp
use mod_wannier
use modxcifc
use mod_util
implicit none
! arguments
integer, intent(in) :: fout
! local variables
integer n,ispn,vl(3),n1,i,j,j1,itp,ir,nwtloc,iloc,nrloc,irloc,nrntp,nrntploc
real(8) t1,t2,vrc(3),pos1(3),pos2(3)
real(8) sic_epot_h,sic_epot_xc
complex(8) z1,wanval(nspinor)
real(8), allocatable :: wanprop(:,:)
real(8), allocatable :: f1tp(:),f2lm(:) 
complex(8), allocatable :: ovlp(:)
complex(8), allocatable :: om(:,:)
complex(8), allocatable :: wantp(:,:,:)
complex(8), allocatable :: wantp_tmp(:,:,:)
complex(8), external :: zdotc
!
allocate(wantp(s_ntp,s_nr,nspinor))
allocate(wantp_tmp(s_ntp,s_nr,nspinor))
allocate(wanprop(nwanprop,sic_wantran%nwan))
allocate(f1tp(s_nr),f2lm(s_nr))
! check spherical and radial grids by integrating sphere volume
if (wproc) then
  write(fout,*)
  write(fout,'(80("="))')
  write(fout,'("generating Wannier functions on a spherical grid")')
  write(fout,'(80("="))')
  write(fout,*)
  t1=0.d0
  t2=fourpi*s_rmax**3/3.d0
  do ir=1,s_nr
    do itp=1,s_ntp
      t1=t1+s_tpw(itp)*s_rw(ir)/t2
    enddo
  enddo
  write(fout,'("spherical grid error : ",G18.10)')abs(1.d0-t1)
  call flushifc(fout)
endif

call timer_start(t_sic_wan,reset=.true.)
call timer_reset(t_sic_wan_gen)
call timer_reset(t_sic_wan_rms)
call timer_reset(t_sic_wan_pot)
call timer_reset(t_sic_wvprod)
s_wlm=zzero
s_wvlm=zzero
wanprop=0.d0
!nrloc=mpi_grid_map(s_nr,dim_k)
nrntp=s_nr*s_ntp
nrntploc=mpi_grid_map(nrntp,dim_k)

! loop over Wannier functions
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
! generate WF on a spherical mesh
  wantp_tmp=zzero
  call timer_start(t_sic_wan_gen)
  call elk_load_wann_unk(n)
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(j1,ir,itp,vrc,wanval)
  do i=1,nrntploc
    j1=mpi_grid_map(nrntp,dim_k,loc=i)
    ir=int((j1-1)/s_ntp)+1
    itp=mod(j1-1,s_ntp)+1
    vrc(:)=s_x(:,itp)*s_r(ir)+m_wanpos(:)
    call get_wanval(vrc,wanval)
    wantp_tmp(itp,ir,:)=wanval(:)
  enddo !i
!$OMP END PARALLEL DO
  call mpi_grid_reduce(wantp_tmp(1,1,1),s_ntp*s_nr*nspinor,wantp(1,1,1),dims=(/dim_k/),all=.true.)
  if (wproc) then
    write(fout,'("n : ",I4, "   image part (average, total) : ",2G18.10)') n, &
      &sum(abs(dimag(wantp)))/s_ntp/s_nr/nspinor, &
      &sum(abs(dimag(wantp)))
  endif
! convert to spherical harmonics
  call sic_genwanlm(fout,n,wantp,s_wlm(1,1,1,j))
! compute norm
  f1tp=0.d0
  f2lm=0.d0
  do ispn=1,nspinor
    do ir=1,s_nr_min
      do itp=1,s_ntp
        f1tp(ir)=f1tp(ir)+abs(wantp(itp,ir,ispn))**2*s_tpw(itp)
      enddo
      f2lm(ir)=f2lm(ir)+abs(zdotc(lmmaxwan,s_wlm(1,ir,ispn,j),1,s_wlm(1,ir,ispn,j),1))
    enddo
  enddo
  wanprop(wp_normtp,j)=rintegrate(s_nr_min,s_r,f1tp)
  wanprop(wp_normlm,j)=rintegrate(s_nr_min,s_r,f2lm)
  call timer_stop(t_sic_wan_gen)
! generate potential
  if (sic_apply(n).eq.3) then
    call timer_start(t_sic_wan_pot)
    call sic_genpot(n,s_wlm(1,1,1,j),s_wvlm(1,1,1,j),wanprop(1,j))
    call timer_stop(t_sic_wan_pot)
  endif
enddo !j
deallocate(wantp,wantp_tmp,f1tp,f2lm)

!#ifdef _MAD_
!call madness_gen_hpot
!call madness_gen_xc
!#endif

! compute overlap integrals 
allocate(ovlp(sic_wantran%nwt))
ovlp=zzero
t1=0.d0
call timer_start(t_sic_wan_ovl,reset=.true.)
nwtloc=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/))
do iloc=1,nwtloc
  i=mpi_grid_map2(sic_wantran%nwt,dims=(/dim_k,dim2/),loc=iloc)
  n=sic_wantran%iwt(1,i)
  j=sic_wantran%idxiwan(n)
  n1=sic_wantran%iwt(2,i)
  j1=sic_wantran%idxiwan(n1)
  vl(:)=sic_wantran%iwt(3:5,i)
  pos1(:)=wanpos(:,n)
  pos2(:)=wanpos(:,n1)+vl(1)*avec(:,1)+vl(2)*avec(:,2)+vl(3)*avec(:,3)
  if (sum(abs(pos1-pos2)).lt.1d-10) then
    ovlp(i)=s_spinor_dotp_lm(s_wlm(1,1,1,j),s_wlm(1,1,1,j1))
  endif
  z1=ovlp(i)
  if (n.eq.n1.and.all(vl.eq.0)) then
    z1=z1-zone
  endif
  t1=max(t1,abs(z1))
enddo
call mpi_grid_reduce(ovlp(1),sic_wantran%nwt)
call mpi_grid_reduce(t1,op=op_max)
call timer_stop(t_sic_wan_ovl)
! print overlap matrix for T=0
if (wproc) then
  allocate(om(sic_wantran%nwan,sic_wantran%nwan))
  om=zzero
  do i=1,sic_wantran%nwt
    n=sic_wantran%iwt(1,i)
    j=sic_wantran%idxiwan(n)
    n1=sic_wantran%iwt(2,i)
    j1=sic_wantran%idxiwan(n1)
    vl(:)=sic_wantran%iwt(3:5,i)
    if (all(vl.eq.0)) then
      om(j,j1)=ovlp(i)
    endif
  enddo
  write(fout,*)
  write(fout,'("overlap matrix for T=0")')
  do j=1,sic_wantran%nwan
    write(fout,'(255F12.6)')(abs(om(j,j1)),j1=1,sic_wantran%nwan)
  enddo
  deallocate(om)
endif
! compute energies
! note: here Hartree potential has a positive sign and XC potential 
!  has a negative sign
sic_energy_kin=0.d0
sic_epot_h=0.d0
sic_epot_xc=0.d0
do j=1,sic_wantran%nwan
  sic_energy_kin=sic_energy_kin+wanprop(wp_vsic,j)
  sic_epot_h=sic_epot_h+0.5d0*wanprop(wp_vha,j)
  sic_epot_xc=sic_epot_xc+wanprop(wp_exc,j)
enddo
sic_energy_pot=sic_epot_h+sic_epot_xc
! total energy engytot=engytot+sic_energy_tot
sic_energy_tot=sic_energy_kin-sic_energy_pot
! print some info
if (wproc) then
  write(fout,*)
  write(fout,'("   n |     normlm      normtp      spread&
    &     normrho       <x>         <y>         <z>")')
  write(fout,'(80("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(151,'(I4," | ",7(F10.6,2X))')n,wanprop(wp_normlm,j),&
      &wanprop(wp_normtp,j),wanprop(wp_spread,j),wanprop(wp_normrho,j),&
      &wanprop(wp_spread_x,j),wanprop(wp_spread_y,j),wanprop(wp_spread_z,j)
  enddo
  write(fout,'(80("-"))')
  write(fout,'("maximum deviation from norm : ",F12.6)')t1
  write(fout,'("total quadratic spread : ",F12.6," [a.u.]^2")')&
    &sum(wanprop(wp_spread,:))
  write(fout,*)
  write(fout,'("   n | ",5X,"V_n^{H}     V_n^{XC}          V_n     E_n^{XC}&
    &        Ex           Ec        E_n^{SIC}")')
  write(fout,'(98("-"))')
  do j=1,sic_wantran%nwan
    n=sic_wantran%iwan(j)
    write(fout,'(I4," | ",7(F12.6,1X))')n,wanprop(wp_vha,j),&
      &wanprop(wp_vxc,j),wanprop(wp_vsic,j),wanprop(wp_exc,j),&
      &wanprop(wp_ex,j),wanprop(wp_ec,j),(0.5d0*wanprop(wp_vha,j)+wanprop(wp_exc,j))
  enddo
  write(fout,'(98("-"))')
  write(fout,'("SIC total energy contribution      : ",G18.10,&
    &"  ! kinetic - potential ")')sic_energy_tot
  write(fout,'("SIC kinetic energy contribution    : ",G18.10,&
    &"  ! sum of V_n" )')sic_energy_kin
  write(fout,'("SIC potential energy contribution  : ",G18.10,&
    &"  ! sum of Hartree and XC terms")')sic_energy_pot
  write(fout,'("  Hartree                          : ",G18.10,&
    &"  ! sum of V_n^{H}/2")')sic_epot_h
  write(fout,'("  XC                               : ",G18.10,&
    &"  ! sum of E_n^{XC}")')sic_epot_xc 
  call flushifc(fout)
endif
deallocate(ovlp)
! save some info
xml_info%wan_spread=0.d0
do j=1,sic_wantran%nwan
  n=sic_wantran%iwan(j)
  xml_info%wan_spread(n)=wanprop(wp_spread,j)
enddo
xml_info%wan_tot_spread=sum(wanprop(wp_spread,:))
xml_info%sic_energy_tot=sic_energy_tot
xml_info%sic_energy_pot=sic_energy_pot
xml_info%sic_energy_kin=sic_energy_kin
deallocate(wanprop)
call timer_stop(t_sic_wan)
if (wproc) then
  write(fout,*)
  write(fout,'("total time          : ",F8.3," sec.")')timer_get_value(t_sic_wan)
  write(fout,'("  generation of WFs : ",F8.3," sec.")')timer_get_value(t_sic_wan_gen)
  write(fout,'("  potential         : ",F8.3," sec.")')timer_get_value(t_sic_wan_pot)
  write(fout,'("   (W*V)            : ",F8.3," sec.")')timer_get_value(t_sic_wvprod)
  write(fout,'("  overlap           : ",F8.3," sec.")')timer_get_value(t_sic_wan_ovl)
  call flushifc(fout)
endif

return
end
