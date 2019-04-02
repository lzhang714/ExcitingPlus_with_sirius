
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writeengy(fnum)
use modmain
use modldapu
use mod_sic
implicit none
! arguments
integer, intent(in) :: fnum
write(fnum,*)
write(fnum,'("Energies :")')
write(fnum,'(" Fermi",T30,": ",G22.12)') efermi
write(fnum,'(" sum of eigenvalues",T30,": ",G22.12)') evalsum
write(fnum,'(" electron kinetic",T30,": ",G22.12)') engykn
write(fnum,'("   core kinetic",T30,": ",G22.12)') engykncr
write(fnum,'("   valence kinetic",T30,": ",G22.12)') engykn-engykncr
write(fnum,'(" Coulomb",T30,": ",G22.12)') engycl
write(fnum,'(" Coulomb potential",T30,": ",G22.12)') engyvcl
write(fnum,'(" nuclear-nuclear",T30,": ",G22.12)') engynn
write(fnum,'(" electron-nuclear",T30,": ",G22.12)') engyen
write(fnum,'(" Hartree",T30,": ",G22.12)') engyhar
write(fnum,'(" Madelung",T30,": ",G22.12)') engymad
write(fnum,'(" xc potential",T30,": ",G22.12)') engyvxc
if (spinpol) then
  write(fnum,'(" xc effective B-field",T30,": ",G22.12)') engybxc
  write(fnum,'(" external B-field",T30,": ",G22.12)') engybext
end if
write(fnum,'(" exchange",T30,": ",G22.12)') engyx
write(fnum,'(" correlation",T30,": ",G22.12)') engyc
if (ldapu.ne.0) then
  write(fnum,'(" LDA+U",T30,": ",G22.12)') engylu
end if
if (stype.eq.3) then
  write(fnum,'(" electron entropic",T30,": ",G22.12)') engyts
end if
write(fnum,'(" total energy",T30,": ",G22.12)') engytot
if (sic) then
  write(fnum,'("   E0",T30,": ",G22.12)') engytot0
  write(fnum,'("   sic_energy_pot",T30,": ",G22.12)') sic_energy_pot  
  write(fnum,'("   sic_energy_kin",T30,": ",G22.12)') sic_energy_kin  
endif
if (spinpol) then
  write(fnum,'(" (external B-field energy excluded from total)")')
end if
call flushifc(fnum)
return
end subroutine

