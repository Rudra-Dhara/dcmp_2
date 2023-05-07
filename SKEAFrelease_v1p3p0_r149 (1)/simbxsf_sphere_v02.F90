program simbxsf

!----------------------------------------------------------
! Patrick's BXSF quadratic dispersion simulation generator
!----------------------------------------------------------

!--- 1. Note: this program assumes a cubic reciprocal unit cell (RUC) and quadratic dispersion, spherical Fermi surface

IMPLICIT NONE

! Numeric parameter declarations
double precision, parameter :: pi = 3.1415926535897932D0
double precision, parameter :: convau2ang = 0.529177209D0 ! Note: this is the conversion from atomic units (Bohr) to angstroms, because BXSF file reciprocal lattice vectors are given in (atomic units)^-1
double precision, parameter :: convfsarea2kt = 10.47576797D0 ! Note: this multiplicative constant is h_bar/(2*pi*elementarycharge) in units so that F ends up in kT (for fsarea given in angstroms^-2)
double precision, parameter :: convfsdade2mstar = 0.089135845D0 ! Note: this multiplicative constant is (h_bar^2)/(2*pi*(m_(e))) in units so that m* is dimensionless (for fsdAde given in (angstroms^-2)*Ryd^-1)
double precision, parameter :: piconvfsdade2mstar = pi*convfsdade2mstar ! Note: this multiplicative constant is (h_bar^2)/(2*(m_(e))) in the same units as convfsdade2mstar

! Variable declarations
integer,save::npoints
double precision,save::rlatx1
double precision,save::rlatx2
double precision,save::rlatx3
double precision,save::rlaty1
double precision,save::rlaty2
double precision,save::rlaty3
double precision,save::rlatz1
double precision,save::rlatz2
double precision,save::rlatz3
integer,save::linepos
integer,save::kxstep
integer,save::kystep
integer,save::kzstep
double precision,save::pctdone

double precision,save::realspacelatt
double precision,save::ksep
double precision,save::kfmax
double precision,save::freqmax
double precision,save::freqmin
double precision,save::freq
double precision,save::mstar
double precision,save::extarea
double precision,save::ef
double precision,dimension(6),save::e
double precision,save::kxdummy
double precision,save::kydummy
double precision,save::kzdummy
double precision,save::kx
double precision,save::ky
double precision,save::kz
double precision,save::eoffset
double precision,save::doscontrib

write(*,'('' '')')
write(*,'(''   ----------------- '')')
write(*,'(''    SimBXSF program '')')
write(*,'(''   ----------------- '')')
write(*,'('' '')')

!--- 2. Ask user for real space lattice constant (in Angstroms), realspacelatt
write(*,'('' What is the real space simple cubic lattice constant (in units of Angstroms)? '')')
read(*,*)realspacelatt

!--- 3. Ask user for # of k-points per RUC side (10-100), npoints
write(*,'('' How many k-points per RUC side? [10 - 100] '')')
read(*,*)npoints
if (npoints<10) then
 npoints = 10
end if
if (npoints>100) then
 npoints = 100
end if

!--- 4. Calculate spacing between points, ksep (includng the factor of 2*Pi)
ksep = 2.0D0*pi/(realspacelatt*dble(npoints - 1))

!--- 5. Calculate maximum safe kf to keep FS inside RUC, kfmax (includng the factor of 2*Pi)
kfmax = (pi/realspacelatt) - (2.0D0*ksep)

!--- 6. Calculate maximum safe freq, freqmax (and minimum safe freq to include at least 1 point, freqmin)
freqmax = convfsarea2kt*pi*kfmax*kfmax
freqmin = convfsarea2kt*2.0D0*ksep*ksep

!--- 7. Ask user for frequency (0-freqmax), freq
write(*,'('' What is the dHvA frequency? ['',f10.4,'' - '',f10.4 ,'' kT] '')')freqmin,freqmax
read(*,*)freq
if (freq<freqmin) then
 freq = freqmin
end if
if (freq>freqmax) then
 freq = freqmax
end if

!--- 8. Ask user for effective mass in units of m_e, mstar
write(*,'('' What is the effective mass (in units of m_e)? '')')
read(*,*)mstar


!--- 9. Calculate extremal area, extarea
extarea = freq/convfsarea2kt

!--- 9.b) Ask user for energy offset, in Ryd, eoffset
write(*,'('' What is energy offset (in units of Ryd)? '')')
read(*,*)eoffset

!--- 10. Calculate Fermi energy in Ryd, ef
ef = eoffset + (convfsdade2mstar*extarea/mstar)

!--create reciprocal lattice vectors (in 1/a.u., without the factor of 2*Pi) for writing to file
rlatx1=convau2ang/realspacelatt
rlatx2=0.0D0
rlatx3=0.0D0
rlaty1=0.0D0
rlaty2=convau2ang/realspacelatt
rlaty3=0.0D0
rlatz1=0.0D0
rlatz2=0.0D0
rlatz3=convau2ang/realspacelatt

!--- 11. Write file header stuff
open(unit=13,file='simulated.bxsf')
rewind 13
write(13,'('' BEGIN_INFO'')')
write(13,'(''   Fermi Energy: '',f12.6)')ef
write(13,'('' END_INFO'')')
write(13,'('' BEGIN_BLOCK_BANDGRID_3D'')')
write(13,'('' band_energies'')')
write(13,'('' BANDGRID_3D_BANDS'')')
write(13,'('' 1'')')
write(13,'(1x,I3,1x,I3,1x,I3)')npoints,npoints,npoints
write(13,'(''      0.00000000      0.00000000      0.00000000'')')
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')rlatx1,rlatx2,rlatx3
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')rlaty1,rlaty2,rlaty3
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')rlatz1,rlatz2,rlatz3
write(13,'('' BAND:    1'')')

linepos=0

write(*,'('' '')')
write(*,'('' Loop started'')')

write(*,'('' '')')
write(*,'('' Generating simulated.bxsf:   0.0 %''$)')

do kxstep=1,npoints
 do kystep=1,npoints
  do kzstep=1,npoints
   kx=((dble(kxstep)-1)/(dble(npoints)-1))*(2.0D0*pi/realspacelatt)
   ky=((dble(kystep)-1)/(dble(npoints)-1))*(2.0D0*pi/realspacelatt)
   kz=((dble(kzstep)-1)/(dble(npoints)-1))*(2.0D0*pi/realspacelatt)

   linepos=linepos+1

!--- 12. Calculate the energy of each point, in Ryd, and write to file
   kxdummy=kx-(pi/realspacelatt)
   kydummy=ky-(pi/realspacelatt)
   kzdummy=kz-(pi/realspacelatt)
   e(linepos)=eoffset+(piconvfsdade2mstar*((kxdummy*kxdummy)+(kydummy*kydummy)+(kzdummy*kzdummy))/mstar)

   if (linepos==6) then
    write(13,'(1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6)')e(1),e(2),e(3),e(4),e(5),e(6)
    linepos=0
   end if

  end do
 end do

 pctdone=(dble(kxstep)/dble(npoints))*100.0D0
 write(*,FMT='(A1,A1,A1,A1,A1,A1,A1$)')char(8),char(8),char(8),char(8),char(8),char(8),char(8)
 write(*,FMT='(F5.1,'' %''$)')pctdone

end do

if (linepos==5) then
 write(13,'(1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6)')e(1),e(2),e(3),e(4),e(5)
 linepos=0
end if

if (linepos==4) then
 write(13,'(1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6)')e(1),e(2),e(3),e(4)
 linepos=0
end if

if (linepos==3) then
 write(13,'(1x,ES13.6,1x,ES13.6,1x,ES13.6)')e(1),e(2),e(3)
 linepos=0
end if

if (linepos==2) then
 write(13,'(1x,ES13.6,1x,ES13.6)')e(1),e(2)
 linepos=0
end if

if (linepos==1) then
 write(13,'(1x,ES13.6)')e(1)
 linepos=0
end if

write(13,'('' END_BANDGRID_3D'')')
write(13,'('' END_BLOCK_BANDGRID_3D'')')

endfile 13
close(13,status='keep')

doscontrib=(2.0D0*pi*mstar/piconvfsdade2mstar)*sqrt(mstar*(ef-eoffset)/piconvfsdade2mstar)

write(*,'('' '')')
write(*,'('' '')')
write(*,'('' DOS contribution of this band = '',f10.6,'' Angstroms^-3 Ryd^-1'')')doscontrib
write(*,'('' '')')

stop
end
