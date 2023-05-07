program simbxsf

!----------------------------------------------------------
! Patrick's BXSF quadratic dispersion simulation generator
!----------------------------------------------------------

!--- 1. Note: this program assumes a cubic reciprocal unit cell (RUC) and quadratic dispersion, spheroidal Fermi surface

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
double precision,save::kshiftx
double precision,save::kshifty
double precision,save::kshiftz
double precision,save::ksep
double precision,save::kfmax
double precision,save::freqmax
double precision,save::freqmin
double precision,save::freqphi0
double precision,save::freqphi90
double precision,save::mstarphi0
double precision,save::mstarphi90
double precision,save::extareaphi0
double precision,save::extareaphi90
double precision,save::ef
double precision,dimension(6),save::e
double precision,save::kxdummy
double precision,save::kydummy
double precision,save::kzdummy
double precision,save::kx
double precision,save::ky
double precision,save::kz
double precision,save::k2nd
double precision,save::k2ndlast
logical,save::boundary
double precision,save::kx1
double precision,save::ky1
double precision,save::kz1
double precision,save::kx2
double precision,save::ky2
double precision,save::kz2
double precision,save::kxdummy1
double precision,save::kydummy1
double precision,save::kzdummy1
double precision,save::kxdummy2
double precision,save::kydummy2
double precision,save::kzdummy2
double precision,save::e1
double precision,save::e2

write(*,'('' '')')
write(*,'(''   ----------------- '')')
write(*,'(''    SimBXSF program '')')
write(*,'(''   ----------------- '')')
write(*,'('' '')')

!--- 2. Ask user for real space lattice constant (in Angstroms), realspacelatt
write(*,'('' What is the real space simple cubic lattice constant (in units of Angstroms)? '')')
read(*,*)realspacelatt

!--- 2b. Ask user for fractional amount of RUC side length to shift the FS centre in x direction, kshiftx
write(*,'('' How much to shift FS centre in x (fraction of RUC side length, e.g. 0.2)? '')')
read(*,*)kshiftx
if (kshiftx<0.0D0) then
 kshiftx = 0.0D0
end if
if (kshiftx>=1.0D0) then
 kshiftx = 0.0D0
end if

!--- 2c. Ask user for fractional amount of RUC side length to shift the FS centre in y direction, kshifty
write(*,'('' How much to shift FS centre in y (fraction of RUC side length, e.g. 0.1)? '')')
read(*,*)kshifty
if (kshifty<0.0D0) then
 kshifty = 0.0D0
end if
if (kshifty>=1.0D0) then
 kshifty = 0.0D0
end if

!--- 2d. Ask user for fractional amount of RUC side length to shift the FS centre in z direction, kshiftz
write(*,'('' How much to shift FS centre in z (fraction of RUC side length, e.g. 0.05)? '')')
read(*,*)kshiftz
if (kshiftz<0.0D0) then
 kshiftz = 0.0D0
end if
if (kshiftz>=1.0D0) then
 kshiftz = 0.0D0
end if

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

!--- 7a. Ask user for frequency for phi=0 (0-freqmax), freqphi0
write(*,'('' What is the dHvA frequency for phi=0? ['',f8.4,'' - '',f8.4 ,'' kT] '')')freqmin,freqmax
read(*,*)freqphi0
if (freqphi0<freqmin) then
 freqphi0 = freqmin
end if
if (freqphi0>freqmax) then
 freqphi0 = freqmax
end if

!--- 7b. Ask user for frequency for phi=90 (0-freqmax), freqphi90
write(*,'('' What is the dHvA frequency for phi=90? ['',f8.4,'' - '',f8.4 ,'' kT] '')')freqmin,freqmax
read(*,*)freqphi90
if (freqphi90<freqmin) then
 freqphi90 = freqmin
end if
if (freqphi90>freqmax) then
 freqphi90 = freqmax
end if

!--- 8. Ask user for effective mass in units of m_e for phi=0, mstarphi0
write(*,'('' What is the effective mass for phi=0 (in units of m_e)? '')')
read(*,*)mstarphi0


!--- 9a. Calculate extremal area for phi=0, extareaphi0
extareaphi0 = freqphi0/convfsarea2kt

!--- 9b. Calculate extremal area for phi=90, extareaphi90
extareaphi90 = freqphi90/convfsarea2kt

!--- 9c. Calculate and display the effective mass in units of m_e for phi=90, mstarphi90
mstarphi90 = mstarphi0*(extareaphi90*extareaphi90)/(extareaphi0*extareaphi0)
write(*,'('' Therefore effective mass for phi=90 (units of m_e) is '',f8.4)')mstarphi90


!--- 10. Calculate Fermi energy in Ryd, ef
ef = convfsdade2mstar*extareaphi0/mstarphi0

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

k2nd=((dble(2)-1.0D0)/(dble(npoints)-1.0D0))*(2.0D0*pi/realspacelatt)
k2ndlast=((dble(npoints-1)-1.0D0)/(dble(npoints)-1.0D0))*(2.0D0*pi/realspacelatt)

do kxstep=1,npoints
 do kystep=1,npoints
  do kzstep=1,npoints

   kx=((dble(kxstep)-1.0D0)/(dble(npoints)-1.0D0))*(2.0D0*pi/realspacelatt)
   ky=((dble(kystep)-1.0D0)/(dble(npoints)-1.0D0))*(2.0D0*pi/realspacelatt)
   kz=((dble(kzstep)-1.0D0)/(dble(npoints)-1.0D0))*(2.0D0*pi/realspacelatt)

   linepos=linepos+1

   boundary=.FALSE.

!--- 12. Calculate the energy of each point, in Ryd, and write to file
   if ((kxstep==1).OR.(kxstep==npoints)) then
    boundary=.TRUE.
    kx1=k2nd
    kx2=k2ndlast
   else
    kx1=kx
    kx2=kx
   end if
   if ((kystep==1).OR.(kystep==npoints)) then
    boundary=.TRUE.
    ky1=k2nd
    ky2=k2ndlast
   else
    ky1=ky
    ky2=ky
   end if
   if ((kzstep==1).OR.(kzstep==npoints)) then
    boundary=.TRUE.
    kz1=k2nd
    kz2=k2ndlast
   else
    kz1=kz
    kz2=kz
   end if
   if (boundary==.TRUE.) then
    kxdummy1=kx1-(pi/realspacelatt)-(kshiftx*2.0D0*pi/realspacelatt)
    kydummy1=ky1-(pi/realspacelatt)-(kshifty*2.0D0*pi/realspacelatt)
    kzdummy1=kz1-(pi/realspacelatt)-(kshiftz*2.0D0*pi/realspacelatt)
    kxdummy2=kx2-(pi/realspacelatt)-(kshiftx*2.0D0*pi/realspacelatt)
    kydummy2=ky2-(pi/realspacelatt)-(kshifty*2.0D0*pi/realspacelatt)
    kzdummy2=kz2-(pi/realspacelatt)-(kshiftz*2.0D0*pi/realspacelatt)
    e1=piconvfsdade2mstar*((((kxdummy1*kxdummy1)+(kydummy1*kydummy1))/mstarphi0)+((kzdummy1*kzdummy1)/mstarphi90))
    e2=piconvfsdade2mstar*((((kxdummy2*kxdummy2)+(kydummy2*kydummy2))/mstarphi0)+((kzdummy2*kzdummy2)/mstarphi90))
    e(linepos)=(e1+e2)/2.0D0
   else
    kxdummy=kx-(pi/realspacelatt)-(kshiftx*2.0D0*pi/realspacelatt)
    kydummy=ky-(pi/realspacelatt)-(kshifty*2.0D0*pi/realspacelatt)
    kzdummy=kz-(pi/realspacelatt)-(kshiftz*2.0D0*pi/realspacelatt)
    e(linepos)=piconvfsdade2mstar*((((kxdummy*kxdummy)+(kydummy*kydummy))/mstarphi0)+((kzdummy*kzdummy)/mstarphi90))
   end if

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

write(*,'('' '')')
write(*,'('' '')')

stop
end
