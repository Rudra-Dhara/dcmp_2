PROGRAM bxsfconverter

IMPLICIT NONE

! Hard-coded parameters
integer*8, parameter :: maxnuminput = 100 ! Maximum number of points per side in the input file RUC; don't change this unless your input BXSF file has nx, ny or nz > 99
integer*8, parameter :: maxinputarraysize = maxnuminput*maxnuminput*maxnuminput

! Numeric parameter declarations
double precision, parameter :: pi = 3.1415926535897932D0

! Variable declarations
character(len=50),save::filenamein
character(len=50),save::filenameout
character(len=100),save::tempstringread
double precision,save::fermienergyin
double precision,save::fermienergyout
double precision,save::numlistedbands
integer*8,save::nxin
integer*8,save::nyin
integer*8,save::nzin
integer*8,save::nxout
integer*8,save::nyout
integer*8,save::nzout
integer*8,save::numkpoints
integer*8,save::numread
logical,save::keepreading
double precision,save::patricklreciplato1
double precision,save::patricklreciplato2
double precision,save::patricklreciplato3
double precision,save::patricklreciplatx1in
double precision,save::patricklreciplatx2in
double precision,save::patricklreciplatx3in
double precision,save::patricklreciplaty1in
double precision,save::patricklreciplaty2in
double precision,save::patricklreciplaty3in
double precision,save::patricklreciplatz1in
double precision,save::patricklreciplatz2in
double precision,save::patricklreciplatz3in
double precision,save::patricklreciplatx1out
double precision,save::patricklreciplatx2out
double precision,save::patricklreciplatx3out
double precision,save::patricklreciplaty1out
double precision,save::patricklreciplaty2out
double precision,save::patricklreciplaty3out
double precision,save::patricklreciplatz1out
double precision,save::patricklreciplatz2out
double precision,save::patricklreciplatz3out
double precision,dimension(maxinputarraysize),save::kreadarray
double precision,dimension(maxnuminput,maxnuminput,maxnuminput),save::masterkarrayin
double precision,dimension(maxnuminput,maxnuminput,maxnuminput),save::masterkarrayout
integer*8::i
integer*8::j
integer*8::k
integer*8::iin
integer*8::jin
integer*8::kin
double precision,dimension(6),save::e
integer,save::linepos
double precision,save::pctdone
logical,save::detperiodicgrid
character(len=1),save::periodicgrid
character(len=1),save::periodicsure
character(len=1),save::twopifactor
character(len=1),save::hartrees
double precision,save::hartmult
character(len=1),save::corminexp
character(len=1),save::knowminexp
double precision,save::maxenergy
double precision,save::minexp
character(len=13),save::maxenergychar
double precision,save::currentmaxexp
double precision,save::divminexp
character(len=1),save::signswitch

! Shows title and warnings
write(*,'('' '')')
write(*,'('' This program is intended to help correct the problems in BXSF files generated'')')
write(*,'('' by the ELK and exciting DFT packages, so that those BXSF files can give'')')
write(*,'('' accurate numerical results when analyzed by the SKEAF program.'')')
write(*,'('' '')')
write(*,'('' ***WARNING: This program can fix problems, but it can also cause them if you'')')
write(*,'(''             pick the wrong parameters! Such BXSF issues can lead to SKEAF'')')
write(*,'(''             results that are incorrect, sometimes in subtle ways!***'')')
write(*,'('' '')')
write(*,'('' Please see the README-forSKEAF.txt file for important instructions and notes!'')')

! This section reads the input BXSF FS file
write(*,'('' '')')
write(*,'('' What is the filename of the input BXSF file? [50 chars. max] '')')
read(*,*)filenamein

write(*,'('' '')')
write(*,'('' Reading BXSF file: '',A50)')filenamein

  open(15,file=filenamein,status='old')

  keepreading=.TRUE.
  do while (keepreading==.TRUE.)
   read(15,fmt='(A100)')tempstringread
   tempstringread=TRIM(tempstringread)
   if ((INDEX(tempstringread,'Fermi Energy')/=0).OR.(INDEX(tempstringread,'fermi energy')/=0)&
&.OR.(INDEX(tempstringread,'FERMI ENERGY')/=0).OR.(INDEX(tempstringread,'Fermi energy')/=0)) then
    tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
    read(tempstringread,*)fermienergyin
    keepreading=.FALSE.
   end if
  end do
  keepreading=.TRUE.

  do i=1,4
   read(15,*)tempstringread
  end do

  read(15,fmt='(A100)')tempstringread
  tempstringread=TRIM(tempstringread)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread,*)numlistedbands

  if (numlistedbands/=1.0D0) then
   write(*,'('' '')')

   write(*,'('' ***BXSF ERROR: Header indicates that file contains multiple bands!***'')')
   write(*,'('' (Need to separate each band into its own BXSF file, usually by'')')
   write(*,'(''   opening in XCrysDen and re-saving individual bands from there.)'')')
   write(*,'('' '')')
   STOP
  end if

  read(15,fmt='(A100)')tempstringread
  tempstringread=TRIM(tempstringread)
  tempstringread=tempstringread(SCAN(tempstringread,'0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)nxin
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)nyin
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'0123456789'):)
  read(tempstringread,*)nzin

  if (((nxin+1)>maxnuminput).OR.((nyin+1)>maxnuminput).OR.((nzin+1)>maxnuminput)) then
   write(*,'('' '')')
   write(*,'('' ***BXSF ERROR: too many k-points for current array bounds!***'')')
   write(*,'('' (Please increase the value of the maxnuminput parameter in the'')')
   write(*,'(''  source code for this program and try compiling and running again.)'')')
   write(*,'('' '')')
   STOP
  end if

  read(15,fmt='(A100)')tempstringread
  tempstringread=TRIM(tempstringread)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplato1
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplato2
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread,*)patricklreciplato3

  if ((patricklreciplato1/=0.0D0).OR.(patricklreciplato2/=0.0D0).OR.(patricklreciplato3/=0.0D0)) then
   write(*,'('' '')')
   write(*,'('' ***BXSF ERROR: reciprocal lattice vector origin is not at (0, 0, 0)!***'')')
   write(*,'('' '')')
   STOP
  end if

  read(15,fmt='(A100)')tempstringread
  tempstringread=TRIM(tempstringread)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplatx1in
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplatx2in
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread,*)patricklreciplatx3in

  read(15,fmt='(A100)')tempstringread
  tempstringread=TRIM(tempstringread)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplaty1in
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplaty2in
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread,*)patricklreciplaty3in
 
  read(15,fmt='(A100)')tempstringread
  tempstringread=TRIM(tempstringread)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplatz1in
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread(:SCAN(tempstringread,' ')-1),*)patricklreciplatz2in
  tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
  tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
  read(tempstringread,*)patricklreciplatz3in

  read(15,*)tempstringread

  numkpoints=nxin*nyin*nzin
  numread=0
  keepreading=.TRUE.
  do while (keepreading==.TRUE.)
   read(15,fmt='(A100)')tempstringread
   tempstringread=TRIM(tempstringread)
   if ((INDEX(tempstringread,'BAND:')/=0).OR.(INDEX(tempstringread,'Band:')/=0)&
&.OR.(INDEX(tempstringread,'band:')/=0).OR.(INDEX(tempstringread,'prod:')/=0)&
&.OR.(INDEX(tempstringread,'Prod:')/=0).OR.(INDEX(tempstringread,'PROD:')/=0)) then
    write(*,'('' '')')
    write(*,'('' ***BXSF ERROR: File contains multiple bands!***'')')
    write(*,'('' (Need to separate each band into its own BXSF file, usually by'')')
    write(*,'(''   opening in XCrysDen and re-saving individual bands from there.)'')')
    write(*,'('' '')')
    STOP
   end if
   if ((INDEX(tempstringread,'END')==0).AND.(INDEX(tempstringread,'End')==0)&
&.AND.(INDEX(tempstringread,'end')==0)) then
    tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
    do while (SCAN(tempstringread,'-0123456789')/=0)
     numread = numread + 1
     read(tempstringread(:SCAN(tempstringread,' ')-1),*)kreadarray(numread)
     tempstringread=tempstringread(SCAN(tempstringread,' ')+1:)
     tempstringread=tempstringread(SCAN(tempstringread,'-0123456789'):)
    end do
   else
    keepreading=.FALSE.
   end if
  end do

  do i=1,nxin
   do j=1,nyin
    do k=1,nzin
     masterkarrayin(i,j,k)=kreadarray(k+(j-1)*nzin+(i-1)*nzin*nyin)
    end do
   end do
  end do
  close(15,status='keep')

  if (masterkarrayin(nxin,nyin,nzin)/=masterkarrayin(1,1,1)) then
   detperiodicgrid=.TRUE.
  else
   detperiodicgrid=.FALSE.
  end if

if (numread/=numkpoints) then
 write(*,'('' '')')
 write(*,'('' ***BXSF ERROR: Number of energies read from file differs from number in file header!***'')')
 write(*,'('' '')')
 STOP
end if

!***************************************** end of read-in section, start of parameter asking section

periodicgrid='u'
periodicsure='u'
twopifactor='u'
hartrees='u'
corminexp='u'
knowminexp='u'
minexp=0.1
signswitch='u'

do while (((periodicgrid/='y').AND.(periodicgrid/='n')).OR.(periodicsure/='y'))
 write(*,'('' '')')
 write(*,'('' Are energies defined on a Periodic Grid (rather than General Grid)? [y or n] '')')
 write(*,'('' (Usually y for BXSFs made with ELK v1.4.22 or exciting helium-3.)'')')
 read(*,*)periodicgrid
 if ((periodicgrid=='y').AND.(detperiodicgrid==.FALSE.)) then
  write(*,'('' '')')
  write(*,'('' ARE YOU SURE? [y or n] '')')
  write(*,'('' (Rough check indicates that it might not be a Periodic Grid.)'')')
  read(*,*)periodicsure
 else if ((periodicgrid=='n').AND.(detperiodicgrid==.TRUE.)) then
  write(*,'('' '')')
  write(*,'('' ARE YOU SURE? [y or n] '')')
  write(*,'('' (Rough check indicates that it is probably a Periodic Grid.)'')')
  read(*,*)periodicsure
 else
  periodicsure='y'
 end if
end do

do while ((twopifactor/='y').AND.(twopifactor/='n'))
 write(*,'('' '')')
 write(*,'('' Do reciprocal lattice vectors contain the factor of 2*Pi'',&
&'' (e.g. 2*Pi/lattconst rather than 1/lattconst)? [y or n] '')')
 write(*,'('' (Usually y for BXSFs made with ELK v1.4.22 or exciting helium-3.)'')')
 read(*,*)twopifactor
end do

do while ((hartrees/='y').AND.(hartrees/='n'))
 write(*,'('' '')')
 write(*,'('' Are energies specified in Hartree units (rather than Rydberg units)? [y or n] '')')
 write(*,'('' (Usually y for BXSFs made with ELK v1.4.22 or exciting helium-3.)'')')
 read(*,*)hartrees
end do

do while ((corminexp/='y').AND.(corminexp/='n'))
 write(*,'('' '')')
 write(*,'('' Has the original minimum exponent been divided out of the energies? [y or n] '')')
 write(*,'('' (Usually n for BXSFs made with ELK v1.4.22; usually y for BXSFs made with exciting helium-3.)'')')
 read(*,*)corminexp
end do

if (corminexp=='y') then
 do while ((knowminexp/='y').AND.(knowminexp/='n'))
  write(*,'('' '')')
  write(*,'('' Do you know the original minimum exponent that was divided out of the energies? [y or n] '')')
  read(*,*)knowminexp
 end do
 if (knowminexp=='y') then
  do while (minexp/=dnint(minexp))
   write(*,'('' '')')
   write(*,'('' What is the original minimum exponent (minexp in exciting fermisurf.f90)? [integer, usually negative] '')')
   read(*,*)minexp
   if (minexp/=dnint(minexp)) then
    write(*,'('' The exponent must be an integer or zero.'')')
   end if
  end do
 else
  write(*,'('' '')')
  write(*,'('' Attempting to guess minexp, by assuming that the original maximum exponent was 0.'')')
  write(*,'('' (This is a very poor assumption... but it seems about right for copper...)'')')
  maxenergy=masterkarrayin(1,1,1)
  do i=1,nxin
   do j=1,nyin
    do k=1,nzin
     if (masterkarrayin(i,j,k)>maxenergy) then
      maxenergy=masterkarrayin(i,j,k)
     end if
    end do
   end do
  end do
  write(maxenergychar,'(ES13.6)')maxenergy
  read(maxenergychar(SCAN(maxenergychar,'eE')+1:),*)currentmaxexp
  minexp=-1.0D0*(currentmaxexp+0.0D0)
  write(*,'('' '')')
  write(*,'('' Guess: minexp = '',I4)')idnint(minexp)
  minexp=dnint(minexp)
 end if
end if

do while ((signswitch/='y').AND.(signswitch/='n'))
 write(*,'('' '')')
 write(*,'('' Has the sign (+ve/-ve) of the energies been switched? [y or n] '')')
 write(*,'('' (Usually n for BXSFs made with ELK v1.4.22; usually y for BXSFs made with exciting helium-3.)'')')
 read(*,*)signswitch
end do

filenameout=filenamein
do while (filenameout==filenamein)
 write(*,'('' '')')
 write(*,'('' What is the filename of the output BXSF file? [50 chars. max; different from input filename] '')')
 read(*,*)filenameout
end do

!***************************************** end of parameter asking section, start of adjustment section

if (twopifactor=='y') then
 patricklreciplatx1out=patricklreciplatx1in/(2.0D0*pi)
 patricklreciplatx2out=patricklreciplatx2in/(2.0D0*pi)
 patricklreciplatx3out=patricklreciplatx3in/(2.0D0*pi)
 patricklreciplaty1out=patricklreciplaty1in/(2.0D0*pi)
 patricklreciplaty2out=patricklreciplaty2in/(2.0D0*pi)
 patricklreciplaty3out=patricklreciplaty3in/(2.0D0*pi)
 patricklreciplatz1out=patricklreciplatz1in/(2.0D0*pi)
 patricklreciplatz2out=patricklreciplatz2in/(2.0D0*pi)
 patricklreciplatz3out=patricklreciplatz3in/(2.0D0*pi)
else
 patricklreciplatx1out=patricklreciplatx1in
 patricklreciplatx2out=patricklreciplatx2in
 patricklreciplatx3out=patricklreciplatx3in
 patricklreciplaty1out=patricklreciplaty1in
 patricklreciplaty2out=patricklreciplaty2in
 patricklreciplaty3out=patricklreciplaty3in
 patricklreciplatz1out=patricklreciplatz1in
 patricklreciplatz2out=patricklreciplatz2in
 patricklreciplatz3out=patricklreciplatz3in
end if

if (hartrees=='y') then
 hartmult=2.0D0
else
 hartmult=1.0D0
end if

if (corminexp=='y') then
 divminexp=(10.0D0**(-1.0*minexp))
else
 divminexp=1.0D0
end if

if (signswitch=='y') then
 hartmult=(-1.0D0)*hartmult
end if

fermienergyout=hartmult*fermienergyin/divminexp

if (periodicgrid=='y') then
 nxout=nxin+1
 nyout=nyin+1
 nzout=nzin+1
 do i=1,nxout
  do j=1,nyout
   do k=1,nzout
    if (i==nxout) then
     iin=1
    else
     iin=i
    end if
    if (j==nyout) then
     jin=1
    else
     jin=j
    end if
    if (k==nzout) then
     kin=1
    else
     kin=k
    end if
    masterkarrayout(i,j,k)=hartmult*masterkarrayin(iin,jin,kin)/divminexp
   end do
  end do
 end do
else
 nxout=nxin
 nyout=nyin
 nzout=nzin
 do i=1,nxout
  do j=1,nyout
   do k=1,nzout
    masterkarrayout(i,j,k)=hartmult*masterkarrayin(i,j,k)/divminexp
   end do
  end do
 end do
end if

!***************************************** end of adjustment section, start of output write section

open(unit=13,file=filenameout)
rewind 13
write(13,'('' BEGIN_INFO'')')
write(13,'(''   Fermi Energy: '',f12.6)')fermienergyout
write(13,'('' END_INFO'')')
write(13,'('' BEGIN_BLOCK_BANDGRID_3D'')')
write(13,'('' band_energies'')')
write(13,'('' BANDGRID_3D_BANDS'')')
write(13,'('' 1'')')
write(13,'(1x,I3,1x,I3,1x,I3)')nxout,nyout,nzout
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')patricklreciplato1,patricklreciplato2,patricklreciplato3
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')patricklreciplatx1out,patricklreciplatx2out,patricklreciplatx3out
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')patricklreciplaty1out,patricklreciplaty2out,patricklreciplaty3out
write(13,'(1x,f15.8,1x,f15.8,1x,f15.8)')patricklreciplatz1out,patricklreciplatz2out,patricklreciplatz3out
write(13,'('' BAND:    1'')')

linepos=0

write(*,'('' '')')
write(*,'('' Generating '',A50,'':   0.0 %''$)')filenameout

do i=1,nxout
 do j=1,nyout
  do k=1,nzout
   linepos=linepos+1
   e(linepos)=masterkarrayout(i,j,k)
   if (linepos==6) then
    write(13,'(1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6,1x,ES13.6)')e(1),e(2),e(3),e(4),e(5),e(6)
    linepos=0
   end if
  end do
 end do

 pctdone=(dble(i)/dble(nxout))*100.0D0
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
write(*,'('' REMINDER: If you are using exciting helium-3, the effective masses and DOS are probably still wrong.'')')
write(*,'(''            (See README-forSKEAF.txt file for details.)'')')
write(*,'('' '')')

STOP

END PROGRAM bxsfconverter
