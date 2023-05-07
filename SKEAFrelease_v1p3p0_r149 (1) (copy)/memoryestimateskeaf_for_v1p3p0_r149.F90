program memoryestimateskeaf

IMPLICIT NONE

integer, parameter :: int8byte = SELECTED_INT_KIND(18) ! for defining 8-byte integers

integer(KIND=int8byte) :: maxnuminput ! Maximum n that can be read in from the input BXSF file
integer(KIND=int8byte) :: maxnumint ! Maximum number of points per single cell side (a full supercell side will have ~4x this many)
integer(KIND=int8byte) :: numthreads ! Number of OpenMP threads for parallel execution
double precision :: freemem ! total free memory (available RAM + swap) in the system
double precision :: skeafmem ! total memory required for SKEAF variables and arrays

! start the program

write(*,'('' '')')
write(*,'('' What is the value of the SKEAF maxnuminput parameter? [default is 100] '')')
write(*,'('' NOTE: maxnuminput is not the same as maxnumint! Just leave maxnuminput '', &
&''at 100 unless you will be reading BXSF files with more than 100 points per side'')')
read(*,*)maxnuminput

write(*,'('' '')')
write(*,'('' How many OpenMP parallel threads are you using? [type 1 if not using OpenMP] '')')
read(*,*)numthreads
if (numthreads<1) then
 numthreads=1
end if

write(*,'('' '')')
write(*,'('' How much total free memory do you have? [available RAM + swap, in gigabytes] '')')
read(*,*)freemem

skeafmem = 0
maxnumint = 0

 do while (skeafmem < freemem)
  maxnumint = maxnumint + 1
  skeafmem = skeafmemory(maxnumint,numthreads)
 end do
 maxnumint = maxnumint - 1
 skeafmem = skeafmemory(maxnumint,numthreads)
 write(*,'('' '')')
 write(*,'('' Free gigabytes of memory in system = '',F11.8)')freemem
 write(*,'('' Estimated gigabytes required by SKEAF = '',F11.8)')skeafmem
 write(*,'('' Value of maxnumint = '',I4)')maxnumint
 write(*,'('' '')')

stop

contains

 double precision function skeafmemory(tmpnumint,tmpnumthreads)
  integer(KIND=int8byte), intent(in)::tmpnumint,tmpnumthreads

  integer(KIND=int8byte) :: scarraydimension
  integer(KIND=int8byte) :: slicedimension
  integer(KIND=int8byte) :: measpardim
  integer(KIND=int8byte) :: measpardim2
  integer(KIND=int8byte) :: Nfloatglobal ! number of global floating point (8 byte double) variables / array elements
  integer(KIND=int8byte) :: Ninteglobal ! number of global integer (8 byte) variables / array elements
  integer(KIND=int8byte) :: Ncharaglobal ! number of global character (1 byte) variables / array elements
  integer(KIND=int8byte) :: Nlogicglobal ! number of global logical (4 byte) variables / array elements
  integer(KIND=int8byte) :: Nfloatthread ! number of floating point (8 byte double) variables / array elements in parallel threads
  integer(KIND=int8byte) :: Nintethread ! number of integer (8 byte) variables / array elements in parallel threads
  integer(KIND=int8byte) :: Ncharathread ! number of character (1 byte) variables / array elements in parallel threads
  integer(KIND=int8byte) :: Nlogicthread ! number of logical (4 byte) variables / array elements in parallel threads
  integer(KIND=int8byte) :: Nfloat ! total number of floating point (8 byte double) variables / array elements
  integer(KIND=int8byte) :: Ninte ! total number of integer (8 byte) variables / array elements
  integer(KIND=int8byte) :: Nchara ! total number of character (1 byte) variables / array elements
  integer(KIND=int8byte) :: Nlogic ! total number of logical (4 byte) variables / array elements

  scarraydimension = (4*tmpnumint)+2
  slicedimension = (4*tmpnumint)*(4*tmpnumint)
  measpardim = (30*tmpnumint)+300
  measpardim2 = 4*measpardim

  Nfloatglobal = (2*maxnuminput*maxnuminput*maxnuminput)+(12*scarraydimension*measpardim)+&
&(13*scarraydimension*measpardim2)+(84*measpardim2)+(2*slicedimension)+178
  Ninteglobal = (3*scarraydimension)+(1*scarraydimension*measpardim)+&
&(2*scarraydimension*measpardim2)+(21*measpardim2)+100
  Ncharaglobal = (1*scarraydimension*scarraydimension*scarraydimension)+(1*scarraydimension)+1412
  Nlogicglobal = (1*scarraydimension*measpardim)+(3*scarraydimension*measpardim2)+71

  Nfloatthread = tmpnumthreads*((3*scarraydimension*scarraydimension)+(9*slicedimension)+77)
  Nintethread = tmpnumthreads*(66)
  Ncharathread = tmpnumthreads*((1*scarraydimension)+1)
  Nlogicthread = tmpnumthreads*((1*scarraydimension*scarraydimension)+6)

  Nfloat = Nfloatglobal + Nfloatthread
  Ninte = Ninteglobal + Nintethread
  Nchara = Ncharaglobal + Ncharathread
  Nlogic = Nlogicglobal + Nlogicthread

  skeafmemory = (dble((Nfloat*8)+(Ninte*8)+(Nchara*1)+(Nlogic*4))/dble(1024*1024*1024))

 end function skeafmemory

end
