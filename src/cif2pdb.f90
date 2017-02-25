PROGRAM CIFTOPDB
 IMPLICIT NONE
 INTEGER            :: IERR,i,j,k,l
 INTEGER            :: n_atoms = 0
 INTEGER            :: n_T = 0
 REAL               :: cell_0(6),rv(3,3),vr(3,3),q,r
 REAL               :: crystal_r(3),cart_r(3)
 CHARACTER (LEN=80) :: line
 CHARACTER (LEN=5 ) :: label(2)
 CHARACTER (LEN=10) :: space = " P1      "
 CHARACTER (LEN=30) :: charac
 integer            :: n_atoms_type,model
 integer            :: id(100,10000) = 0
 real               :: coordinate_cart(100,3,10000) = 0.0
 integer            :: histogram_atom_type(100)
 character (len=5)  :: atoms_type_max_label(100)
 atoms_type_max_label = 'Xx'
 n_atoms_type         = 0
 !open(101,file='input',IOSTAT=IERR,STATUS='OLD',ACTION='READ')
 !if(ierr/=0)then
 model=1
 !else
 ! read(101,*) model
 !end if
 !close(101)
!
 i=0
 OPEN(111,FILE='input.cif',STATUS='OLD',ACTION='READ',IOSTAT=IERR )
 WRITE(6,'(A6,I4)')'MODEL ',model
 j=0
 fileopen: IF( IERR == 0) THEN
  read_: DO
    READ (111,'(A)',IOSTAT=IERR) line
    IF (line(1:13)=='_cell_length_') THEN
      i=i+1
      READ(line,*) charac,cell_0(i)
    ENDIF
    IF (line(1:12)=='_cell_angle_') THEN
      j=j+1
      READ(line,*) charac,cell_0(3+j)
    ENDIF
    IF (line(1:15)=='_symmetry_space') THEN
      READ(line,'(35x,a10)') space
    ENDIF
    IF (line(1:)=='_atom_site_fract_z') EXIT read_
  END DO read_
 ENDIF fileopen
 CALL cell(rv,vr,cell_0)
 i=0
 q=0.0
 r=0.0
 read_atoms: DO
  i=i+1
  read (111,'(A)',iostat=IERR) line
  IF (IERR/=0) exit read_atoms
  READ(line,*) label(2),crystal_r(1),crystal_r(2),crystal_r(3)
  if(label(1)(1:2)=="sh") cycle read_atoms
  label(1)=label(2)
  n_atoms = n_atoms + 1
  do j=1,3
   crystal_r(j)=MOD(crystal_r(j)+100.0,1.0)
   cart_r(j) = rv(j,1)*crystal_r(1) + rv(j,2)*crystal_r(2) + rv(j,3)*crystal_r(3)
  end do
  if(n_atoms_type==0)then
   n_atoms_type=1
   atoms_type_max_label(1)=label(2)
   histogram_atom_type(1)=1
   id(1,1)=i
   do j=1,3
    coordinate_cart(1,j,i)=cart_r(j)
   end do
   cycle read_atoms
  end if
  do j=1,n_atoms_type
   if(label(2)==atoms_type_max_label(j))then
    histogram_atom_type(j)=histogram_atom_type(j)+1
    do k=1,3
     coordinate_cart(j,k,i)=cart_r(k)
    end do
    id(j,histogram_atom_type(j))=i
    cycle read_atoms
   end if
  end do
  n_atoms_type=n_atoms_type+1
  atoms_type_max_label(n_atoms_type)=label(2)
  histogram_atom_type(n_atoms_type)=1
  do k=1,3
   coordinate_cart(n_atoms_type,k,i)=cart_r(k)
  end do
  id(n_atoms_type,1)=i
 END DO read_atoms 
 write(6,'(a6,1x,a11,1x,i5)') 'REMARK','Atoms type:', n_atoms_type
 write(6,'(a6,1x,999(a5,1x))')'REMARK',(atoms_type_max_label(j),j=1,n_atoms_type)
 write(6,'(a6,1x,999(i5,1x))')'REMARK',(histogram_atom_type(j),j=1,n_atoms_type)
 WRITE(6,'(A6,3f9.3,3f7.2)') &
  'CRYST1',cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6)
 i=0
 do k=1,n_atoms_type
   do j=1,histogram_atom_type(k)
    i=i+1
    label(1)=atoms_type_max_label(k)
    label(2)=atoms_type_max_label(k)
    WRITE(6,'(a6,i5,1x,2(a4,1x),i4,4x,3f8.3,2f6.2,10x,a4)') &
    'ATOM  ',i,label(1),'MOL ',0,(coordinate_cart(k,l,id(k,j)),l=1,3),0.0,0.0,label(1)
   end do
 end do
 WRITE(6,'(A6)')'ENDMDL'
END PROGRAM
!
subroutine cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: pi,DEGTORAD
 pi = ACOS(-1.0)
 DEGTORAD=pi/180.0
 IF(cell_0(4) == 90.0) THEN
   cosa = 0.0
 ELSE
   ALP=cell_0(4)*degtorad
   COSA=cos(ALP)
 ENDIF
 IF(cell_0(5) == 90.0) THEN
   cosb = 0.0
 ELSE
   bet = cell_0(5)*degtorad
   cosb = cos(bet)
 ENDIF
 IF(cell_0(6) == 90.0) then
   sing = 1.0
   cosg = 0.0
 ELSE
   gam = cell_0(6)*degtorad
   sing = sin(gam)
   cosg = cos(gam)
 ENDIF
 rv(1,1) = cell_0(1)
 rv(1,2) = cell_0(2)*cosg
 rv(1,3) = cell_0(3)*cosb
 rv(2,1) = 0.0
 rv(2,2) = cell_0(2)*sing
 rv(2,3) = cell_0(3)*(cosa - cosb*cosg)/sing
 rv(3,1) = 0.0
 rv(3,2) = 0.0
 rv(3,3) = sqrt( cell_0(3)*cell_0(3) - rv(1,3)*rv(1,3) - rv(2,3)*rv(2,3)) 
 call inverse(rv,vr,3)
 RETURN
end subroutine cell
!
subroutine inverse(a,c,n)
!============================================================
! Inverse matrix
! Method: Based on Doolittle LU factorization for Ax=b
! Alex G. December 2009
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - dimension
! output ...
! c(n,n) - inverse matrix of A
! comments ...
! the original matrix a(n,n) will be destroyed 
! during the calculation
!===========================================================
 implicit none 
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
!
! step 0: initialization for matrices L and U and b
! Fortran 90/95 aloows such operations on matrices
 L=0.0
 U=0.0
 b=0.0
!
! step 1: forward elimination
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
!
! Step 2: prepare L and U matrices 
! L matrix is a matrix of the elimination coefficient
! + the diagonal elements are 1.0
 do i=1,n
  L(i,i) = 1.0
 end do
! U matrix is the upper triangular part of A
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
!
! Step 3: compute columns of the inverse matrix C
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
! Step 3a: Solve Ld=b using the forward substitution
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
! Step 3b: Solve Ux=d using the back substitution
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
! Step 3c: fill the solutions x(n) into column k of C
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 return
end subroutine inverse
