program zif_cif2gin
 use iso_fortran_env
 implicit none
 integer             :: i,j,k,ierr,nn
 integer             :: num_args
 integer             :: n_atoms = 0
 real                :: cell_0(1:6) = 0.0
 real                :: sumcharge = 0.0
 real                :: rv(3,3),vr(3,3)
 real                :: pressure    = 1.0e5
 real                :: temperature = 298.0
 real                :: mass_sorbate = 131.293 ! Xe,  Ar:40.0
 real,parameter      :: k_B = 8.617332478e-5
 real                :: beta = 0.0
 real                :: chemical_potential=0.0
 character(len=100)  :: line
 character(len=20)   :: spam,simulation_type="single"
 character(len=80)   :: string_stop_head= "_atom_site_charge" !"_atom_site_charge"
 character(len=50)   :: forcefield = "WHCJ"
 character(len=100)  :: CIFFilename=" "
 character(len=100)  :: filename=" "
 ! allocatable
 real,allocatable             :: xcrystal(:,:),xcartes(:,:),charge(:)
 real,allocatable             :: DistanceMatrix(:,:)
 character(len=100), dimension(:), allocatable :: args
 character(len=4),allocatable :: label(:),newlabel(:)
 num_args = command_argument_count()
 allocate(args(num_args))
 do i = 1, num_args
  call get_command_argument(i,args(i))
 end do
 write(6,'(a,1x,i2)')'Arguments:',command_argument_count()
 write(6,'(a)')( args(i),i=1,num_args)
 if(num_args==0) then
  call print_help()
  stop
 end if
 do i=1,num_args
  select case(args(i))
   case ('-h','--help')
    call print_help()
   case ('-c','--cif')
    CIFFilename=args(i+1)
    filename=CIFFilename(1:Clen_trim(CIFFilename)-4)
    write(6,'(a)') filename
   case ('-T','--temperature')
    read(args(i+1),*) temperature
   case ('-GCMC','--adsorption')
    read(args(i+1),*) pressure
    simulation_type="montecarlo"
   case ('-fff','--forcefield')
    read(args(i+1),*) forcefield
  end select
 end do
 beta = 1.0/(k_B*temperature)
 chemical_potential=chem(pressure,temperature,mass_sorbate)
 write(6,*)'Temperature:',temperature
 write(6,*)'Chemical potential:', chemical_potential
 open(100,file=CIFfilename,status='old',iostat=ierr)
 if(ierr/=0) stop 'Error opening CIF file'
 read_cif: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_cif
  if(line(1:14)=="_cell_length_a")then
   read(line,*)spam,cell_0(1)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_b")then
   read(line,*)spam,cell_0(2)
   cycle read_cif
  end if
  if(line(1:14)=="_cell_length_c")then
   read(line,*)spam,cell_0(3)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_alpha")then
   read(line,*)spam,cell_0(4)
   cycle read_cif
  end if
  if(line(1:16)=="_cell_angle_beta")then
   read(line,*)spam,cell_0(5)
   cycle read_cif
  end if
  if(line(1:17)=="_cell_angle_gamma")then
   read(line,*)spam,cell_0(6)
   cycle read_cif
  end if
  if(line(1:)==string_stop_head) exit read_cif
 end do read_cif
 call cell(rv,vr,cell_0)
 n_atoms=0
 read_natoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit read_natoms
  n_atoms=n_atoms+1
 end do read_natoms
!
 allocate(xcrystal(3,n_atoms),xcartes(3,n_atoms))
 allocate(charge(n_atoms))
 allocate(DistanceMatrix(n_atoms,n_atoms))
 allocate(label(n_atoms),newlabel(n_atoms))
!
 rewind(100)
 write(6,*)'Atoms:',n_atoms
 do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  if(line(1:)==string_stop_head) exit
 end do
 i=0
 read_atoms: do
  read(100,'(a)',iostat=ierr) line
  if(ierr/=0) exit
  i=i+1
  read(line,*)label(i),( xcrystal(j,i),j=1,3),charge(i)
  do j=1,3
   xcrystal(j,i)=MOD(xcrystal(j,i)+100.0,1.0)
   xcartes(j,i) =rv(j,1)*xcrystal(1,i) + rv(j,2)*xcrystal(2,i)+rv(j,3)*xcrystal(3,i)
  end do
  write(6,'(a4,1x,3(f14.12,1x),3(f14.7,1x))',ADVANCE='yes') &
   label(i),( xcrystal(j,i),j=1,3),( xcartes(j,i),j=1,3)
 end do read_atoms
 call make_dist_matrix(n_atoms,cell_0,rv,vr,xcrystal,DistanceMatrix)
 write(6,*)'=========='
 close(100)
! Analysis:
 do i=1,n_atoms
  if(label(i)=="Zn  ")then
   charge(i)=0.6918
   newlabel(i)="Zn  "
  else if(label(i)=="N   ".or.label(i)=="N1  ")then
   charge(i)=-0.3879
   newlabel(i)="N   "
  else if(label(i)=="C   ".or.label(i)=="C1  ".or.label(i)=="C2  ")then
   !we analise the neigbourhod of the C-atom
   nn=0
   do j=1,n_atoms
    if(i==j) cycle
    if(label(j)=="N   ".or.label(j)=="N1  ")then 
     !write(6,*)label(i),label(j),DistanceMatrix(i,j)
     if(DistanceMatrix(i,j)<1.5.and.DistanceMatrix(i,j)>=1.2) nn=nn+1
    end if
   end do
   if (nn==1) then
    charge(i)=-0.0839
    newlabel(i)="C2  "
   else if (nn==2) then
    charge(i)=0.260
    !0.4291-0.4526+3*0.1308-0.1128
    newlabel(i)="C4  "
   else
    write(6,*)i,label(i),nn
    STOP 'Error with C assignation'
   end if
  else if(label(i)=="H   ".or.label(i)=="H1  ".or.label(i)=="H2  ")then
   charge(i)=0.1128
   newlabel(i)="H2  "
  else
   newlabel(i)=label(i)
   charge(i)=0.0
   write(6,'(a)') "Adsorbate found!:",newlabel(i)
  end if
  !write(6,'(a4,1x,3(f14.12,1x))')newlabel(i),( xcrystal(j,i),j=1,3)
 end do
 write(6,*)'Total charge:', SUM(charge)
 do i=1,n_atoms
  charge(i)=charge(i)-SUM(charge)/n_atoms
 end do
 write(6,*)'Total charge (recalibration):', SUM(charge)
 call write_gin(cell_0,xcrystal,n_atoms,newlabel,charge,filename,simulation_type,pressure,beta,forcefield)
!
 deallocate(xcrystal)
 deallocate(xcartes)
 deallocate(DistanceMatrix)
 deallocate(charge)
 deallocate(label)
 deallocate(newlabel)
 !stop 'Done'
 contains
 PURE INTEGER FUNCTION Clen(s)      ! returns same result as LEN unless:
 CHARACTER(*),INTENT(IN) :: s       ! last non-blank char is null
 INTEGER :: i
 Clen = LEN(s)
 i = LEN_TRIM(s)
 IF (s(i:i) == CHAR(0)) Clen = i-1  ! len of C string
 END FUNCTION Clen
 PURE INTEGER FUNCTION Clen_trim(s) ! returns same result as LEN_TRIM unless:
 CHARACTER(*),INTENT(IN) :: s       ! last char non-blank is null, if true:
 INTEGER :: i                       ! then len of C string is returned, note:
                                    ! Ctrim is only user of this function
 i = LEN_TRIM(s) ; Clen_trim = i
 IF (s(i:i) == CHAR(0)) Clen_trim = Clen(s)   ! len of C string
 END FUNCTION Clen_trim
!
 real(8) function chemh(pressure,temperature,mass)
  implicit none
  real,intent(in) :: pressure,temperature,mass
  real(8)            :: fugacity_coefficient = 1.0
  real(8)            :: fugacity
  real(8)            :: lambda,n
  real(8),parameter  :: h_planck = 6.6260695729e-34 ! J s
  real(8),parameter  :: PI       = acos(-1.0)
  real(8),parameter  :: k_B_JK   = 1.380648813e-23  ! J / K
  real(8),parameter  :: uma2kg   = 1.660538921e-27  ! kg 
  fugacity  = pressure*fugacity_coefficient
  n         = k_B_JK*temperature/fugacity
  lambda    = (2*PI*mass*uma2kg*k_B_JK*temperature/(h_planck*h_planck))**(-3.0/2.0)
  chemh     = -log(n/lambda)/beta
  return
 end function chemh
 real function chem(pressure,temperature,mass) ! eV
  !ideal gas
  implicit none
  real,intent(in) :: pressure,temperature,mass
  real,parameter  :: k_B_JK   = 1.380648813e-23 ! J/K
  real,parameter  :: k_b      = 8.617332478e-5  ! eV/K
  real(16)        :: n_Q,n_ideal
  n_Q     = 1e30*mass**(3/2.0)*1.0**(3/2.0)
  n_ideal = pressure/(k_b_jk*temperature)
  chem    = k_b*temperature*log(n_ideal/n_Q)
  return
 end function chem
 SUBROUTINE cell(rv,vr,cell_0)
 implicit none
 integer :: i,j
 real, intent(in)  :: cell_0(6)
 real, intent(out) :: rv(3,3),vr(3,3)
 real, parameter   :: pi = ACOS(-1.0)
 real :: alp,bet
 real :: cosa,cosb,cosg
 real :: gam,sing
 real :: DEGTORAD
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
 WRITE(6,'(a)') 'Cell:'
 WRITE(6,'(6F14.7)')( cell_0(j), j=1,6 )
 WRITE(6,'(a)')'Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 WRITE(6,'(a)')'Inverse Linear Transformation Operator:'
 DO i=1,3
  WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
  implicit none
  real,intent(out)   :: cell_0(6)
  real,intent(in)    :: rv(3,3)
  integer            :: i,j
  real               :: temp(6)
  REAL               :: radtodeg
  REAL, PARAMETER    :: pi=ACOS(-1.0) 
  radtodeg=180.0/PI
  do i = 1,3
    temp(i) = 0.0
    do j = 1,3
      temp(i) = temp(i) + rv(j,i)**2
    enddo
    temp(i) = sqrt(temp(i))
  enddo
  cell_0(1) = abs(temp(1))
  cell_0(2) = abs(temp(2))
  cell_0(3) = abs(temp(3))
  do i = 1,3
    temp(3+i) = 0.0
  enddo
  do j = 1,3
    temp(4) = temp(4) + rv(j,2)*rv(j,3)
    temp(5) = temp(5) + rv(j,1)*rv(j,3)
    temp(6) = temp(6) + rv(j,1)*rv(j,2)
  enddo
  temp(4) = temp(4)/(temp(2)*temp(3))
  temp(5) = temp(5)/(temp(1)*temp(3))
  temp(6) = temp(6)/(temp(1)*temp(2))
  cell_0(4) = radtodeg*acos(temp(4))
  cell_0(5) = radtodeg*acos(temp(5))
  cell_0(6) = radtodeg*acos(temp(6))
  DO i=4,6
     if (abs(cell_0(i) - 90.0 ).lt.0.00001) cell_0(i) = 90.0
     if (abs(cell_0(i) - 120.0).lt.0.00001) cell_0(i) = 120.0
  ENDDO
  RETURN
 END SUBROUTINE uncell
!
 SUBROUTINE inverse(a,c,n)
 implicit none
 integer n
 real a(n,n), c(n,n)
 real L(n,n), U(n,n), b(n), d(n), x(n)
 real coeff
 integer i, j, k
 L=0.0
 U=0.0
 b=0.0
 do k=1, n-1
   do i=k+1,n
      coeff=a(i,k)/a(k,k)
      L(i,k) = coeff
      do j=k+1,n
         a(i,j) = a(i,j)-coeff*a(k,j)
      end do
   end do
 end do
 do i=1,n
  L(i,i) = 1.0
 end do
 do j=1,n
  do i=1,j
    U(i,j) = a(i,j)
  end do
 end do
 do k=1,n
  b(k)=1.0
  d(1) = b(1)
  do i=2,n
    d(i)=b(i)
    do j=1,i-1
      d(i) = d(i) - L(i,j)*d(j)
    end do
  end do
  x(n)=d(n)/U(n,n)
  do i = n-1,1,-1
    x(i) = d(i)
    do j=n,i+1,-1
      x(i)=x(i)-U(i,j)*x(j)
    end do
    x(i) = x(i)/u(i,i)
  end do
  do i=1,n
    c(i,k) = x(i)
  end do
  b(k)=0.0
 end do
 RETURN
 END SUBROUTINE inverse
!
subroutine make_dist_matrix(n,cell_0,rv,vr,x,dist_matrix)
 implicit none
 integer,intent(in) :: n
 real,intent(in)    :: cell_0(6),rv(3,3),vr(3,3),x(3,n)
 real,intent(out)   :: dist_matrix(n,n)
 integer            :: i,j,k
 real               :: r1(3),r2(3),s
 DO i=1,n
    dist_matrix(i,i)=0.0
    DO j=i+1,n
       forall ( k=1:3 )
        r1(k)=x(k,i)
        r2(k)=x(k,j)
       end forall
       call make_distances(cell_0,r1,r2,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
 return
end subroutine make_dist_matrix
!
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
  k=0
  do l=-1,1
   do m=-1,1
      do n=-1,1
         k = k + 1
         ouratom(1) = r1(1)
         ouratom(2) = r1(2)
         ouratom(3) = r1(3)
         atom(1) = r2(1) + l
         atom(2) = r2(2) + m
         atom(3) = r2(3) + n
         d_image(k) = distance(atom,ouratom,rv)
         forall ( i=1:3)
           image(i,k) = atom(i)
         end forall
     enddo
   enddo
  enddo
  dist=MINVAL(d_image)
  RETURN
 END SUBROUTINE
!
 REAL FUNCTION DISTANCE(atom,ouratom,rv)
  IMPLICIT NONE
  INTEGER :: j
  REAL :: atom(3),ouratom(3),per_atom(3),dist(3),o_atom(3),o_ouratom(3)
  REAL :: rv(3,3)
  FORALL ( j=1:3 )
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*atom(1) + rv(j,2)*atom(2) + rv(j,3)*atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
  END FORALL
  DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
 END FUNCTION
!
 SUBROUTINE write_gin(cell_0,xcryst,n_atoms,labels,charge,filename,simulation_type,pressure,beta,fff)
  IMPLICIT NONE
  INTEGER,          intent(in) :: n_atoms
  character (len=4),intent(in) :: labels(n_atoms)
  REAL,             intent(in) :: xcryst(3,n_atoms),cell_0(1:6),charge(n_atoms)
  integer                      :: i,j
  character(len=100)           :: filename
  character(len=104)           :: GINFilename
  character(len=7)             :: CIF="_Xe.cif"
  character(len=4)             :: GIN=".gin"
  character(len=4)             :: RES=".res"
  character(len=20),intent(in) :: simulation_type
  character(len=50),intent(in) :: fff
  real             ,intent(in) :: pressure,beta
  GINFilename=CIFFilename(1:Clen_trim(CIFFilename)-4)//GIN
  OPEN(999,file=GINFilename)
  select case(simulation_type)
   case('montecarlo')
    write(999,'(a)')'montecarlo conv molq qok'
    write(999,'(a,1x,f20.10,1x,f20.10,1x,f20.10)')&
          '# (pressure,temperature,mass)',pressure,temperature,mass_sorbate
    write(999,'(a,1x,f20.10,1x,f20.10,1x,f20.10)')'# (\mu,\muQ,\beta)',&
          chem(pressure,temperature,mass_sorbate),chemh(pressure,temperature,mass_sorbate),beta
   case default
    write(999,'(A)')'single conv molq qok'
  end select
  WRITE(999,'(A)')'cell'
  WRITE(999,'(6(f9.5,1x))') (cell_0(j) , j=1,6)
  WRITE(999,'(A)')'fractional'
  DO i=1,n_atoms
   WRITE(999,'(a4,1x,3(f14.7,1x),1x,f14.7)') labels(i),xcryst(1,i),xcryst(2,i),xcryst(3,i),charge(i)
  END DO
  WRITE(999,'(A)')'species 10'
  select case(fff)
  case('HZJ')
   WRITE(999,'(A)')'Zn     core  1.0' 
   WRITE(999,'(A)')'N      core -0.5'
   WRITE(999,'(A)')'C2     core -0.1'
   WRITE(999,'(A)')'C4     core  0.4'
   WRITE(999,'(A)')'H2     core  0.1'
  case default
  WRITE(999,'(A)')'Zn     core'
  WRITE(999,'(A)')'N      core'
  WRITE(999,'(A)')'C2     core'
  WRITE(999,'(A)')'C4     core'
  WRITE(999,'(A)')'H2     core'
  end select
  write(999,'(a)')'C5     core 0.0'
  write(999,'(a)')'He     core 0.0'
  write(999,'(a)')'Ar     core 0.0'
  write(999,'(a)')'Kr     core 0.0'
  write(999,'(a)')'Xe     core 0.0'
  select case(simulation_type)
   case('montecarlo')
   write(999,*)'# montecarlo'
   write(999,*)'temperature',temperature,'K'
   write(999,*)'mcchemicalpotential', chemh(pressure,temperature,mass_sorbate)
   write(999,*)'mcmaxdisplacement    0.5'
   write(999,*)'mccreate             0.333333'
   write(999,*)'mcdestroy            0.333333'
   write(999,*)'mcmove               0.333333'
   write(999,*)'mctrial           2000'
   write(999,*)'gcmcmolecule      1'
   write(999,*)'Xe core           0.0 0.0 0.0'
   write(999,*)'#'
   case default
   write(999,*)'# single'
  end select
! WHCJ:
  select case(fff)
   case default ! WHCJ
   write(999,'(a)')'# Xuanjun Wu, Jin Huang, Weiquan Cai, and Mietek Jaroniecb'
   write(999,'(a)')'# RSC Adv., 2014,4, 16503-16511'
   write(999,'(a)')'# DOI: 10.1039/C4RA00664J'
   write(999,'(a)')'lenn epsilon geometric 12 6'
   write(999,'(a)')'C2 core C2 core 0.0 12.0'
   write(999,'(a)')'C4 core C4 core 0.0 12.0'
   write(999,'(a)')'H2 core H2 core 0.0 12.0'
   write(999,'(a)')'N  core N  core 0.0 12.0'
   write(999,'(a)')'Zn core Zn core 0.0 12.0'
   write(999,'(a)')'C5 core C5 core 0.0 12.0'
   write(999,'(a)')'He core He core 0.0 12.0'
   write(999,'(a)')'Ar core Ar core 0.0 12.0'
   write(999,'(a)')'Kr core Kr core 0.0 12.0'
   write(999,'(a)')'Xe core Xe core 0.0 12.0'
   write(999,'(a)')'epsilon/sigma kcal'
   write(999,'(a)')'Zn core 0.0787   2.462'
   write(999,'(a)')'N  core 0.0438 3.261'
   write(999,'(a)')'C2 core 0.0667 3.431'
   write(999,'(a)')'C4 core 0.0667 3.431'
   write(999,'(a)')'H2 core 0.0279 2.571'
   write(999,'(a)')'C5 core 0.294104 3.73'
   write(999,'(a)')'Ar core 0.238066 3.34'
   write(999,'(a)')'Xe core 0.439169 4.1'
   write(999,'(a)')'He core 0.02166 2.64'
   write(999,'(a)')'Kr core 0.330669 3.636'
   WRITE(999,'(A)')'harmonic intra bond'
   WRITE(999,'(A)')'C2 core N core  42.323376 1.339 '
   WRITE(999,'(A)')'C4 core N core  42.323376 1.339 '
   WRITE(999,'(A)')'C2 core C2 core 44.925224 1.346 '
   WRITE(999,'(A)')'H2 core C2 core 31.82926  0.929 '
   WRITE(999,'(A)')'H2 core C4 core 31.82926  0.929 '
   WRITE(999,'(A)')'Zn core N  core 7.458628  1.987 '
   WRITE(999,'(A)')'three kcal intra bond'
   WRITE(999,'(A)')'C4 core N  core N  core  35.0     112.17 '
   WRITE(999,'(A)')'C4 core N  core H2 core  35.0     123.89 '
   WRITE(999,'(A)')'C2 core C2 core N  core  35.0     108.67'
   WRITE(999,'(A)')'C2 core C2 core H2 core  25.0     125.67'
   WRITE(999,'(A)')'C2 core N  core H2 core  25.0     125.66'
   WRITE(999,'(A)')'N  core C4 core C2 core  35.0     105.24'
   WRITE(999,'(A)')'N  core C4 core Zn core  25.0     127.50'
   WRITE(999,'(A)')'N  core C2 core Zn core  17.5     128.00'
   WRITE(999,'(A)')'Zn core N  core N  core   5.25    109.47'
   write(999,'(a)')'torsion kcal bond'
   write(999,'(a)')'C2 N  C4 N  4.800 +2.0 180.0 '
   write(999,'(a)')'C2 N  C4 H2 4.150 +2.0 180.0 '
   write(999,'(a)')'C4 N  C2 C2 4.800 +2.0 180.0 '
   write(999,'(a)')'N  C2 C2 N  4.0   +2.0 180.0 '
   write(999,'(a)')'N  C2 C2 H2 4.0   +2.0 180.0 '
   write(999,'(a)')'H2 C2 C2 H2 4.0   +2.0 180.0 '
   write(999,'(a)')'Zn N  C4 N  0.100 +2.0 180.0 '
   write(999,'(a)')'Zn N  C4 H2 0.100 +2.0 180.0 '
   write(999,'(a)')'Zn N  C2 C2 0.100 +2.0 180.0 '
   write(999,'(a)')'N  Zn N  C2 0.174 +3   0.0   '
   write(999,'(a)')'N  Zn N  C4 0.174 +3   0.0   '
   write(999,'(a)')'torsion kcal intra'
   write(999,'(a)')'N  H2 C4 N  1.100 +2   180.0 2.5 2.5 2.5 2.5'
   write(999,'(a)')'C2 H2 C2 N  1.000 +2   180.0 2.5 2.5 2.5 2.5'
  case ('HZJ')
   write(999,'(a)')'# J. Chem. Phys. 136, 244703, 2012, Z. Hu, L. Zhang, and J. Jiang'
   write(999,'(a)')'# doi: 10.1063/1.4729314'
   write(999,'(a)')'lenn epsilon geometric 12 6'
   write(999,'(a)')'C2 core C2 core 0.0 14.0'
   write(999,'(a)')'C4 core C4 core 0.0 14.0'
   write(999,'(a)')'H2 core H2 core 0.0 14.0'
   write(999,'(a)')'N  core N  core 0.0 14.0'
   write(999,'(a)')'Zn core Zn core 0.0 14.0'
   write(999,'(a)')'C5 core C5 core 0.0 14.0'
   write(999,'(a)')'He core He core 0.0 14.0'
   write(999,'(a)')'Ar core Ar core 0.0 14.0'
   write(999,'(a)')'Kr core Kr core 0.0 14.0'
   write(999,'(a)')'Xe core Xe core 0.0 14.0'
   write(999,'(a)')'epsilon/sigma'
   write(999,'(a)')'Zn core  0.0005420554  1.960              # Zn'                
   write(999,'(a)')'N  core  0.0073721607  3.250              # N '           
   write(999,'(a)')'C2 core  0.0037290924  3.400              # C4'             
   write(999,'(a)')'C4 core  0.0037290924  3.400              # C1'            
   write(999,'(a)')'H2 core  0.0006508811  2.421              # H2'            
   write(999,'(a)')'C5 core  0.0000000104  1.0                # C5'            
   write(999,'(a)')'Ar core  0.0000000104  1.0                # He'
   write(999,'(a)')'Xe core  0.0000000104  1.0                # Ar'            
   write(999,'(a)')'He core  0.0000000104  1.0                # Kr'            
   write(999,'(a)')'Kr core  0.0190443133  4.100              # Xe'
!  bonds: 
   WRITE(999,'(A)')'harmonic intra bond'
   WRITE(999,'(A,1x,f14.7,1x,f14.7)')'C2 core N  core',kjmolnmnm2evaa(3.4309*100000),1.371
   WRITE(999,'(A,1x,f14.7,1x,f14.7)')'C4 core N  core',kjmolnmnm2evaa(4.0836*100000),1.340
   WRITE(999,'(A,1x,f14.7,1x,f14.7)')'C2 core C2 core',kjmolnmnm2evaa(4.3346*100000),1.346
   WRITE(999,'(A,1x,f14.7,1x,f14.7)')'H2 core C2 core',kjmolnmnm2evaa(3.0710*100000),0.9290
   WRITE(999,'(A,1x,f14.7,1x,f14.7)')'H2 core C4 core',kjmolnmnm2evaa(2.8451*100000),0.9600
   WRITE(999,'(A,1x,f14.7,1x,f14.7)')'Zn core N  core',kjmolnmnm2evaa(0.7000*100000),1.987
   WRITE(999,'(A)')'three intra bond kjmol'
   WRITE(999,'(A)')'C4 core N  core N  core  5.8576 112.17'
   WRITE(999,'(A)')'C4 core N  core H2 core  5.8576 123.89 '
   WRITE(999,'(A)')'C2 core C2 core N  core  5.8576 108.67'
   WRITE(999,'(A)')'C2 core C2 core H2 core  4.1840 125.67'
   WRITE(999,'(A)')'C2 core N  core H2 core  4.1840 125.66'
   WRITE(999,'(A)')'N  core C4 core C2 core  5.8576 105.24'
   WRITE(999,'(A)')'N  core C4 core Zn core  1.7000 128.35'
   WRITE(999,'(A)')'N  core C2 core Zn core  1.8000 126.40'
   WRITE(999,'(A)')'Zn core N  core N  core  1.0000 109.47'
   write(999,'(a)')'torsion bond kjmol'
   write(999,'(a)')'C2 N  C4 N  20.0832 +2.0 180.0'
   write(999,'(a)')'C2 N  C4 H2 20.0832 +2.0 180.0'
   write(999,'(a)')'C4 N  C2 C2 20.0832 +2.0 180.0'
   write(999,'(a)')'N  C2 C2 N  16.7360 +2.0 180.0'
   write(999,'(a)')'N  C2 C2 H2 16.7360 +2.0 180.0'
   write(999,'(a)')'H2 C2 C2 H2 16.7360 +2.0 180.0'
   write(999,'(a)')'Zn N  C4 N   0.5000 +2.0 180.0'
   write(999,'(a)')'Zn N  C4 H2  0.5000 +2.0 180.0'
   write(999,'(a)')'Zn N  C2 C2  0.5000 +2.0 180.0'
   write(999,'(a)')'Zn N  C2 H2  0.5000 +2.0 180.0'
   write(999,'(a)')'torsion intra kjmol'
   write(999,'(a)')'N  Zn H2 C2  0.4000 +2.0 180.0 2.5 3.5 1.5 3.5'
   write(999,'(a)')'H2 N  N  C4  4.6024 +2.0 180.0 3.5 3.5 2.0 2.0'
   write(999,'(a)')'C2 N  C2 H2  4.6024 +2.0 180.0 2.0 2.5 1.5 3.5'
  end select
  WRITE(999,'(A)')'nobond Zn C2'
  WRITE(999,'(A)')'nobond Zn C4'
  WRITE(999,'(A)')'nobond Zn H2'
  WRITE(999,'(A)')'nobond N  H2'
  WRITE(999,'(A)')'nobond C2 C4'
  write(999,'(a)')'nobond Xe Zn'
  write(999,'(a)')'nobond Xe C2' 
  write(999,'(a)')'nobond Xe C4'
  write(999,'(a)')'nobond Xe H2'
  write(999,'(a)')'nobond Xe N '
  write(999,'(a)')'nobond Xe Xe'
  GINFilename=CIFFilename(1:Clen_trim(CIFFilename)-4)//RES
  WRITE(999,'(a,1x,a)')'dump every 1 ',GINFilename
  WRITE(999,'(a,1x,a)')'output lammps ',CIFFilename(1:Clen_trim(CIFFilename)-4)
  if(simulation_type=="montecarlo")then
   GINFilename=CIFFilename(1:Clen_trim(CIFFilename)-4)//CIF
   write(999,'(a,1x,a)')'output cif',GINFilename
  end if
  CLOSE(999)
  RETURN
 END SUBROUTINE write_gin
 real function kjmol2ev(x)
  implicit none
  real,intent(in) :: x
  kjmol2ev=x*(1036.427/100000.0)
  return
 end  function kjmol2ev
 real function kjmolnmnm2evAA(x)
  real, intent(in) :: x
  kjmolnmnm2evaa = kjmol2ev(x)/100.0
  return
 end function  kjmolnmnm2evaa
 subroutine print_help()
    print '(a)', '  -h, --help  print usage information and exit'
    print '(a)', '  -c, --cif   CIF File input'
 end subroutine print_help
end program zif_cif2gin

