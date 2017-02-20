PROGRAM farthrest_nodes_subtitutions
 IMPLICIT NONE
 INTEGER            :: n_atoms = 0
 INTEGER            :: n_Si = 0
 INTEGER            :: n_Al = 0
 INTEGER            :: walker_position,walker_drunked
 REAL               :: r1  = 1.2
 REAL               :: r2  = 2.0
 REAL,    PARAMETER :: infinite = 9999999.999999
 REAL,    PARAMETER :: solera   = 1e-10
 REAL               :: MC_cycles 
 INTEGER            :: IERR,SEED
 INTEGER            :: spacegroup = 1
 INTEGER            :: i,j,k,l,m
 REAL               :: atom(3),ouratom(3),pot_dist
 REAL               :: cell_0(6),cost,p,q,r,s,rv(3,3),vr(3,3)
 REAL,              ALLOCATABLE :: xcryst(:,:),xinbox(:,:),dist_matrix(:,:)
 CHARACTER (LEN=2), ALLOCATABLE :: label(:)
 INTEGER,           ALLOCATABLE :: adj(:,:),pivots(:),id(:)
 CHARACTER (LEN=1), ALLOCATABLE :: adj_char(:,:)
 CHARACTER (LEN=6)  :: potential,atomc
 CHARACTER (LEN=80) :: line,string
 CHARACTER (LEN=4)  :: mol
 LOGICAL            :: FLAG      = .true.
! {{
 OPEN(100,FILE='INPUT',STATUS='OLD',ACTION='READ',IOSTAT=IERR )
 fileopen: IF( IERR == 0) THEN
  DO
   READ (100,'(A)',IOSTAT=IERR) line
   IF( IERR /= 0 ) EXIT
   IF (line(1:5)=='MODEL') THEN
    READ (line,*) atomc,cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6),spacegroup
   ELSE IF(line(1:5)=='PROGR') THEN
      READ(line,*) atomc, n_Al, MC_cycles, potential
   ELSE IF(line(1:5)=='ATOMS') THEN
      n_atoms = 0
   ELSE
      n_atoms = n_atoms +1
   ENDIF
  ENDDO
 ENDIF fileopen
 REWIND(100)
! {{ primer output }}
 CALL get_seed(SEED)
 WRITE(6,'(a,2x,i10)')'Random Seed:',SEED
 WRITE(6,'(a,2x,i4)')'Number of atoms:',n_atoms
 WRITE(6,'(a,2x,i3)')'Space Group:',spacegroup
 CALL cell(rv,vr,cell_0)
! {{ alicatamos variables }}
 ALLOCATE(xcryst(0:3,1:n_atoms)       ,STAT=IERR)
 ALLOCATE(xinbox(1:3,1:n_atoms)       ,STAT=IERR)
 ALLOCATE(label(1:n_atoms)            ,STAT=IERR)
 ALLOCATE(pivots(1:n_Al)              ,STAT=IERR)
 ALLOCATE(id(1:n_atoms)               ,STAT=IERR)
 ALLOCATE(dist_matrix(n_atoms,n_atoms),STAT=IERR)
 IF(IERR/=0) STOP '[ERROR] variables "xcryst xinbox label" sin alicatar en memoria.'
! {{ leemos las coordenadas de los atomos }}
 DO i=1,3
    READ (100,'(A)') line
 ENDDO
 read_coor: DO i=1,n_atoms
    READ (100,'(A)',IOSTAT=IERR) line
    IF( IERR /= 0 ) EXIT read_coor
    READ(line,*) label(i),xcryst(1,i),xcryst(2,i),xcryst(3,i)
    FORALL ( j=1:3 ) 
      xinbox(j,i) = rv(j,1)*xcryst(1,i) + rv(j,2)*xcryst(2,i) + rv(j,3)*xcryst(3,i)
    END FORALL
 ENDDO read_coor
 CLOSE(100)
! {{ creamos el grafo de atomos }}
 make_graph: IF( FLAG.EQV..true. ) THEN
    ALLOCATE(adj(1:n_atoms,1:n_atoms)      ,STAT=IERR)
    IF(IERR/=0) STOP '[ERROR] variable adj sin alicatar en memoria.'
    ALLOCATE(adj_char(1:n_atoms,1:n_atoms) ,STAT=IERR)
    IF(IERR/=0) STOP '[ERROR] variable adj_char sin alicatar en memoria.'
    pairs: DO i=1,n_atoms
      DO j=i,n_atoms
      adj(i,j)=0.0
      adj_CHAR(i,j)=' '
      IF(j/=i)THEN
       FORALL ( k=1:3 )
          atom(k)= xcryst(k,j)
          ouratom(k) = xcryst(k,i)
       END FORALL
       call make_distances(cell_0,ouratom,atom,rv,r)
       IF(r>=r1.AND.r<r2)THEN
         adj(i,j)=1.0
         adj_CHAR(I,J)='@'
       ENDIF
      ENDIF
      adj(j,i)=adj(i,j)
      adj_CHAR(J,I)=adj_CHAR(I,J)
     ENDDO
    ENDDO pairs
    conectivity: DO i=1,n_atoms ! conectividad de los ATOM  ! k=4 y k=2    en zeolitas.
     k=0
     DO j=1,n_atoms
        k=k+adj(i,j)
     ENDDO
     xcryst(0,i)=k
    ENDDO conectivity
    WRITE(6,'(A)')'[loops] Connectivity between nodes [Degree]:'
    WRITE(6,'(1000(I1))')(int(xcryst(0,k)),k=1,n_atoms)
 END IF make_graph
! make dist_matrix
 DO i=1,n_atoms
    dist_matrix(i,i)=0.0
    DO j=i+1,n_atoms
       forall ( k=1:3 )
        ouratom(k)=xcryst(k,i)
        atom(k)   =xcryst(k,j)
       end forall
       call make_distances(cell_0,ouratom,atom,rv,s)
       dist_matrix(i,j)=s
       dist_matrix(j,i)=dist_matrix(i,j)
    END DO
 END DO
! {{ sustituimos los aluminios
!    {{ el primer Al es aleatorio
 pivots(1)=INT(r4_uniform(1.0,n_atoms+1.0,SEED))
 DO WHILE ( label(pivots(1))/='Si' ) 
  pivots(1)=INT(r4_uniform(1.0,n_atoms+1.0,SEED))
 ENDDO
 label(pivots(1))='Al'
 q=3.40
 WRITE(6,'(a)')'Repulsive substitution with cost-function:'
 WRITE(6,'(a,f5.2,a)')'cost = 1/(s - p0),  p0 =',q,' 10^-10 m, if s >= p0'
 WRITE(6,'(a,f14.3,a)')'cost = ',infinite,' if s <  p0'
 WRITE(6,'(a)')'     Crystalographic positions                  Cost     Label'
 WRITE(6,'(3(f14.6,1x),f14.5)')xcryst(1,pivots(1)),xcryst(2,pivots(1)),xcryst(3,pivots(1)),0.0
! }}
! {{ colocamos los siguienets Al
 walker_position = pivots(1)
 walker_drunked = 0
 choose_Al: DO k=2,n_Al
  m = SEED
  r = infinite
  jump: DO walker_position=1,n_atoms
     IF (label(walker_position)=='Si') THEN
       pot_dist=0.0
       p=energy(n_atoms,n_Al,xcryst,label,dist_matrix,cell_0,rv,potential)
       DO j=1,k-1
          walker_drunked = pivots(j)
          !FORALL ( i=1:3 )
          !  ouratom(i)=xcryst(i,walker_position)
          !  atom(i)   =xcryst(i,pivots(j))
          !END FORALL
          !call make_distances(cell_0,ouratom,atom,rv,s)
          s = dist_matrix(walker_position,walker_drunked)
          q = repulsive_potential_Lowenstein(s,3.40,infinite,potential)
          pot_dist = pot_dist + q/real(k)
       END DO
       pot_dist = pot_dist + p
       IF ( pot_dist < r ) THEN
          r = pot_dist
          m = walker_position
       END IF
     END IF
  END DO jump
  IF (m==SEED) THEN
     m=INT(r4_uniform(1.0,n_atoms+1.0,SEED))
     DO WHILE ( label(m)/='Si' )
       m=INT(r4_uniform(1.0,n_atoms+1.0,SEED))
     ENDDO
  END IF
  pivots(k) = m
  label(pivots(k))='Al'
  WRITE(6,'(3(f14.6,1x),f14.5,1x,i3,1x,a2)')xcryst(1,m),xcryst(2,m),xcryst(3,m),pot_dist,m,label(m)
 END DO choose_Al                       ! }}i
 q=infinite
 k=0
 WRITE(6,'(a)')'=============================================='
 WRITE(6,'(a)')'MonteCarlo: Metropolis with identity changes'
 WRITE(6,'(a,f10.5)')'Si > Al, Al > Si. STOP when occurrence < ',solera
 WRITE(6,'(a)')'=============================================='
 WRITE(6,'(a,5x,a,5x,a)')'Cost_function','Occurrence %','Scan-Cost'
 !DO WHILE ( q >= solera .and. k<=MC_cycles )
 DO WHILE ( k<=MC_cycles )
    ! Metropolis:
    ! identity changes: Si - Al
    CALL MonteCarlo(potential,n_atoms,dist_matrix,n_Al,SEED,xcryst,label,cell_0,rv,q,cost)
    k=k+1
 END DO
 CALL write_gin(cell_0,xcryst,n_atoms,label)
 WRITE(6,'(a)')'=============================================='
 WRITE(6,'(a,i3,a)')'geometrical properties of ',n_Al,' Al:'
 q=0.0
 r=0.0
 m=0
 do i=1,n_atoms
   if(label(i)=='Al')then
    do j=1,n_atoms
     if(i/=j)then
      if(label(j)=='Al')then
       !forall ( k=1:3 )
       !    atom(k)    = xcryst(k,i)
       !    ouratom(k) = xcryst(k,j)
       !end forall 
       !call make_distances(cell_0,atom,ouratom,rv,s)
       s = dist_matrix(i,j)
       m=m+1
       r=r+s
       q=q+s*s
      end if
     endif
    end do
   end if
 end do
 WRITE(6,*)'distancia media y dispersion'
 WRITE(6,*)r/real(m),sqrt(q/real(m)-(r*r)/real(m*m))
 pot_dist=r/real(m)
 q=0.0
 r=0.0
 m=0
 do i=1,n_atoms
   if(label(i)=='Al')then
    p=infinite
    do j=1,n_atoms
     if(i/=j)then
      if(label(j)=='Al')then
       !forall ( k=1:3 )
       !    atom(k)    = xcryst(k,i)
       !    ouratom(k) = xcryst(k,j)
       !end forall
       !call make_distances(cell_0,atom,ouratom,rv,s)
       s=dist_matrix(i,j)
       if( s <= p ) p=s
      end if
     endif
    end do
    m=m+1
    r=r+p
    q=q+p*p
   endif
 end do
 WRITE(6,*)'distancia minima media y dispersion'
 WRITE(6,*)r/real(m),sqrt(q/real(m)-(r*r)/real(m*m))
 DEALLOCATE(adj,adj_char,xcryst,xinbox) 
 WRITE(6,*)'STOP',pot_dist,r/real(m),cost
 CONTAINS
!
 REAL FUNCTION repulsive_potential_Lowenstein(r,p0,p1,p2)
  IMPLICIT NONE
  REAL              :: r,p0,p1
  CHARACTER (LEN=6) :: p2
  if(p2/='coulom'.and.p2/='london') p2='zerooo'
  IF (r <  p0) repulsive_potential_Lowenstein = p1
  if (r >= p0.and.p2=='coulom') repulsive_potential_Lowenstein = (r-p0)**(-1)
  if (r >= p0.and.p2=='london') repulsive_potential_Lowenstein = (r-p0)**(-6)
  if (r >= p0.and.p2=='zerooo') repulsive_potential_Lowenstein = 0.00
 END FUNCTION repulsive_potential_Lowenstein
!
 SUBROUTINE MonteCarlo(potential,n_atoms,dist_matrix,n_Al,SEED,xcryst,label,cell_0,rv,exito,coste)
  IMPLICIT NONE
  INTEGER, intent(in) :: n_atoms,n_Al,SEED
  REAL,    intent(in) :: xcryst(0:3,1:n_atoms),dist_matrix(n_atoms,n_atoms)
  REAL,    intent(in) :: cell_0(1:6),rv(1:3,1:3)
  CHARACTER (LEN=2)   :: label(1:n_atoms)
  INTEGER             :: i,j,k
  REAL                :: energia(1:2) = 0.0
  REAL                :: eta = 0.0
  REAL                :: temperature = 0.00001
  REAL                :: delta = 0.0
  REAL,    PARAMETER  :: infinite = 9999999.999999
  REAL,    intent(out):: exito,coste
  CHARACTER (LEN=6),intent(in):: potential
!
  exito = 0.0
  MC_step: DO k=0,n_atoms-1
    IF(k==0)then
      energia(1) = energy(n_atoms,n_Al,xcryst,label,dist_matrix,cell_0,rv,potential)
      eta = energia(1)
    END IF
   !energia(2) = infinite + 1.0
   !DO WHILE ( energia(2) >= infinite )
     i = INT(R4_UNIFORM(1.0,real(n_atoms)+1.0,SEED))
     DO WHILE ( label(i) /= 'Si' )
       i = INT(R4_UNIFORM(1.0,real(n_atoms)+1.0,SEED))
     END DO
     j = INT(R4_UNIFORM(1.0,real(n_atoms)+1.0,SEED))
     DO WHILE ( label(j) /= 'Al' )
       j = INT(R4_UNIFORM(1.0,real(n_atoms)+1.0,SEED))
     END DO
! identity_change
     label(i) = 'Al'
     label(j) = 'Si'
     energia(2) = energy(n_atoms,n_Al,xcryst,label,dist_matrix,cell_0,rv,potential)
   !END DO
    IF ( energia(2) > energia(1) ) THEN
       eta = R4_UNIFORM(0.0,1.0,SEED)
       IF ( eta > exp(-(energia(2)-energia(1))/temperature) ) THEN
       !   ! se rechaza el movimento
          label(i) = 'Si'
          label(j) = 'Al'
       ELSE
          exito = exito + 1.0
          energia(1) = energia(2)
          ! se aprueba aun siendo la energia mayor
       ENDIF
    ELSE
       energia(1) = energia(2)
       exito = exito + 1.0
    END IF
  END DO MC_step
  exito   = exito/REAL(n_atoms*n_atoms)
  coste = energia(1)
  eta = coste - eta
! coste, ocurrencia, exito
  WRITE(6,'(f14.4,2x,2f14.6)')coste,exito,eta
  RETURN
 END SUBROUTINE MonteCarlo
!
 REAL FUNCTION energy(n,m,x,lab,matrix,cell_0,rv,potential)
  IMPLICIT NONE
  INTEGER, intent(in) :: n,m
  REAL,    intent(in) :: cell_0(1:6),rv(1:3,1:3),x(0:3,1:n),matrix(n,n)
  CHARACTER (LEN=2)   :: lab(1:n,1:2)
  REAL                :: atom(3),ouratom(3),s
  REAL,    PARAMETER  :: infinite = 9999999.999999
  INTEGER             :: l = 0
  CHARACTER (LEN=6),intent(in):: potential
!
  energy = 0.0
  l=0
  DO i =1, n
    IF(label(i)=='Al')THEN
     DO j =i+1,n
       IF(label(j)=='Al')THEN
        l=l+1
        !forall ( k=1:3)
        !  atom(k)    = x(k,i)
        !  ouratom(k) = x(k,j)
        !end forall
        !call make_distances(cell_0,ouratom,atom,rv,s)
        s = matrix(i,j)
        energy = energy + repulsive_potential_Lowenstein(s,3.4,infinite,potential)
       ENDIF
     END DO
    END IF
  END DO
  energy = energy/real(l)
 END FUNCTION energy
!
 SUBROUTINE cell(rv,vr,cell_0)
! ======================
! GULP
! ======================
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
 WRITE(6,'(a)')'Box:'
 DO i=1,3
    WRITE(6,'(F14.7,F14.7,F14.7)')( rv(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 WRITE(6,'(a)')'Inverse Box:'
 DO i=1,3
    WRITE(6,'(F14.7,F14.7,F14.7)')( vr(i,j), j=1,3 )
 ENDDO
 WRITE(6,'(a)')'----------------------------------------'
 RETURN
 END SUBROUTINE cell
!
 SUBROUTINE uncell(rv,cell_0)
! ======================
! GULP
! ======================
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
!==========================================================
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
 SUBROUTINE make_distances(cell_0,r2,r1,rv,dist)
! {{ prepara para el calculo de distancias en una celda triclinica }}
 IMPLICIT NONE
 REAL,    intent(in)  :: r1(3),r2(3),rv(3,3),cell_0(6)    ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist
 REAL                 :: d_image(1:27),image(3,27)        ! array de distancias
! REAL                 :: distance
 INTEGER              :: k,l,m,n,o,i,j                    ! variables mudas
 REAL                 :: atom(3),ouratom(3)               ! coordenadas preparadas
! {{ calculamos la matriz de distancias }}
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
 subroutine get_seed(seed)
  implicit none
  integer day,hour,i4_huge,seed,milli,minute,month,second,year
  parameter (i4_huge=2147483647)
  double precision temp
  character*(10) time
  character*(8) date
  call date_and_time (date,time)
  read (date,'(i4,i2,i2)')year,month,day
  read (time,'(i2,i2,i2,1x,i3)')hour,minute,second,milli
  temp=0.0D+00
  temp=temp+dble(month-1)/11.0D+00
  temp=temp+dble(day-1)/30.0D+00
  temp=temp+dble(hour)/23.0D+00
  temp=temp+dble(minute)/59.0D+00
  temp=temp+dble(second)/59.0D+00
  temp=temp+dble(milli)/999.0D+00
  temp=temp/6.0D+00
!  Force 0 < TEMP <= 1.
  DO
    IF(temp<=0.0D+00 )then
       temp=temp+1.0D+00
       CYCLE
    ELSE
       EXIT
    ENDIF
  ENDDO
  DO
    IF(1.0D+00<temp)then
       temp=temp-1.0D+00
       CYCLE
    else
       EXIT
    ENDIF
  ENDDO
  seed=int(dble(i4_huge)*temp)
  if(seed==0)       seed = 1
  if(seed==i4_huge) seed = seed-1
  RETURN
 END subroutine get_seed

 REAL function r4_uniform(b1,b2,seed)
  implicit none
  real b1,b2
  integer i4_huge,k,seed
  parameter (i4_huge=2147483647)
  if(seed==0)then
       write(*,'(b1)')' '
       write(*,'(b1)')'R4_UNIFORM - Fatal error!'
       write(*,'(b1)')'Input value of SEED = 0.'
      stop '[ERROR] r4_uniform'
  end if
  k=seed/127773
  seed=16807*(seed-k*17773)-k*2836
  if(seed<0) then
    seed=seed+i4_huge
  endif
  r4_uniform=b1+(b2-b1)*real(dble(seed)* 4.656612875D-10)
  RETURN
 end function r4_uniform
!
 SUBROUTINE write_grafo(A,N,U)
  integer, intent(in) :: N,A(N,N),U
  integer :: X,Y
  nodos: DO X=1,N
   DO Y=X+1,N
      IF(A(X,Y)>0.0) THEN
        IF(X>Y)THEN
         WRITE(U,'(I5,1x,I5,1x,A)') X,Y,'+'
        ELSE
         WRITE(U,'(I5,1x,I5,1x,A)') Y,X,'+'
        ENDIF
      ENDIF
   ENDDO
  ENDDO nodos
  RETURN
 END SUBROUTINE write_grafo
!
 SUBROUTINE ESCRITURA_GRAFO(A,N,U)
!--------------------------------------------
! ADAPTAMOS LA MATRIX DE ADYACENCIA PARA QUE EL PROGRAMA DE
! REPRESENTACION DE GRAFOS LO ENTIENDA. LA SALIDA LA TIENEN UN ARCHIVO
! LLAMADO: "GRAFO.DOT"
!--------------------------------------------
  INTEGER N,X,Y,A(N,N),U
  OPEN(U,file='grafo.dot')
  write(U,*)'graph G {'
  WRITE(U,*)'node[label="",height=0.1,width=0.1,fontsize=1]'
  do X=1,N
    do Y=X+1,N
       if(A(X,Y)>0.0)then
         write(U,'(A,I4,A,I4,A)')'   ',X,' -- ',Y,' ;'
       endif
    enddo
  enddo
  write(U,*)'}'
  close(u)
  RETURN
 END SUBROUTINE
!
 SUBROUTINE write_gin(cell_0,xcryst,n_atoms,label)
  ! Catlow potential
  IMPLICIT NONE
  INTEGER,          intent(in) :: n_atoms
  CHARACTER (LEN=2),intent(in) :: label(1:n_atoms)
  REAL,             intent(in) :: xcryst(0:3,1:n_atoms),cell_0(1:6)
  OPEN(999,file='out.gin')
  WRITE(999,'(A)')'single conv nosym qok'
  WRITE(999,'(A)')'cell'
  WRITE(999,'(6f9.5)') (cell_0(j) , j=1,6)
  WRITE(999,'(A)')'fractional'
  do i=1,n_atoms
    WRITE(999,'(a,1x,3f10.5)')label(i),xcryst(1,i),xcryst(2,i),xcryst(3,i)
  enddo
  WRITE(999,'(A)')'species 6'
  WRITE(999,'(A)')'Na     core    1.000000'
  WRITE(999,'(A)')'Sr     core    2.000000'
  WRITE(999,'(A)')'Al     core    3.000000'
  WRITE(999,'(A)')'Si     core    4.000000'
  WRITE(999,'(A)')'O2     core    0.869020'
  WRITE(999,'(A)')'O2     shel   -2.869020'
  WRITE(999,'(A)')'buck Na core O2 shel  1226.84000  0.306500  0.000000 0.00 10.00'
  WRITE(999,'(A)')'buck Na core Na core  7895.40000  0.170900  0.000000 0.00 10.00'
  WRITE(999,'(A)')'buck Sr core O2 shel  1952.39000  0.336850 19.22000  0.00 10.00'
  WRITE(999,'(A)')'buck Al core O2 shel  1460.30000  0.299120  0.000000 0.00 10.00'
  WRITE(999,'(A)')'buck Si core O2 shel  1283.90700  0.320520 10.66158  0.00 10.00'
  WRITE(999,'(A)')'buck O2 shel O2 shel 22764.0000   0.149000 27.87900  0.00 12.00'
  WRITE(999,'(A)')'spring'
  WRITE(999,'(A)')'O2 74.920000'
  WRITE(999,'(A)')"three regular Si core O2 shel O2 shel 2.0972 109.47000 &"
  WRITE(999,'(A)')'  0.000  1.800  0.000  1.800  0.000  3.200'
  WRITE(999,'(A)')"three regular Al core O2 shel O2 shel 2.0972 109.47000 &"
  WRITE(999,'(A)')'  0.000  1.800  0.000  1.800  0.000  3.200'
  WRITE(999,'(A)')'cutp     16.00000    1.00000'
  WRITE(999,'(A)')'switch_min rfo   gnorm     0.001000'
  WRITE(999,'(A)')'dump every 1     out.res'
  WRITE(999,'(A)')'output cif out.cif'
  CLOSE(999)
  RETURN
 END SUBROUTINE write_gin
END PROGRAM farthrest_nodes_subtitutions
