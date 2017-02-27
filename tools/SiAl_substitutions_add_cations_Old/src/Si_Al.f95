PROGRAM Si_Al
! {{ programa principal }}
 IMPLICIT NONE
! Modified:
 INTEGER :: i,j,k,CONT,CONTO,IERR
 INTEGER :: n_Si        ! numero de posiciones Si + Al
 INTEGER :: n_Al        ! numero de Al
 INTEGER :: n_Sr,n_Na
 INTEGER :: n_atoms
 INTEGER :: n_lines     ! numero de lines del fichero input.cif 
 INTEGER, PARAMETER :: HISTMAX=100,SPACE=1
 INTEGER, ALLOCATABLE :: AL(:)
 INTEGER :: SEED        ! random.f95
 REAL :: R4_UNIFORM     ! random.f95
 REAL :: D1,D2,D3,MU,LD,RV(3,3)
 REAL :: X,Y,Z,Q,P,HIST(1:HISTMAX),M,SUMA,RHO,volume,area,xi,yi,zi
 REAL :: atom(3),ouratom(3),distance,cell_0(6)
 REAL, ALLOCATABLE :: A(:,:),DIST(:,:)
 REAL, ALLOCATABLE :: xfrac(:,:),xinbox(:,:)
 CHARACTER (LEN=4) :: ATOML,MOL
 CHARACTER (LEN=1), ALLOCATABLE :: A_CHAR(:,:)
 CHARACTER (LEN=3) :: ELEM,ELEM0
 CHARACTER (LEN=2), ALLOCATABLE :: LABEL(:,:)
 LOGICAL :: SUBTITUTION,CORE_SHELL
 OPEN(1,FILE="input.cif")
 OPEN(2,FILE="output")
 OPEN(22,FILE="aluminios.txt")
 OPEN(4,FILE='cations.txt')
! {{ lectura de input }}
 READ(*,*)n_Al,n_Sr,n_Na
 READ(*,*)SUBTITUTION,CORE_SHELL
 READ(*,*)n_lines,n_Si
 n_atoms=n_lines+n_Sr+n_Na 
 READ(*,*)SEED
 READ(*,*)cell_0(1),cell_0(2),cell_0(3),cell_0(4),cell_0(5),cell_0(6)
! {{ alicatamos las variables en memoria }}
 ALLOCATE(Al(1:n_Al),STAT=IERR)
 ALLOCATE(xfrac(0:3,1:n_atoms),STAT=IERR)
! ALLOCATE(xinbox(1:n_Si,0:3),STAT=IERR)
 ALLOCATE(LABEL(1:n_atoms,2),STAT=IERR)
 ALLOCATE(A(n_Si,n_Si),STAT=IERR)
 ALLOCATE(A_CHAR(n_Si,n_Si),STAT=IERR)
 ALLOCATE(DIST(n_Si,n_Si),STAT=IERR)
 IF(IERR/=0) STOP '[ERROR] variables sin alicatar en memoria.'
! {{ RamdomSeed, GET_SEED esta en random.f95
 IF (seed==0) CALL GET_SEED(seed)
 WRITE(*,*)'[RandomSeed] seed: ',seed
! }}
! {{ CELL esta en cell.f95, VOLUME y AREA estan en volume.f95
! calcula la matriz de cambio para pasar de variables fraccionarias
! a coordenadas reales
 CALL CELL(rv,cell_0)
 WRITE(*,*)'Volume = ',VOLUME(rv)
 WRITE(*,*)'Area = ',AREA(rv)
! }}
 WRITE(*,*)'[composition] atoms:',n_atoms
 WRITE(*,*)'Posiciones tetraedricas: ',n_Si
 WRITE(*,*)'ratio Si/Al: ',(n_Si-n_Al)/REAL(n_Al) 
 j  = 0  ! contador de tetraedros
 k = 0   ! contador de oxigenos
 read_cif: DO i=1,n_lines
   READ(1,*)elem,ELEM0,X,Y,Z,Q
   IF(ELEM0=='O'.or.ELEM=='O')THEN
     k=k+1
     xfrac(1,n_Si+k) = mod( x + 100.0, 1.0)
     xfrac(2,n_Si+k) = mod( y + 100.0, 1.0)
     xfrac(3,n_Si+k) = mod( z + 100.0, 1.0)
!     WRITE(2,'(A,1X,A,1X,F14.7,F14.7,F14.7)') 'O2 ','core',xfrac(1,n_Si+k),xfrac(2,n_Si+k),xfrac(3,n_Si+k)
!     IF(CORE_SHELL) WRITE(2,'(A,1X,A,1X,F14.7,F14.7,F14.7)') 'O2 ','shel',xfrac(1,n_Si+k),xfrac(2,n_Si+k),xfrac(3,n_Si+k)
!     COOR_O(k,1)=xfrac(1,n_Si+k)
!     COOR_O(k,2)=xfrac(2,n_Si+k)
!     COOR_O(k,3)=xfrac(3,n_Si+k)
     label(n_Si+k,1)=elem
     label(n_Si+k,2)=elem0
   ELSE IF(ELEM0=='Si'.or.elem0=='Al'.or.elem0=='Ge') THEN
     j=j+1
     xfrac(1,j) = mod( x + 100.0, 1.0)
     xfrac(2,j) = mod( y + 100.0, 1.0)
     xfrac(3,j) = mod( z + 100.0, 1.0)
!     COOR(j,1)=xfrac(1,j)
!     COOR(j,2)=xfrac(2,j)
!     COOR(j,3)=xfrac(3,j)
     LABEL(j,1)=ELEM
     LABEL(j,2)=ELEM0
   ENDIF
 ENDDO read_cif
! reiniciamos las variables mudas i,j,k
 distancias: DO i=1,n_Si
   DO j=i,n_Si
     A(i,j)=0.0
     A_CHAR(i,j)=' '
     IF(j/=i)THEN
       DO k=1,3
          atom(k)= xfrac(k,j)
          ouratom(k) = xfrac(k,i)
       ENDDO
       call make_distances(ouratom,atom,rv,distance) 
       DIST(i,j)=distance
       IF(dist(i,j)>=3.0.AND.dist(i,j)<4.0)THEN 
         a(i,j)=1.0
         A_CHAR(I,J)='@'
       ENDIF
     ENDIF
     DIST(j,i)=DIST(i,j)
     a(j,i)=a(i,j)
     A_CHAR(J,I)=A_CHAR(I,J)
   ENDDO
 ENDDO distancias
 PRINT*,'GRAFO:'
 CALL write_grafo(A,n_Si,200)
 CALL ESCRITURA_GRAFO(a,n_Si,500)
 DO I=1,n_Si
   WRITE(*,'(100(1X,A))')(A_CHAR(I,J),J=1,n_Si)
 ENDDO
! {{ el siguiente bloque revida la conectividad de cada tetraedro, debe ser 4
! independiente del grupo espacial SPACE
 conectivity: DO I=1,n_Si ! conectividad de los ATOM
   k=0                   ! k=4 en zeolitas
   DO j=1,n_Si
     k=k+a(i,j)
   ENDDO
   xfrac(0,i)=k
 ENDDO conectivity
 PRINT*,'Conectividad entre nodos (4 for zeolite):'
 WRITE(*,'(100(1X,I0))')(INT(xfrac(0,k)),k=1,n_Si)
! }}
! {{ ahora realizamos las sustituciones de Si por Al
! se realiza si la variable SUBSTITUTIONS=.TRUE.
 PRINT*,'# SI > AL '
 P=R4_UNIFORM(1.0,REAL(n_Si)+1.0,SEED)
 CONT=0
 subtitutions: IF(SUBTITUTION)THEN 
  choose_Al: DO WHILE (CONT<N_AL)     ! elijo N_AL Si y los cambio por Al
   Q=R4_UNIFORM(1.0,REAL(n_Si)+1.0,SEED)! 1) numero random que no se repita (Q/=P)
   K=INT(Q)                             ! 2) que no se sustituya un aluminio ya sustituido
   J=0                                  ! 3) que se respete Lowenstein (J=0)
!   {{ subbloque de revision de la regla de Lowenstein (J/=0)
! tambien se añade la posibilidad de que el generador !
! de numeros aleatorios haya fallado                 (P==Q)
   low: DO I=1,n_Si
     IF(label(k,1)=='Al') J=J+1
     IF(A(K,I)==1.0.AND.LABEL(I,1)=='Al') J=J+1
   ENDDO low
   IF(P==Q.OR.J/=0) CYCLE choose_Al
!   }}
   CONT=CONT+1
   P=Q
   label(K,1)='Al'
   label(K,2)='Al'
   AL(CONT)=K
  ENDDO choose_Al
! }}
! {{ el siguiente bloque calcula a cuantos alumnios esta conectado cada tetraedro:
  N_Al_Si: DO K=1,n_Si
    xfrac(0,k)=0.0
    inner: DO I=1,n_Si
      IF(A(K,I)==1.0.AND.LABEL(I,1)=='Al') xfrac(0,k)=xfrac(0,k)+1.0
    ENDDO inner
  ENDDO N_Al_Si
  PRINT*,'Number of Al by Si:'
  WRITE(*,'(100(1X,I0))')(INT(xfrac(0,k)),k=1,n_Si)
! }}
! {{ calculo la distancia media de atomos de Al en la cela unidad <Al-Al> 
  SUMA=0.0
  K=0
  IF (n_Al>0) THEN
   averages: DO I=1,N_AL
    DO J=I+1,N_AL
      IF(DIST(AL(I),AL(J))>0.0)THEN
        SUMA=SUMA+DIST(AL(I),AL(J))
        K=K+1
      ENDIF
    ENDDO
   ENDDO averages
  ENDIF
  MU=SUMA/REAL(k)
  WRITE(*,*)'[Aluminiums]:'
  WRITE(*,'(100(1X,I0))')(AL(K),K=1,N_AL)
  WRITE(*,*)'<Al Al>:',MU,'1e-10 [m]'
 ENDIF subtitutions
! }}
 CALL ADD_CATIONS(space,n_Sr,n_Na,n_atoms,xfrac,label,SEED,RV)
!  CALL escritura_cif(COOR,COOR_O,COOR_CAT,n_Si,n_lines,N_AL,N_SR,N_NA,LABEL,cell_0,rv)
 call uncell(cell_0,rv)
 call GULP(xfrac,n_atoms,label,cell_0,CORE_SHELL)
! }}
! {{ cerramos las salidas y desalicatamos la memoria
 CLOSE(1)
 CLOSE(2)
 CLOSE(22)
 CLOSE(4)
 CLOSE(5)
 deallocate(xfrac)
 DEALLOCATE(LABEL)
 DEALLOCATE(A,A_CHAR)
 DEALLOCATE(DIST)
! }}
 STOP "[subtitutions done]"
END PROGRAM
