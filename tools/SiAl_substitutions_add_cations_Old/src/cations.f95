SUBROUTINE ADD_CATIONS(SPACE,n_Sr,n_Na,n_atoms,xfrac,label,SEED,RV)
 IMPLICIT NONE
 integer, intent(in) :: n_Sr,n_Na,n_atoms,SEED,SPACE
 INTEGER :: I,J,K,O
 INTEGER :: cont,cont0
!    
 REAL :: xfrac(0:3,1:n_atoms)
 REAL :: R4_UNIFORM,r,xvirtual(0:100,0:3)
 REAL :: atom(3),ouratom(3),rv(3,3),distance
 LOGICAL    :: OVERLOAD
 CHARACTER (LEN=2) :: LABEL(1:n_atoms,2)
 PRINT*,'OVERLOAD?'
 PRINT*,'========='
 cations: DO i=1,N_SR+N_NA ! numero de cationes anadidos
 cont0=1
 xvirtual(0,0)=56024.0 ! lejísimo
 virtuales_x: do cont=1,10
   OVERLOAD=.true. ! para que entre en el bucle
   k=0
   IF (i<=N_SR) label(n_atoms-n_Sr-n_Na+i,1)='Sr'
   IF (i<=N_SR) label(n_atoms-n_Sr-n_Na+i,2)='Sr'
   IF (i >N_SR) label(n_atoms-n_Sr-n_Na+i,1)='Na'
   IF (i >N_SR) label(n_atoms-n_Sr-n_Na+i,2)='Na'
   WRITE(*,*)LABEL(n_atoms-n_Sr-n_Na+i,2),i,':'
   over: DO WHILE (OVERLOAD.EQV..true.)
     OVERLOAD=.false. ! se activa si hay overload
     k=k+1
     r=R4_UNIFORM(0.0,1.0,SEED)
     IF (K>=100) THEN
       xvirtual(cont,1)=R4_UNIFORM(0.0,1.0,SEED)
       xvirtual(cont,2)=R4_UNIFORM(0.0,1.0,SEED)
       xvirtual(cont,3)=R4_UNIFORM(0.0,1.0,SEED)
     else
     space_group: SELECT CASE (SPACE)
       CASE (229)
         IF(LABEL(n_atoms-n_Sr-n_Na+i,2)=='Na') THEN
           position_Na_im3m: SELECT CASE (int(100.0*r))
           CASE (0:453) ! 0.453 es la occupancy
               xvirtual(cont,1)=MOD(0.3145+0.5*INT(R4_UNIFORM(-2.0,2.0,SEED))+100.0,1.0)
               xvirtual(cont,2)=MOD(0.3145+0.5*INT(R4_UNIFORM(-2.0,2.0,SEED))+100.0,1.0)
               xvirtual(cont,3)=MOD(0.3145+0.5*INT(R4_UNIFORM(-2.0,2.0,SEED))+100.0,1.0)
           CASE DEFAULT
               xvirtual(cont,1)=R4_UNIFORM(0.0,1.0,SEED)
               xvirtual(cont,2)=R4_UNIFORM(0.0,1.0,SEED)
               xvirtual(cont,3)=R4_UNIFORM(0.0,1.0,SEED)
           END SELECT position_Na_im3m
         ELSE
           xvirtual(cont,1)=R4_UNIFORM(0.0,1.0,SEED)
           xvirtual(cont,2)=R4_UNIFORM(0.0,1.0,SEED)
           xvirtual(cont,3)=R4_UNIFORM(0.0,1.0,SEED)
         ENDIF
       CASE (217)
         IF(LABEL(n_atoms-n_Sr-n_Na+i,2)=='Na') THEN
           position_Na_i43m: SELECT CASE (int(100.0*r))
           CASE (0:284)
              xvirtual(cont,1)=MOD(0.2989+0.5*INT(R4_UNIFORM(-21.0,21.0,SEED))+100.0,1.0)
              xvirtual(cont,2)=MOD(0.2989+0.5*INT(R4_UNIFORM(-21.0,21.0,SEED))+100.0,1.0)
              xvirtual(cont,3)=MOD(0.2989+0.5*INT(R4_UNIFORM(-21.0,21.0,SEED))+100.0,1.0)
           CASE (285:637)
              xvirtual(cont,1)=MOD(0.0226+0.5*INT(R4_UNIFORM(-21.0,21.0,SEED))+100.0,1.0)
              xvirtual(cont,2)=MOD(0.0226+0.5*INT(R4_UNIFORM(-21.0,21.0,SEED))+100.0,1.0)
              xvirtual(cont,3)=MOD(0.4224+0.5*INT(R4_UNIFORM(-21.0,21.0,SEED))+100.0,1.0)
           CASE DEFAULT
              xvirtual(cont,1)=R4_UNIFORM(0.0,1.0,SEED)
              xvirtual(cont,2)=R4_UNIFORM(0.0,1.0,SEED)
              xvirtual(cont,3)=R4_UNIFORM(0.0,1.0,SEED)
           END SELECT position_Na_i43m
         ELSE
           xvirtual(cont,1)=R4_UNIFORM(0.0,1.0,SEED)
           xvirtual(cont,2)=R4_UNIFORM(0.0,1.0,SEED)
           xvirtual(cont,3)=R4_UNIFORM(0.0,1.0,SEED)
         ENDIF
       CASE DEFAULT
         xvirtual(cont,1)=R4_UNIFORM(0.0,1.0,SEED)
         xvirtual(cont,2)=R4_UNIFORM(0.0,1.0,SEED)
         xvirtual(cont,3)=R4_UNIFORM(0.0,1.0,SEED)  
       END SELECT space_group
     endif
     DO o=1,3
        ouratom(o)=xvirtual(cont,o)
     ENDDO
     overload_search: DO j=1,n_atoms-n_Sr-n_Na+i
       IF(j/=i)THEN
         DO o=1,3
           atom(o)=xfrac(o,j)
         ENDDO
         call make_distances(ouratom,atom,rv,distance)
         IF( distance < 2.2 )THEN
          OVERLOAD=.true.
          WRITE(*,*)label(j,2),j,k,distance
          CYCLE overload_search
         ENDIF        
         xvirtual(cont,0) = distance
         if(xvirtual(cont,0)<=xvirtual(0,0)) cont0=cont
        ENDIF
     ENDDO overload_search
   ENDDO over
  ENDDO virtuales_x
  xfrac(1,n_atoms-n_Sr-n_Na+i)=xvirtual(cont0,1)
  xfrac(2,n_atoms-n_Sr-n_Na+i)=xvirtual(cont0,2)
  xfrac(3,n_atoms-n_Sr-n_Na+i)=xvirtual(cont0,3)
 ENDDO cations
 RETURN
END SUBROUTINE
