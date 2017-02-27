SUBROUTINE make_distances(r1,r2,rv,dist)
 implicit none
 REAL,    intent(in)  :: r1(3),r2(3)                       ! coordenadas y matriz de cambio
 REAL,    intent(out) :: dist,rv(1:3,1:3)                  ! matriz de distancias N X N
 REAL                 :: d_image(1:27),distance            ! array de distancias
 INTEGER              :: k,l,m,n                           ! variables mudas
 REAL                 :: atom(3),ouratom(3)                ! coordenadas preparadas
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
          d_image(k)=distance(atom,ouratom,rv)
!          WRITE(*,'(I2,1X,I2,1X,I2,2X,F14.7)')l,m,n,d_image(k)
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
 PBC: DO j=1,3
    IF(atom(j)-ouratom(j)>0.5)THEN
      per_atom(j) = atom(j) - 1.0
    ELSE IF(ouratom(j)-atom(j)>0.5)THEN
        per_atom(j) = 1.0 + atom(j)
    ELSE
        per_atom(j) = atom(j)
    ENDIF
 ENDDO PBC
 frac2real: DO j=1,3
   o_ouratom(j) = rv(j,1)*ouratom(1)  + rv(j,2)*ouratom(2)  + rv(j,3)*ouratom(3)
   o_atom(j)    = rv(j,1)*per_atom(1) + rv(j,2)*per_atom(2) + rv(j,3)*per_atom(3)
   dist(j) = o_ouratom(j) - o_atom(j)
 ENDDO frac2real
 DISTANCE = sqrt(dist(1)*dist(1) + dist(2)*dist(2) + dist(3)*dist(3))
END FUNCTION
