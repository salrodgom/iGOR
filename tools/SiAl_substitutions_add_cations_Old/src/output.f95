SUBROUTINE write_grafo(A,N,U)
! escribe grafo.red
 integer, intent(in) :: N,A(N,N),U
 integer :: X,Y
 OPEN(U,FILE="grafo.red")
 nodos: DO X=1,N
   DO Y=1,N
      IF((A(X,Y)>0.0)) WRITE(U,'(I3,1x,I3,1x,A)') X,Y,'+'
   ENDDO
 ENDDO nodos
 WRITE(*,*)'[loops] Escrito grafo.red en formato: i j +'
 CLOSE(U)
 RETURN 
end subroutine write_grafo

SUBROUTINE ESCRITURA_GRAFO(A,N,U)
!--------------------------------------------
! ADAPTAMOS LA MATRIS DE ADYACENCIA PARA QUE EL PROGRAMA DE
! REPRESENTACION DE GRAFOS LO ENTIENDA. LA SALIDA LA TIENEN UN ARCHIVO
! LLAMADO: "GRAFO.DOT"
!--------------------------------------------
 INTEGER N,X,Y,A(N,N),U
 OPEN(U,FILE='grafo.dot')
 write(U,*)'graph G {'
 WRITE(U,*)'node[label="",height=0.1,width=0.1,fontsize=1]'
 do X=1,N
    do Y=1,N
       if((Y>X).and.(A(X,Y)>0.0))then
         WRITE(U,'(A,I3,A,I3,A)')'   ',X,' -- ',Y,' ;'
       endif
    enddo
 enddo
 WRITE(*,*)'[graph] Escrito grafo.dot en formato: graph G { i -- j ; }'
 write(U,*)'}'
 CLOSE(U)
 RETURN
end SUBROUTINE

subroutine GULP(xfrac,n,label,cell_0,CORE_SHELL)
 implicit none
 integer, intent(in)              :: n
 real, intent(in)                 :: xfrac(0:3,1:n),cell_0(6)
 CHARACTER (LEN=2), intent(in)    :: label(1:n,2)
 LOGICAL, intent(in)              :: CORE_SHELL
 integer :: k
!
 OPEN(1,FILE='gin')
 WRITE(1,'(A)')'opti conp phon freq noden'
 WRITE(1,'(A)')'title'
 WRITE(1,'(A,1X,I3)')' Al/Si substitutions',n
 WRITE(1,'(A)')'end'
 WRITE(1,'(A)')'cell'
 WRITE(1,'(6(F14.7))')(cell_0(k), k=1,6)
 WRITE(1,'(A)')'frac'
 write_: do k=1,n
  WRITE(1,'(A,1X,A,1X,F14.7,F14.7,F14.7)')LABEL(k,1),'core',xfrac(1,k),xfrac(2,k),xfrac(3,k)
  IF(CORE_SHELL.and.label(k,2)=='O') THEN
   WRITE(1,'(A,1X,A,1X,F14.7,F14.7,F14.7)')LABEL(K,1),'shel',xfrac(1,k),xfrac(2,k),xfrac(3,k)
  ENDIF
 enddo write_
 close(1)
 RETURN
end subroutine

SUBROUTINE escritura_cif(COOR,COOR_O,COOR_CAT,NSI,NLINES,N_AL,N_SR,N_NA,LABEL,cell_0,rv)
! Escribe el CIF P1.cif
 IMPLICIT NONE
 REAL :: COOR(1:NSI,0:3),COOR_O(1:NLINES-NSI,0:3),COOR_CAT(1:N_SR+N_NA,0:3)
 REAL :: volume,cell_0(6),rv(3)
 INTEGER :: I,NSI,NLINES,N_AL,N_SR,N_NA,U
 CHARACTER(2) :: LABEL(1:NSI,2)
 U=1000
 OPEN(U,FILE="P1.cif")
 WRITE(U,'(A)')'data_subtitutions'
 WRITE(U,'(A)')'_audit_creation_method IGOR'
 WRITE(U,'(A)')"_audit_author_name 'S.R.G.Balestra'"
 WRITE(U,'(A,F14.7)')'_cell_length_a',cell_0(1)
 WRITE(U,'(A,F14.7)')'_cell_length_b',cell_0(2)
 WRITE(U,'(A,F14.7)')'_cell_length_c',cell_0(3)
 WRITE(U,'(A,F14.7)')'_cell_angle_alpha',cell_0(4)
 WRITE(U,'(A,F14.7)')'_cell_angle_beta',cell_0(5)
 WRITE(U,'(A,F14.7)')'_cell_angle_gamma',cell_0(6)
 WRITE(U,'(A,F14.7)')'_cell_volume',volume(rv)
 WRITE(U,'(A)')'_symmetry_cell_setting cubic'
 WRITE(U,'(A)')"_symmetry_space_group_name_Hall 'P 1'"
 WRITE(U,'(A)')"_symmetry_space_group_name_H-M 'P 1'"
 WRITE(U,'(A)')'_symmetry_Int_Tables_number 1'
 WRITE(U,'(A)')"_symmetry_equiv_pos_as_xyz 'x,y,z'"
 WRITE(U,'(A)')'loop_'
 WRITE(U,'(A)')'_atom_site_label'
 WRITE(U,'(A)')'_atom_site_fract_x'
 WRITE(U,'(A)')'_atom_site_fract_y'
 WRITE(U,'(A)')'_atom_site_fract_z' 
 tetraedros: DO I=1,NSI
   WRITE(U,'(A,F14.7,F14.7,F14.7)')LABEL(I,1),COOR(I,1),COOR(I,2),COOR(I,3)
 END DO tetraedros
 oxigens: DO I=1,NLINES-NSI
   WRITE(U,'(A,F14.7,F14.7,F14.7)')'O2',COOR_O(I,1),COOR_O(I,2),COOR_O(I,3)
 END DO oxigens
 cations: DO I=1,N_SR+N_NA
   IF(I<=N_SR) WRITE(U,'(A,F14.7,F14.7,F14.7)')'Sr',COOR_CAT(I,1),COOR_CAT(I,2),COOR_CAT(I,3)
   IF(I>N_SR)  WRITE(U,'(A,F14.7,F14.7,F14.7)')'Na',COOR_CAT(I,1),COOR_CAT(I,2),COOR_CAT(I,3)
 ENDDO cations
 CLOSE(U)
 RETURN
END SUBROUTINE escritura_cif
