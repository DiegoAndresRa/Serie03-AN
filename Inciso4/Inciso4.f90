PROGRAM ResolucionLU !Resolución mediante descomposición LU con pivoteo y escalamiento
    IMPLICIT NONE
     INTEGER,PARAMETER:: n=5 !orden del sistema de ecuaciones
     INTEGER:: o(1:n),i,j,er
     CHARACTER (LEN=5):: TAMANO,INDICE
     REAL(KIND=8)::a(1:n,1:n),b(1:n),x(1:n)
     !Declaración de la matriz [A] a la cual se obtendrá su descomposición [L][U]  
     a(1,1)=6.1d0;a(1,2)=-2.4d0;a(1,3)=23.3d0;a(1,4)=-16.4d0;a(1,5)=-8.9d0
     a(2,1)=-14.2d0;a(2,2)=-31.6d0;a(2,3)=-5.8d0;a(2,4)=9.6d0;a(2,5)=23.1d0
     a(3,1)=10.5d0;a(3,2)=46.1d0;a(3,3)=-19.6d0;a(3,4)=-8.8d0;a(3,5)=-41.2d0
     a(4,1)=37.3d0;a(4,2)=-14.2d0;a(4,3)=62.d0;a(4,4)=14.7d0;a(4,5)=-9.6d0
     a(5,1)=0.8d0;a(5,2)=17.7d0;a(5,3)=-47.5d0;a(5,4)=-50.2d0;a(5,5)=29.8d0

     !Declaración del vector {B}
     b(1)=121.7d0;b(2)=-87.7d0;b(3)=10.8d0;b(4)=61.3d0;b(5)=-27.8d0
     !impresión de la matriz original
     WRITE(*,*) '      Matriz [A] '
     WRITE(TAMANO,'(I3)')n
     DO i=1,n
         WRITE(*,'('//TAMANO//'(f12.7,2x))') (a(i,j),j=1,n)  !se usa DO IMPLICITO para escribir la matriz
     ENDDO
     !Obtención de la matriz que contiene cifradas [L] y [U]
     CALL DESCOMPOSICION_LU(a,n,o,er)
     !Obtención mediante sustitución adelante y hacia atrás del vector solución {X}
     IF (er.NE.-1)THEN
         CALL SUSTITUCION(a,o,n,b,x)
     ENDIF
         !Impresión de resultados para {X}
     WRITE(*,*) 
     WRITE(*,*) '      Solución {X}'
     DO i=1,n
        WRITE(*,'(A2,I2,A3,f10.4)') 'x(',i,')= ',x(i)
     ENDDO    
 ENDPROGRAM
 
 
 !*****************************************************************************
 !**************Subrutina para descomposicion LU, en la que se llama a********* 
 !************************** subrutina PIVOTEO*********************************
 SUBROUTINE DESCOMPOSICION_LU(a,n,o,er)
     IMPLICIT NONE
     REAL (KIND=8),PARAMETER::tol=10.d0**(-8)
     REAL(KIND=8)::a(1:n,1:n),s(1:n),factor
     INTEGER:: o(1:n),i,j,k,n,er
     
     DO i=1,n
         o(i)=i    !este vector lleva la cuenta el número de pivoteos.
         s(i)=DABS(a(i,i))
         DO j=2,n      !Se determina cual es el elemento máximo de cada renglón en [A], es decir,s(i)
             IF (DABS(a(i,j)).GT.s(i)) THEN 
                 s(i)=DABS(a(i,j))
             ENDIF
         ENDDO
     ENDDO
     
     DO k=1,n-1,1
         CALL PIVOTEO(a,o,s,n,k)    !Se reordena la matriz [A] usando pivoteo parcial
         IF (DABS(a(o(k),k)/s(o(k))).LT.tol) THEN    !Determina si la matriz [A] es singular
             er=-1
             WRITE(*,*) a(o(k),k)/s(o(k))
             WRITE(*,*) 'Matriz singular'
             GOTO 100    !ESCAPA DEL CICLO Y CALCULOS SI [A] ES SINGULAR 
         ENDIF
         
             DO i=k+1,n
                 factor=a(o(i),k)/a(o(k),k)
                 a(o(i),k)=factor
                     DO j=k+1,n
                         a(o(i),j)=a(o(i),j)-factor*a(o(k),j)
                     ENDDO
             ENDDO              
     ENDDO
 100 IF (DABS(a(o(n),n)/s(o(n))).LT.tol) THEN    !Determina si la matriz [A] es singular
             er=-1
             WRITE(*,*) a(o(n),k)/s(o(n))     
     ENDIF    
  ENDSUBROUTINE DESCOMPOSICION_LU
 
 
 !*****************************************************************
 !***********Subrutina para realizar pivoteo parcial***************
 SUBROUTINE PIVOTEO(a,o,s,n,k)
     IMPLICIT NONE
     REAL(KIND=8)::a(1:n,1:n),s(1:n),factor
     REAL(KIND=8) a_mayor,titere
     INTEGER:: o(1:n),i,j,k,n,p
     p=k     !Posición del pivote en la columna "k"
     a_mayor=DABS(a(o(k),k)/s(o(k))) !elemento pivote normalizado
     
     DO i=k+1,n
         titere=DABS(a(o(i),k))/s(o(i))
         IF (titere.GT.a_mayor) THEN
             a_mayor=titere
             p=i
         ENDIF
     ENDDO
     !En el siguiente proceso solo se intercambian las posiciones de o(p) y o(k), eso evita tener 
     !que intercambiar todos los elementos a(p,j) por a(k,j)  (j=1,2,...,n)
     titere=o(p)
     o(p)=o(k)
     o(k)=titere
 ENDSUBROUTINE PIVOTEO
 !*****************************************************************************
 !Subrutina que realiza sustitución hacia adelante y hacia atrás para encontrar
 !**********resolver  [L]{D}={B}, y posteriormente resolver [U]{X}={D}*********
 SUBROUTINE SUSTITUCION(a,o,n,b,x)
     IMPLICIT NONE
     REAL(KIND=8)::a(1:n,1:n),b(1:n),x(1:n),sum
     INTEGER:: o(1:n),i,j,k,n
     !SUSTITUCIÓN HACIA ADELANTE
     DO i=2,n
         sum=b(o(i))
         DO j=1,(i-1),1
             sum=sum-a(o(i),j)*b(o(j))
         ENDDO
         b(o(i))=sum
     ENDDO
     !SUSTITUCIÓN HACIA ATRÁS
     x(n)=b(o(n))/a(o(n),n)
     DO i=n-1,1,-1
         sum=0.d0
         DO j=i+1,n
             sum=sum+a(o(i),j)*x(j)
         ENDDO
         x(i)=(b(o(i))-sum)/a(o(i),i)
     ENDDO
 ENDSUBROUTINE SUSTITUCION
 
 