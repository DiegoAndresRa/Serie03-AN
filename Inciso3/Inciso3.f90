PROGRAM ResolucionLU !Resolución mediante descomposición LU con pivoteo y escalamiento
    IMPLICIT NONE
     INTEGER,PARAMETER:: n=4 !orden del sistema de ecuaciones
     INTEGER:: o(1:n),i,j,er
     CHARACTER (LEN=4):: TAMANO,INDICE
     REAL(KIND=8)::a(1:n,1:n),b(1:n),x(1:n)
     !Declaración de la matriz [A] a la cual se obtendrá su descomposición [L][U]  
     a(1,1)=23.42d0;a(1,2)=-16.89d0;a(1,3)=57.31d0;a(1,4)=82.6d0
     a(2,1)=-14.77d0;a(2,2)=-38.29d0;a(2,3)=92.36d0;a(2,4)=-4.36
     a(3,1)=-77.21d0;a(3,2)=71.26d0;a(3,3)=-16.55d0;a(3,4)=43.09d0
     a(4,1)=91.82d0;a(4,2)=81.43d0;a(4,3)=33.94d0;a(4,4)=57.22d0
     
     !Declaración del vector {B}
     b(1)=2158.36d0;b(2)=-1123.02d0;b(3)=3248.71d0;b(4)=235.25
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
 
 