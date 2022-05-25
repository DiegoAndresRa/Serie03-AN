PROGRAM GAUSS_PIVOTEOPARCIAL_ESCALAMIENTO
    IMPLICIT NONE
    INTEGER, PARAMETER:: n=4                   !orden del sistema de ecuaciones lineales
    REAL (KIND=8)::a(1:n,1:n),b(1:n),x(1:n),s(1:n)  !Asignación de las respectivas dimensiones de A=a(i,j);b(i),x(i)
    INTEGER er,i,j
    CHARACTER (LEN=4) CONTADOR
   CHARACTER (Len=5) LONGITUD_RENGLON
    !Declaración de la matriz de coeficientes A=(a_ij)
    a(1,1)=0.d0;a(1,2)=2.d0;a(1,3)=-1.d0;a(1,4)=-4.d0
    a(2,1)=1.d0;a(2,2)=-1.d0;a(2,3)=5.d0;a(2,4)=2.d0
    a(3,1)=3.d0;a(3,2)=3.d0;a(3,3)=-7.d0;a(3,4)=-1.d0
    a(4,1)=-1.d0;a(4,2)=-2.d0;a(4,3)=3.d0;a(4,4)=0.d0
    !Declaración de vector de coeficientes b=b(i)
    b(1)=2.d0;b(2)=-4.d0;b(3)=4.d0;b(4)=-7.d0

    !Se convierte el "n" númerico en texto y se almacena en Longitud_renglon para
    !poder mandar a pantalla con el formato establecido

    WRITE(LONGITUD_RENGLON,'(I3)') n

    !Solo para visualizar los datos de entrada
    WRITE(*,*)'Matriz A'
    DO i=1,n,1
        WRITE(*,'('//LONGITUD_RENGLON//'(f13.8,2x))') (a(i,j),j=1,n)   !DO IMPLICITO
    END DO
    WRITE(*,*)'Vector de coeficientes  B'
     DO i=1,n,1
        WRITE(*,'(f13.8)') b(i)  !DO IMPLICITO
    END DO

    !Eliminación Gaussiana del sistema de ecuaciones a la forma triangular
    CALL Gauss(a,b,s,n)
    !Obtención del vector solución x(i) a partir de sustitución hacia atrás
    CALL sustitucion_atras(a,b,x,n)
    !IMPRESIÓN DE RESULTADOS
    WRITE(*,*)'Soluciones  X=x(i)'
    DO i=1,n
        WRITE(*,'(2x,A,I1,A2,f21.16)') 'x',i,'= ',x(i)
        
    END DO
ENDPROGRAM

!********************Subrutina para elminación Gaussiana con pivoteo y escalamiento
SUBROUTINE Gauss(a,b,s,n)
    IMPLICIT NONE
    REAL (KIND=8),PARAMETER:: tol=10.D0**(-12)   !Tolerancia
    REAL (KIND=8)::a(1:n,1:n),b(1:n),x(1:n),s(1:n)
    INTEGER i,j,n,er
    !s(j)  almacena los elementos los "n" elementos a(i,1), es decir, la primera columna (véase [1]),
    !para posteriormente compararlos con los a(i,j) (j=2,3,4,...,n) del renglón "i"
    !y determinar cual es el mayor, es decir, toma el elemento máximo de cada renglón (véase el ciclo [2])

    !"er" es una variable entera que tiene la función de advertir si algún elemento de la diagonal (pivote)
    !es menor a una tolerancia, es decir, a(i,i)<tol, donde "tol" es un valor muy pequeño y positivo
    !que hace las veces de "cero".
    er=0
    DO i=1,n
        s(i)=DABS(a(i,1)) ![1]
        DO j=2,n  ![2]
            IF (DABS(a(i,j)).GT.s(i)) THEN
                s(i)=DABS(a(i,j))
            ENDIF
        ENDDO
    ENDDO
        CALL Eliminacion(a,b,s,n,er,tol) !Subrutina que realiza la eliminación Gaussiana
                                         !(triangular superior)
    IF (er.NE.-1) THEN                      !Se realiza la sustitución hacia atrás si
        CALL sustitucion_atras(a,b,x,n)
    ENDIF
ENDSUBROUTINE

!***********************************************************
SUBROUTINE Eliminacion(a,b,s,n,er,tol)
    IMPLICIT NONE
    REAL (KIND=8)::a(1:n,1:n),b(1:n),s(1:n)
    REAL (KIND=8):: tol,pivote
    INTEGER k,i,j,n,er
    DO k=1,n-1,1    !desplaza las operaciones al k-ésimo renglón, esto se nota para i=k+1, que realmente i[2,n]
        CALL PIVOTEO_PARCIAL_ESCALAMIENTO(k,a,b,s,n) !pivoteo realizado según sus elementos normalizados
        IF (DABS(a(k,k)/s(k)).LT. tol) THEN !Prueba si el elemento a(k,k) NORMALIZADO por s(k) es cercano a CERO
            er=-1
            GOTO 100!termina el ciclo si detecta un elemento de la diagonal (pivote NORMALIZADO) cercano a CERO
        END IF
        DO i=k+1,n   !Se realiza la eliminación Gaussiana volviendo la matriz de coeficientes en
                     !una triangular superior
            pivote=a(i,k)/a(k,k)
                DO j=k+1,n,1
                    a(i,j)=a(i,j)-pivote*a(k,j)
                ENDDO
                b(i)=b(i)-pivote*b(k)
        ENDDO
   ENDDO
100 IF (DABS(a(n,n)/s(n)).LT.tol) THEN !
            er=-1
    END IF
ENDSUBROUTINE

!*****************************Subrutina para pivoteo parcial con Escalamiento
SUBROUTINE PIVOTEO_PARCIAL_ESCALAMIENTO(k,a,b,s,n) !pivoteo realizado según sus elementos normalizados
    IMPLICIT NONE
    REAL (KIND=8):: a(1:n,1:n),b(1:n),s(1:n)
    REAL (KIND=8) a_mayor,titere
    INTEGER ii,jj,k,n,p
    p=k                     !Posición del pivote en la columna k
    a_mayor=ABS(a(k,k)/s(k))     !Elemento pivote NORMALIZADO en la columna k

!Este ciclo compara el pivote NORMALIZADO a(k,k)/s(k) con los elementos también NORMALIZADOS debajo de él
!y determina cual de todos es el mayor
    DO ii=k+1,n
        titere=ABS(a(ii,k)/s(ii))
        IF (titere .GT. a_mayor) THEN
            a_mayor=titere
            p=ii
        ENDIF
    ENDDO
 !La siguiente condición intercambia los renglones respecto al elemento más grande
 !almacenado al final del ciclo anterior (junto con su posición), si es que el pivote original
 !no es el de mayor magnitud
   IF (p .NE. k) THEN   !intercambio del renglón "k" con el renglón "j"

        DO jj=k,n   !justo en este ciclo se intercambian cada uno de los elementos
                    !del renglón a(k,jj) con el  mayor a(p,jj)
            titere=a(p,jj)
            a(p,jj)=a(k,jj)
            a(k,jj)=titere
        END DO
        titere=b(p) !Se intercambian b(p) y b(k) de posición
        b(p)=b(k)
        b(k)=titere
        !En este momento se intercambian los elementos s(k) por s(p)
        !(recuerde que s(j) (j=1,2,...,n) es la magnitud
        !el elemento más grande de cada renglón
        titere=s(p)
        s(p)=s(k)
        s(k)=titere
    END IF
ENDSUBROUTINE
!********************Subrutina para sustitución hacia atrás
SUBROUTINE sustitucion_atras(a,b,x,n)
    IMPLICIT NONE
    REAL (KIND=8)::a(1:n,1:n),b(1:n),x(1:n)
    REAL (KIND=8) sum
    INTEGER i,j,n
    x(n)=b(n)/a(n,n)
    DO i=n-1,1,-1
        sum=0.d0
        DO j=i+1,n,1
            sum=sum+a(i,j)*x(j)
        ENDDO
        x(i)=(b(i)-sum)/a(i,i)
    ENDDO
ENDSUBROUTINE

