program GAUSS_SIMPLE
    implicit none
    integer, parameter :: n=4 !orden del sistema de ecuaciones
    real(kind=8) :: a(1:n,1:n), b(1:n),x(1:n)
    integer i,j
    !Declaración de la matriz de coeficientes A=(a_ij)
    a(1,1)=0.d0;a(1,2)=2.d0;a(1,3)=-1.d0;a(1,4)=-4.d0
    a(2,1)=1.d0;a(2,2)=-1.d0;a(2,3)=5.d0;a(2,4)=2.d0
    a(3,1)=3.d0;a(3,2)=3.d0;a(3,3)=-7.d0;a(3,4)=-1.d0
    a(4,1)=-1.d0;a(4,2)=-2.d0;a(4,3)=3.d0;a(4,4)=0.d0
    !Declaración de vector de coeficientes b=b(i)
    b(1)=2.d0;b(2)=-4.d0;b(3)=4.d0;b(4)=-7.d0

    !Solo para visualizar los datos de entrada
    write(*,*) 'Matriz A'
    do i=1,n,1
        write(*,'(3(f20.16,2x))') (a(i,j),j=1,2) !Do implicito
    end do 

    write(*,*)'Vetro de coeficientes B'
    
    do i=1,n,1
        write(*,'(f13.8,2x)') b(i)
    end do 

    !Eliminación Gaussiana simple del sistma de ecuaciones a la forma triangular
    call eliminacion_gaussiana(a,b,n)
    !Obtención del vector solución X=x(i) a partir de la sustitución hacia atras
    CALL sustitucion_atras(a,b,x,n)
    !Solo para visualizar los datos de entrada

    WRITE(*,*)'Soluciones  X=x(i)'
     DO i=1,n,1
        WRITE(*,'(2x,A,I1,A2,f13.8)') 'x',i,'= ',x(i)  !DO IMPLICITO
    END DO
end program GAUSS_SIMPLE


!****** Eliminación Gaussiana ******!
SUBROUTINE eliminacion_gaussiana(a,b,n)
    IMPLICIT NONE 
    REAL (KIND=8):: a(1:n,1:n),b(1:n)
    REAL (KIND=8):: pivote
    INTEGER k,i,j,n
    DO k=1,n-1,1 !desplaza las operaciones al k-ésimo renglón, esto se nota para i=k+1, que realmente i[2,n]
        DO i=k+1,n
        pivote=a(i,k)/a(k,k)
            DO j=k+1,n,1 !
                a(i,j)=a(i,j)-pivote*a(k,j)
            END DO
        b(i)=b(i)-pivote*b(k)
        END DO
    END DO
END SUBROUTINE

!********************Subrutina para sustitución hacia atrás
SUBROUTINE sustitucion_atras(a,b,x,n)
    IMPLICIT NONE
    REAL (KIND=8)::a(1:n,1:n),b(1:n),x(1:n)
    REAL (KIND=8) sum
    INTEGER i,j,n
    x(n)=b(n)/a(n,n)
    DO i=n-1,1,-1
        sum=b(i)
        DO j=i+1,n,1
            sum=sum-a(i,j)*x(j)
        END DO
        x(i)=sum/a(i,i)
    END DO
END SUBROUTINE

