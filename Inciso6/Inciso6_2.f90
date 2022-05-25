PROGRAM Simson3_8
    IMPLICIT NONE
    REAL(KIND=8) a,b,deltaX,f,fx0,fxn,xi,suma1,suma2,suma3,integracion 
    INTEGER n,i
    a=1.d0 !Límite inferior
    b=0.6d0 !Límite superior
    n=90 !Número de divisiones, debe ser múltiplo de 3
    
    deltaX=(b-a)/DBLE(n) !Incremento
    
    fx0=f(a) 
    fxn=f(b)
    
    suma1=0.d0 
    DO i=1,(n-2),3
        xi=a+deltaX*DBLE(i)
        suma1=suma1+f(xi) 
    END DO

    suma2=0.d0 
    DO i=2,(n-1),3
        xi=a+deltaX*DBLE(i)
        suma2=suma2+f(xi) 
    END DO

    suma3=0.d0 
    DO i=3,(n-3),3
        xi=a+deltaX*DBLE(i)
        suma3=suma3+f(xi) 
    END DO
    !Resultado de la integración numérica 
    Integracion=(3.d0/8.d0)*deltaX*(fx0+3.d0*suma1+3.d0*suma2+2.d0*suma3+fxn) 
    WRITE(*,'(A9,2x,f17.10)') 'Integral=',integracion
    WRITE(*,'(A2,2x,I6)') 'n=',n
    END PROGRAM
    !******f(x) 
    FUNCTION f(x)
        IMPLICIT NONE
        REAL (KIND=8) x,f 
        f=1.d0/(dsqrt(4-x**2))
    END FUNCTION