program SIMSON1_3
    implicit none
    real (kind=8) a,b,f,deltaX,suma1,suma2,fxn,fx0,xi,integral
    integer n,i

    a=1.d0 !límite inferior 
    b=2.5d0 !límite superior
    
    n =150 ! División del dominio
    deltaX=(b-a)/DBLE(n) 

    fx0=f(a)
    fxn=f(b)

    suma1 = 0.d0
    suma2 = 0.d0
    !SUMA PARES
    do i=1, (n-1),2
        xi = a + deltaX*dble(i)
        suma1 = suma1 + f(xi)
    end do

    !SUMA IMPARES
    do i=2, (n-2),2
        xi = a + deltaX*dble(i)
        suma2 = suma2 + f(xi)
    end do

    ! se calcula la integral segun el método simpson 1/3
    integral = (deltaX/3.d0)*(fx0+4.d0*suma1+2.d0*suma2+fxn)
    write(*,*) 'Integral = ',integral
    write(*,*) 'n = ',n
end program

function f(x)
    implicit none
    real (kind=8) f,x

    f=((x**2)+8.d0)**(1.d0/3.d0)
end function 