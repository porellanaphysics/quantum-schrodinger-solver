program schr
    implicit none

    integer*8 :: N, n_ciclos, t, j, unit_Inicial, unit_resultado, unit_norma, k, tiempo, unit_potencial, unit_correlacion
    real*8 :: lambda, k_0, pi, s_barra, A, omega, periodo
    real*8, allocatable :: V(:), periodos(:)
    complex*16, allocatable :: onda(:), ondaInicial(:), gammaMenos(:), A_0(:), alfa(:), beta(:), b(:), psi(:)
    real*8, allocatable :: correlaciones(:)
    integer*8 :: numPeriodos, tiempoPeriodo(2), condTiempoPeriodo(2)
    real*8 :: numeroLambda(4)
    integer*8 :: numeroN(4)

    ! Factores de entrada
    N = 500
    tiempo = 4000
    !n_ciclos = N/4.d0
    lambda = 0.001d0
    
    ! Otros valores que se calculan a posterior
    unit_Inicial = 20
    unit_resultado = 21
    unit_norma = 22
    unit_potencial = 23



    

    pi = 4.d0*atan(1.d0)
    k_0 =2* pi / 8
    print *, k_0
    s_barra = 1/(4*k_0) ! Este puede cambiar por lo que se necesite
    A = 1.d0


    ! Voy a calcular el periodo del oscilador
    omega = sqrt(lambda)
    periodo = 2*pi/omega
    
    ! Aquí defino el tamaño de los distintos vectores
    allocate(onda(N), ondaInicial(N), gammaMenos(N-1), A_0(N), alfa(N-1), beta(N-1), b(N), psi(N), V(N))
    numPeriodos = 5
    allocate(periodos(numPeriodos), correlaciones(tiempo))





    ! Aquí empieza el codigo
    open(unit=unit_Inicial, file="inicial.dat", status="unknown", action="write")
    open(unit=unit_resultado, file="onda_resultado.dat", status="unknown", action="write")
    open(unit=unit_norma, file="norma.dat", status="unknown", action="write")
    open(unit=unit_potencial, file="potencial.dat", status="unknown", action="write")
    open(unit=unit_correlacion, file="correlacion.dat", status="unknown", action="write")


    
    call construirArmonico
    call construirFuncionOnda
    call construirAlfaGammaA
    
    do k = 1, N
        write(unit_Inicial,*) abs(onda(k))
        write(unit_potencial,*) V(k)
    end do
    
    write(unit_resultado,*) (abs(onda(k))**2, k=1, N)
    write(unit_norma,*) calcularNorma()
    
    do t = 1, tiempo
        call calcularOnda
        write(unit_resultado,*) (abs(onda(k))**2, k=1, N)
        write(unit_norma,*) calcularNorma()
        write(unit_correlacion,*) correlacion()
    end do
    


    print *, "Hola"
    
    close(unit_Inicial)
    close(unit_resultado)
    close(unit_norma)
    close(unit_potencial)
    close(unit_correlacion)
    
    contains

        ! Esta subrutina sirve para definir el potencial a principio del programa
        subroutine construirPotencial
            integer*8 :: j, inicio, fin
            real*8 :: V_valor

            V_valor = lambda*(k_0**2)
            inicio = 2*N/5 ! Se cambia si se quiere cambiar el valor del pontecial
            fin = 3*N/5

            do j=1, inicio
                V(j) = 0
            end do
            do j=inicio, fin
                V(j) = V_valor
            end do
            do j=fin, N
                V(j) = 0
            end do

        end subroutine construirPotencial
        
        subroutine construirArmonico
            integer*8 :: j, x

            do j = 1, N
                x = j - N/2
                V(j) = lambda/2 * (x**2)
            end do
            
            
        end subroutine construirArmonico
        
        ! Esta rutina servirá para construir el primer potencial
        subroutine construirFuncionOnda
            integer*8 :: j
            
            onda(1) = 0
            onda(N) = 0
            do j = 2, N-1
                onda(j) = exp(cmplx(0, k_0*(j - N/2))) * exp((-8.d0) * (4.d0 * real(j - N/2, kind=8))**2 / (real(N, kind=8)**2))
            end do
            
            onda = onda/sqrt(calcularNorma())

            do j = 1, N
                ondaInicial(j) = onda(j)
            end do
        end subroutine construirFuncionOnda

        subroutine construirAlfaGammaA
            integer*8 :: j
            complex*16 :: valorComplejo
            
            valorComplejo = cmplx(0, 2.d0/s_barra)
            
            do j=1, N
                A_0(j) = -2.d0 + valorComplejo - V(j)
            end do
            
            alfa(N-1) = 0
            do j = N-1, 1, -1
                gammaMenos(j) = A_0(j) + A*alfa(j)
                alfa(j-1) = -A/gammaMenos(j)
            end do

        end subroutine construirAlfaGammaA

        subroutine calcularOnda
            integer*8 :: j
            
            do j = 1, N
                b(j) = cmplx(0, 4/s_barra) * onda(j)
            end do

            beta(N-1) = 0
            do j = N-1, 2, -1
                beta(j-1) = (b(j) - A*beta(j)) / gammaMenos(j)
            end do
            
            psi(N) = 0.d0
            psi(1) = 0.d0
            do j = 1, N-1
                psi(j+1) = alfa(j)*psi(j) + beta(j)
            end do
            
            do j = 1, N
                onda(j) = psi(j) - onda(j)
            end do

        end subroutine calcularOnda
        
        function calcularNorma() result(norma)
            integer*8 :: j
            real*8 :: norma
            
            norma = 0.d0
            do j = 1, N
                norma = norma + abs(onda(j)) ** 2
            end do

        end function calcularNorma

        function correlacion()
            integer*8 :: j
            real*8 :: correlacion

            correlacion = 0
            do j = 1, N
                correlacion = correlacion + abs(conjg(ondaInicial(j))*onda(j))
            end do

            
        end function correlacion
    


        
end program schr