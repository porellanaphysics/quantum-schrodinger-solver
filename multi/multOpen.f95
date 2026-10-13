program schr
    use omp_lib
    implicit none

    integer*8 :: N, n_ciclos, t, j, k, pasos, numeroN(10), l, m
    real*8 :: lambda, k_0, pi, s_barra, A, omega, periodo, numeroCiclos(4), tiempo
    real*8, allocatable :: V(:)
    complex*16, allocatable :: onda(:), ondaInicial(:), gammaMenos(:), A_0(:), alfa(:), beta(:), b(:), psi(:)
    real*8, allocatable :: correlaciones(:)
    integer*8 :: numPeriodos, tiempoPeriodo(2), condTiempoPeriodo(2)
    character(len=60) :: nombreArchivo

    ! Inicializar OpenMP
    integer :: num_threads
    num_threads = omp_get_max_threads()
    print *, "Usando", num_threads, "hilos OpenMP"

    numeroN = (/500, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000/)
    numeroCiclos = (/4, 6, 8, 10/)

    lambda = 0.001d0
    pi = 4.d0*atan(1.d0)
    pasos = 3000
    A = 1.d0

    do l=1, size(numeroN)
        do m = 1, size(numeroCiclos)
            N = numeroN(l)
            k_0 = 2*pi / numeroCiclos(m)
            s_barra = 1/(4*k_0) 
    
            ! Verifica y desasigna los arreglos si ya están asignado
            ! Asigna memoria para los arreglos dinámicos
            allocate(onda(N), ondaInicial(N), gammaMenos(N-1), A_0(N), alfa(N-1), beta(N-1), b(N), psi(N), V(N))

            call construirArmonico
            call construirFuncionOnda
            call construirAlfaGammaA

            ! Abrir archivo de salida para este caso
            write(nombreArchivo, '(A,I0,A,I0,A)') 'resultados_N', N, '_ciclos', nint(numeroCiclos(m)), '.dat'
            open(unit=20, file=trim(nombreArchivo), status='replace')
    
            do t = 1, pasos
                call calcularOnda
                tiempo = s_barra*t
                write(20,*) tiempo, correlacion()
            end do
            
            ! Desasigna los arreglos al final de cada iteración
            deallocate(onda, ondaInicial, gammaMenos, A_0, alfa, beta, b, psi, V)
            close(20)
        end do
    end do

contains

    subroutine construirArmonico
        integer*8 :: j, x
        
        !$OMP PARALLEL DO PRIVATE(x)
        do j = 1, N
            x = j - N/2
            V(j) = lambda/2 * (x**2)
        end do
        !$OMP END PARALLEL DO
        
    end subroutine construirArmonico
    
    subroutine construirFuncionOnda
        integer*8 :: j
        
        onda(1) = 0
        onda(N) = 0
        
        !$OMP PARALLEL DO
        do j = 2, N-1
            onda(j) = exp(cmplx(0, k_0*(j - N/2))) * exp((-8.d0) * (4.d0 * real(j - N/2, kind=8))**2 / (real(N, kind=8)**2))
        end do
        !$OMP END PARALLEL DO
        
        onda = onda/sqrt(calcularNorma())

        !$OMP PARALLEL DO
        do j = 1, N
            ondaInicial(j) = onda(j)
        end do
        !$OMP END PARALLEL DO
    end subroutine construirFuncionOnda

    subroutine construirAlfaGammaA
        integer*8 :: j
        complex*16 :: valorComplejo
        
        valorComplejo = cmplx(0, 2.d0/s_barra)
        
        !$OMP PARALLEL DO
        do j=1, N
            A_0(j) = -2.d0 + valorComplejo - V(j)
        end do
        !$OMP END PARALLEL DO
        
        alfa(N-1) = 0
        ! Este bucle no se puede paralelizar fácilmente por dependencias
        do j = N-1, 2, -1
            gammaMenos(j) = A_0(j) + A*alfa(j)
            alfa(j-1) = -A/gammaMenos(j)
        end do
    end subroutine construirAlfaGammaA

    subroutine calcularOnda
        integer*8 :: j
        
        !$OMP PARALLEL DO
        do j = 1, N
            b(j) = cmplx(0, 4/s_barra) * onda(j)
        end do
        !$OMP END PARALLEL DO

        beta(N-1) = 0.d0
        ! Este bucle no se puede paralelizar fácilmente por dependencias
        do j = N-1, 2, -1
            beta(j-1) = (b(j) - A*beta(j)) / gammaMenos(j)
        end do
        
        psi(N) = 0.d0
        psi(1) = 0.d0
        ! Este bucle no se puede paralelizar fácilmente por dependencias
        do j = 1, N-1
            psi(j+1) = alfa(j)*psi(j) + beta(j)
        end do
        
        !$OMP PARALLEL DO
        do j = 1, N
            onda(j) = psi(j) - onda(j)
        end do
        !$OMP END PARALLEL DO
    end subroutine calcularOnda
    
    function calcularNorma() result(norma)
        integer*8 :: j
        real*8 :: norma, suma_local
        
        norma = 0.d0
        suma_local = 0.d0
        
        !$OMP PARALLEL PRIVATE(j) REDUCTION(+:suma_local)
        !$OMP DO
        do j = 1, N
            suma_local = suma_local + abs(onda(j)) ** 2
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
        norma = suma_local
    end function calcularNorma

    function correlacion()
        integer*8 :: j
        real*8 :: correlacion, suma_local
        
        correlacion = 0.d0
        suma_local = 0.d0
        2
        !$OMP PARALLEL PRIVATE(j) REDUCTION(+:suma_local)
        !$OMP DO
        do j = 1, N
            suma_local = suma_local + abs(onda(j) * conjg(ondaInicial(j))) ** 2
        end do
        !$OMP END DO
        !$OMP END PARALLEL
        
        correlacion = suma_local / calcularNorma()
    end function correlacion
    
end program schr